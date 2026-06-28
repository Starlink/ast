/*
 * gen_transform_oracle: scan the transform fixture corpus, sample input
 * points, transform them, and write two text oracle files used by
 * check_transform_oracle as a tolerance-based regression reference.
 *
 * Usage:
 *   gen_transform_oracle <root> <simplify_oracle_out> <headers_oracle_out>
 *
 * <root> is the ast_tester source directory.  This program is the only
 * component that scans the filesystem; the checker is driven entirely by
 * the files written here.
 */
#include "ast.h"
#include "transform_oracle.h"

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double **alloc_cols( int ncol, int nrow ) {
    double **c = malloc( sizeof(double *) * (size_t) ncol );
    for ( int a = 0; a < ncol; a++ ) c[a] = malloc( sizeof(double) * (size_t) nrow );
    return c;
}

static void free_cols( double **c, int ncol ) {
    if ( !c ) return;
    for ( int a = 0; a < ncol; a++ ) free( c[a] );
    free( c );
}

/* Return 1 if <root>/<relpath> can be opened for reading. */
static int file_exists( const char *root, const char *relpath ) {
    char path[1024];
    snprintf( path, sizeof path, "%s/%s", root, relpath );
    FILE *fp = fopen( path, "r" );
    if ( fp ) { fclose( fp ); return 1; }
    return 0;
}

/* Return 1 if any line of <root>/<relpath> contains needle. */
static int file_contains( const char *root, const char *relpath,
                          const char *needle ) {
    char path[1024];
    snprintf( path, sizeof path, "%s/%s", root, relpath );
    FILE *fp = fopen( path, "r" );
    if ( !fp ) return 0;
    int found = 0;
    char line[1024];
    while ( fgets( line, (int) sizeof line, fp ) ) {
        if ( strstr( line, needle ) ) { found = 1; break; }
    }
    fclose( fp );
    return found;
}

/* Fill lo[]/hi[] for the naxis input axes of a FITS-header fixture from
   its NAXISj keywords; axes with no NAXISj (degenerate or extra WCS axes)
   fall back to a small [1,2] range. */
static void head_bounds( const char *root, const char *relpath,
                         int naxis, double *lo, double *hi ) {
    char path[1024];
    snprintf( path, sizeof path, "%s/%s", root, relpath );
    AstFitsChan *fc = astFitsChan( NULL, NULL, " " );
    FILE *fp = fopen( path, "r" );
    if ( fp ) {
        char line[256];
        while ( fgets( line, (int) sizeof line, fp ) ) {
            size_t n = strlen( line );
            while ( n > 0 && ( line[n-1] == '\n' || line[n-1] == '\r' ) )
                line[--n] = '\0';
            astPutFits( fc, line, 0 );
        }
        fclose( fp );
    }
    for ( int a = 0; a < naxis; a++ ) {
        char key[16]; int val = 0;
        snprintf( key, sizeof key, "NAXIS%d", a + 1 );
        lo[a] = 1.0;
        if ( astGetFitsI( fc, key, &val ) && val >= 1 ) hi[a] = (double) val;
        else hi[a] = 2.0;
    }
    fc = astAnnul( fc );
    if ( !astOK ) astClearStatus;
}

/* Fill lo[]/hi[] for naxis input axes.  FITS headers use their pixel grid
   (NAXISj); native dumps have no declared domain, so use a fixed symmetric
   range.  Sampling outside a map's valid domain merely yields AST__BAD,
   which the golden comparison records and matches faithfully. */
static void axis_bounds( const char *root, const char *relpath, int is_head,
                         int naxis, double *lo, double *hi ) {
    if ( is_head ) {
        head_bounds( root, relpath, naxis, lo, hi );
    } else {
        for ( int a = 0; a < naxis; a++ ) { lo[a] = -1000.0; hi[a] = 1000.0; }
    }
}

/* Transform pre-sampled inputs `in` (in_n x np) through map in the given
   direction and write a section.  The transform runs first: if it raises
   an AST error the status is cleared and no section is written.  On
   success `out` (out_n x np) holds the outputs.  Returns 1 if written. */
static int transform_write( FILE *fp, const char *relpath, AstMapping *map,
                            int forward, int in_n, int out_n, int np,
                            double **in, double **out ) {
    astTranP( map, np, in_n, (const double **) in, forward, out_n, out );
    if ( !astOK ) {
        astClearStatus;
        fprintf( stderr, "skip %s dir=%s (transform raised an error)\n",
                 relpath, forward ? "forward" : "inverse" );
        return 0;
    }
    fprintf( fp, "[%s  nin=%d nout=%d dir=%s]\n",
             relpath, in_n, out_n, forward ? "forward" : "inverse" );
    char buf[64];
    for ( int p = 0; p < np; p++ ) {
        for ( int a = 0; a < in_n;  a++ ) {
            oracle_format_double( buf, sizeof buf, in[a][p] );
            fprintf( fp, "%s%s", a ? " " : "  ", buf );
        }
        fprintf( fp, "  " );
        for ( int a = 0; a < out_n; a++ ) {
            oracle_format_double( buf, sizeof buf, out[a][p] );
            fprintf( fp, "%s%s", a ? " " : "", buf );
        }
        fprintf( fp, "\n" );
    }
    fprintf( fp, "\n" );
    return 1;
}

/* Emit a fixture's golden sections: a forward section (if the forward
   transform is defined) and an inverse section (if the inverse is
   defined).  The inverse section's inputs are the forward outputs when a
   forward exists, so they are genuine output-space coordinates (the
   inverse's natural domain); otherwise they are sampled directly.  Both
   are recorded as golden, so neither needs the inputs to be "valid" -- a
   change in either kernel's output is detected regardless.  Returns the
   number of sections written. */
static int emit_fixture( FILE *fp, const char *root, const char *relpath ) {
    int wrote = 0;
    /* The simplify fixtures' "simplifyidentity" IntraMap is a private test
       extension that cannot be loaded without registration; skip any
       fixture that uses an IntraMap. */
    if ( file_contains( root, relpath, "IntraMap" ) ) return 0;

    astBegin;
    AstMapping *map = oracle_load_mapping( root, relpath );
    if ( map ) {
        int nin     = astGetI( map, "Nin" );
        int nout    = astGetI( map, "Nout" );
        int is_head = ( strstr( relpath, ".head" ) != NULL );
        int np      = oracle_sample_axis_count();
        int has_fwd = astGetI( map, "TranForward" );
        int has_inv = astGetI( map, "TranInverse" );

        double **fwd_out = NULL;   /* forward outputs, reused as inverse inputs */

        if ( has_fwd ) {
            double *lo = malloc( sizeof(double) * (size_t) nin );
            double *hi = malloc( sizeof(double) * (size_t) nin );
            axis_bounds( root, relpath, is_head, nin, lo, hi );
            double **in  = alloc_cols( nin, np );
            double **out = alloc_cols( nout, np );
            oracle_sample_points( nin, lo, hi, np, in );
            if ( transform_write( fp, relpath, map, 1, nin, nout, np, in, out ) ) {
                wrote++;
                fwd_out = out;            /* keep for inverse inputs */
            } else {
                free_cols( out, nout );
            }
            free_cols( in, nin );
            free( lo ); free( hi );
        }

        if ( has_inv ) {
            /* inverse takes nout-dim input -> nin-dim output */
            double **inv_in;
            int own_inv_in;
            if ( fwd_out ) { inv_in = fwd_out; own_inv_in = 0; }
            else {
                double *lo = malloc( sizeof(double) * (size_t) nout );
                double *hi = malloc( sizeof(double) * (size_t) nout );
                axis_bounds( root, relpath, is_head, nout, lo, hi );
                inv_in = alloc_cols( nout, np );
                oracle_sample_points( nout, lo, hi, np, inv_in );
                free( lo ); free( hi );
                own_inv_in = 1;
            }
            double **inv_out = alloc_cols( nin, np );
            if ( transform_write( fp, relpath, map, 0, nout, nin, np,
                                  inv_in, inv_out ) )
                wrote++;
            free_cols( inv_out, nin );
            if ( own_inv_in ) free_cols( inv_in, nout );
        }

        free_cols( fwd_out, nout );
    }
    if ( !astOK ) astClearStatus;
    astEnd;
    return wrote;
}

/* Comparator for qsort of strings, for deterministic file order. */
static int cmp_str( const void *a, const void *b ) {
    return strcmp( *(const char *const *) a, *(const char *const *) b );
}

/* Collect entries of dir matching suffix into a sorted relpath list.
   prefix is prepended to each name (e.g. "simplify_fixtures/" or ""). */
static char **scan_dir( const char *root, const char *subdir,
                        const char *prefix, const char *suffix, int *count ) {
    char path[1024];
    snprintf( path, sizeof path, "%s%s%s", root, subdir[0] ? "/" : "", subdir );
    DIR *d = opendir( path );
    *count = 0;
    if ( !d ) return NULL;
    size_t cap = 64, n = 0;
    char **list = malloc( sizeof(char *) * cap );
    struct dirent *e;
    size_t slen = strlen( suffix );
    while ( ( e = readdir( d ) ) ) {
        size_t nlen = strlen( e->d_name );
        if ( nlen <= slen || strcmp( e->d_name + nlen - slen, suffix ) != 0 )
            continue;
        if ( n == cap ) { cap *= 2; list = realloc( list, sizeof(char *) * cap ); }
        char rel[1024];
        snprintf( rel, sizeof rel, "%s%s", prefix, e->d_name );
        list[n++] = strdup( rel );
    }
    closedir( d );
    qsort( list, n, sizeof(char *), cmp_str );
    *count = (int) n;
    return list;
}

int main( int argc, char *argv[] ) {
    int status_value = 0;
    int *status = &status_value;
    if ( argc < 4 ) {
        fprintf( stderr,
            "Usage: gen_transform_oracle <root> <simplify_out> <headers_out>\n" );
        return 1;
    }
    const char *root         = argv[1];
    const char *simplify_out = argv[2];
    const char *headers_out  = argv[3];
    astWatch( status );

    const char *hdr =
        "# transform oracle - regenerate with gen_transform_oracle\n"
        "# tol: rtol=1e-12 atol=1e-12  equiv_rtol=1e-9 equiv_atol=1e-9\n\n";

    /* Native-dump corpus: .map and .simp interleaved by stem so paired
       sections sit together in the file. */
    FILE *fs = fopen( simplify_out, "w" );
    if ( !fs ) { fprintf( stderr, "cannot write %s\n", simplify_out ); return 1; }
    fputs( hdr, fs );
    int nmap = 0;
    char **maps = scan_dir( root, "simplify_fixtures", "simplify_fixtures/",
                            ".map", &nmap );
    int nwritten = 0;
    for ( int i = 0; i < nmap; i++ ) {
        char simp[1024];
        size_t L = strlen( maps[i] );
        snprintf( simp, sizeof simp, "%.*s.simp", (int)(L - 4), maps[i] );
        nwritten += emit_fixture( fs, root, maps[i] );
        /* "neg_*" fixtures (and any other already-simple map) have no .simp
           partner; only emit one when the file exists. */
        if ( file_exists( root, simp ) ) nwritten += emit_fixture( fs, root, simp );
        free( maps[i] );
    }
    free( maps );
    fclose( fs );

    /* FITS-header corpus. */
    FILE *fh = fopen( headers_out, "w" );
    if ( !fh ) { fprintf( stderr, "cannot write %s\n", headers_out ); return 1; }
    fputs( hdr, fh );
    int nhead = 0;
    char **heads = scan_dir( root, "", "", ".head", &nhead );
    for ( int i = 0; i < nhead; i++ ) {
        nwritten += emit_fixture( fh, root, heads[i] );
        free( heads[i] );
    }
    free( heads );
    fclose( fh );

    printf( "gen_transform_oracle: wrote %d sections\n", nwritten );
    return astOK ? 0 : 1;
}
