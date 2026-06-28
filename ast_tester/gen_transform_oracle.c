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

/* Pick the sampling range for one axis of one fixture.  FITS-header
   mappings sample the pixel grid 1..N (approximated by a fixed generous
   pixel range here, since NAXISn is not exposed on the bare Mapping);
   native-dump mappings use a fixed symmetric range. */
static void axis_range( int is_head, double *lo, double *hi ) {
    if ( is_head ) { *lo = 1.0;     *hi = 2000.0;  }    /* pixel-ish */
    else           { *lo = -1000.0; *hi =  1000.0; }
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

/* Transform npoint points through map in the given direction and write a
   section to fp.  in_n / out_n are the arities for that direction.  The
   transform runs first: if it raises an AST error (e.g. a map that cannot
   evaluate the sampled domain) the status is cleared and no section is
   written.  Returns 1 if a section was written. */
static int write_section( FILE *fp, const char *relpath, AstMapping *map,
                          int forward, int in_n, int out_n, int is_head ) {
    int np = oracle_sample_axis_count();

    double **in  = malloc( sizeof(double *) * (size_t) in_n );
    double **out = malloc( sizeof(double *) * (size_t) out_n );
    double *lo = malloc( sizeof(double) * (size_t) in_n );
    double *hi = malloc( sizeof(double) * (size_t) in_n );
    for ( int a = 0; a < in_n;  a++ ) { in[a]  = malloc(sizeof(double)*(size_t)np);
                                        axis_range(is_head,&lo[a],&hi[a]); }
    for ( int a = 0; a < out_n; a++ )   out[a] = malloc(sizeof(double)*(size_t)np);

    oracle_sample_points( in_n, lo, hi, np, in );
    astTranP( map, np, in_n, (const double **) in, forward, out_n, out );

    int wrote = 0;
    if ( astOK ) {
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
        wrote = 1;
    } else {
        astClearStatus;
        fprintf( stderr, "skip %s (transform raised an error)\n", relpath );
    }

    for ( int a = 0; a < in_n;  a++ ) free( in[a] );
    for ( int a = 0; a < out_n; a++ ) free( out[a] );
    free( in ); free( out ); free( lo ); free( hi );
    return wrote;
}

/* Emit one fixture (forward if defined, else inverse-only). Returns 1 if
   a section was written. */
static int emit_fixture( FILE *fp, const char *root, const char *relpath ) {
    int wrote = 0;
    /* The simplify fixtures' "simplifyidentity" IntraMap is a private test
       extension that cannot be loaded without registration; skip any
       fixture that uses an IntraMap. */
    if ( file_contains( root, relpath, "IntraMap" ) ) return 0;

    astBegin;
    AstMapping *map = oracle_load_mapping( root, relpath );
    if ( map ) {
        int nin  = astGetI( map, "Nin" );
        int nout = astGetI( map, "Nout" );
        int is_head = ( strstr( relpath, ".head" ) != NULL );
        if ( astGetI( map, "TranForward" ) ) {
            wrote = write_section( fp, relpath, map, 1, nin, nout, is_head );
        } else if ( astGetI( map, "TranInverse" ) ) {
            wrote = write_section( fp, relpath, map, 0, nout, nin, is_head );
        }
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
