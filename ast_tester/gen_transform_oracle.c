/*
 * gen_transform_oracle: scan the transform fixture corpus, sample input
 * points, transform them, and write two text oracle files used by
 * check_transform_oracle as a tolerance-based regression reference.
 *
 * Usage:
 *   gen_transform_oracle <root> <simplify_oracle_out> <headers_oracle_out>
 *
 * <root> is the ast_tester/fixtures directory.  Paths written in section
 * headers are relative to this root.  This program is the only
 * component that scans the filesystem; the checker is driven entirely by
 * the files written here.
 *
 * SAMPLING DOMAINS
 *
 * A FITS header is sampled over its NAXIS pixel grid, and a FrameSet dump
 * over the grid of the header it was read from (see axis_bounds).  A native
 * dump under simplify/ carries no domain, so by default its inputs take
 * either sign with magnitudes spread logarithmically from 1e-3 to 1e3: an
 * angle-based Mapping is then sampled within its principal range as well as
 * far outside it, and a linear one over many scales.  A fixture that is
 * defined, or is the identity, only over some other range is given that
 * range in <root>/oracle/transform_oracle_domains.txt (see load_domains).
 * Its .simp partner is always sampled at the same points as its .map, since
 * the checker compares the two row by row.
 *
 * PLATFORM SENSITIVITY OF THE RECORDED VALUES
 *
 * The oracle files record whatever the generating platform computes, and
 * two effects make that platform-dependent:
 *
 * 1. Fused multiply-add contraction.  Compilers may contract a*b + c into
 *    a single fused instruction that skips rounding the product.  Arm64
 *    has FMA in its base ISA and compilers contract by default; baseline
 *    x86-64 has no FMA instruction, so the same expression rounds twice
 *    (an x86 build with -march=haswell or later behaves like arm64, so
 *    this is an instruction-set property, not an architecture split).
 *    The difference is normally ~1 ulp and vanishes inside the golden
 *    tolerance, but at a mathematical singularity it changes the answer
 *    in kind: e.g. UnitNormMap's inverse computes in*norm + centre, and
 *    for a fixture whose inverse chain reconstructs another UnitNormMap's
 *    centre, separate rounding lands on the centre exactly (forward
 *    output AST__BAD) while FMA preserves a ~1e-16 residual whose
 *    direction is pure rounding noise (forward output an arbitrary unit
 *    vector).  No tolerance can bridge BAD vs a value, so such rows can
 *    never be recorded stably.  Per-fixture "golden off" lines in the
 *    overrides file (see transform_oracle_overrides.txt) are the escape
 *    hatch, unless the build forces consistent floating point math
 *    everywhere (-ffp-contract=off for libast AND these test programs
 *    - the Halton sampling below contracts too), which would pin every
 *    pure-arithmetic value bit-for-bit across platforms at the cost of
 *    diverging from builds without the flag.
 *
 * 2. libm differences.  Transcendental functions (sin, cos, atan2, ...)
 *    are not correctly-rounded and differ between C libraries (glibc vs
 *    Apple libm) by up to ~5e-9 relative near singular points.  This is
 *    why the golden tolerance is 1e-7 rather than a few ulp; see
 *    transform_oracle.h.  No compiler flag removes this class.
 */
#include "ast.h"
#include "transform_oracle.h"

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

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

/* A fixture missing from an oracle silently stops being tested, so abort
   rather than carry on with a short manifest. */
static void *xalloc( void *p, size_t size, const char *what ) {
    void *q = realloc( p, size );
    if ( !q ) {
        fprintf( stderr, "gen_transform_oracle: out of memory for %s\n", what );
        exit( 1 );
    }
    return q;
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

/* Sampling-domain kind for a fixture, chosen from its file extension. */
enum { DOM_NATIVE, DOM_HEAD, DOM_FRAMESET };

static int path_kind( const char *relpath ) {
    if ( strstr( relpath, ".head" ) ) return DOM_HEAD;
    if ( strstr( relpath, ".ast"  ) ) return DOM_FRAMESET;
    return DOM_NATIVE;
}

/* Magnitudes between which a native dump with no declared domain is
   sampled; see oracle_sample_points_log. */
#define NATIVE_MIN_ABS 1e-3
#define NATIVE_MAX_ABS 1e3

/* A declared sampling domain for a simplify fixture: one [lo, hi] range for
   every input axis (naxis == 1) or one per input axis. */
typedef struct {
    char relpath[1024];         /* the .map, relative to root */
    int naxis;
    double *lo, *hi;
    int used;                   /* matched a scanned fixture */
} Domain;

static void domain_fail( const char *path, int lineno, const char *why ) {
    fprintf( stderr, "gen_transform_oracle: %s:%d: %s\n", path, lineno, why );
    exit( 1 );
}

/* Ranges that a domain entry can give by name, for inputs whose valid
   values are fixed by what they are rather than by the fixture.  A
   longitude covers a full turn; the checker compares angles modulo 2*pi,
   so the two ends are the same position.  A latitude stops short of the
   poles: longitude is undefined at a pole and ill-conditioned close to
   one, so a Mapping and its simplification may legitimately return
   different longitudes there for the same position. */
static const struct {
    const char *name;
    double lo, hi;
} named_ranges[] = {
    { "lon", -3.14159265358979323846, 3.14159265358979323846 },
    { "lat", -1.5, 1.5 },
};

/* Read <root>/oracle/transform_oracle_domains.txt.  Each line that is not
   blank or a '#' comment is

     <relpath> <range>                   one range for every input axis
     <relpath> <range1> ... <rangeN>     one range per input axis

   where <relpath> names a simplify .map fixture and each <range> is either
   two numbers, <lo> <hi>, or the name of one of named_ranges.  The fixture
   and its .simp partner are then sampled uniformly over the ranges instead
   of over the default logarithmic spread.  An axis whose <lo> equals its
   <hi> is held at that value.  The file must exist, since without it every
   declared domain would silently revert to the default; a malformed line, a
   duplicate, or an entry that names no fixture stops the program. */
static Domain *load_domains( const char *root, int *ndom ) {
    char path[1024];
    snprintf( path, sizeof path, "%s/oracle/transform_oracle_domains.txt", root );
    FILE *fp = fopen( path, "r" );
    if ( !fp ) {
        fprintf( stderr, "gen_transform_oracle: cannot read %s\n", path );
        exit( 1 );
    }
    Domain *doms = NULL;
    int n = 0, cap = 0, lineno = 0;
    char line[4096];
    while ( fgets( line, (int) sizeof line, fp ) ) {
        lineno++;
        char *p = line;
        while ( *p == ' ' || *p == '\t' ) p++;
        if ( *p == '#' || *p == '\n' || *p == '\0' ) continue;

        char rel[1024];
        int used;
        if ( sscanf( p, "%1023s%n", rel, &used ) != 1 )
            domain_fail( path, lineno, "no fixture path" );
        p += used;
        double v[64];
        int nv = 0;
        int half = 0;       /* a <lo> has been read and its <hi> is due */
        for ( ;; ) {
            char *end;
            while ( *p == ' ' || *p == '\t' ) p++;
            if ( *p == '\n' || *p == '\0' || *p == '#' ) break;
            if ( nv + 2 > (int)( sizeof v / sizeof v[0] ) )
                domain_fail( path, lineno, "too many ranges" );
            double x = strtod( p, &end );
            if ( end != p ) {
                v[nv++] = x;
                half = !half;
                p = end;
                continue;
            }
            char name[32];
            int len = 0;
            while ( p[len] && p[len] != ' ' && p[len] != '\t' &&
                    p[len] != '\n' && p[len] != '#' ) len++;
            if ( half || len >= (int) sizeof name )
                domain_fail( path, lineno, "value is not a number" );
            memcpy( name, p, (size_t) len );
            name[len] = '\0';
            size_t k;
            for ( k = 0; k < sizeof named_ranges / sizeof named_ranges[0]; k++ )
                if ( !strcmp( name, named_ranges[k].name ) ) break;
            if ( k == sizeof named_ranges / sizeof named_ranges[0] )
                domain_fail( path, lineno, "unknown range name" );
            v[nv++] = named_ranges[k].lo;
            v[nv++] = named_ranges[k].hi;
            p += len;
        }
        if ( nv == 0 || half )
            domain_fail( path, lineno,
                         "need a range for one axis or each axis" );
        for ( int i = 0; i < nv; i += 2 )
            if ( !( v[i] <= v[i+1] ) )
                domain_fail( path, lineno, "<lo> is above <hi>" );
        size_t rl = strlen( rel );
        if ( strncmp( rel, "simplify/", 9 ) || rl < 4 ||
             strcmp( rel + rl - 4, ".map" ) )
            domain_fail( path, lineno,
                         "a domain can only be declared for a simplify/*.map" );
        for ( int d = 0; d < n; d++ )
            if ( !strcmp( doms[d].relpath, rel ) )
                domain_fail( path, lineno, "fixture already has a domain" );

        if ( n == cap ) {
            cap = cap ? 2 * cap : 32;
            doms = xalloc( doms, sizeof(Domain) * (size_t) cap, "domain list" );
        }
        Domain *d = &doms[n++];
        snprintf( d->relpath, sizeof d->relpath, "%s", rel );
        d->naxis = nv / 2;
        d->lo = xalloc( NULL, sizeof(double) * (size_t) d->naxis, "domain" );
        d->hi = xalloc( NULL, sizeof(double) * (size_t) d->naxis, "domain" );
        for ( int a = 0; a < d->naxis; a++ ) {
            d->lo[a] = v[2*a];
            d->hi[a] = v[2*a+1];
        }
        d->used = 0;
    }
    fclose( fp );
    *ndom = n;
    return doms;
}

static Domain *find_domain( Domain *doms, int ndom, const char *relpath ) {
    for ( int d = 0; d < ndom; d++ )
        if ( !strcmp( doms[d].relpath, relpath ) ) return &doms[d];
    return NULL;
}

/* Fill lo[]/hi[] for naxis input axes of a FITS header or FrameSet dump.
   FITS headers use their pixel grid (NAXISj).  FrameSet dumps (.ast) have
   a GRID base frame but carry no NAXIS, so they take the grid of the header
   they were read from where it can be found, or a default positive pixel
   range.  Sampling outside a map's valid domain merely yields AST__BAD,
   which the golden comparison records and matches faithfully. */
static void axis_bounds( const char *root, const char *relpath, int kind,
                         int naxis, double *lo, double *hi ) {
    if ( kind == DOM_HEAD ) {
        head_bounds( root, relpath, naxis, lo, hi );
    } else {
/* A fixture under wcsconv/framesets/ is the dump of the object read from one
   wcsconv input, and its name is that input's filename with each "." replaced by
   "-".  Recover the input and take its pixel grid, so the mapping is sampled over
   the domain it was built for.

   Sampling wider is not harmless for every mapping.  A SIP distortion is a
   polynomial fitted over the image, and its inverse is iterative: outside the
   image the iteration need not converge, so a default 1..1000 range over a
   256-pixel image reports round-trip errors of 100 pixels and more that say
   nothing about the library.  The same mappings round-trip cleanly when sampled
   over their own NAXIS grid, which is what the .head sections already do. */
        char origin[ 1024 ];
        const char *base = strrchr( relpath, '/' );
        base = base ? base + 1 : relpath;
        if ( !strncmp( relpath, "wcsconv/framesets/", 18 ) &&
             snprintf( origin, sizeof origin, "wcsconv/inputs/%s", base )
             < (int) sizeof origin ) {
            char *dot = strrchr( origin, '.' );      /* drop the ".ast" */
            if ( dot ) *dot = 0;
            char *dash = strrchr( origin, '-' );     /* the "." that became "-" */
            if ( dash ) *dash = '.';
        } else {
            origin[ 0 ] = 0;
        }

        if ( origin[ 0 ] && strstr( origin, ".head" ) &&
             file_exists( root, origin ) ) {
            head_bounds( root, origin, naxis, lo, hi );
        } else {
            for ( int a = 0; a < naxis; a++ ) { lo[a] = 1.0; hi[a] = 1000.0; }
        }
    }
}

/* Fill in (naxis x np) with the inputs at which to sample a fixture in one
   direction.  A declared domain, given only for the direction it describes,
   is sampled uniformly; any other native dump takes the default logarithmic
   spread, and a FITS header or FrameSet dump its pixel grid. */
static void sample_inputs( const char *root, const char *relpath, int kind,
                           const Domain *dom, int naxis, int np,
                           double **in ) {
    if ( kind == DOM_NATIVE && !dom ) {
        oracle_sample_points_log( naxis, NATIVE_MIN_ABS, NATIVE_MAX_ABS, np, in );
        return;
    }
    double *lo = xalloc( NULL, sizeof(double) * (size_t) naxis, "bounds" );
    double *hi = xalloc( NULL, sizeof(double) * (size_t) naxis, "bounds" );
    if ( dom ) {
        for ( int a = 0; a < naxis; a++ ) {
            int d = dom->naxis == 1 ? 0 : a;
            lo[a] = dom->lo[d];
            hi[a] = dom->hi[d];
        }
    } else {
        axis_bounds( root, relpath, kind, naxis, lo, hi );
    }
    oracle_sample_points( naxis, lo, hi, np, in );
    free( lo ); free( hi );
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
    /* Right-align every value in a field wide enough for the largest %.17g
       rendering (e.g. "-1.7976931348623157e+308") so the columns of a
       section line up for human readers. */
    enum { COLW = 24 };
    char buf[64];
    for ( int p = 0; p < np; p++ ) {
        for ( int a = 0; a < in_n;  a++ ) {
            oracle_format_double( buf, sizeof buf, in[a][p] );
            fprintf( fp, "%s%*s", a ? " " : "  ", COLW, buf );
        }
        fprintf( fp, "  " );
        for ( int a = 0; a < out_n; a++ ) {
            oracle_format_double( buf, sizeof buf, out[a][p] );
            fprintf( fp, "%s%*s", a ? " " : "", COLW, buf );
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
   change in either kernel's output is detected regardless.

   dom is the declared domain of the fixture's input axes, or NULL.  A .simp
   is passed the domain of its .map so that the two are sampled at the same
   points.  Returns the number of sections written. */
static int emit_fixture( FILE *fp, const char *root, const char *relpath,
                         const Domain *dom ) {
    int wrote = 0;
    astBegin;
    AstMapping *map = oracle_load_mapping( root, relpath, NULL, NULL );
    if ( map ) {
        int nin     = astGetI( map, "Nin" );
        int nout    = astGetI( map, "Nout" );
        int kind    = path_kind( relpath );
        int np      = oracle_sample_axis_count();
        int has_fwd = astGetI( map, "TranForward" );
        int has_inv = astGetI( map, "TranInverse" );

        /* A domain describes the forward inputs, so it is an error for a
           fixture that has none or has a different number of them. */
        if ( dom && ( !has_fwd || ( dom->naxis != 1 && dom->naxis != nin ) ) ) {
            fprintf( stderr, "gen_transform_oracle: the domain declared for "
                     "%s does not describe the forward inputs of %s "
                     "(%d ranges, Nin=%d, TranForward=%d)\n", dom->relpath,
                     relpath, dom->naxis, nin, has_fwd );
            exit( 1 );
        }

        double **fwd_out = NULL;   /* forward outputs, reused as inverse inputs */

        if ( has_fwd ) {
            double **in  = alloc_cols( nin, np );
            double **out = alloc_cols( nout, np );
            sample_inputs( root, relpath, kind, dom, nin, np, in );
            if ( transform_write( fp, relpath, map, 1, nin, nout, np, in, out ) ) {
                wrote++;
                fwd_out = out;            /* keep for inverse inputs */
            } else {
                free_cols( out, nout );
            }
            free_cols( in, nin );
        }

        if ( has_inv ) {
            /* inverse takes nout-dim input -> nin-dim output */
            double **inv_in;
            int own_inv_in;
            if ( fwd_out ) { inv_in = fwd_out; own_inv_in = 0; }
            else {
                inv_in = alloc_cols( nout, np );
                sample_inputs( root, relpath, kind, NULL, nout, np, inv_in );
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

/* Recursively collect files below subdir whose names end in suffix.  Every
   returned name is the full path relative to root. */
static void scan_tree_into( const char *root, const char *subdir,
                            const char *suffix, char ***list, size_t *n,
                            size_t *cap ) {
    char path[1024];
    if ( snprintf( path, sizeof path, "%s%s%s",
                   root, subdir[0] ? "/" : "", subdir )
         >= (int) sizeof path ) {
        fprintf( stderr, "gen_transform_oracle: path too long, not scanned: "
                         "%s/%s\n", root, subdir );
        return;
    }
    DIR *d = opendir( path );
    if ( !d ) return;
    struct dirent *e;
    size_t slen = strlen( suffix );
    while ( ( e = readdir( d ) ) ) {
        if ( strcmp( e->d_name, "." ) == 0 ||
             strcmp( e->d_name, ".." ) == 0 ) continue;

        char rel[1024];
        char child[1024];
        /* Truncation would leave a path that no longer names the fixture, so
           report it instead of letting the stat() below quietly drop it. */
        if ( snprintf( rel, sizeof rel, "%s%s%s", subdir,
                       subdir[0] ? "/" : "", e->d_name )
             >= (int) sizeof rel ||
             snprintf( child, sizeof child, "%s/%s", root, rel )
             >= (int) sizeof child ) {
            fprintf( stderr, "gen_transform_oracle: path too long, skipped: "
                             "%s/%s%s%s\n", root, subdir,
                             subdir[0] ? "/" : "", e->d_name );
            continue;
        }
        struct stat st;
        if ( stat( child, &st ) != 0 ) continue;
        if ( S_ISDIR( st.st_mode ) ) {
            scan_tree_into( root, rel, suffix, list, n, cap );
            continue;
        }

        size_t nlen = strlen( e->d_name );
        if ( nlen <= slen || strcmp( e->d_name + nlen - slen, suffix ) != 0 )
            continue;
        if ( *n == *cap ) {
            *cap *= 2;
            *list = xalloc( *list, sizeof(char *) * *cap, "fixture list" );
        }
        char *copy = xalloc( NULL, strlen( rel ) + 1, "fixture path" );
        strcpy( copy, rel );
        (*list)[(*n)++] = copy;
    }
    closedir( d );
}

static char **scan_tree( const char *root, const char *subdir,
                         const char *suffix, int *count ) {
    size_t cap = 64, n = 0;
    char **list = xalloc( NULL, sizeof(char *) * cap, "fixture list" );
    scan_tree_into( root, subdir, suffix, &list, &n, &cap );
    qsort( list, n, sizeof(char *), cmp_str );
    *count = (int) n;
    return list;
}

int main( int argc, char *argv[] ) {
    int status_value = 0;
    int *status = &status_value;
    if ( argc < 5 ) {
        fprintf( stderr, "Usage: gen_transform_oracle <root> <simplify_out> "
                         "<headers_out> <framesets_out>\n" );
        return 1;
    }
    const char *root          = argv[1];
    const char *simplify_out  = argv[2];
    const char *headers_out   = argv[3];
    const char *framesets_out = argv[4];
    astWatch( status );

    oracle_register_intramaps( status );

    int ndom = 0;
    Domain *doms = load_domains( root, &ndom );

    /* The tolerance line is informational (the checker takes its tolerances
       from the transform_oracle.h macros); format it from those same macros
       so it cannot drift from what the checker actually applies. */
    char hdr[256];
    snprintf( hdr, sizeof hdr,
        "# transform oracle - regenerate with gen_transform_oracle\n"
        "# tol: rtol=%g atol=%g  equiv_rtol=%g equiv_atol=%g"
        "  rtrip_rtol=%g rtrip_atol=%g\n\n",
        ORACLE_DEF_RTOL, ORACLE_DEF_ATOL,
        ORACLE_DEF_EQUIV_RTOL, ORACLE_DEF_EQUIV_ATOL,
        ORACLE_DEF_RTRIP_RTOL, ORACLE_DEF_RTRIP_ATOL );

    /* Native-dump corpus: .map and .simp interleaved by stem so paired
       sections sit together in the file. */
    FILE *fs = fopen( simplify_out, "w" );
    if ( !fs ) { fprintf( stderr, "cannot write %s\n", simplify_out ); return 1; }
    fputs( hdr, fs );
    int nmap = 0;
    char **maps = scan_tree( root, "simplify", ".map", &nmap );
    int nwritten = 0;
    for ( int i = 0; i < nmap; i++ ) {
        char simp[1024];
        size_t L = strlen( maps[i] );
        snprintf( simp, sizeof simp, "%.*s.simp", (int)(L - 4), maps[i] );
        Domain *dom = find_domain( doms, ndom, maps[i] );
        if ( dom ) dom->used = 1;
        nwritten += emit_fixture( fs, root, maps[i], dom );
        /* "neg_*" fixtures (and any other already-simple map) have no .simp
           partner; only emit one when the file exists. */
        if ( file_exists( root, simp ) )
            nwritten += emit_fixture( fs, root, simp, dom );
        free( maps[i] );
    }
    free( maps );
    fclose( fs );

    /* A domain that names no fixture has outlived it, or names it wrongly. */
    int stale = 0;
    for ( int d = 0; d < ndom; d++ ) {
        if ( !doms[d].used ) {
            fprintf( stderr, "gen_transform_oracle: a domain is declared for "
                     "%s, which is not a fixture\n", doms[d].relpath );
            stale++;
        }
        free( doms[d].lo );
        free( doms[d].hi );
    }
    free( doms );
    if ( stale ) return 1;

    /* FITS-header corpus. */
    FILE *fh = fopen( headers_out, "w" );
    if ( !fh ) { fprintf( stderr, "cannot write %s\n", headers_out ); return 1; }
    fputs( hdr, fh );
    int nhead = 0;
    char **heads = scan_tree( root, "", ".head", &nhead );
    for ( int i = 0; i < nhead; i++ ) {
        nwritten += emit_fixture( fh, root, heads[i], NULL );
        free( heads[i] );
    }
    free( heads );
    fclose( fh );

    /* FrameSet-dump corpus (.ast).  Non-FrameSet .ast dumps (KeyMap, STC,
       etc.) are not Mappings and are skipped by the loader. */
    FILE *ff = fopen( framesets_out, "w" );
    if ( !ff ) { fprintf( stderr, "cannot write %s\n", framesets_out ); return 1; }
    fputs( hdr, ff );
    int nast = 0;
    char **asts = scan_tree( root, "", ".ast", &nast );
    for ( int i = 0; i < nast; i++ ) {
        nwritten += emit_fixture( ff, root, asts[i], NULL );
        free( asts[i] );
    }
    free( asts );
    fclose( ff );

    printf( "gen_transform_oracle: wrote %d sections\n", nwritten );
    return astOK ? 0 : 1;
}
