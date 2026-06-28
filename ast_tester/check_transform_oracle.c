/*
 * check_transform_oracle: read a transform oracle file (written by
 * gen_transform_oracle), reload each named fixture, re-transform the
 * recorded inputs, and compare outputs within tolerance.  The oracle
 * file is the sole manifest; this program does no corpus scanning.
 *
 * Usage:
 *   check_transform_oracle <root> <oracle_file> [<rtrip_overrides>]
 *   check_transform_oracle --selftest
 */
#include "ast.h"
#include "transform_oracle.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    char relpath[1024];
    int nin, nout, forward;
    int npoint;
    double **in;
    double **out;
} Section;

static double **alloc_cols( int ncol, int nrow ) {
    double **c = malloc( sizeof(double *) * (size_t) ncol );
    for ( int a = 0; a < ncol; a++ ) c[a] = malloc( sizeof(double) * (size_t) nrow );
    return c;
}

static int oracle_read_file( const char *path, Section **out_secs, int *out_n ) {
    FILE *fp = fopen( path, "r" );
    if ( !fp ) { fprintf( stderr, "cannot open oracle %s\n", path ); return 1; }

    size_t cap = 64, n = 0;
    Section *secs = malloc( sizeof(Section) * cap );
    char line[8192];
    Section *cur = NULL;
    int row = 0, rowcap = 0;

    while ( fgets( line, (int) sizeof line, fp ) ) {
        if ( line[0] == '#' || line[0] == '\n' ) continue;
        if ( line[0] == '[' ) {
            if ( n == cap ) { cap *= 2; secs = realloc( secs, sizeof(Section)*cap ); }
            cur = &secs[n++];
            char dir[16] = "forward";
            /* header: [<relpath>  nin=N nout=M dir=word] */
            if ( sscanf( line, "[%1023s nin=%d nout=%d dir=%15[^]]",
                         cur->relpath, &cur->nin, &cur->nout, dir ) < 3 ) {
                fprintf( stderr, "bad section header: %s", line );
                fclose( fp ); return 1;
            }
            cur->forward = ( strncmp( dir, "inverse", 7 ) != 0 );
            cur->npoint = 0;
            rowcap = oracle_sample_axis_count() + 8;
            cur->in  = alloc_cols( cur->nin,  rowcap );
            cur->out = alloc_cols( cur->nout, rowcap );
            row = 0;
            continue;
        }
        if ( !cur ) continue;
        /* a data row: nin inputs then nout outputs */
        char *p = line, *tok;
        int col = 0, ok = 1;
        if ( row >= rowcap ) {
            rowcap *= 2;
            for ( int a = 0; a < cur->nin;  a++ )
                cur->in[a]  = realloc( cur->in[a],  sizeof(double)*(size_t)rowcap );
            for ( int a = 0; a < cur->nout; a++ )
                cur->out[a] = realloc( cur->out[a], sizeof(double)*(size_t)rowcap );
        }
        while ( ( tok = strtok( p, " \t\r\n" ) ) ) {
            p = NULL;
            int good = 0;
            double v = oracle_parse_double( tok, &good );
            if ( !good ) { ok = 0; break; }
            if ( col < cur->nin ) cur->in[col][row] = v;
            else if ( col < cur->nin + cur->nout ) cur->out[col - cur->nin][row] = v;
            col++;
        }
        if ( ok && col == cur->nin + cur->nout ) { row++; cur->npoint = row; }
    }
    fclose( fp );
    *out_secs = secs; *out_n = (int) n;
    return 0;
}

static int selftest( void ) {
    const char *tmp = "/tmp/oracle_selftest.txt";
    FILE *fp = fopen( tmp, "w" );
    fputs( "# header\n\n"
           "[foo.map  nin=2 nout=1 dir=forward]\n"
           "  1 2  3\n"
           "  4 5  BAD\n\n", fp );
    fclose( fp );
    Section *s = NULL; int n = 0;
    if ( oracle_read_file( tmp, &s, &n ) || n != 1 ) {
        fprintf( stderr, "selftest: parse failed\n" ); return 1;
    }
    int ok = ( strcmp( s[0].relpath, "foo.map" ) == 0 && s[0].nin == 2 &&
               s[0].nout == 1 && s[0].forward == 1 && s[0].npoint == 2 &&
               s[0].in[0][0] == 1 && s[0].in[1][1] == 5 &&
               s[0].out[0][0] == 3 && s[0].out[0][1] == AST__BAD );
    printf( "selftest: %s\n", ok ? "ok" : "FAIL" );
    return ok ? 0 : 1;
}

int main( int argc, char *argv[] ) {
    int status_value = 0; int *status = &status_value;
    astWatch( status );
    if ( argc == 2 && strcmp( argv[1], "--selftest" ) == 0 ) return selftest();
    fprintf( stderr, "checking not yet implemented\n" );
    return 2;   /* replaced in Task 6 */
}
