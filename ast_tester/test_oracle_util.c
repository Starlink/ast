/* Unit tests for transform_oracle_util.c pure helpers. */
#include "ast.h"
#include "transform_oracle.h"
#include <stdio.h>
#include <math.h>

static int failures = 0;
#define CHECK(cond) do { if(!(cond)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #cond); failures++; } } while(0)

static void test_within_tol(void) {
    /* exact equality */
    CHECK(oracle_within_tol(1.0, 1.0, 0.0, 0.0) == 1);
    /* a few ULP apart, within default rtol */
    CHECK(oracle_within_tol(1.0 + 4e-16, 1.0, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 1);
    /* gross difference rejected */
    CHECK(oracle_within_tol(1.0001, 1.0, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 0);
    /* absolute floor near zero */
    CHECK(oracle_within_tol(5e-13, 0.0, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 1);
    /* BAD matches BAD only */
    CHECK(oracle_within_tol(AST__BAD, AST__BAD, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 1);
    CHECK(oracle_within_tol(AST__BAD, 1.0, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 0);
    CHECK(oracle_within_tol(1.0, AST__BAD, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 0);
    /* NaN never matches */
    CHECK(oracle_within_tol(NAN, 1.0, 1e9, 1e9) == 0);
}

int main(void) {
    test_within_tol();
    if (failures) { fprintf(stderr, "%d failure(s)\n", failures); return 1; }
    printf("test_oracle_util: all passed\n");
    return 0;
}
