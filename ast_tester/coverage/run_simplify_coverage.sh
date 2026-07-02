#!/usr/bin/env bash
# Regenerate the simplify coverage gap ledger.
# Requires gcovr on PATH (e.g. `source ~/pyenv/bin/activate`).
set -euo pipefail
cd "$(dirname "$0")/../.."          # repo root
ROOT="$PWD"
BUILD="build-cov"
COVDIR="ast_tester/coverage"

if ! command -v gcovr >/dev/null 2>&1; then
    echo "ERROR: gcovr not found. Activate your gcovr env, e.g.:" >&2
    echo "  source ~/pyenv/bin/activate" >&2
    exit 1
fi

# In-scope file filter (32 MapMerge files + mapping.c).
FILTER='src/(mapping|box|cmpmap|dssmap|grismmap|interval|intramap|lutmap|mathmap|matrixmap|normmap|nullregion|pcdmap|permmap|pointlist|polymap|prism|ratemap|selectormap|shiftmap|slamap|specmap|sphmap|splinemap|switchmap|timemap|tranmap|unitmap|unitnormmap|wcsmap|winmap|xphmap|zoommap)\.c$'

cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_C_FLAGS="--coverage -O0" \
      -DCMAKE_EXE_LINKER_FLAGS="--coverage" \
      -DCMAKE_SHARED_LINKER_FLAGS="--coverage" >/dev/null
cmake --build "$BUILD" >/dev/null

# The gcov tool must match the compiler that produced the .gcno data. For
# clang that is the `llvm-cov gcov` shipped beside the compiler; using a
# mismatched (e.g. Apple xcrun) llvm-cov yields version-mismatch errors.
if [ -z "${GCOV:-}" ]; then
    CC_PATH=$(sed -n 's/^CMAKE_C_COMPILER:[^=]*=//p' "$BUILD/CMakeCache.txt")
    case "$(basename "$CC_PATH")" in
        clang*|cc)
            LLVM_COV="$(dirname "$CC_PATH")/llvm-cov"
            if [ -x "$LLVM_COV" ]; then GCOV="$LLVM_COV gcov"
            elif [ "$(uname)" = "Darwin" ]; then GCOV="xcrun llvm-cov gcov"
            else GCOV="llvm-cov gcov"; fi
            ;;
        *) GCOV="gcov" ;;
    esac
fi

run_gcovr() {  # $1 = output json
    gcovr --root "$ROOT" "$BUILD" \
          --gcov-executable "$GCOV" \
          --filter "$FILTER" --json --output "$1" --delete >/dev/null
}

echo ">> simplify-only pass"
ctest --test-dir "$BUILD" -R '^simplify_' --output-on-failure >/dev/null
run_gcovr "$COVDIR/simplify.json"

echo ">> full-suite pass"
ctest --test-dir "$BUILD" >/dev/null || true   # coverage even if some non-simplify test fails
run_gcovr "$COVDIR/full.json"

python3 "$COVDIR/simplify_coverage_gaps.py" \
    --simplify "$COVDIR/simplify.json" \
    --full "$COVDIR/full.json" \
    --ledger "$COVDIR/simplify_coverage_gaps.md"
