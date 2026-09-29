# Performance benchmarks

Small stand-alone programs for measuring the cost of specific AST code
paths. They are not tests and nothing runs them automatically. CMake
builds them with the rest of the tree; the autotools build compiles them
under `make check`, so that `make distcheck` keeps them building.

## `bench_lock.c`: object locking

Measures the cost of an `astUnlock`/`astLock` pair on a FITS-WCS FrameSet
(a TAN projection read through a FitsChan). Locking a FrameSet also locks
every Object it contains, so each `astLock` runs the Object class's lock
management, including `ChangeThreadVtab`, once per component Object.

Before timing, the program creates one object of each of about 25 classes,
so that the calling thread's list of known vtabs is about as long as in a
real application. `ChangeThreadVtab` searches that list.

### Building

Time it in an optimised build without sanitizers, with thread support
(the default when pthreads are available):

```sh
cmake -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel --target bench_lock
```

With autotools, `make check` builds `perf/bench_lock`.

### Running

```sh
build-rel/perf/bench_lock [N]      # N unlock/lock pairs, default 200000
```

It prints the mean time per unlock/lock pair.

To see where the time goes, `perf` is the lighter option where it is
permitted. Otherwise use callgrind with a smaller N:

```sh
valgrind --tool=callgrind --callgrind-out-file=cg.out ./bench_lock 20000
callgrind_annotate --inclusive=yes cg.out
```

### Results

Measured 2026-09-29 on Linux x86-64, GCC, Release build, five runs each.
They record the effect of changing how `ChangeThreadVtab` (`src/object.c`)
finds the calling thread's vtab for an Object's class:

| `ChangeThreadVtab` lookup | ns per unlock/lock pair |
|---|---|
| `strcmp` of class names over the thread's known vtabs | 2350-2550 |
| Comparison of class identifiers (`vtab->top_id->check`) | 1470-1590 |
| Class identifiers, plus an early return when the Object already uses one of the calling thread's vtabs | 990-1060 |

With class-name comparison, callgrind showed about 300 `strcmp` calls per
`astLock` of this FrameSet, and `ChangeThreadVtab` accounted for about 52%
of all instructions executed.

The benchmark measures only the case where the Object already uses the
calling thread's vtabs, which is what the early return skips. The middle
row was measured with the early return disabled. It slightly overstates
the cost of the identifier comparison, since without the early return an
Object that already uses the right vtab also takes and releases a
reference on the calling thread's global data.
