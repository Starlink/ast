# IsSimple tag placement depends on aliasing history, not tree shape

**Status:** analysis complete, no fix written. For team discussion.
**Date:** 9 August 2026
**Affects:** `src/mapping.c`, `src/mapping.h`, and every class that probes a Mapping's Invert flag.

## Summary

`AST__ISSIMPLE_FLAG` is a single bit carrying two unrelated meanings: "do not simplify this again" and "write `IsSimp = 1` when dumping".
`astSetInvert` clears the bit unconditionally, including when the value being set is the value already stored.
Routines with a documented read-only contract — `astEqual` among them — flip a component's Invert flag to line it up for inspection and then restore the *value*, but not the bit that the flip incidentally cleared.

The consequence is that which node in a serialized Mapping tree carries `IsSimp = 1` is decided by which nodes happened to be aliased into some earlier probe, not by the shape of the tree or by which parts of it have actually been simplified.
Two Mappings that `astEqual` reports as equal can serialize differently.

## Motivation: the observable defect

All three demonstrations below use only the public C API against an unmodified `build` tree.

### 1. Setting Invert to the value it already holds clears the tag

```c
AstMapping *s = astSimplify( astCmpMap( astZoomMap( 1, 2.0, "" ),
                                        astZoomMap( 1, 3.0, "" ), 1, "" ) );
int inv = astGetI( s, "Invert" );
astSetI( s, "Invert", inv );          /* a no-op by any reasonable reading */
```

```
1. astSetInvert to the SAME value: IsSimp 1 -> 0  (Invert still 0)
```

Nothing about the Mapping changed, yet it now claims not to have been simplified — so `astSimplify` will walk it again, and a dump of it will disagree with a dump taken a moment earlier.

### 2. A read-only comparison silently drops an interior tag

`astEqual` is documented as a comparison.
`cmpmap.c`'s `Equal` implements it by decomposing both operands with `astMapList`, which returns reference-counted pointers to the **live** components rather than copies, and lining them up with `astSetInvert` (`src/cmpmap.c:385-397`):

```c
/* Temporarily set the Mapping Invert flags to the required values,
   saving the original values so that they can be re-instated later.*/
               this_inv = astGetInvert( this_map_list[ i ] );
               astSetInvert( this_map_list[ i ], this_invert_list[ i ] );
               that_inv = astGetInvert( that_map_list[ i ] );
               astSetInvert( that_map_list[ i ], that_invert_list[ i ] );

/* Compare the two component Mappings for equality. */
               result = astEqual( this_map_list[ i ], that_map_list[ i ] );

/* Re-instate the original Invert flags. */
               astSetInvert( this_map_list[ i ], this_inv );
               astSetInvert( that_map_list[ i ], that_inv );
```

The comment says exactly what the routine intends, and the Invert *value* is indeed reinstated.
The `IsSimple` bit is not — and cannot be, because the routine has no idea it was cleared.

Comparing a two-node tree against an unrelated throwaway Mapping:

```
2. after astEqual (read-only):     IsSimp 2 -> 1
```

The tag on the interior `ZoomMap` is gone.
The outer `CmpMap` keeps its own, because it was not one of the components handed back by `astMapList`.

### 3. Equal Mappings serialize differently

Build two structurally identical trees the same way, pass only the first to `astEqual`, then dump both:

```
two structurally identical Mappings dump identically: NO
astEqual( t, u ) still reports them equal:            yes
```

The diff is one line, on the interior node:

```
  Begin CmpMap
  Nin = 2
  IsSimp = 1
  IsA Mapping
  MapA =
  Begin ZoomMap
  Nin = 2
- IsSimp = 1          <-- present in the untouched Mapping, absent in the probed one
  IsA Mapping
  Zoom = 6
  End ZoomMap
```

This is the property worth fixing.
Tag placement should be derivable from first principles — "this node was returned by `astSimplify`, so it is tagged" — rather than depending on the internal order in which some unrelated routine happened to clone, flip and restore a component.

## Mechanism

`src/mapping.h:410` defines the bit, one of four in `Mapping.flags`:

```c
#define AST__ISSIMPLE_FLAG 1  /* Mapping has been simplified */
```

Two consumers read it, for two different purposes:

| Purpose | Site |
|---|---|
| **Simplify guard** — short-circuits `astSimplify_` | `src/mapping.c:24744`, `if( !astGetIsSimple( this ) && !astDoNotSimplify( this ) )` |
| **Dump tag** — writes the `IsSimp` card | `src/mapping.c:23763`, read back by `astLoadMapping` at `:24117` |

and it is raised in two places: `astSimplify_` on its result (`:24763`) and the loader from the `IsSimp` card (`:24117`).

The clear happens here (`src/mapping.c:23260`):

```c
astMAKE_SET(Mapping,Invert,int,invert,(astClearIsSimple(this),(value!=0)))
```

`astClearIsSimple` fires before the assignment and without comparing the incoming value against the stored one, which is why demonstration 1 works.
The corresponding `astMAKE_CLEAR` (`:23257`) does **not** touch the flag, so `astClearInvert` is already tag-preserving — a detail that matters below.

There is also a fifth site: `GetAttrib` exposes `IsSimple` as a **public, read-only attribute** (`src/mapping.c:1555`), documented at `:23307-23345`.
Its documentation describes only the guard meaning:

> This attribute indicates whether a Mapping has been simplified by the astSimplify method. If the IsSimple value is non-zero, then the Mapping has been simplified and so there is nothing to be gained by simplifying it again.

## Analysis

An instrumented build (both the setter and clearer wrapped to log object address, class, before/after value, whether the value changed, whether `IsSimple` was raised on entry, and a resolved backtrace; object and memory caching short-circuited so addresses are not recycled) was run over the full 1039-test suite.
The instrumentation was behaviour-neutral — 1039/1039, matching baseline — and produced 584,161 events, of which **25,644 `astSetInvert` calls cleared an already-raised `IsSimple` bit**.

Those 25,644 split by the calling routine's contract:

| Group | Clears | Routines | Files |
|---|---:|---:|---:|
| **Probe — routine's contract is read-only. This is the defect.** | **9,458** | **38** | **17** |
| Genuine mutation — routine changes Invert deliberately | 16,186 | 10 | 8 |

The mutation group is defensible under the reading that a flipped Mapping is a different Mapping and may simplify differently; it is dominated by `frame.c:3203` (`ConvertX`, 15,724 clears).

Probe-group clears by file:

```
cmpmap.c 4150 · winmap.c 1897 · frameset.c 1003 · wcsmap.c 705 · matrixmap.c 571
fitschan.c 544 · ratemap.c 198 · tranmap.c 149 · switchmap.c 95 · polymap.c 43
cmpregion.c 32 · mapping.c 28 · zoommap.c 15 · moc.c 10 · region.c 10
specfluxframe.c 6 · yamlchan.c 2
```

Heaviest individual sites: `cmpmap.c:389` `Equal` (1215), `frameset.c:5444` `GetMapping` (997), `cmpmap.c:393` `Equal` (680), `cmpmap.c:3681`/`:3705` `Simplify` (519/507), `cmpmap.c:1856` `MapMerge` (463), `winmap.c:347` `CanSwap` (453), `wcsmap.c:1021` `CanSwap` (447).

A runtime classification was also computed (7,513 clears where the value did not change, 12,677 restored by the same routine, 5,454 kept) but is not the authority: runtime pairing cannot distinguish "this routine restored the value" from "the same address was flipped again later for an unrelated reason", and `ConvertX` alone inflates the "restored" column by 11,038.
The static, contract-based classification is what the table above reports.

### The `astClearInvert` exemption does not exist

Because `astMAKE_CLEAR` is tag-preserving, probes that restore Invert via `astClearInvert` would need no change.
Nine probe routines do use one — `cmpmap.c` `CombineMaps`/`MapMerge`/`Simplify`, `frameset.c` `CombineMaps`/`Simplify`, `fitschan.c` `AnalysePoly`, `polymap.c` `MergeShift`, `yamlchan.c` `IsPolyMap`, `zoommap.c` `MapMerge` — 20 calls out of 26 tree-wide.

But **not one restores exclusively that way.**
Every one uses the standard two-branch shape (`src/cmpmap.c:3711-3716`):

```c
if ( !set ) {
   astClearInvert( map_list[ i ] );
} else {
   astSetInvert( map_list[ i ], invert );
}
```

Both branches restore the value; only the `astSetInvert` branch cleared the tag on the way in; neither restores it.
So no site can be skipped on these grounds, and the anticipated exemption is empty.

## Scale of the two candidate fixes

### Option A — split the flag

Keep `AST__ISSIMPLE_FLAG` as the simplify guard that `astSetInvert` clears, since a Mapping whose Invert genuinely changed may simplify differently.
Add a second bit set only by `astSimplify_` and the loader, read only by `Dump` and the attribute getter, and never touched by `astSetInvert`.

**Files: 2.** `src/mapping.h` (one `#define` plus three macros, alongside the four flags already there) and `src/mapping.c` (five sites).

```c
/* src/mapping.h, beside the existing four flags (1, 2, 4, 8) */
#define AST__SIMPLIFIED_FLAG 16   /* Mapping is the result of a simplification */

#define astSetSimplified(this) \
((void)(this&&(((AstMapping*)this)->flags|=AST__SIMPLIFIED_FLAG)))
#define astClearSimplified(this) \
((void)(this&&(((AstMapping*)this)->flags&=~AST__SIMPLIFIED_FLAG)))
#define astSimplified(this) \
(this&&((((AstMapping*)this)->flags&AST__SIMPLIFIED_FLAG)!=0))
```

The five `mapping.c` edits:

| Site | Change |
|---|---|
| `:23260` | unchanged — `astSetInvert` still clears the guard, and only the guard |
| `:24763` | `astSetIsSimple( result )` also raises the new bit |
| `:24117` | loader raises the new bit from the `IsSimp` card |
| `:23763` | `Dump` reads the new bit |
| `:1555` / `:23347` | the public `IsSimple` attribute reads one of the two — see the open question below |

Every probe in the library stays exactly as written.
`astSimplify_`'s guard keeps its current behaviour.

Two details the implementation has to respect:

- `Mapping.flags` is a `char` (`src/mapping.h:461`) holding four bits so far, so bit 16 is available, but the field is close enough to full that a fifth bit is worth a comment rather than being added silently.
- `astGetIsSimple` is a **virtual** accessor, and `frame.c` overrides it to return 0 unconditionally (`src/frame.c:5590-5621`, because a Frame is mutable and so can never be guaranteed simple).
  `Dump` currently reads the bit through that virtual accessor, which is what keeps `IsSimp` out of Frame dumps — confirmed: a bare `astFrame` dumps no `IsSimp` card and reports `IsSimple = 0`.
  A raw `astSimplified(this)` macro read in `Dump` would lose that suppression and start emitting the card for Frames, so the new bit needs its own virtual getter, or `Dump` must keep gating on the existing one.

### Option B — make probes non-mutating

Save and restore `IsSimple` alongside `Invert` at each probe site.

**Files: 28. Routines: 98. Brackets: 211** — 157 save/restore brackets (`x = astGetInvert(M)` … restore) and 54 `astInvert` … read … `astInvert` brackets.
Those routines contain **283 of the library's 307 `astSetInvert` calls**, so the edit touches most of the library's Invert traffic.

Each of the 211 brackets needs the shape below, which does not currently exist anywhere in the tree because there is no public or protected accessor for reading the bit back other than `astGetIsSimple`:

```c
/* Save. */
   this_inv  = astGetInvert( map );
   this_simp = astGetIsSimple( map );          /* new */
   astSetInvert( map, wanted );

   ... inspect ...

/* Restore. */
   astSetInvert( map, this_inv );
   if( this_simp ) astSetIsSimple( map );      /* new */
```

Two properties of this option are worth weighing:

**60 of the 98 routines produce no signal in a full 1039-test run.**
They were never observed clearing a raised bit, so an error in them ships silently.
Among them: `fitschan.c` `WcsNative`, `GrismSpecWcs`, `TabMapping`, `CelestialAxes`, `NearestPix`, `AnalysePoly`; `tranmap.c` `Equal`/`MapSplit`/`Rate`/`Transform`; all four `pcdmap.c` probes; `matrixmap.c` `MatWin`/`MatWin2`/`MatZoom`/`MatPermSwap`; ten of the eleven `yamlchan.c` probes; `plot.c` `Boundary`/`ToggleLogLin`.

**Per-site judgement is required, and it is genuinely subtle.**
Several routines look like probes but are safe because they flip a private copy — `fitschan.c` `SplitMap2` (334 clears, all no-ops, on an `astCopy` at `:33186`/`:33199`/`:33250`) and `frameset.c` `GetMapping` (997 clears, on `astCopy( path[0] )` at `:5443`).
Worse, `cmpmap.c` `CombineMaps` decides per argument, within one routine, based on whether the two arguments alias (`src/cmpmap.c:615-621`):

```c
/* If both Mappings are actually the same but we need different Invert
   flag values to be set, then this can only be achieved by making a
   copy. Note if this is necessary. */
   copy = ( ( mapping1 == mapping2 ) && ( invert1 != invert2 ) );

/* Clone the first Mapping pointer. Do likewise for the second but
   make a copy instead if necessary. */
   map1 = astClone( mapping1 );
   map2 = copy ? astCopy( mapping2 ) : astClone( mapping2 );
```

`map1` is always a clone and so always shares the caller's object; `map2` is a clone or a copy depending on runtime aliasing.
A mechanical save/restore pass adds dead restores in the safe cases and still has to get the aliasing analysis right in the rest.

### Side by side

| | Option A: split the flag | Option B: non-mutating probes |
|---|---|---|
| Files changed | 2 | 28 |
| Routines changed | 1 (`astSimplify_`) + loader + `Dump` | 98 |
| Edit sites | ~5 | 211 brackets |
| Sites with no test coverage | 0 | 60 of 98 routines |
| Per-site judgement needed | no | yes, including per-argument within a routine |
| Fixes the cause or the symptom | cause — the bit stops meaning two things | symptom — the ambiguity survives for the next caller |

## Proposed outcome

**Recommendation: Option A, split the flag.**

The reasoning, in the order it mattered:

1. **The edit surfaces are two orders of magnitude apart, and one of them is largely untestable.** Two files of heavily-exercised code against 211 brackets in 28 files, 60 routines of which the suite cannot check.

2. **The two meanings want different lifetimes.**
   Clearing the *guard* when Invert genuinely changes is correct and should be kept; clearing the *dump tag* on a transient probe flip never is.
   One bit cannot express both, and that is precisely why the defect exists.
   Option B patches 211 consumers while leaving the flag ambiguous for the next person who touches Invert.

3. **Containment and a checkable invariant.** The bit is touched in five places tree-wide, so after the split the rule "only `astSimplify_` and the loader set the dump bit" is locally verifiable by `grep`, permanently.

4. **The exemption that would have shrunk Option B is empty**, as shown above, so its 98 routines cannot be reduced by excluding already-safe sites.

Fixture churn is not an argument against this.
Roughly 304 committed fixtures carry `IsSimp` cards and some will move, but that movement *is* the fix becoming visible: afterwards, tag placement follows from which nodes `astSimplify` returned, and can be reasoned about from the tree alone rather than reconstructed from the order in which some routine cloned and flipped a component.

### Open question for the team

The public `IsSimple` attribute (`src/mapping.c:1555`, documented at `:23307`) must read one of the two bits.
Its documentation describes the guard meaning — "there is nothing to be gained by simplifying it again" — which argues for keeping it on the guard.
Against that, reporting the dump bit keeps `astGetC( map, "IsSimple" )` consistent with the `IsSimp` card in a dump of the same object, which is what a user comparing the two would expect.

The attribute is read-only and has no in-tree consumers outside `mapping.c` itself, so either choice is safe to implement; it is an API-semantics decision rather than a mechanical one.
Whichever is chosen, the attribute's documentation should be updated to say which of the two it reports.

### Deliberately out of scope

Whether `astMAKE_SET` should fire `astClearIsSimple` at all when the value is unchanged (demonstration 1) is a separate question.
Under the split it stops being a correctness issue and becomes a performance one — a redundant re-simplification — and can be addressed independently by making the clear conditional on an actual change.

## Reproducing

The instrumented build was a throwaway `git worktree` under `/tmp`, since removed; reproducing it means re-applying the instrumentation described above.
The three demonstrations need only the public API and an ordinary `build` tree.
On this host they link with:

```
gcc -std=c11 -I build -I . demo.c -o demo -Wl,-rpath,$PWD/build -Lbuild -last \
    -Wl,--whole-archive build/libast_err.a build/libast_grf_2.0.a \
    build/libast_grf_3.2.a build/libast_grf_5.6.a build/libast_grf3d.a \
    -Wl,--no-whole-archive -lyaml -lpthread -lm
```

The `-Wl,-rpath` must precede the conda RPATH or the probe silently resolves `libast.so` from an installed conda copy instead of the tree, and the graphics stubs need whole-archive linkage or the link fails on `astGLine` and friends.
