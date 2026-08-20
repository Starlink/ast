# Design: C library defect fixes in KeyMap, FitsChan and Mapping

Date: 2026-08-08

## Purpose

A review of `keymap.c`, `fitschan.c`, `cmpmap.c` and `mapping.c` identified eighteen defects.
Each was verified against the source before this design was written.
A nineteenth, incorrect rounding of negative values, was found while specifying the fix for one of them and turned out to affect twenty-two files rather than one.
They fall into four groups: the rounding correction, unambiguous defects with local fixes, determinism problems that change serialized output, and design-level issues that change behavior on real data.

This document specifies the fix for each and the order in which they land.

## Scope

All four groups are in scope.
One defect, `MapSz`, receives no code change; the reasoning is recorded below so the decision is not revisited.
One defect, the `IsSimple` tag placement, begins with an investigation whose result selects between two candidate fixes.

## Structure

Work lands on `u/timj/misc-bugs` as one commit per defect.
Each commit updates the history section of the prologue in the file it touches, as required by `CLAUDE.md`.
Commits are ordered so that fixes which change no output land before fixes which change serialized bytes, and the design-level group lands last.

Group R lands first, as a single sweep, because it touches lines that Group A then modifies again.

## Group R: incorrect rounding of negative values

### R1. `(int)( x + 0.5 )` rounds every negative value toward zero by one

The idiom `(int)( x + 0.5 )` is used throughout the library to round a floating-point value to the nearest integer.
It is correct only for non-negative values.
A C cast from floating-point to integer truncates toward zero, so for a negative `x` the expression evaluates `ceil( x + 0.5 )` rather than the intended nearest-integer.

Measured, compiled at `-O0`:

| value | `(int)( x + 0.5 )` | `(int)round( x )` |
|---|---|---|
| 2.3 | 2 | 2 |
| 2.5 | 3 | 3 |
| 2.7 | 3 | 3 |
| 2.0 | 2 | 2 |
| -0.6 | 0 | -1 |
| -2.0 | -1 | -2 |
| -2.3 | -1 | -2 |
| -2.5 | -2 | -3 |
| -2.7 | -2 | -3 |

Every negative input is wrong, including exact integers: a value of `-2.0` yields `-1`.
As one consequence, `astMapGet0I` on a KeyMap `Double` entry holding `-2.0` returns `-1` today.

Commit `e9515162` established both the fix and the reasoning for it, replacing the idiom with C99 `round()` in `cmpframe.c`, `moc.c`, `plot.c`, `pointlist.c`, `region.c`, `switchmap.c` and `timeframe.c`.
That change converted some occurrences in each file rather than all of them, and did not cover the whole library.
Seventy cast-rounding sites remain across sixteen files: `fitschan.c` (13), `keymap.c` (12), `moc.c` (12), the `grf_*` and `grf3d_*` graphics stubs (12 between them), `yamlchan.c` (6), `mathmap.c` (4), `matrixmap.c`, `pcdmap.c`, `wcsmap.c` and `winmap.c` (2 each), and `axis.c`, `mapping.c` and `polymap.c` (1 each).

Fix: replace each remaining `(type)( expr + 0.5 )` with `(type)round( expr )`, completing `e9515162`.

#### Which sites can actually receive a negative value

The sweep is applied to every site, including those that cannot change behavior.

Converting the provably non-negative sites is worthwhile on its own terms.
It leaves one rounding idiom in the library instead of two, so a reader does not have to work out whether a given `+ 0.5` was chosen deliberately or merely inherited, nor re-derive the sign argument at each site to satisfy themselves it cannot misbehave.
The remaining occurrences of `+ 0.5` are then all genuine arithmetic, which makes the next sweep of this kind unnecessary.

Only the groups identified below can change behavior.

Most sites are provably non-negative and the substitution is a no-op:

- `PermGet` in `matrixmap.c`, `pcdmap.c`, `wcsmap.c` and `winmap.c` (eight sites). A PermMap permutation entry may legitimately be negative, as a reference into the constants array, but `PermGet` handles that case in a separate branch that assigns `-( nc + 1 )` directly. The cast is reached only for a copied axis index, which is non-negative.
- `axis.c:987`, where the operand is `absgap / b` with `absgap` an absolute value and `b` the power of ten below it, so the quotient lies in `[1,10)`.
- `mapping.c:2423`, `0.3f * nmax` for a count.
- `mathmap.c` argument counts and variable indices, taken from a compiled opcode stream.
- `polymap.c:2048` and `yamlchan.c` polynomial powers.
- `fitschan.c` SIP polynomial powers, `naxis`, `nwcs`, the `-TAB` axis number and the interpolation code.
- The `grf_*` and `grf3d_*` graphics attribute sites: color index, font size and line width.

One site is reachable with a negative value and is easily testable through the public API: `fitschan.c:14322`, described next.

The `moc.c` sites are reachable with negative values, since grid coordinates go negative for a region extending beyond the array, but constructing that state requires assembling a Moc, a Region and a FrameSet.
They are converted with the rest of the sweep and left to the existing Moc tests rather than given a dedicated regression test.

#### R1a. A negative grism interference order is corrupted or rejected

`GrismSpecWcs` sets the interference order of a GrismMap from FITS keyword `PVi_1` with `astSetGrismM( gmap, ( pv != AST__BAD ) ? (int) ( pv + 0.5 ) : 0 )` (`fitschan.c:14322`).
`GrismM` is documented as "the interference order" and is a plain signed `int` (`grismmap.c:1923-1927`) with no non-negativity constraint.
Negative diffraction orders are physically routine.

Measured against the library built from this tree, reading a `WAVE-GRI` header and writing it back through a FITS-WCS FitsChan:

| `PV1_1` in header | stored `GrismM` | result |
|---|---|---|
| 0.0 | 0 | `astRead` fails |
| 1.0 | 1 | round-trips as 1.0 |
| -1.0 | 0 | `astRead` fails, identically to an order of 0 |
| -2.0 | -1 | round-trips as **-1.0** |

An order of -1, the commonest grism order, is rounded to 0, which is degenerate, so `astRead` rejects an otherwise valid header with "Unknown algorithm code or unusable parameter values".
An order of -2 is silently stored as -1, producing a WCS with the wrong dispersion rather than a diagnostic.

This is the only site outside `keymap.c` whose repair is worth a dedicated regression test, and it is the strongest evidence that the sweep is a correctness fix rather than a tidy-up.

This is safe to apply without auditing each site for whether its input can be negative.
For any non-negative value within the range of the destination type the two expressions are identical, including at ties: the largest such value is below 2^31, far below 2^53, so `expr + 0.5` is computed exactly and its truncation equals `round`.
The substitution is therefore a strict no-op wherever the input cannot be negative, and a correction wherever it can.

The sweep matches only the cast-rounding idiom.
Other uses of `+ 0.5` in the library are legitimate arithmetic, such as FITS pixel-center conventions and midpoint calculations, and are not touched.

`round()` requires `<math.h>` and linking against the maths library, both of which are already in place following `e9515162`.

This correction is independent of the range-clamping fix in A4, which addresses undefined behavior in the same casts.
Rounding lands first, and A4 then adds clamping on top of the corrected expressions.

## Group A: unambiguous defects

All ten are in `src/keymap.c`.

### A1. `MapIterate` never terminates when `SortBy` is set

`AddToSortedList` (`keymap.c:802`) builds a circular doubly linked list.
The empty case links an entry to itself (`keymap.c:864-866`) and the insert-at-front and insert-at-back cases wrap the ends together (`keymap.c:872-885`), so `snext` is never NULL for an entry in the list.
`MapKey` accounts for this and breaks when the walk returns to `this->first` (`keymap.c:3526-3533`).
`MapIterate`'s sorted branch does not (`keymap.c:8703-8706`).
Its only termination signal is `key == NULL`, and `key` is read from an entry that is never NULL, so a caller's loop never ends.

Fix: apply the guard `MapKey` already uses.

```c
if( entry ) {
   key = entry->key;
   this->iter_entry = ( entry->snext == this->first ) ? NULL : entry->snext;
}
```

The `SORTBY_NONE` branch (`keymap.c:8660-8688`) is correct as written and is not touched, because it walks the per-bucket `next` chains where NULL genuinely marks the end.

### A2. `astMapRename` on a locked KeyMap loses the entry

`MapRename` (`keymap.c:7608`) removes the entry from the table before it checks `MapLocked`.
`AddTableEntry`, the only call that returns the entry to the table, is reached only `if( astOK )` (`keymap.c:7735`).
The `AST__BADKEY` branch above it sets the error status, so on that path the entry falls through to `FreeMapEntry` and the KeyMap is left holding neither the old key nor the new one.

Fix: hoist the `MapLocked` check above the first `RemoveTableEntry`.
Whether the rename would introduce a new key can be determined non-destructively with `SearchTableEntry( this, itab_new, newkey )`.
Reporting the error before any mutation avoids having to reverse a partially applied rename, which would otherwise require restoring both `entry->key` and `entry->hash`.

### A3. Zero-length vectors produce a scalar that dumps pointer bytes

`AstMapEntry.nel` is documented as "0 => scalar, >0 => array" (`keymap.h:116`), so the representation has no encoding for a vector of length zero.
`MAKE_MAPPUT1` (`keymap.c:4894`) applies no check to `size`, because `CHECK_I` and its numeric siblings are empty macros (`keymap.c:4996-5003`).
A call with `size == 0` therefore stores `nel = 0` and a value pointer that `astMalloc_`'s `if ( size > (size_t) 0 )` gate (`memory.c:2931`) leaves NULL.

`DumpEntry` then reads `entry->nel`, takes the scalar branch, and aliases the `Entry1<X>` structure as an `Entry0<X>`, writing the NULL pointer field reinterpreted as the payload type.
The Object branch guards against this (`keymap.c:2895`).
The String branch does not: `((Entry0C *)entry)->value` is passed to `astWriteString` as a `const char *` (`keymap.c:2882`), which dereferences NULL.

Fix: reject `size < 1` in `MAKE_MAPPUT1` with an error, so the malformed entry cannot be created.
Separately, add to `DumpEntry`'s String branch the same NULL guard the Object branch has, as defense in depth for any entry that predates the check.

This is the one change in this design that alters a public API contract rather than repairing an internal defect.
It is preferred to silently degrading the call to a scalar, because a caller passing zero has a bug that should be reported rather than absorbed.

### A4. Out-of-range floating-point to narrow integer conversion is undefined

C11 6.3.1.4p1 makes conversion of a floating value to an integer type that cannot represent it undefined behavior.
The affected arms are the `Double` cases at `keymap.c:1952-1965`, the `Float` cases at `keymap.c:2011-2024`, and the `%lf` fallbacks in the `String` arms, none of which range-check before casting.

This is a separate defect from R1, which corrects the same casts' rounding.
R1 lands first, so this fix is applied to expressions already using `round()`.

Fix: add static clamping helpers to `keymap.c` and route every affected cast through them.

```c
static int DtoI( double dval ) {
   double r = round( dval );
   return ( r >= (double) INT_MAX ) ? INT_MAX :
          ( r <= (double) INT_MIN ) ? INT_MIN : (int) r;
}
```

with siblings for `short int`, `unsigned char` and `int64_t`.

`(double) INT_MAX` and `(double) INT_MIN` are exactly representable, so those comparisons are safe as written.
`(double) INT64_MAX` is not; the 64-bit helper compares against the literal `9223372036854775808.0` instead.

The `unsigned char` helper clamps to `0` and `UCHAR_MAX`, so a negative value saturates at zero rather than converting to a large unsigned value.

Only inputs whose behavior is already undefined change.
Every in-range value converts to the same integer that R1 leaves it converting to.

### A5. Numeric-string overflow through `%d` is undefined, and glibc wraps

C11 7.21.6.2p10 makes a `sscanf` numeric conversion undefined when the value is outside the range of the destination.
`ConvertValue`'s `String` arms scan with `astSscanf( cval, " %d %n", &ival, &nc )` and accept the result when `nval == 1 && nc >= strlen( cval )` (`keymap.c:2060-2073`, with the same shape in the `SInt`, `Byte` and `KInt` arms at `keymap.c:2076-2121`).
glibc does not fail on overflow; it wraps the accumulated value and reports success, so the acceptance test passes and the existing `%lf` fallback is never reached.

Fix: replace the primary scan with `strtol` (or `strtoll` for the 64-bit arm) and an `errno == ERANGE` check, falling through to the existing `%lf` attempt when the range check fails.
`errno`-based range detection is already used for this purpose in `axis.c`, `channel.c`, `fitschan.c`, `mathmap.c`, `memory.c` and `object.c`, so it is consistent with the codebase.

```c
char *end;
errno = 0;
long lval = strtol( cval, &end, 10 );
int consumed = ( end != cval );
while( isspace( (int) *end ) ) end++;
if( consumed && errno != ERANGE && *end == '\0' &&
    lval >= INT_MIN && lval <= INT_MAX ) {
   if( out ) *( (int *) out ) = (int) lval;
} else {
   /* existing "%lf" attempt */
```

`consumed` must be captured immediately after the `strtol` call, before the trailing-whitespace skip.
For an all-whitespace input such as `"   "`, `strtol` consumes no digits and leaves `end == cval`, but the skip loop would then advance `end` to the terminating NUL before the test ran, so a later test would spuriously pass and return zero for an input the current code correctly refuses.

Behavior changes only for digit strings whose magnitude exceeds the destination type.
Empty, whitespace-only, trailing-garbage and in-range inputs all keep their current results.

### A6. `LoadKeyMap` injects a spurious `KyCas` card

`KeyCase` is the one KeyMap attribute whose unset sentinel is `-1` (`keymap.c:10160-10161`).
`KeyError` (`keymap.c:10201`), `MapLocked` (`keymap.c:10242`) and `SortBy` (`keymap.c:10306`) all use `-INT_MAX`.
`LoadKeyMap` reads `KyCas` with a default of `-INT_MAX` (`keymap.c:11001`), which is not equal to `-1`, so `TestKeyCase` reports an absent attribute as set and the next dump writes it.
No freshly written KeyMap matches its own re-dump.

Fix: change the default at `keymap.c:11001` from `-INT_MAX` to `-1`.

The alternative, changing the sentinel to `-INT_MAX` for consistency with the other three, is rejected because it touches the public-facing default-value expression in `astMAKE_GET`.
The load-path default is the smaller change with the same effect.

### A7. `MpLck` is applied before any entry is loaded

`astLoadKeyMap_` reads and applies every attribute, `MapLocked` included, before the entry-reading loop starts (`keymap.c:11008-11036`).
Each `MapPut0<X>` and `MapPut1<X>` the loop then calls reports `AST__BADKEY` for a key not already present when `MapLocked` is set, through the shared guard in `MAKE_MAPPUT0` (`keymap.c:4736`) and `MAKE_MAPPUT1` (`keymap.c:4971`).
A KeyMap is empty when `astLoadObject` returns, so the first entry always fails: dumping a locked, non-empty KeyMap and reading it back is impossible.

Fix: move the `MpLck` read and its `SetMapLocked` call below the entry loop.

`SetMapLocked` (`keymap.c:9381`) does not error on a non-empty KeyMap, so this is safe.
It also propagates the flag into nested KeyMaps held as entries (`keymap.c:9426`), which the current placement cannot do because those entries do not exist yet when it runs.
Moving the call therefore repairs a second latent defect.

### A8. `astMapGetElem<X>` reports success on an undefined entry

`MAKE_MAPGETELEM` sets `result = 1` as soon as `SearchTableEntry` finds the entry (`keymap.c:6808`).
The `AST__UNDEFTYPE` arm sets `raw = NULL` and touches nothing else (`keymap.c:6886`).
The bounds check fires first, so an out-of-range index still reports `AST__MPVIN` (`keymap.c:6901`), but for an in-range index the conversion is guarded by `} else if( raw )` (`keymap.c:6910`), so `ConvertValue` is skipped and nothing after that point touches `result`.
The function returns 1 without writing the caller's `value` output.
`MapGetElemC` has the identical shape (`keymap.c:7015`, `7093-7094`, `7116`).

A caller who trusts the non-zero return reads whatever was already in its own buffer.
This is inconsistent with the rest of the family: `MAKE_MAPGET0` has the equivalent `if( !raw ) { result = 0; }` guard (`keymap.c:5525`) and `MAKE_MAPGET1` has it at `keymap.c:6014`.

Fix: set `result = 0` when `raw` is NULL and the index is in range, in both `MAKE_MAPGETELEM` and `MapGetElemC`.

### A9. `LoadKeyMap` resets the member counter before every put

`astLoadKeyMap_` reads a per-entry `mem%d` card and assigns it to `new->member_count` before adding each entry (`keymap.c:11067-11068`).
`DumpEntry` (`keymap.c:2700-2930`) writes `Key%d`, `Com%d`, `Typ%d`, `Nel%d` and `Val%d`/`V%d_%d`, and nothing named `Mem%d`.
`astReadInt`'s default of 0 therefore fires on every iteration.
`AddTableEntry` assigns `entry->member = (this->member_count)++` (`keymap.c:686`), so the counter is reset to 0 before each put and every loaded entry receives `member = 0`.

This is not merely a read with no matching write.
It actively destroys the per-entry member numbering that `SortBy=AgeUp` and `SortBy=AgeDown` order by.

Fix: delete the read.
Entries then receive sequential member numbers in dump order, which is the natural reconstruction.
`MemCnt` continues to be read once after the loop and restores the total.

### A10. `astMapCopyEntry` normalizes the key with the source KeyMap's `KeyCase`

`MapCopyEntry` converts the caller's key using the source KeyMap's `KeyCase` (`keymap.c:4350`) and then uses the same converted string for both the source lookup and the destination lookup.
When the two KeyMaps have different `KeyCase` settings, the destination lookup runs under the source's casing rule and can miss a key that its own rule would match.

Fix: call `ConvertKey` once per KeyMap, each with its own `KeyCase`, and use the appropriate string for each lookup.

## Group B: determinism

### B1. KeyMap dumps never converge

`Dump` walks the hash table bucket by bucket and follows each bucket's `next` chain from the head (`keymap.c:10516-10525`), ignoring `SortBy` entirely.
`LoadKeyMap` re-puts the entries in the order it reads them, and `AddTableEntry` head-inserts into the bucket (`keymap.c:675`).
Any two keys sharing a bucket therefore swap places on every write/read round trip, and the dump never stabilizes.
This affects shipped fixtures: `ast_tester/stcschan-test1-doc3-props.ast` has `SIZE` and `FRAME` colliding, and `TIME` and `REFPOS` colliding.

Fix: emit entries from `Dump` in `KeyCmp` order rather than bucket order.
Because `Dump` already ignores `SortBy`, this changes no documented behavior, and alphabetical output is idempotent under any insertion order.
The hash table itself and the `MapSz` card that records its size are unaffected.

### B2. `EncodeFloat` injects a stray space into wide values

`EncodeFloat` (`fitschan.c:9754`) formats into a right-justified 20-column field, so exponents arrive with C's mandatory two digits.
It then removes a redundant leading zero (`fitschan.c:9858-9872`) with a length-preserving shift: rather than closing the gap it slides everything left of the zero one place right and writes a space into the vacated first character.

The shift is the right-justification, so it is correct whenever there is padding to consume.
When the two-digit-exponent form already exceeds the field width there is no padding, the space survives, and the value is written starting one column further right than every other value in the header, running one character past column 30.

Fix: consume a leading pad space when one exists, and otherwise close the gap with `memmove`.

This matches the idiom already used twenty lines below in the same function, where the `".0"` insertion guards its shift on `buf[ 0 ] == ' ' && buf[ 1 ] == ' ' ` (`fitschan.c:9913`).
Narrow values are unaffected; only the malformed wide cards change.

### B3. `WcsNative` stamps a vestigial `Invert` on a shared PermMap

When the celestial axes are not already at slots 0 and 1, `WcsNative` builds one `AstPermMap` and uses it twice (`fitschan.c:38396-38409`): once forward for stage 1, then inverted for stage 3.
`astCmpMap` clones rather than copies, since `CombineMaps` (`cmpmap.c:616`) copies only when the same pointer is supplied for both components of one CmpMap, which is not the case here.
The stage-1 CmpMap therefore holds the object that `astInvert` subsequently flips.

The stage-1 CmpMap records `invert1 = 0` at construction and dumps no `InvA`, yet its child PermMap dumps `Invert = 1`, a state it acquired after being encapsulated.
The flag never reaches the mathematics, because `Transform` forces the parent's recorded `invert1` onto the child before transforming (`cmpmap.c:3073`), but it is serialized, and a reader that trusted it would need the parent's `InvA` to cancel it.

Fix: use a copy of the PermMap for stage 3, so neither CmpMap's child carries a flag set after encapsulation.

### B4. `MapSz` records collision history: no change

`mapsize` changes only in `DoubleTableSize` (`keymap.c:2578`), which `AddTableEntry` calls when a bucket exceeds `MAX_ENTRIES_PER_TABLE_ENTRY` (`keymap.c:708`).
Nothing reduces it: `RemoveTableEntry` decrements `nentry[itab]` (`keymap.c:9036`) and stops.
`Dump` writes the current value unconditionally (`keymap.c:10510`).
Two KeyMaps with identical contents can therefore disagree on `MapSz` if one of them grew and shrank.

No change is made.
Emitting a content-derived size would change `MapSz`'s meaning from "current table size" to "table size implied by these keys", which is a deliberate semantic change rather than a defect repair.
This paragraph exists so the decision is not revisited.

## Group C: design-level

### C1. `IsSimple` tag placement depends on shared-pointer aliasing

`IsSimple` is one bit of `Mapping.flags` (`mapping.h:410`) carrying two meanings.
`astSimplify_` raises it on whatever it returns (`mapping.c:24759`) and short-circuits at the top for any Mapping whose bit is already set (`mapping.c:24740`), so the bit means both "dump `IsSimp = 1`" and "do not simplify again".
`astMAKE_SET(Mapping,Invert,...)` clears it as a side effect of every set, even when the stored value does not change (`mapping.c:23256`); `astMAKE_CLEAR` does not.

Probing routines decompose a candidate with `astMapList`, which returns reference-counted pointers to the live component objects rather than copies.
To line components up they call `astSetInvert` and restore the Invert value afterward, but not the `IsSimple` bit that `astSetInvert` incidentally cleared.
`cmpmap.c`'s `Equal` (`cmpmap.c:386-398`), `matrixmap.c`'s `CanSwap` (`416-421`, restored at `582-583`) and `winmap.c`'s `CanSwap` (`102-106`, restored at `201-202`) all follow this shape, and the save-flip-restore idiom recurs across at least a dozen further Mapping subclasses.

Which interior node ends up tagged in a dump is therefore decided by which nodes were aliased into a probe, not by the tree's shape.
The tag is advisory and does not affect any transform, so the consequence is confined to serialized output, but it makes that output a function of allocation and reference-counting history.

Two candidate fixes exist and they differ by an order of magnitude in footprint.

The first splits the flag: keep `AST__ISSIMPLE_FLAG` as the simplify guard that `astSetInvert` clears, and add a separate bit that only `astSimplify_` sets and only `Dump` reads.
This is confined to `mapping.c` and `mapping.h`, and leaves every probe as written.

The second makes probes non-mutating: save and restore `IsSimple` alongside `Invert` at each site.
This is truer to each routine's read-only contract, but touches a dozen or more files and each site needs individual judgment about whether it is a probe or a genuine mutation.

This design does not choose between them.
An investigation runs first: an instrumented `libast` logs every `astSetInvert` that clears `IsSimple`, together with its call site, and the full `ast_tester` suite is run against it.
The resulting distribution of clears selects the fix, and the choice is reported before any fix is written.

### C2. `CDjjjiii` cards are consumed under encodings that discard their values

`SpecTrans` reads each element of the legacy six-digit `CDjjjiii` matrix with `GetValue2` unconditionally, but stores the result as `PCj_i` plus `CDELTj = 1.0` only when `encoding == FITSIRAF_ENCODING` (`fitschan.c:31274-31300`).
`GetValue2` marks the card used when it finds it in the second FitsChan (`fitschan.c:18376`).
Under any other encoding the values are read, discarded, and the cards still fail to reach the output, so the rotation is silently lost rather than either honored or passed through.

Fix: move the `GetValue2` read inside the encoding test, so the cards survive under encodings that do not use them.

### C3. Duplicate-keyword consumption is inconsistent between the two read paths

`WcsFcRead` (`fitschan.c:36233`) is a card sweep: it clears the card index, walks every card, pattern-matches each keyword name against a template list, and marks every match.
Every copy of a repeated keyword is consumed.

`SpecTrans` (`fitschan.c:30880`) reads through `GetValue2`, which delegates to `GetValue`, `SearchCard` and `FindKeyCard`, a first-match forward search (`fitschan.c:18185`).
It marks exactly one card however many copies exist.
The keywords `SpecTrans` alone handles have no `WcsFcRead` branch, so they are only ever marked once.

For a header holding two HDUs' worth of WCS cards, a single `astRead` strips both copies of `CRPIXi`, `CRVALi`, `CTYPEi` and `EQUINOX` but leaves the second copy of `CDi_j`, `RADECSYS` and `DATE-OBS`.
What survives is a partial, self-inconsistent WCS description, and which keywords survive is decided by which read path happens to handle them.

A related asymmetry: `WcsFcRead`'s `MJD-OBS` branch sets `mark = 0` deliberately so the card survives a read, but when `DATE-OBS` is also present, `SpecTrans`'s translation-9 block (`fitschan.c:31703`) reaches `MJD-OBS` through `GetValue2` with marking enabled and consumes it after all.

Fix: bring `SpecTrans` in line with `WcsFcRead`.
Consumption marks every copy of a repeated keyword, and the deliberate `mark = 0` exception for `MJD-OBS` is respected on both paths.
`WcsFcRead`'s sweep is treated as the intended semantics because it is the more considered of the two: it carries explicit per-keyword exceptions, whereas `SpecTrans`'s single-copy marking is a consequence of the lookup helper it happens to use.

## Verification

New C tests are added under `ast_tester/` following the existing conversion patterns in `CLAUDE.md`, for the defects whose repair is observable through the public API:

- R1, that retrieving a negative `Double` or `Float` KeyMap entry as an integer rounds to nearest rather than toward zero (specified in full below);
- R1a, that a FITS grism header declaring a negative interference order reads back with that order intact;
- A1, that iteration over a KeyMap with `SortBy` set terminates;
- A2, that a rename refused by `MapLocked` leaves the entry intact under its original key;
- A3, that dumping a KeyMap does not crash and that a zero-length put is refused;
- A4 and A5, that out-of-range conversions produce clamped values rather than wrapped ones;
- A8, that `astMapGetElem<X>` reports failure for an undefined entry, matching `astMapGet0<X>`.

Every one of these tests is written and confirmed failing against the unfixed library before the corresponding fix is committed.
A test that passes in both states pins nothing.

### The R1 test

`ast_tester/testkeymap.c` already exercises this conversion, but only through a positive value: `1999.9` read back as `2000` at `testkeymap.c:559`.
Positive values round identically under both forms, which is why the defect survived the existing suite.
The new checks are added alongside that block, storing negative values and reading them back through `astMapGet0I`, `astMapGet0S` and `astMapGet0K`.

Values, with the result each produces before and after the fix, measured at `-O0`:

| stored | pre-fix | post-fix |
|---|---|---|
| `-1999.9` (Double) | -1999 | -2000 |
| `-1999.9f` (Float) | -1999 | -2000 |
| `-2.0` | -1 | -2 |
| `-2.5` | -2 | -3 |
| `-0.6` | 0 | -1 |
| `-0.4` | 0 | 0 |

`-2.0` is the most important of these: an exact integer that the current code does not return unchanged.
`-0.4` is a control that changes under neither form, included so the test would catch an over-correction that rounded away from zero unconditionally.

Failures report through the existing `stopit` helper, following the convention already used in the file.

These checks must be confirmed to fail against the unfixed library before the R1 fix is committed.
A test that passes both before and after demonstrates nothing, which is precisely the position the existing `1999.9` check is in.

The `Byte` arm is deliberately excluded from this test.
A negative value converted to `unsigned char` is the range-clamping concern of A4, not a rounding question, and it is tested there.

### The R1a test

Added to the FitsChan tests rather than `testkeymap.c`, since it exercises the FITS-WCS read path.

The test builds a minimal `WAVE-GRI` spectral header, sets `PV1_1` to a negative interference order, reads it with `astRead` and writes the result back through a FitsChan with `Encoding=FITS-WCS`, then compares the returned `PV1_1` with the value supplied.

Two cases carry the defect:

- `PV1_1 = -2.0` currently round-trips as `-1.0`, and must round-trip as `-2.0`;
- `PV1_1 = -1.0` currently causes `astRead` to fail outright, and must read successfully and round-trip as `-1.0`.

A positive order is included as a control, since positive orders round identically under both forms and must be unaffected.

The other grism parameters are held at values for which the header is physically usable, because `astRead` legitimately rejects some combinations of order and grating parameters.
The control case establishes that the chosen parameters are usable, so a failure in the negative cases cannot be attributed to the parameter choice.

A6, A7, A9 and B1 are covered by existing checkdump round-trips, which is where their effect appears.
B2, B3, C2 and C3 are verified by diffing `ast_tester` reference output.

Reference `.ast` and FITS fixtures are expected to change for B1, B2, B3 and C3.
Those diffs are reviewed before the corresponding commit is made, so that a fixture update is never a silent consequence of a source change.

R1 is expected to change no fixture, because the substitution is a no-op for non-negative inputs.
Any fixture it does change identifies a site that was rounding a negative value incorrectly, which is a finding worth reporting rather than a cost of the sweep.
Such diffs are reported rather than absorbed.

Every commit is gated on the sanitizer build staying clean, per `CLAUDE.md`:

```shell
cmake -B build-dev -DCMAKE_BUILD_TYPE=Debug \
  -DAST_ENABLE_WARNINGS=ON -DAST_ENABLE_SANITIZERS=ON
cmake --build build-dev
ctest --test-dir build-dev --output-on-failure
```

No change may introduce a new compiler warning.

## Out of scope

`CheckCircle`'s `AST__KYCIR` guard (`keymap.c:951`) was reviewed and is correct as written.
It is listed here only to record that it was examined and needs no change.
