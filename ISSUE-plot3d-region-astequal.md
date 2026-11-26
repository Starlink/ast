# Plot3D over a Region is not astEqual to its own round trip

Found 2026-11-26 while writing the Plot3D logging tests (testplot3dlog).

A Plot3D whose supplied Frame is a 3D Region (a Box) is written and read
back unchanged (the two Native dumps are byte-identical), yet `astEqual`
reports the original and the reloaded Plot3D as different.
A 2D Plot over a 2D Box, a FrameSet over the Box, and a Plot3D over a plain
3D Frame all compare equal after the same round trip.

```c
AstFrame *f3 = astFrame( 3, "Domain=PIXEL" );
double lbnd[3] = { 20.0, 20.0, 300.0 }, ubnd[3] = { 120.0, 80.0, 1200.0 };
AstBox *box3 = astBox( f3, 1, lbnd, ubnd, NULL, " " );
float gbox[6] = { 0, 0, 0, 1, 1, 1 };
double bbox[6] = { 0.5, 0.5, 0.5, 155.5, 107.5, 1640.5 };
AstPlot3D *p = astPlot3D( box3, gbox, bbox, " " );
char *s1 = astToString( p );
AstObject *copy = astFromString( s1 );
char *s2 = astToString( copy );
/* strcmp( s1, s2 ) == 0, but astEqual( p, copy ) == 0 */
```

| Object | astEqual after round trip | dumps equal |
|---|---|---|
| Box3 | 1 | 1 |
| FrameSet(Box3) | 1 | 1 |
| Plot(Box2) | 1 | 1 |
| Plot3D(Frame3) | 1 | 1 |
| Plot3D(Box3) | 0 | 1 |
| Plot3D(FrameSet(Box3)) | 0 | 1 |

Plot3D's `Equal` (src/plot3d.c) compares the three component Plots with
`astEqual` after the parent FrameSet comparison passes, so one of
`plotxy`, `plotxz` or `plotyz` compares unequal to its reloaded copy even
though it dumps identically.
The component Plots are built by `SplitFrameSet` from the Region, so the
2D Regions picked out of the 3D Box are the likely place to look.

Not fixed here: the plan for the logging tests leaves `plot3d.c` alone.
The logging fixtures therefore build their Plot3D from a plain 3D Frame
(`fixtures/plot3d/inputs/frame3d.ast`) and pass the Box to the `outline`
action separately.

## Analysis (2026-11-26, after the logging tests landed)

The defect is not in `plot3d.c`.
Descending through the dump with the protected interface shows that in a
Plot3D over a Box the XY component Plot compares equal, while the XZ and YZ
component Plots do not: their third Frame, a CmpFrame that `SplitFrameSet`
builds from a one-axis Box picked out of the input Region and a dummy
one-axis Box, is unequal to its reloaded copy even though each of its two
component Boxes compares equal and its dump is byte-identical.

The CmpFrame class has no `Equal` of its own, so Frame's `Equal` applies.
Frame's `Equal` calls `astConvert` on the two Frames, asks the resulting
two-Frame FrameSet for its base-to-current Mapping, and accepts the Frames
as equal only if that Mapping simplifies to a UnitMap.

The Mapping does not simplify by design, not by accident.
A Frame used as a Mapping is a unit Mapping, except that a Region (and so
a CmpFrame containing one) masks: positions outside the Region come out
bad (`src/region.c`, class description; `src/frame.c:116`).
`FrameSet::Span`, which `astGetMapping` uses to build the path between two
Frames (`src/frameset.c:10404-10412` and `10450-10492`), therefore inserts
every Frame along the path whose `astIsUnitFrame` is zero into the Mapping
it returns, so that transforming through a FrameSet honours the masking.
`Region::IsUnitFrame` always returns 0 (`src/region.c:5065`) and
`CmpFrame::IsUnitFrame` is the AND of its components (`src/cmpframe.c:5046`).
The conversion between the two CmpFrames of Boxes is thus
`CmpMap( CmpMap( CmpFrame_a, UnitMap ), CmpFrame_b )`, and nothing in
`astSimplify` turns a Region-as-Mapping into a UnitMap because it is not
one; `astRemoveRegions` exists precisely to strip that masking on request,
and applying it to this Mapping does give a UnitMap.

So the design choice is deliberate and documented; the missing piece is
on the `Equal` side.
Frame's `Equal` predates the idea of a non-unit Frame and silently assumes
that identical Frames convert through a UnitMap.
The Region class already works around this with its own `Equal`
(`src/region.c:2746`), comparing points, both encapsulated Frames, the
encapsulated Mapping, Negated and Closed directly.
CmpFrame never received the matching override, so any CmpFrame with a
Region component compares unequal to everything, including itself after a
copy.
The public API reproduces it without any Plot3D:

```c
AstBox *ba = astBox( astFrame( 1, "Domain=PIXEL" ), 1, &l1, &u1, NULL, " " );
AstBox *bb = astBox( astFrame( 1, "Domain=PIXEL" ), 1, &l2, &u2, NULL, " " );
AstCmpFrame *cf1 = astCmpFrame( ba, bb, " " );
AstCmpFrame *cf2 = astCmpFrame( ba, bb, " " );
astEqual( cf1, cf2 );   /* 0, though built from the same objects */
```

| CmpFrame | astEqual after round trip |
|---|---|
| CmpFrame(Frame, Frame) | 1 |
| CmpFrame(Box, Box) | 0 |
| CmpFrame(Box2D, Frame) | 0 |
| CmpFrame(Frame, Box) | 0 |

Region has its own `Equal` (points, both Frames and the Mapping of its
FrameSet, Negated, Closed), which is why a lone Box and a Plot over a Box
compare fine; only a Region nested inside a CmpFrame goes through the
`astConvert` route.

### Size of the fix

Small.
Give CmpFrame an `Equal` override in `src/cmpframe.c`, modelled on the one
in `src/cmpmap.c`: same class, `astEqual` on `frame1` and `frame2`, and
the same axis permutation (`perm`).
About 40 lines plus a History entry, and a C test (`testcmpframe.c` has the
right home) asserting that two CmpFrames built from identical Regions, and
a CmpFrame of Regions and its round trip, are `astEqual`.
Once that lands, a Plot3D over a Region round-trips equal and the
`norm`, `markers`, `outline` and `svg_behind` fixture rows can go back to
building their Plot3D over `plot3d/inputs/box3d.ast` if wanted (their logs
would change, since the current Frame would be the Box).
Two alternatives are worse.
Applying `astRemoveRegions` inside Frame's `Equal` before the UnitMap test
would make a CmpFrame over the Box 0 to 1 equal to one over the Box 0 to 2,
since it discards exactly the part that differs.
Teaching a Region-as-Mapping to simplify to a UnitMap would break the
documented masking behaviour of FrameSets that hold Regions.
