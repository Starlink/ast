# Plot3D astGrid fails for RootCorner ULU and UUU with the testplot3d FrameSet

Found 2026-11-26 while writing the Plot3D logging tests (testplot3dlog).

With the FrameSet in `fixtures/programs/testplot3d/plot3d-test1.ast`
(RA, Dec, radio velocity; graphics box -1 to 1 on each axis; physical box
0.5 to 155.5, 0.5 to 107.5, 0.5 to 1640.5), `astGrid` on a Plot3D succeeds
for six of the eight `RootCorner` values and fails for `ULU` and `UUU`
with

```
astGrid(Plot): No grid curves can be drawn for axis 2.
```

reported with status `AST__INTER`, which `plot.c` reserves for internal
programming errors.
The same eight corners all succeed with the 3D FITS cube
`fixtures/wcsconv/inputs/skyspec3d.head` and with a plain 3D Frame.

Reproduce from `build-dev/ast_tester`:

```
./testplot3dlog --in fixtures/programs/testplot3d/plot3d-test1.ast \
    --box "-1 -1 -1 1 1 1" --bbox "0.5 0.5 0.5 155.5 107.5 1640.5" \
    --attr RootCorner=ULU --do grid
```

Not fixed here: the plan for the logging tests leaves `plot3d.c` alone.
The `rootcorner*` fixture rows use the FITS cube instead.

## Root cause (2026-11-26, from an instrumented tick hand-off)

Plot3D draws the three planes in turn.
The sky plane (XY here) is drawn first, and `Grid` (`src/plot3d.c`, the
two blocks around `astGetDrawnTicks` and `astSetTickValues`) then hands the
positions of the ticks it drew along the shared graphics axis to the next
plane, so that the planes' tick marks line up.
The hand-off assumes those ticks lie on the edge of the plane: it takes the
positions returned by `astGetDrawnTicks`, overwrites graphics axis 1 with
the plane constant, extends the set past both ends of the plane and
transforms it into the next plane's coordinates.

That assumption fails when the 2D Plot labels its axes in the interior.
Plot's `Grid` (`src/plot.c`, the block starting `labelling =
astGetLabelling`) uses exterior labelling only if `EdgeLabels` judges that
enough grid curves cross the chosen edges; otherwise it falls back to
interior labelling, and it cannot try the other edges because Plot3D has
set `Edge(1)` and `Edge(2)` explicitly for the corner.
With `RootCorner=ULU` the X-upper corner puts the Dec labels on the right
edge of the sky plane, which only two Dec curves cross for this FrameSet,
so the Plot labels both axes along interior curves.
The RA ticks are then on a diagonal (raw drawn ticks for ULU:
`(-0.59,0.88) ... (-0.03,-0.87)` against `(-0.86,-1) (-0.23,-1) (0.40,-1)`
for LLU), the extension runs along the wrong axis and produces duplicates,
and the XZ plane receives eight values between pixel 32 and 76 instead of
positions spanning the plane.
Its own grid along the other axis then has no room and `GridLines` reports
"No grid curves can be drawn for axis 2" with `AST__INTER`.

Evidence:

- `RootCorner=ULU,Labelling=interior` succeeds, and so does
  `RootCorner=LLL,Labelling=interior`: when the hand-off happens to span
  the plane the grid draws, but the ticks it passes on are the interior
  positions of the *other* sky axis, so even the "working" interior cases
  hand over the wrong ticks.
- `RootCorner=ULU,ForceExterior=1` and `RootCorner=UUU,ForceExterior=1`
  succeed with exterior labels (Dec labelled twice on the right edge).
  `ForceExterior` is the Plot attribute that makes it attempt exterior
  labels even when it judges them insufficient.
- The FITS cube's sky grid is nearly aligned with the plane, so every edge
  carries enough crossings and exterior labelling never falls back.
- `Direction(1)` plays no part.

### Size of the fix

Small, two parts, in `src/plot3d.c`:

1. In `Grid`, set `ForceExterior=1` on each encapsulated Plot before its
   `astGrid` (and clear it afterwards unless the user set it), since the
   hand-off is only meaningful for ticks on the plane edges.  About ten
   lines.  With it, all eight corners draw for `plot3d-test1.ast`.
2. In the hand-off, check that the drawn tick positions are constant on
   the graphics axis being replaced; if they are not (a Plot that still
   labelled in the interior, or `Labelling=interior` set by the user), skip
   `astSetTickValues` so the next plane computes its own ticks instead of
   receiving positions from an interior curve.  About ten lines, and it
   turns the `AST__INTER` error into a drawable, if unaligned, grid.

Then regenerate the `rootcorner*` references (they can return to
`plot3d-test1.ast`) and add a `Labelling=interior` row, and tell the Rust
port, whose Plot3D copies this hand-off.
