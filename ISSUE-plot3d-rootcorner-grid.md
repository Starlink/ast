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

## Analysis (2026-11-26, after the logging tests landed)

Facts established with the logging module and a temporary print in
`Grid` (`src/plot3d.c`, at the `astSetTickValues` hand-off):

- The base XY plane (the sky plane) draws completely for every corner;
  the error is raised by `astGrid` on the second plane, XZ, for its
  second WCS axis.
- The failing corners are exactly those with the X bit set together with
  the Z bit (`ULU`, `UUU`).  `Direction(1)=0` on the sky axis is not the
  cause: forcing `Direction(1)=1` still fails, and the FITS cube with
  `Direction(1)=0` still succeeds for all eight corners.
- The tick values handed from the XY plane to the XZ plane differ in
  kind.  For `LLU` the XZ plane receives 5 major values spanning
  -30.5 to 186.5, the extended pixel positions along the shared X
  graphics axis.  For `ULU` it receives 8 values between 32.1 and 75.7,
  which do not span the plane, so the XZ plane's own grid along the
  other axis has no room and `GridLines` reports zero curve length.
  With the FITS cube the `ULU` hand-off is sane (-19.5 to 120.5).

`Grid` finds the ticks to hand over with `astGetDrawnTicks( baseplot,
base_wax2d, 1 )`, where `base_wax2d` is chosen from the base Plot's
`Edge(1)` attribute, then overwrites graphics axis 1 of every tick
position with the plane constant `gcon` on the assumption that the ticks
vary along graphics axis 0.
For `ULU` the labelled edges chosen for the X-upper corner put the sky
axis ticks on an edge where that assumption no longer holds for this
FrameSet, so the positions are corrupted before they are transformed into
the XZ plane's WCS.
Why the test1 FrameSet and the FITS cube differ here (both are RA, Dec
plus a spectral axis) is not yet pinned down; the DSBSpecFrame velocity
axis and the FrameSet's reversed sky axis are the visible differences.

### Size of the fix

Medium: a focused debugging session of a few hours.
The change itself is likely small, in the two tick hand-off blocks of
`Grid` (`src/plot3d.c` around lines 3470 to 3550 and 3575 to 3655): decide
which graphics axis the drawn ticks vary along from the tick positions
themselves (as `ExtendTicks` already does) instead of inferring it from
`Edge(1)`, and overwrite the other axis with `gcon`.
It needs the eight `rootcorner*` rows run against `plot3d-test1.ast` as
well as the FITS cube, since the cube never showed the problem.
