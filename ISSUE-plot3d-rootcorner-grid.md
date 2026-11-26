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
