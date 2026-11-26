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
