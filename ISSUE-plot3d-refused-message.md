# Plot3D's refused-method errors name the method as "ast##name"

Found 2026-11-26 while writing the Plot3D logging tests (testplot3dlog).

The `METHOD1` and `METHOD2` macros near the end of `src/plot3d.c` build the
error for a Plot method a Plot3D does not support with

```c
astError( AST__INTER, "ast##name(%s): The ast##name " \
          "method cannot be used with a %s (programming error).", ...
```

The `##` operator does not paste inside a string literal, so the message
reads literally `ast##name(Plot3D): The ast##name method cannot be used
with a Plot3D (programming error).` for `astBoundingBox`, `astClip`,
`astCurve`, `astGenCurve` and `astGridLine`.
Only `astPolyCurve` and `astSetTickValues`, written out by hand, name
themselves.
The fix is a stringising helper (`#name` through a two-level macro) in
place of the `##` inside the literal.

Not fixed here: the plan for the logging tests leaves `plot3d.c` alone.
`fixtures/plot3d/expected/refused.log` records the message as it is today
and will need regenerating when this is fixed.

## Analysis (2026-11-26)

Trivial: about ten lines.
Replace the `##` inside the string literal with a stringising macro, for
example

```c
#define STR2(x) #x
#define STR(x) STR2(x)
#define METHOD1(name) \
static void name(ARGLIST,int *status){ \
   if( !astOK ) return; \
   astError( AST__INTER, "ast" STR(name) "(%s): The ast" STR(name) " " \
             "method cannot be used with a %s (programming error).", status, \
             astGetClass( this ), astGetClass( this ) ); \
}
```

and the same in `METHOD2`, with a History entry.
Then regenerate `fixtures/plot3d/expected/refused.log`, which pins the
message text, and tell the Rust port.
