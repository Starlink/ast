# JsonChan Phase 1 (Skeleton + Write Path) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `JsonChan` AST class that can *write* (serialize) any AST object to lossless JSON-native text, leaving the read path and schema generator to later plans.

**Architecture:** `JsonChan` is a new direct subclass of `Channel`, a sibling of `XmlChan`/`YamlChan`, created by templating the existing `yamlchan.{c,h}` boilerplate. It overrides the `Channel` `WriteXxx` virtual callbacks to build an `AstKeyMap` node tree, then emits that tree as JSON text through the channel sink. The ~70 class `Dump`/`Load` functions and the native `Channel` are untouched.

**Tech Stack:** C99 (AST source), C11 (serialization tests), CMake build, AST object system (`astBegin`/`astEnd`, reference counting, `astOK` status), `AstKeyMap` as the intermediate representation.

**Design reference:** `docs/superpowers/specs/2026-06-27-jsonchan-design.md`

## Global Constraints

- Build and test with the CMake dev configuration: `cmake -B build-dev -DCMAKE_BUILD_TYPE=Debug -DAST_ENABLE_WARNINGS=ON -DAST_ENABLE_SANITIZERS=ON` then `cmake --build build-dev`.
- Changes must add **no new compiler warnings** and must pass the sanitizer (ASan/UBSan) tests.
- AST C source is **C99, no C++ features**.
- `testjsonchan.c` must compile as **C11** (the project requires C11 for serialization tests because `AST__DBL_DIG = 18` under C11 vs 17 under C99).
- Tests use the **public C API** (`ast.h`) only, with `astWatch(&status)` + `astOK` for error checking (no Starlink EMS).
- Use `astMalloc`/`astFree`, `astBegin`/`astEnd`, and AST reference counting.
- When modifying any file under `src/`, **update the prologue History section** with the change.
- No external runtime dependencies in Phase 1 (the write path emits text directly; cJSON is a Phase 2 concern).
- The `WHOLE_ARCHIVE` linking pattern in `ast_tester/CMakeLists.txt` is required for tests to link satellite libraries; new tests must follow it.
- Function naming: `astFoo` (public macro), `Foo`/`astFoo_` (internal). Public class predicate `astIsAJsonChan`.
- JSON format rules from the spec: discriminator field is `"$type"` (class name); `AST__BAD` and any non-finite double serialize as JSON `null`; doubles formatted with AST's own `%.*g` at `AST__DBL_DIG`; integers as JSON numbers (no booleans); child objects as nested JSON objects keyed by their slot name; vectors as JSON arrays. Document-root object additionally carries `"$schema"`, `"$schemaVersion"`, `"$minReadVersion"`; nested objects carry only `"$type"` + attributes.
- `AST__SCHEMA_VERSION` starts at `"1.0.0"`; `$minReadVersion` is `1`.

---

## File Structure

- `src/jsonchan.h` — public + protected interface for the JsonChan class (templated from `yamlchan.h`). Also defines `AST__SCHEMA_VERSION` and `AST__JSON_MINREADVER`.
- `src/jsonchan.c` — class implementation: boilerplate (vtab, init, load) templated from `yamlchan.c`, plus the hand-written write path (KeyMap build + JSON emission).
- `src/fjsonchan.c` — Fortran 77 interface (templated from `fyamlchan.c`).
- `ast_tester/testjsonchan.c` — C11 test program for the write path.
- `CMakeLists.txt` — add header to `AST_H_FILES` and the install header list, add `src/jsonchan.c` to the source list, add `src/fjsonchan.c` to the Fortran source list.
- `src/loader.c` — add `LOAD(JsonChan)`.
- `ast_err.msg` — add the `BDJSN` error code.
- `ast_tester/CMakeLists.txt` — register `testjsonchan` as a test.

---

## Task 1: JsonChan class skeleton that builds and constructs

Produce a minimal but complete `Channel` subclass (C + Fortran) that compiles into `libast`, with a working `astJsonChan` constructor and `astIsAJsonChan` predicate. No serialization behaviour yet — `astWrite` may fall back to the parent for now.

**Files:**
- Create: `src/jsonchan.h`
- Create: `src/jsonchan.c`
- Create: `src/fjsonchan.c`
- Create: `ast_tester/testjsonchan.c`
- Modify: `CMakeLists.txt` (lines ~440 `AST_H_FILES`, ~445 install headers, ~551 sources, ~701 Fortran sources)
- Modify: `src/loader.c:189` (the `LOAD(...)` block)
- Modify: `ast_err.msg` (insert `BDJSN` alphabetically near the other `BD*` codes)
- Modify: `ast_tester/CMakeLists.txt`

**Interfaces:**
- Produces:
  - `AstJsonChan *astJsonChan( const char *(* source)( void ), void (* sink)( const char * ), const char *options, ... )` — public constructor (macro `astJsonChan`).
  - `int astIsAJsonChan( AstJsonChan *this )` — class predicate.
  - `AstJsonChan *astInitJsonChan_( void *mem, size_t size, int init, AstJsonChanVtab *vtab, const char *name, const char *(* source)( void ), char *(* source_wrap)( const char *(*)( void ), int * ), void (* sink)( const char * ), void (* sink_wrap)( void (*)( const char * ), const char *, int * ), int * )` — initialiser (signature mirrors `astInitYamlChan_`).
  - `AstJsonChan *astLoadJsonChan_( void *mem, size_t size, AstJsonChanVtab *vtab, const char *name, AstChannel *channel, int * )` — loader.
  - Struct `AstJsonChan` and vtab `AstJsonChanVtab` (mirroring the YamlChan structs, minus YAML-specific members).
  - Header guard symbol `JSONCHAN_INCLUDED`.

- [ ] **Step 1: Template the header from yamlchan.h**

Copy `src/yamlchan.h` to `src/jsonchan.h`. Then apply these exact transformations:
- Replace every `YamlChan` → `JsonChan`, `YAMLCHAN` → `JSONCHAN`, `yamlchan` → `jsonchan`, `Yaml` → `Json`.
- Delete the `#include "keymap.h"` line is **kept** (we need KeyMap); delete the `#if defined( YAML ) #include <yaml.h> #endif` block.
- Delete any YAML/ASDF-specific attribute declarations and accessors (`VerboseRead`, `PreserveName`, `YamlEncoding` and their `astGet*/astSet*/astClear*/astTest*` prototypes and vtab slots). The Phase-1 JsonChan adds no new attributes; it relies only on inherited `Channel` attributes.
- Update the prologue: Purpose "Define the interface to the JsonChan class", Description mentions JSON, add a `History:` entry dated today with author and "Original version."
- After the `#endif` guard but inside the protected (`astCLASS`) section is not needed; instead, add near the top of the public macro section:

```c
/* Schema version reported in JSON output written by a JsonChan. */
#define AST__SCHEMA_VERSION "1.0.0"
#define AST__JSON_MINREADVER 1
```

- [ ] **Step 2: Template the implementation from yamlchan.c**

Copy `src/yamlchan.c` to `src/jsonchan.c`. Apply the same name substitutions as Step 1. Then reduce it to a compiling skeleton:
- Delete all YAML/ASDF-specific code: every `ReadYAML*`, `WriteAsdf*`, `ReadValues`, `ReadNative`, `IsAsdf*`, `StartAsdf*`, `AddR2D`, the `#include <yaml.h>` and any `#if defined(YAML)` blocks, and the YAML-specific attribute machinery (`VerboseRead`/`PreserveName`/`YamlEncoding`).
- Keep and rename the structural boilerplate: includes, `astMAKE_*` for the class, `astInitJsonChanVtab`, `astInitJsonChan_`, `astLoadJsonChan_`, `astJsonChan_` (the C constructor body), the `Delete`/`Copy`/`Dump` for the class itself, and the global vtab handling.
- For the overridable `Write`/`Read` methods: in this task, do **not** override them yet — leave the parent `Channel` implementations in place (i.e. the vtab does not reassign `WriteBegin`/`WriteDouble`/etc.). This makes the skeleton compile and construct without serialization behaviour.
- Update the prologue Purpose ("I/O Channel that uses JSON to represent Objects"), Description, and add today's `History:` entry.

- [ ] **Step 3: Template the Fortran interface**

Copy `src/fyamlchan.c` to `src/fjsonchan.c`, apply the same name substitutions, delete any YAML-attribute-specific wrappers, and update the prologue + History. The result must expose `AST_JSONCHAN`.

- [ ] **Step 4: Wire the build — headers**

In `CMakeLists.txt`, add `src/jsonchan.h` to both header lists. After line 441 (`src/yamlchan.h`) in `AST_H_FILES`:

```cmake
    src/channel.h
    src/yamlchan.h
    src/jsonchan.h
```

and add the same `src/jsonchan.h` line to the install-headers list near line 445 (next to `src/xmlchan.h`).

- [ ] **Step 5: Wire the build — sources**

In `CMakeLists.txt`, after line 551 (`src/yamlchan.c`) add:

```cmake
    src/yamlchan.c
    src/jsonchan.c
```

and after the Fortran source `src/fyamlchan.c` (line ~701) add `src/fjsonchan.c`:

```cmake
        src/fyamlchan.c
        src/fjsonchan.c
```

- [ ] **Step 6: Register the loader**

In `src/loader.c`, add `LOAD(JsonChan);` in alphabetical position (between `LOAD(IntraMap);` at line 149 and `LOAD(KeyMap);` at line 150):

```c
   LOAD(IntraMap);
   LOAD(JsonChan);
   LOAD(KeyMap);
```

Update the `src/loader.c` prologue History with today's entry.

- [ ] **Step 7: Add the error code**

In `ast_err.msg`, add a JSON error code among the `BD*` entries (after `BDFTS`):

```
BDJSN           <bad JSON input or output>
```

- [ ] **Step 8: Write the construction test**

Create `ast_tester/testjsonchan.c`:

```c
/* Test the JsonChan class (Phase 1: construction + write path).
 *
 * Differences from any Fortran original: none (new C test). Uses the
 * public C API and astWatch/astOK error checking. Compiled as C11 so
 * AST__DBL_DIG == 18 for serialization precision.
 */
#include "ast.h"
#include <stdio.h>
#include <string.h>

int main( void ) {
   int status = 0;
   astWatch( &status );

   astBegin;

   AstJsonChan *chan = astJsonChan( NULL, NULL, " " );
   if ( !astOK ) { printf( "astJsonChan failed\n" ); return 1; }
   if ( !astIsAJsonChan( chan ) ) { printf( "not a JsonChan\n" ); return 1; }

   char buf[ 50 ];
   astGetC( chan, "Class" ); /* smoke: attribute access works */
   (void) buf;

   astEnd;

   if ( !astOK ) { printf( "status set at end\n" ); return 1; }
   printf( "All JsonChan construction tests passed\n" );
   return 0;
}
```

- [ ] **Step 9: Register the test**

In `ast_tester/CMakeLists.txt`, add `ast_add_test(testjsonchan)` alongside the other channel tests (e.g. near `testyamlchan`/`testxmlchan`). Match the surrounding pattern exactly, including any C11 standard property the other serialization tests set. If serialization tests set `C_STANDARD 11`, do the same for `testjsonchan`.

- [ ] **Step 10: Build and run**

Run:
```bash
cmake -B build-dev -DCMAKE_BUILD_TYPE=Debug -DAST_ENABLE_WARNINGS=ON -DAST_ENABLE_SANITIZERS=ON
cmake --build build-dev 2>&1 | tee /tmp/jsonbuild.log
ctest --test-dir build-dev -R testjsonchan --output-on-failure
```
Expected: build succeeds with no new warnings referencing `jsonchan`; `testjsonchan` passes printing "All JsonChan construction tests passed".

- [ ] **Step 11: Commit**

```bash
git add src/jsonchan.h src/jsonchan.c src/fjsonchan.c ast_tester/testjsonchan.c \
        CMakeLists.txt src/loader.c ast_err.msg ast_tester/CMakeLists.txt
git commit -m "feat(jsonchan): add JsonChan class skeleton that builds and constructs"
```

---

## Task 2: Write a flat scalar object to JSON

Override the scalar `WriteXxx` callbacks to accumulate into an `AstKeyMap`, and on `WriteEnd` of the outermost object emit the KeyMap as JSON text through the sink. Handle `WriteBegin` (sets `$type`), `WriteInt`, `WriteInt64`, `WriteDouble`, `WriteString`, and `WriteEnd`.

**Files:**
- Modify: `src/jsonchan.c`
- Modify: `ast_tester/testjsonchan.c`

**Interfaces:**
- Consumes: `AstKeyMap` API (`astKeyMap`, `astMapPut0D`, `astMapPut0I`, `astMapPut0C`, `astMapGet0*`, `astMapKey`, `astMapSize`, `astMapType`), and the inherited `Channel` sink (`astPutNextText`).
- Produces (file-internal, static): `static void WriteBegin( AstChannel *, const char *class, const char *comment, int * )`, `static void WriteDouble( AstChannel *, const char *name, int set, int helpful, double value, const char *comment, int * )`, `static void WriteInt(...)`, `static void WriteInt64(...)`, `static void WriteString(...)`, `static void WriteEnd( AstChannel *, const char *class, int * )`, plus a helper `static void EmitJson( AstJsonChan *this, AstKeyMap *km, int *status )` and `static char *FormatDouble( double value, char *buf, int *status )`.

- [ ] **Step 1: Write the failing test**

Add to `ast_tester/testjsonchan.c` a helper sink that captures output into a global buffer, and a test that writes a `ZoomMap`:

```c
static char jsonbuf[ 8192 ];
static void Sink( const char *line ) {
   strncat( jsonbuf, line, sizeof( jsonbuf ) - strlen( jsonbuf ) - 1 );
   strncat( jsonbuf, "\n", sizeof( jsonbuf ) - strlen( jsonbuf ) - 1 );
}

static int test_flat_scalar( void ) {
   jsonbuf[ 0 ] = '\0';
   AstJsonChan *chan = astJsonChan( NULL, Sink, " " );
   AstZoomMap *zm = astZoomMap( 2, 5.0, " " );
   astWrite( chan, zm );
   if ( !astOK ) return 1;
   if ( !strstr( jsonbuf, "\"$type\"" ) ) { printf("no $type\n"); return 1; }
   if ( !strstr( jsonbuf, "ZoomMap" ) )   { printf("no class\n"); return 1; }
   if ( !strstr( jsonbuf, "\"Zoom\"" ) )  { printf("no Zoom\n"); return 1; }
   if ( jsonbuf[ 0 ] != '{' )             { printf("not an object\n"); return 1; }
   return 0;
}
```

Call `test_flat_scalar()` from `main` (after the construction test) and return 1 on non-zero.

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build-dev && ctest --test-dir build-dev -R testjsonchan --output-on-failure`
Expected: FAIL — the skeleton does not yet emit JSON, so `jsonbuf` is empty (`not an object`).

- [ ] **Step 3: Implement the KeyMap-building write callbacks**

In `src/jsonchan.c`, add a per-instance current-KeyMap stack to the `AstJsonChan` struct (a `AstKeyMap **km_stack`, `int km_depth`) in `jsonchan.h`, initialised to NULL/0 in `astInitJsonChan_`. Implement the callbacks (register them in `astInitJsonChanVtab` by overriding the parent slots):

```c
static void WriteBegin( AstChannel *this_channel, const char *class,
                        const char *comment, int *status ) {
   AstJsonChan *this = (AstJsonChan *) this_channel;
   AstKeyMap *km = astKeyMap( " " );
   astMapPut0C( km, "$type", class, NULL );
   this->km_stack = astGrow( this->km_stack, this->km_depth + 1,
                             sizeof( AstKeyMap * ) );
   if ( astOK ) this->km_stack[ this->km_depth++ ] = km;
}

static void WriteDouble( AstChannel *this_channel, const char *name, int set,
                         int helpful, double value, const char *comment,
                         int *status ) {
   AstJsonChan *this = (AstJsonChan *) this_channel;
   if ( !set ) return;                       /* match native: only set attrs */
   AstKeyMap *km = this->km_stack[ this->km_depth - 1 ];
   if ( value == AST__BAD || !astISFINITE( value ) ) {
      astMapPutU( km, name, NULL );          /* null sentinel for bad */
   } else {
      astMapPut0D( km, name, value, NULL );
   }
}

static void WriteInt( AstChannel *this_channel, const char *name, int set,
                      int helpful, int value, const char *comment,
                      int *status ) {
   AstJsonChan *this = (AstJsonChan *) this_channel;
   if ( !set ) return;
   astMapPut0I( this->km_stack[ this->km_depth - 1 ], name, value, NULL );
}

/* WriteInt64 identical but astMapPut0K / int64_t; WriteString uses astMapPut0C. */

static void WriteEnd( AstChannel *this_channel, const char *class, int *status ) {
   AstJsonChan *this = (AstJsonChan *) this_channel;
   AstKeyMap *km = this->km_stack[ --this->km_depth ];
   if ( this->km_depth == 0 ) {              /* outermost object: emit */
      EmitJson( this, km, status );
   }
   km = astAnnul( km );
}
```

Notes: `astMapPutU` stores a valueless (undefined) entry; the emitter renders it as JSON `null`. `astISFINITE` and `AST__BAD` come from `pointset.h` (already included transitively via `channel.h`).

- [ ] **Step 4: Implement the JSON emitter and double formatter**

```c
static char *FormatDouble( double value, char *buf, int *status ) {
   (void) sprintf( buf, "%.*g", AST__DBL_DIG, value );
   return buf;
}

static void EmitJson( AstJsonChan *this, AstKeyMap *km, int *status ) {
   AstChannel *chan = (AstChannel *) this;
   int n = astMapSize( km );
   astPutNextText( chan, "{" );
   for ( int i = 0; i < n && astOK; i++ ) {
      const char *key = astMapKey( km, i );
      char line[ 256 ], dbuf[ 64 ];
      int type = astMapType( km, key );
      if ( type == AST__UNDEFTYPE ) {
         (void) sprintf( line, "  \"%s\": null%s", key, i < n-1 ? "," : "" );
      } else if ( type == AST__DOUBLETYPE ) {
         double d; astMapGet0D( km, key, &d );
         (void) sprintf( line, "  \"%s\": %s%s", key,
                         FormatDouble( d, dbuf, status ), i < n-1 ? "," : "" );
      } else if ( type == AST__INTTYPE ) {
         int v; astMapGet0I( km, key, &v );
         (void) sprintf( line, "  \"%s\": %d%s", key, v, i < n-1 ? "," : "" );
      } else { /* AST__STRINGTYPE */
         const char *s; astMapGet0C( km, key, &s );
         (void) sprintf( line, "  \"%s\": \"%s\"%s", key, s, i < n-1 ? "," : "" );
      }
      astPutNextText( chan, line );
   }
   astPutNextText( chan, "}" );
}
```

(Nested objects and vectors are added in Task 3; string escaping is hardened in Task 4. `AST__UNDEFTYPE`, `AST__DOUBLETYPE`, etc. are from `keymap.h`.)

- [ ] **Step 5: Run the test to verify it passes**

Run: `cmake --build build-dev && ctest --test-dir build-dev -R testjsonchan --output-on-failure`
Expected: PASS.

- [ ] **Step 6: Verify the output is valid JSON**

Run (manual sanity, not a committed test yet):
```bash
ctest --test-dir build-dev -R testjsonchan -V 2>/dev/null | sed -n '/^{/,/^}/p' | python3 -m json.tool
```
Expected: `python3 -m json.tool` parses it without error.

- [ ] **Step 7: Commit**

```bash
git add src/jsonchan.c src/jsonchan.h ast_tester/testjsonchan.c
git commit -m "feat(jsonchan): emit flat scalar objects as JSON"
```

---

## Task 3: Nested child objects and vector attributes

Extend the write path to handle `WriteObject` (nested AST objects → nested JSON objects keyed by slot name) and vector-valued attributes (`WriteDouble`/`WriteInt`/`WriteString` with an index > 0, used by classes such as `WinMap` for `Ina`/`Inb`) → JSON arrays.

**Files:**
- Modify: `src/jsonchan.c`
- Modify: `ast_tester/testjsonchan.c`

**Interfaces:**
- Consumes: `AstKeyMap` nested-map support (`astMapPut0A` to store a child KeyMap; vector puts `astMapPut1D`/`astMapPut1I`/`astMapPut1C`).
- Produces: `static void WriteObject( AstChannel *, const char *name, int set, int helpful, AstObject *value, const char *comment, int * )`; updated `EmitJson` recursion for nested KeyMaps and arrays.

- [ ] **Step 1: Write the failing test**

Add to `ast_tester/testjsonchan.c`:

```c
static int test_nested_and_vectors( void ) {
   jsonbuf[ 0 ] = '\0';
   AstJsonChan *chan = astJsonChan( NULL, Sink, " " );
   AstZoomMap *zm = astZoomMap( 2, 5.0, " " );
   double ina[2] = { 0.0, 1.0 }, inb[2] = { 1.0, 2.0 };
   AstWinMap *wm = astWinMap( 2, ina, inb, ina, inb, " " );
   AstCmpMap *cm = astCmpMap( zm, wm, 1, " " );   /* series */
   astWrite( chan, cm );
   if ( !astOK ) return 1;
   if ( !strstr( jsonbuf, "CmpMap" ) )  { printf("no CmpMap\n"); return 1; }
   if ( !strstr( jsonbuf, "ZoomMap" ) ) { printf("no nested ZoomMap\n"); return 1; }
   if ( !strstr( jsonbuf, "[" ) )       { printf("no array\n"); return 1; }
   return 0;
}
```

Call it from `main`. (Use whatever child slot names `CmpMap`'s `Dump` actually uses; the test only checks that the nested class name and an array appear.)

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build-dev && ctest --test-dir build-dev -R testjsonchan --output-on-failure`
Expected: FAIL — `WriteObject` is not overridden, so nested objects are missing.

- [ ] **Step 3: Implement WriteObject**

`astDump(value, this)` drives the child object's own `WriteBegin` ... `WriteEnd` callbacks. Because our `WriteEnd` only *emits* at depth 0 and otherwise hands the popped child KeyMap to `this->pending_child` (see the `WriteEnd` change below), `WriteObject` just needs to trigger the recursion and then attach the resulting child KeyMap to the parent under `name`:

```c
static void WriteObject( AstChannel *this_channel, const char *name, int set,
                         int helpful, AstObject *value, const char *comment,
                         int *status ) {
   AstJsonChan *this = (AstJsonChan *) this_channel;
   if ( !set ) return;
   AstKeyMap *parent = this->km_stack[ this->km_depth - 1 ];
   /* astDump drives WriteBegin (pushes child km) ... WriteEnd (pops it). */
   int depth_before = this->km_depth;
   astDump( value, this_channel );
   /* WriteEnd popped the child into `this->pending_child`; attach it. */
   astMapPut0A( parent, name, this->pending_child, NULL );
   this->pending_child = astAnnul( this->pending_child );
   (void) depth_before;
}
```

Adjust `WriteEnd` so that when it pops a non-outermost object it stores the popped KeyMap in `this->pending_child` (a new struct member) instead of annulling it:

```c
static void WriteEnd( AstChannel *this_channel, const char *class, int *status ) {
   AstJsonChan *this = (AstJsonChan *) this_channel;
   AstKeyMap *km = this->km_stack[ --this->km_depth ];
   if ( this->km_depth == 0 ) {
      EmitJson( this, km, status );
      km = astAnnul( km );
   } else {
      this->pending_child = km;   /* handed to WriteObject; not annulled */
   }
}
```

Add `AstKeyMap *pending_child;` to the struct in `jsonchan.h`, initialised NULL.

- [ ] **Step 4: Implement vector handling in the callbacks**

The `Channel` vector convention: `Dump` writes vector elements by calling `WriteDouble`/`WriteInt`/`WriteString` repeatedly with names like `"Ina1"`, `"Ina2"` (a base name plus a 1-based index). Inspect `src/winmap.c`'s `Dump` to confirm the exact naming, then aggregate consecutive same-base entries into a KeyMap vector. Simplest robust approach: store each scalar as written, and in `EmitJson` leave them as separate keys — but that does not produce arrays. Instead, detect the trailing-integer pattern at write time:

```c
/* Helper: split "Ina3" into base "Ina" and index 3 (1-based). Returns 1 if a
   trailing positive integer was found, else 0. */
static int SplitIndexedName( const char *name, char *base, int *index );
```

In `WriteDouble`/`WriteInt`/`WriteString`, if `SplitIndexedName` succeeds, append to the existing vector entry for `base` using `astMapPut1D`/`astMapPutElem...` at position `index-1`; otherwise store the scalar as before. Confirm the actual vector naming in `winmap.c` before implementing (do not assume — read the `Dump`).

- [ ] **Step 5: Implement array + nested-object emission in EmitJson**

Extend `EmitJson` so that for `astMapType == AST__OBJECTTYPE` it recurses (`astMapGet0A` → `EmitJson` indented), and for vector entries (`astMapLength( km, key ) > 1`) it emits a JSON array of the element type, using `FormatDouble` for doubles and rendering bad/undefined elements as `null`.

- [ ] **Step 6: Run the test to verify it passes**

Run: `cmake --build build-dev && ctest --test-dir build-dev -R testjsonchan --output-on-failure`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/jsonchan.c src/jsonchan.h ast_tester/testjsonchan.c
git commit -m "feat(jsonchan): emit nested objects and vector attributes"
```

---

## Task 4: Double precision, bad values, and string escaping

Lock down lossless double formatting, `null` for `AST__BAD`/non-finite, and JSON string escaping for special characters.

**Files:**
- Modify: `src/jsonchan.c`
- Modify: `ast_tester/testjsonchan.c`

**Interfaces:**
- Produces: `static const char *EscapeJson( const char *in, char **out, int *status )` (returns a heap buffer the caller frees, or writes into a grown buffer).

- [ ] **Step 1: Write the failing tests**

```c
static int test_precision_bad_escape( void ) {
   jsonbuf[ 0 ] = '\0';
   AstJsonChan *chan = astJsonChan( NULL, Sink, " " );
   /* ShiftMap with a bad shift and a high-precision shift. */
   double shifts[2] = { 0.1, AST__BAD };
   AstShiftMap *sm = astShiftMap( 2, shifts, " " );
   astSetC( sm, "Ident", "a\"b\\c" );   /* forces string escaping */
   astWrite( chan, sm );
   if ( !astOK ) return 1;
   if ( !strstr( jsonbuf, "null" ) )          { printf("no null for bad\n"); return 1; }
   if ( !strstr( jsonbuf, "0.1000000000000000" ) &&
        !strstr( jsonbuf, "0.1" ) )           { printf("no precise double\n"); return 1; }
   if ( !strstr( jsonbuf, "a\\\"b\\\\c" ) )   { printf("string not escaped\n"); return 1; }
   return 0;
}
```

(Use the actual ShiftMap attribute/vector names from `src/shiftmap.c`'s `Dump`; adjust the bad-value field name accordingly. The precise-double substring should match `%.*g` at `AST__DBL_DIG=18`.)

- [ ] **Step 2: Run to verify it fails**

Run: `cmake --build build-dev && ctest --test-dir build-dev -R testjsonchan --output-on-failure`
Expected: FAIL on the escaping assertion (raw quotes/backslashes not yet escaped).

- [ ] **Step 3: Implement JSON string escaping**

Add `EscapeJson` that replaces `"` → `\"`, `\` → `\\`, control chars (`\n`, `\t`, `\r`, and `\u00xx` for others) per RFC 8259, growing an output buffer with `astGrow`. Use it for every string value and for the `$type`/key names (keys are class/attribute identifiers so will not normally need escaping, but values can).

- [ ] **Step 4: Confirm double formatting**

Verify `FormatDouble` uses `AST__DBL_DIG` (18 under the C11 test build). Add an assertion in the test that re-reads the printed value with `strtod` and compares bit-exactly:

```c
   /* find the substring after "Shift1": and strtod it back */
```

(Implement the substring extraction inline; assert `parsed == 0.1` exactly.)

- [ ] **Step 5: Run to verify it passes**

Run: `cmake --build build-dev && ctest --test-dir build-dev -R testjsonchan --output-on-failure`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/jsonchan.c ast_tester/testjsonchan.c
git commit -m "feat(jsonchan): precise doubles, null bad values, string escaping"
```

---

## Task 5: Document-root version fields

Add `"$schema"`, `"$schemaVersion"`, and `"$minReadVersion"` to the outermost object only, sourced from `AST__SCHEMA_VERSION` / `AST__JSON_MINREADVER`.

**Files:**
- Modify: `src/jsonchan.c`
- Modify: `ast_tester/testjsonchan.c`

**Interfaces:**
- Consumes: `AST__SCHEMA_VERSION`, `AST__JSON_MINREADVER` from `jsonchan.h`.

- [ ] **Step 1: Write the failing test**

```c
static int test_root_version_fields( void ) {
   jsonbuf[ 0 ] = '\0';
   AstJsonChan *chan = astJsonChan( NULL, Sink, " " );
   AstZoomMap *zm = astZoomMap( 2, 5.0, " " );
   AstCmpMap *cm = astCmpMap( zm, astUnitMap( 2, " " ), 1, " " );
   astWrite( chan, cm );
   if ( !astOK ) return 1;
   if ( !strstr( jsonbuf, "\"$schemaVersion\": \"1.0.0\"" ) ) { printf("no version\n"); return 1; }
   if ( !strstr( jsonbuf, "\"$minReadVersion\": 1" ) )        { printf("no minread\n"); return 1; }
   /* Nested objects must NOT carry version fields: there should be exactly
      one occurrence of "$schemaVersion" in the whole document. */
   const char *p = strstr( jsonbuf, "$schemaVersion" );
   if ( p && strstr( p + 1, "$schemaVersion" ) ) { printf("version repeated in nested\n"); return 1; }
   return 0;
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cmake --build build-dev && ctest --test-dir build-dev -R testjsonchan --output-on-failure`
Expected: FAIL — version fields not emitted.

- [ ] **Step 3: Inject root fields in EmitJson**

`EmitJson` already only runs for the outermost object (called from `WriteEnd` at depth 0). Emit the three extra members immediately after the opening `{` and the `$type` line:

```c
   astPutNextText( chan, "{" );
   /* $type is the first stored key; render it, then root metadata. */
   ...
   char vline[ 128 ];
   (void) sprintf( vline, "  \"$schema\": \"https://www.starlink.ac.uk/ast/json/v%d/envelope.schema.json\",",
                   AST__JSON_MINREADVER );
   astPutNextText( chan, vline );
   (void) sprintf( vline, "  \"$schemaVersion\": \"%s\",", AST__SCHEMA_VERSION );
   astPutNextText( chan, vline );
   (void) sprintf( vline, "  \"$minReadVersion\": %d,", AST__JSON_MINREADVER );
   astPutNextText( chan, vline );
```

Ensure comma handling stays valid JSON (root metadata is always followed by at least the object's attributes; if an object could be empty apart from `$type`, guard the trailing comma). Keep the nested-object emission path (depth > 0) free of these fields.

- [ ] **Step 4: Run to verify it passes**

Run: `cmake --build build-dev && ctest --test-dir build-dev -R testjsonchan --output-on-failure`
Expected: PASS.

- [ ] **Step 5: Validate full document is valid JSON**

Run:
```bash
ctest --test-dir build-dev -R testjsonchan -V 2>/dev/null | sed -n '/^{/,/^}$/p' | python3 -m json.tool
```
Expected: parses cleanly, showing the root version fields and no version fields on nested objects.

- [ ] **Step 6: Commit**

```bash
git add src/jsonchan.c ast_tester/testjsonchan.c
git commit -m "feat(jsonchan): add document-root schema version fields"
```

---

## Self-Review notes (for the implementer)

- **`set` flag semantics:** native `Dump` calls `WriteXxx` for every attribute but passes `set=0` for unset ones; Phase 1 skips unset attributes (`if (!set) return;`) to mirror the native default output. The schema-harvesting mode (Phase 3) will instead record *all* calls regardless of `set`.
- **Vector naming:** Step 4 of Task 3 must be implemented against the *actual* `Dump` naming in `winmap.c`/`shiftmap.c`; read those before coding rather than assuming `BaseN`.
- **`astDump` recursion in `WriteObject`:** verify against `channel.c`'s own `WriteObject` how the native channel recurses into a child object, and mirror that call sequence so Begin/IsA/End nesting is consistent.
- **Comma/escaping correctness:** every `EmitJson` path must produce strictly valid JSON; the Task 5 `python3 -m json.tool` check is the gate.
- This plan covers spec sections: Architecture/placement, Shared object→KeyMap converter (write direction), JSON format (discriminator, callback-type mapping, example), numeric precision + bad values, document-root fields, and `AST__SCHEMA_VERSION`. **Deferred to Plan 2:** read path, cJSON, full round-trip + golden fixtures. **Deferred to Plan 3:** schema generator + CI drift/version guard.
