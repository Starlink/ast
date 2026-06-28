# JsonChan: JSON-native serialization for AST

Date: 2026-06-27
Status: Design (approved for planning)

## Summary

Add a new `JsonChan` class that serializes AST objects to and from JSON, as an alternative to the bespoke native text format produced by `Channel`/`astToString`.
The goal is to let non-AST tools consume AST objects using off-the-shelf JSON parsers, while remaining a lossless 1:1 representation of AST's object model (JSON syntax, AST semantics).
The work also publishes a per-class JSON Schema set so the format is self-documenting and validatable.

## Goals

- A `JsonChan` class supporting full lossless round-trip (write and read) of any AST object that the native `Channel` can serialize.
- JSON output readable by any standard JSON parser, with no AST-specific parsing required to read the bytes.
- A published, generated JSON Schema (draft 2020-12) per AST class, kept in sync with the code.
- Maximum reuse of existing AST infrastructure (`KeyMap`, the `Channel` callback contract, vendored-dependency pattern).
- Code factored so that a future *native* `YamlChan` encoding can reuse the same object/KeyMap core.

## Non-goals

- Interoperability with a foreign community schema (e.g. ASDF/gwcs). That is `YamlChan`'s existing job and is out of scope here.
- Changing or refactoring the native `Channel` text format or the per-class `Dump`/`Load` functions.
- Implementing the native `YamlChan` encoding itself (only enabling it structurally).
- JSON Schema validation performed *inside* AST at read time. The schemas are an external deliverable; runtime validation is a possible later phase.

## Background: why this fits AST cleanly

AST's `Channel` already separates *what* a class serializes from *how* a format is rendered.
Each class `Dump()` calls format-agnostic virtual callbacks (`WriteBegin`, `WriteInt`, `WriteInt64`, `WriteDouble`, `WriteString`, `WriteObject`, `WriteEnd`), and each `Load()` calls the matching `ReadInt`/`ReadDouble`/`ReadObject`/`ReadString` accessors.
A channel subclass decides what those callbacks do.

Two existing subclasses already build an intermediate node tree behind those callbacks rather than streaming:

- `XmlChan` builds an `AstXmlElement` DOM and serializes it.
- `YamlChan` uses an `AstKeyMap` as its intermediate representation (object to KeyMap to YAML, and YAML to KeyMap to object).

The native `Channel` deliberately streams (`WriteDouble` does `sprintf` straight to the sink) to avoid materializing large objects and to keep byte-exact output that many tests compare against.

Because the callback contract is already format-agnostic, adding JSON does not require touching the ~70 `Dump`/`Load` functions.
`JsonChan` is therefore a sibling of `XmlChan`/`YamlChan`, modelled on the same node-tree pattern, leaving the native path untouched.

## Architecture

### Placement and files

`JsonChan` is a direct subclass of `Channel`, a sibling of `XmlChan` and `YamlChan`.
New and modified files mirror the existing channel classes:

- `src/jsonchan.c`, `src/jsonchan.h` — class implementation and protected interface.
- `src/fjsonchan.c` — Fortran 77 interface.
- Public constructor `astJsonChan` (C) / `AST_JSONCHAN` (Fortran).
- `src/cjson/` — vendored cJSON (MIT licence; LGPL-compatible) for the read path only.
- `ast_tester/testjsonchan.c` plus `ast_add_test(testjsonchan)`.
- Schema generator program (see "Per-class schemas").
- Build wiring: `makeh`/`ast.h`, `AST_PAR`, new status codes in `ast_err.msg`, and both the CMake and autotools source lists.

### Shared object/KeyMap converter

The reusable, format-independent core is a converter between an AST object and an `AstKeyMap` node tree:

- Write: drive the object's `Dump()` and capture the `WriteXxx` callbacks into a nested `KeyMap`.
- Read: take a nested `KeyMap` and drive the class `Load()` functions to rebuild the object.

This converter is factored as a shared protected facility rather than living inside `jsonchan.c`.
`JsonChan` is then only responsible for KeyMap-to-JSON-text and JSON-text-to-KeyMap.
This is deliberate: `YamlChan` already owns YAML-to-KeyMap and KeyMap-to-YAML machinery, so a future *native* `YamlChan` encoding becomes mostly a matter of pointing that existing code at this shared converter.
Native YAML is not implemented here, but the factoring makes it nearly free later.

### Data flow

Write (`astWrite`):

1. Class `Dump()` fires the `WriteXxx` callbacks.
2. `JsonChan` handlers build an `AstKeyMap` node tree (nested KeyMaps for child objects, vector entries for array-valued items).
3. The KeyMap is walked and emitted as JSON text, using AST's own number formatting.

Read (`astRead`):

1. Source text is parsed by cJSON into a cJSON DOM.
2. The DOM is converted to an `AstKeyMap` tree.
3. Class `Load()` functions query the KeyMap via the `ReadXxx` accessors to rebuild the object.

The KeyMap uses nested KeyMap entries (not stored `AstObject` pointers) for child objects, so the intermediate representation is pure data, matching the `YamlChan` approach.

## JSON format specification

Each AST object maps to a JSON object containing:

- A discriminator field `"$type"` whose value is the AST class name (for example `"ZoomMap"`).
- Each serialized attribute, keyed by its AST attribute name, with values typed by the callback that produced them.
- Nested JSON objects for child AST objects, keyed by the slot name the class used in `WriteObject`.

The discriminator is `"$type"`: collision-proof (no AST attribute name begins with `$`), pairs with `"$schema"`, and maps directly to JSON Schema's `discriminator.propertyName`, so OpenAPI/pydantic/serde-style tooling can consume it.

### Document-level fields (root object only)

The top-level (root) object additionally carries, alongside `"$type"`:

- `"$schema"` — reference to the published envelope schema for the format version.
- `"$schemaVersion"` — the AST-JSON format version the document was written with (semantic version string).
- `"$minReadVersion"` — the minimum reader version required to interpret the document (integer). A reader whose supported version is lower refuses or warns.

These three fields appear once, at the document root. Nested AST objects carry only `"$type"` and their attributes, so there is no per-node version bloat.

### Callback type to JSON type mapping

The format is faithful to the four scalar callback types; there is no inference beyond them.

- `WriteInt` / `WriteInt64` → JSON number (integer).
- `WriteDouble` → JSON number (with AST precision; see below).
- `WriteString` → JSON string.
- `WriteObject` → nested JSON object.
- Vector-valued items → JSON array of the corresponding scalar type.

There are deliberately no JSON booleans: the callback layer has no boolean type, so a flag such as `Series` arrives as an integer and round-trips as the JSON number `1`/`0`, preserving native semantics exactly.

### Example

A `CmpMap` combining a `ZoomMap` and a `WinMap`:

```json
{
  "$type": "CmpMap",
  "$schema": "https://www.starlink.ac.uk/ast/json/v1/envelope.schema.json",
  "$schemaVersion": "1.0.0",
  "$minReadVersion": 1,
  "Nin": 2,
  "Nout": 2,
  "Series": 1,
  "MapA": { "$type": "ZoomMap", "Nin": 2, "Zoom": 5.0 },
  "MapB": { "$type": "WinMap", "Nin": 2, "Ina": [0.0, null], "Inb": [1.0, 1.0] }
}
```

(The `null` in `Ina` above is an `AST__BAD` element; see below.)

### Numeric precision and bad values

- Doubles must round-trip bit-exactly. The write path formats numbers itself using AST's `%.*g` at `AST__DBL_DIG`, the same logic the native `Channel` uses. cJSON's printer is not used for output because it truncates precision. On read, cJSON's `strtod`-based parse recovers the exact value.
- `AST__BAD` is `-(DBL_MAX)` (finite, defined in `src/pointset.h`). The native format writes it as the magic string `"<bad>"`, but for JSON cleanliness `JsonChan` represents it as JSON `null` instead, and reads `null` back as `AST__BAD`. This is a deliberate, format-specific choice: `null` keeps every numeric field a clean `number`-or-`null` union rather than a `number`-or-`string` union, which is far friendlier to schema tooling and codegen (it maps to, for example, `Optional[float]` in pydantic).
- To declare that union once rather than inlining it at every field, the schemas define a reusable type (for example `$defs/AstReal = {"type": ["number", "null"]}`) and `$ref` it wherever a double appears.
- Semantics: an absent key means the attribute is unset (default); a key present with `null` means an explicitly bad value. This matches AST, where `AST__BAD` denotes missing/bad data.
- The write path guards every double with `astISFINITE`/equality-to-`AST__BAD`, so `AST__BAD` and any stray NaN or infinity are all emitted as `null` rather than as invalid JSON. Stock cJSON parses `null` natively on read, so no parser patching is required.

## Versioning and compatibility

The AST-JSON format carries a single version for the format as a whole, declared once at the document root (`$schemaVersion` plus `$minReadVersion`).
There is deliberately no per-class (`$type`) version field.

This mirrors how AST's native format has actually evolved.
The native format has no version number anywhere; compatibility is purely additive (`channel.c` documents this).
When a class gains a new attribute, an older AST reading a newer dump ignores the unknown attribute, and a newer AST reading an older dump supplies the default for the missing one (for example `astReadDouble(channel, "zoom", 0.0)`).
No mapping class branches on a schema version, and there is no version-aware reader anywhere in the library.
In roughly 30 years there has never been a breaking native format version bump.

`JsonChan` follows the same model:

- Per-class evolution is additive. New attributes get defaults on read and are ignored when unknown. The inherited `Channel` `Strict` attribute governs whether unknown keys are silently ignored (default) or raise an error.
- `$minReadVersion` is cheap forward-compatibility insurance for a hypothetical breaking change, but in practice is expected to remain `1`, because a breaking JSON change would imply a breaking native change that AST has historically never made.
- The per-class schema files are versioned collectively under the format version (for example `schemas/v1/`), not individually per class.

## Vendored cJSON

cJSON is vendored under `src/cjson/`, following the existing pattern for `wcslib`, `pal`, `erfa`, and `cminpack`.
It is MIT-licensed and LGPL-compatible.
It is used only for the read path (JSON text to DOM); the write path emits text directly so AST controls number formatting.
cJSON's DOM maps almost 1:1 onto a `KeyMap`, which is why it is preferred over a pure tokenizer such as jsmn.

## Per-class schemas

A dedicated generator program produces one JSON Schema (draft 2020-12) file per AST class, committed to the repository (for example under `schemas/`), plus a generic envelope schema describing the recursive container shape (an object with a `$type`, named attributes, and nested objects).

Mechanism:

- A schema-harvesting mode of `JsonChan` intercepts every `WriteXxx` call regardless of the `set` flag and records the attribute `name`, its callback type, and the human-readable `comment` string that every `Dump()` already passes.
- The generator constructs one representative instance of each class from a small curated factory list and dumps it through the harvesting channel.
- Comments become JSON Schema `description` text, so the schemas are self-documenting.

A CI check regenerates the schemas and diffs them against the committed copies to catch drift when classes change.

Known limitation, stated explicitly rather than hidden: harvesting reliably captures scalar attributes, because `Dump()` calls the scalar `WriteXxx` functions unconditionally (even for unset attributes).
Conditional or purely structural items (vector lengths, child-object slots that appear only in some configurations) are described by the generic envelope schema and by hand-authored annotations where a class needs them.
The schemas are therefore best-effort derived from the code, with targeted manual augmentation.

## Attributes and error handling

- Attributes: reuse the inherited `Channel` attributes (`Indent`, `Full`, `Comment`, `SinkFile`, `SourceFile`, and `Strict`). No new encoding attribute is needed, because JSON carries only the native semantics. The inherited `Strict` attribute governs unknown-key tolerance on read (ignore by default, error when set), exactly as for the native `Channel`. A built-in schema-validation toggle is a possible later phase, not part of v1.
- Error handling: standard AST `astOK` plus `astError` with `AST__` status codes. New status codes are added to `ast_err.msg` for malformed JSON input and for an unrecognized or missing `$type`.

## Testing

`testjsonchan.c` follows the established checkdump test pattern and is compiled as C11, per the project rule that serialization tests need `AST__DBL_DIG = 18`.

Three layers of testing:

1. Semantic round-trip: build a battery of varied objects (ZoomMap, FrameSet with SkyFrame, CmpMap, PolyMap, LutMap, KeyMap), write to JSON, read back, and compare with `astEqual`/dump comparison.
2. Golden JSON fixtures: commit reference `.json` files generated once from a curated subset of the existing fixture corpus, and diff freshly emitted JSON against them to catch unintended format drift.
3. Schema validation: validate emitted JSON against the generated per-class schemas.

Fixture corpus available for round-trip and schema validation:

- `ast_tester/simplify_fixtures/`: 248 `.map` and 165 `.simp` files (the largest corpus, covering most Mapping classes and their simplifications).
- Top-level `ast_tester/`: 12 `.ast`, 3 `.map`, 3 `.simp`, and 132 `.head` files.
- `.head` files are turned into objects via the existing `wcsconverter` path.

The `simplify_fixtures` corpus in particular exercises a wide range of Mapping classes, making it valuable for both round-trip and schema-validation coverage.

Error-path tests cover malformed JSON, unknown `$type`, and bad-value handling.

## Phasing

1. JsonChan skeleton and write path: class boilerplate (vtab, `astInitJsonChan`, `astLoadJsonChan`), Fortran interface, build wiring, the shared object-to-KeyMap converter, and KeyMap-to-JSON emission with AST-precision numbers and `null` bad-value handling.
2. Read path: vendor cJSON, implement JSON-to-KeyMap and the KeyMap-to-object direction, and add the full round-trip and golden-fixture tests.
3. Schema generator: harvesting channel mode, curated factory list, committed per-class schema files, the generic envelope schema, and the CI drift check.

## Risks and mitigations

- Class boilerplate volume: adding a new class touches many build and interface files. Mitigation: copy the `XmlChan`/`YamlChan` structure directly.
- Schema drift and incompleteness: mitigated by generation plus CI diffing, with the conditional/structural limitation documented and hand-augmented.
- Memory: building a full KeyMap materializes the whole object, unlike native streaming. Acceptable because JSON is a sharing/interchange format, not the hot persistence path; noted so it is a conscious trade-off.
- Double precision: the dedicated AST formatting on write and cJSON `strtod` on read together guarantee bit-exact round-trip; covered by the C11 round-trip tests.
