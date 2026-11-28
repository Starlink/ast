# News fragments

This directory holds one small text file for each change that users of AST should hear about.
At release time `towncrier` assembles them into a new section at the top of `ast.news` and deletes them; see `RELEASING.md`.
Do not edit `ast.news` directly for a change you are making: add a fragment here instead.

## Naming

Name the file `<PR number>.<type>`, for example `104.hardening`.
A change with no pull request uses `+<slug>.<type>`, for example `+erfa-2.0.1.behaviour`.
The name is only used to keep fragments distinct; it is not written to the notes.

The type decides which section of the release notes the fragment appears in:

| Type | Section | Use it for |
| --- | --- | --- |
| `breaking` | IMPORTANT - changes that may break existing code | Code that worked before and now reports an error or gives different results |
| `behaviour` | Changes in behaviour | Documented or long-standing behaviour that changes, including numerical results |
| `feature` | New features | New classes, functions, attributes, encodings or build options |
| `hardening` | Hardening against corrupt input | Crashes, memory errors or unbounded work on corrupt or malicious input |
| `performance` | Performance | Changes that make something materially faster or smaller |
| `serialisation` | Serialisation fixes | Reading and writing dumps, FITS headers, MOCs, STC-S, XML or YAML |
| `bugfix` | Bug fixes | Everything else that gave a wrong answer or an error |
| `doc` | Documentation | Changes to SUN/210, SUN/211 or the prologues that users should know about |
| `build` | Build and test changes | Build systems, packaging, test suite |

A `breaking` or `behaviour` fragment means the next release needs at least a minor version bump.

## Writing the text

Write one paragraph of one to three sentences describing the effect a user sees, not how the code changed.
Say what was wrong and what happens now, name the public functions, classes and attributes involved, and give the trigger if it is not obvious.
Use British spelling, as the rest of AST does.
Do not wrap the text by hand; towncrier wraps it to 79 columns.
Do not start with a dash or bullet; towncrier adds it.
Avoid Markdown and LaTeX markup, because the same text is copied into `sun_master.tex` by hand at release time.

A good fragment:

```
astRemoveTables with a blank key now removes every table in the FitsChan.
It used to remove every other table and then report an invalid index.
```

A fragment that is too detailed:

```
With a blank key, astRemoveTables removed the table at each index in turn
while the KeyMap shrank, so it removed every other table and then reported
AST__MPIND. The loop now iterates from the end.
```

## Checking

Preview the section that the current fragments would produce:

```shell
uvx towncrier build --draft --version "$(cat version.number)"
```

Pull requests are checked in CI for the presence of a fragment.
A pull request that changes nothing a user would notice, such as a test-only or housekeeping change, can be labelled `no news` to skip the check.
