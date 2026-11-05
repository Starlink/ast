# Oracles

Captures of what this library does, for a port to be held to without reading the
C.  Four files: three transform oracles, described first, and
`keymap_mapsz.txt`, described at the end.

## Transform oracles

`simplify_fixtures.oracle`, `headers.oracle` and `framesets.oracle` record what
the library's transforms produce for every fixture the generator can drive, so a
port has something to be checked against that does not depend on reading the C.

Regenerate all three together, to temporary files first so the diff can be read
before it is installed:

    cmake --build build --target gen_transform_oracle
    build/ast_tester/gen_transform_oracle ast_tester/fixtures \
        /tmp/simplify.oracle /tmp/headers.oracle /tmp/framesets.oracle

Section headers name paths relative to `ast_tester/fixtures`, so a section
identifies its fixture unambiguously even where two directories hold a file of
the same name — `car1.head` through `car5.head` exist under both
`wcsconv/inputs/` and `programs/joye_car_headers/`.

### Fixtures with no section

A fixture absent from an oracle is silently untested, so every absence is
accounted for here.  An absence for any other reason is a bug.

The generator writes a section whenever a fixture's object has a forward *or* an
inverse transformation.  These nine have neither, so there is nothing to record:

    simplify/cap_matrixmap_04.map
    simplify/cap_specmap_02.map
    simplify/cap_specmap_02.simp
    simplify/neg_ratemap_different_inner.map
    simplify/neg_ratemap_diff_indices.map
    simplify/neg_ratemap_diff_inner.map
    simplify/ratemap_inverse_cancel.map
    simplify/selectormap_inverse_cancel.map
    simplify/selectormap_inverse_cancel_tail.map

Nothing else is skipped.  The generator used to skip any fixture containing an
IntraMap, because such a fixture names a private transformation function that has
to be registered before its dump can be read; it now makes the same four
registrations `simplify.c` does, so those ten fixtures are covered.

### Sampling domain

A fixture's inputs are sampled over the domain it was built for, which for a
FITS header is its NAXIS pixel grid.  A dump under `wcsconv/framesets/` is the
object read from one wcsconv input, and its name is that input's filename with
each `.` replaced by `-`, so the generator recovers the input and takes its grid.

That matters rather than being tidy.  A SIP distortion is a polynomial fitted
over the image and its inverse is iterative: sampled a few times beyond the
image the iteration need not converge, and a default 1..1000 range over a
256-pixel image produced round-trip errors above 100 pixels that said nothing
about the library.  The same mappings round-trip cleanly over their own grid.

### Checks that are disabled, and why

`transform_oracle_overrides.txt` carries them, each with its reason.  Two entries
name the same object twice, once as the header and once as the dump under
`wcsconv/framesets/`, because an override is keyed on the relpath: `tnx-cheb` has
no inverse at all, its TNX Chebyshev distortion building a ChebyMap that cannot
represent its own Jacobian.

A simplify fixture's `.map` and `.simp` are compared only where their recorded outputs agree, so a pair that agrees nowhere would never be compared.
Such a pair fails the check unless the overrides file accepts it with `equiv off`, and an accepted pair that is compared after all fails too, so the list cannot go stale.
The pairs accepted are those whose common domain the -1000 to +1000 sampling misses, and one whose inverse is a fitted approximation; the overrides file says how each was checked instead.

### Values

The recorded values are this library's, on the host that generated them.  Ports
compare against them with a tolerance -- last-bit differences between platforms'
maths libraries are expected, and an angular axis may legitimately come out as 0
where another platform gives 2*pi, or -pi where another gives +pi.

The checker accepts `|got - ref| <= atol + rtol*|ref|`, with the `rtol` and `atol` in each file's `# tol:` line.
`atol` is for values that cancel toward zero, so it is measured against the size of the values on that axis rather than against 1.
On an axis whose largest recorded |value| is below 1, the checker uses `atol` times that value.
Otherwise, including on an axis whose recorded values are all zero, it uses `atol` as written.
Without that, an axis in metres, such as a wavelength of 5e-7, would be allowed an error of 20%, and an axis in joules any error at all.

### Grism values checked against FITS-WCS Paper III

Most sections record what this library does and nothing more.
The grism headers below were also checked against an independent source, so their sections can be relied on as correct values rather than only as a record of current behavior.

`wcsconv/inputs/awav-gra-coude.head`, `awav-gra-hydra.head` and `awav-gra-mars.head` are the KPNO examples printed with Figs. 3-5 of FITS-WCS Paper III (Greisen et al. 2006, A&A 446, 747).
`wave-gra.head` and `awav-gri.head` use the Fig. 5 grism parameters, plus a non-zero epsilon and theta, on an axis in vacuum wavelength for a grism in air and the reverse.
`grism.head` is a `WAVE-GRI` axis for a 300 lines/mm grating.

Their sections, and those of the matching dumps under `wcsconv/framesets/`, were compared with a direct implementation of Paper III eqs. (64)-(83) written independently of this library.
That implementation finds an air wavelength by inverting eq. (64) exactly, as this library does, rather than with the approximation of eq. (67).
The two agree to within 2e-11 pixels for `grism.head` and the three `AWAV-GRA` headers, and 2e-8 pixels for the other two, where the difference is the numerical derivative of the air-vacuum conversion.
The same four headers also agree with WCSLIB 8.6 to within 2e-11 pixels, and the departures of the three `AWAV-GRA` headers from a linear axis match the curves plotted in Paper III.

WCSLIB differs from the other two headers by about 0.003 pixels.
It converts between air and vacuum wavelengths with the Edlén (1953) relation, the IAU standard, whereas Paper III eq. (65), which this library uses, is the IUGG (1999) relation; WCSLIB's `spx.h` explains the choice.

## `keymap_mapsz.txt`

Eighteen rows of `<label> sizeguess=<n> nkey=<n> MapSz=<n>`: the `MapSz` card a
freshly built KeyMap carries for a range of key populations.

`MapSz` is the size of the hash table and the only observable of the table's
growth policy, because `HashFun` and `DoubleTableSize` are both `static` in
`keymap.c`.  A port cannot be checked against those directly, so the populations
here sit either side of each growth threshold instead, which pins the policy from
outside.  The colliding-key rows are the ones that actually drive
`DoubleTableSize` -- forty keys landing in one bucket give `MapSz=128` where
forty spread keys do not -- and `gen_keymap_mapsz.c` chooses them with a copy of
`keymap.c`'s djb2 arithmetic, used only to pick keys and never to check a result.

Regenerate with:

    cmake --build build --target gen_keymap_mapsz
    ./build/ast_tester/gen_keymap_mapsz > ast_tester/fixtures/oracle/keymap_mapsz.txt

Unlike the transform oracles this one is checked, by the `keymap_mapsz` ctest,
which re-runs the generator and requires byte-identical output.  A byte
comparison is right here and a tolerance is not: every value is an integer
produced by integer arithmetic over ASCII keys, with no floating point anywhere,
so there is no platform variation to absorb and a difference is a real change in
the growth policy, in `HashFun`, or in the dump.

The rows are consumed by the Rust port
(`crates/ast-core/tests/test_keymap_serial.rs`, `map_sz_matches_c_oracle`), which
rebuilds each population and compares its own `MapSz`.  Nothing in this suite
reads them -- the ctest regenerates rather than reads -- so the file is here
because the fixture tree is shared with that port and a capture of this library
has to be taken from this library.
