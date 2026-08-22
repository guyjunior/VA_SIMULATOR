# external_tools

The converters ship with this repository, as portable builds — nothing is
installed into the operating system. The code resolves these paths from the
project root, so they work from any folder and on any machine.

```
external_tools/
├── thermo_raw_file_parser/   ThermoRawFileParser.exe + its dependencies
└── msconvert/                msconvert.exe (ProteoWizard portable)
```

There is no environment variable to point somewhere else, on purpose: an
absolute path in a config file is what hides a broken default until someone
changes machine.

## Licences

These binaries are redistributed here and carry terms that are not this
project's. **ThermoRawFileParser** is Apache-2.0 but embeds Thermo Fisher's
`RawFileReader`, under Thermo's own licence. **ProteoWizard** is Apache-2.0 at
its core, while its vendor readers (Waters, Agilent, Sciex, Bruker, Shimadzu,
Thermo) are proprietary libraries under per-vendor agreements. Check both before
publishing the repository.

Upstream, if you need to rebuild the folders:

- **ThermoRawFileParser** — <https://github.com/compomics/ThermoRawFileParser/releases>
  (the self-contained build)
- **ProteoWizard / msconvert** — <https://proteowizard.sourceforge.io/download.html>
  (the *portable* version, not the installer)

## Linux

TRFP has a self-contained Linux build (CoreCLR embedded) — swapping the folder is
enough.

**msconvert has no native Linux build**: the Waters, Agilent and Sciex vendor
readers are Windows DLLs. The options are the ProteoWizard Docker image
(`chambm/pwiz-skyline-i-agree-to-the-vendor-licenses`) or Wine. On Linux
msconvert stops being "an executable in a folder" and becomes a container or a
Wine prefix — an infrastructure decision, not a code one.
