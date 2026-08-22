# examples/

The data the one-command bootstrap uses:

```
raw/
├── Sample.raw        a CQCL injection, standing in for the sample under analysis
├── CQCL.raw          the positive control - the reference
├── CQCL_reinj.raw    the reinjection of that control - the fallback reference
└── CQN.raw           the negative control - the baseline
```

**Every one of them is a control injection. No athlete sample is published here.**
A real urine cannot be redistributed, and a fabricated one would not exercise the
pipeline honestly, so the sample slot is filled by a further CQCL injection. Its
internal header still identifies it as a CQCL - that is not an oversight, it is
the point.

It also makes for a better demonstration: a fortified control judged against a
fortified control lights up for the compounds it was spiked with, so the results
screen has something on it. What it is **not** is a validation of the method - a
control judged against itself says nothing about sensitivity or specificity.

They are Thermo `.raw` files, acquired on the laboratory's LC-HRMS under the TVII
screening method, and they are loaded by:

```bat
venv\Scripts\python -m tools.seed_example
```
