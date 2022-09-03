# 10.1021/acs.jpca.1c03061

[![DOI](https://zenodo.org/badge/448337446.svg)](https://zenodo.org/badge/latestdoi/448337446)

Supporting information (calculation outputs, structures) for https://doi.org/10.1021/acs.jpca.1c03061

The folder structure is method, basis set, then "standard" calculations.  Variations are in further subfolders.

| folder                                                           | description of contents                                                                       |
|------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| [`b3lyp/methane`](b3lyp/methane)                                 | Calculations for establishing spring and angle force coefficients for methane, B3LYP/6-31G(d) |
| [`m05-2x/def2-tzvp`](m05-2x/def2-tzvp)                           | All calculations based on M05-2X/def2-TZVP structures.                                        |
| [`wb97x-d/def2-tzvp`](wb97x-d/def2-tzvp)                         | All calculations based on wB97X-D/def2-TZVP structures, including SAPT.                       |
| [`wb97x-d/def2-tzvp/swap_metals`](wb97x-d/def2-tzvp/swap_metals) | Placing magnesium in the calcium-bound structure and vice-versa                               |

## Mapping to text in SI document

Originally, some of the calculations were combined geometry optimization and vibrational frequency calculations in the same input.  These inputs and outputs have all been split.  (As a result, some of the frequency inputs repeated in their respective outputs have `$molecule read`, where the extracted input has been modified to work standalone.)

- `methane.out` is really the entire contents of `b3lyp/methane`.
- `{apo,ca,mg}_opt_freq_{m,w}.out` are `{m05-2x,wb97x-d}/def2-tzvp/{apo,ca_ct,mg_ct}_{opt,freq}.out`, where the optimization and frequency parts have been split as described above.
- `{ca,mg}_{eda_w,sapt}.out` are `wb97x-d/def2-tzvp/{ca,mg}_ct_{eda,sapt_gas}.out`.
- `{cainmg,mginca}_{eda_w,sapt}.out` are `wb97x-d/def2-tzvp/swap_metals/{cainmg,mginca}_ct_{eda,sapt_gas}.out`.

## Analysis notes

The script to make EDA/SAPT analysis tables is also provided under [`wb97x-d/def2-tzvp/analysis.py`](wb97x-d/def2-tzvp/analysis.py).  It depends on [cclib](https://pypi.org/project/cclib/), [pydantic](https://pypi.org/project/pydantic/), and [pylatex](https://pypi.org/project/PyLaTeX/).
