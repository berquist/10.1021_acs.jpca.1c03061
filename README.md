# 10.1021/acs.jpca.1c03061
Supporting information (calculation outputs, structures) for https://doi.org/10.1021/acs.jpca.1c03061

The folder structure is method, basis set, then "standard" calculations.  Variations are in further subfolders.

| folder                                                           | description of contents                                                                       |
|------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| [`b3lyp/methane`](b3lyp/methane)                                 | Calculations for establishing spring and angle force coefficients for methane, B3LYP/6-31G(d) |
| [`m05-2x/def2-tzvp`](m05-2x/def2-tzvp)                           | All calculations based on M05-2X/def2-TZVP structures.                                        |
| [`wb97x-d/def2-tzvp`](wb97x-d/def2-tzvp)                         | All calculations based on wB97X-D/def2-TZVP structures, including SAPT.                       |
| [`wb97x-d/def2-tzvp/swap_metals`](wb97x-d/def2-tzvp/swap_metals) | Placing magnesium in the calcium-bound structure and vice-versa                               |

## Mapping to text in SI document

Originally, some of the calculations were combined geometry optimization and vibrational frequency calculations in the same input.  These inputs and outputs have all been split.

- `methane.out` is really the entire contents of `b3lyp/methane`.

## Analysis notes

The script to make EDA/SAPT analysis tables is also provided under `wb97x-d/def2-tzvp/analysis.py`.  It depends on [cclib](https://pypi.org/project/cclib/), [pydantic](https://pypi.org/project/pydantic/), and [pylatex](https://pypi.org/project/PyLaTeX/).
