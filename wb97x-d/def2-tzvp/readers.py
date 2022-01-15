import re
from pathlib import Path
from typing import Union

KJ_TO_KCAL = 4.184
RE_EDA_FRAGMENT_ENERGY = re.compile(r"\d{1,}\s*(-\d{1,}\.\d{10})")
DASHES = "--------------------"
re_number_basic_eda_quantities = re.compile(r"=\s(-?\d+\.\d+)")
re_number_str = r"(-?[0-9]*\.[0-9]*)"
re_number = re.compile(re_number_str)

SAPT_HEADERS_MONOMER = [
    "//         Monomer Basis SAPT        //",
]

SAPT_HEADERS_DIMER = [
    "//               SAPT0               //",
    "//          Dimer Basis SAPT         //",
]

SAPT_HEADERS = SAPT_HEADERS_MONOMER + SAPT_HEADERS_DIMER

SAPT_BASES = ("monomer", "dimer")


def read_qchem_eda_v2(filename: Union[str, Path]):
    """ALMO-EDA interaction energy components are in units of kcal/mol, but
    fragment energies are in units of hartrees.
    """
    almo_data = dict()
    fh = open(str(filename))
    line = next(fh)
    while "Results of EDA2" not in line:
        line = next(fh)
    line = next(fh)
    assert line.strip() == "================================"
    line = next(fh)
    assert line.strip() == "Basic EDA Quantities"
    line = next(fh)
    assert line.strip() == DASHES
    line = next(fh)
    assert line.strip() == "Fragment Energies (Ha):"
    line = next(fh)

    frgm_energies = list()
    while line.strip() != DASHES:
        frgm_energies.append(float(RE_EDA_FRAGMENT_ENERGY.search(line).groups()[0]))
        line = next(fh)
    line = next(fh)

    while line.strip() != DASHES:
        if "_prp " in line:
            label = "prp"
        elif "_sol " in line:
            label = "sol"
        elif "_frz(solv)" in line:
            label = "frz_solv"
        elif "_frz " in line:
            label = "frz"
        elif "_pol " in line:
            label = "pol"
        elif "_vct " in line:
            label = "vct"
        elif "_int " in line:
            label = "int"
        else:
            raise RuntimeError(f"Don't know how to handle {line}")
        almo_data[label] = float(re_number_basic_eda_quantities.search(line).groups()[0])
        line = next(fh)

    # TODO this assumes that the preparation term is zero: that we are not doing bonded EDA.
    #
    # keys = almo_data.keys()
    # for key in keys:
    #     value = almo_data[key]
    #     if value < 1.0e-12:
    #         del almo_data[key]
    if "prp" in almo_data:
        del almo_data["prp"]

    gas_headers = {
        "Orthogonal Frozen Decomposition:",
        "Classical Frozen Decomposition:",
    }
    solvated_headers = {
        "Classical Frozen Decomposition (unscreened):",
        "Classical Decomposition Terms with Solvent Contributions:",
        "Orthogonal Frozen Decomposition (unscreened):",
        "Orthogonal Decomposition Terms with Solvent Contribution:",
    }
    headers = gas_headers | solvated_headers

    while line.strip() != "Perturbative CT Analysis:":
        if line.strip() in headers:
            line = next(fh)
            assert set(line.strip()) == {"-"}
            line = next(fh)
            while set(line.strip()) != {"-"}:
                tokens = line.split()
                key = tokens[0][2:]
                value = float(
                    re_number_basic_eda_quantities.search(line).groups()[0]
                )
                if key in almo_data:
                    print(f"{key} already in almo_data")
                    print(f"  old: {almo_data[key]}")
                    print(f"  new: {value}")
                    diff = abs(almo_data[key] - value)
                    if diff > 1.0e-12:
                        raise RuntimeError(diff)
                almo_data[key] = value
                line = next(fh)
        line = next(fh)

    assert line.strip() == "Perturbative CT Analysis:"
    line = next(fh)
    assert line.strip() == DASHES
    line = next(fh)
    almo_data["pct"] = float(line.split()[3])
    line = next(fh)
    almo_data["HO"] = float(line.split()[3])
    line = next(fh)
    assert line.strip() == "---------------"
    # TODO PCT Energy lowering
    # TODO PCT Charge displacement
    fh.close()
    for k in almo_data:
        almo_data[k] /= KJ_TO_KCAL
    return almo_data, frgm_energies


def read_psi4_sapt_section(fi, pos: int = 1, calculation_thresh: float = 1.0e-7):
    """All returned values have units of kcal/mol."""
    sapt_single_basis_data = dict()

    line = ""
    while "SAPT Results" not in line:
        line = next(fi)
    line = next(fi)
    assert "------" in line
    line = next(fi)
    assert "Electrostatics" in line
    val_electrostatics = float(re_number.findall(line)[pos])
    sapt_single_basis_data["el"] = val_electrostatics
    line = next(fi)
    assert "Elst10,r" in line
    line = next(fi)
    assert line.strip() == ""
    line = next(fi)
    assert "Exchange" in line
    val_exchange = float(re_number.findall(line)[pos])
    sapt_single_basis_data["exch"] = val_exchange
    line = next(fi)
    assert "Exch10" in line
    line = next(fi)
    assert "Exch10(S^2)" in line
    line = next(fi)
    assert line.strip() == ""
    line = next(fi)
    assert "Induction" in line
    line = next(fi)
    assert "Ind20,r" in line
    val_induction = float(re_number.findall(line)[pos])
    sapt_single_basis_data["ind"] = val_induction
    line = next(fi)
    assert "Exch-Ind20,r" in line
    val_exchange_induction = float(re_number.findall(line)[pos])
    sapt_single_basis_data["exch-ind"] = val_exchange_induction
    line = next(fi)
    assert "delta HF,r (2)" in line
    val_induction_delta_hf = float(re_number.findall(line)[pos])
    sapt_single_basis_data["ind_HO"] = val_induction_delta_hf
    line = next(fi)
    assert line.strip() == ""
    line = next(fi)
    assert "Dispersion" in line
    line = next(fi)
    assert "Disp20" in line
    val_dispersion = float(re_number.findall(line)[pos])
    sapt_single_basis_data["disp"] = val_dispersion
    line = next(fi)
    assert "Exch-Disp20" in line
    val_exchange_dispersion = float(re_number.findall(line)[pos])
    sapt_single_basis_data["exch-disp"] = val_exchange_dispersion

    while "Total SAPT0" not in line:
        line = next(fi)
    sapt0_total_calculated = float(re_number.findall(line)[pos])

    sapt0_total = (
        val_electrostatics
        + val_exchange
        + val_induction
        + val_exchange_induction
        + val_induction_delta_hf
        + val_dispersion
        + val_exchange_dispersion
    )

    assert abs(sapt0_total - sapt0_total_calculated) < calculation_thresh

    sapt_single_basis_data["total"] = sapt0_total

    return sapt_single_basis_data


def make_file_iterator(filename):
    """Return an iterator over the contents of the given file name."""
    # pylint: disable=C0103
    with open(filename) as f:
        contents = f.read()
    return iter(contents.splitlines())


def read_psi4_sapt0(filename: str, pos: int = 1):
    """SAPT interaction energy components are in units of kcal/mol, but fragment
    energies are in units of hartrees.
    """
    fi = make_file_iterator(filename)
    sapt_data = dict()
    frgm_energies = list()

    # Collect both dimer-centered and monomer-centered SAPT basis
    # data.
    for line in fi:
        if "Final Energy:" in line:
            frgm_energies.append(float(line.split()[3]))
        # Dimer results always come before monomer results.
        if any(sapt_header in line for sapt_header in SAPT_HEADERS_DIMER):
            sapt_data["dimer"] = read_psi4_sapt_section(fi, pos)
        if any(sapt_header in line for sapt_header in SAPT_HEADERS_MONOMER):
            sapt_data["monomer"] = read_psi4_sapt_section(fi, pos)
            break

    # Finally, check to see if a charge transfer (CT) calculation has
    # been performed.
    for line in fi:
        if "SAPT Charge Transfer Analysis" in line:
            line = next(fi)
            assert list(set(line.strip())) == ["-"]
            line = next(fi)
            # Asserts here comparing to induction values that were
            # parsed earlier?
            assert "SAPT Induction (Dimer Basis)" in line
            line = next(fi)
            assert "SAPT Induction (Monomer Basis)" in line
            line = next(fi)
            assert "SAPT Charge Transfer" in line
            ct = float(re_number.findall(line)[pos])
            sapt_data["ct"] = ct

    return sapt_data, frgm_energies
