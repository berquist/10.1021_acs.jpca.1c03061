from copy import deepcopy
from pathlib import Path
from pprint import pprint
from typing import Any, Tuple, Union

from cclib.io import ccread
from cclib.parser.utils import convertor

from pydantic import BaseModel, Field

from pylatex.utils import NoEscape

from readers import read_psi4_sapt0, read_qchem_eda_v2


SINGLEROWFMT = r"{} & {} & {} \\".format
DOUBLEROWFMT = r"{} & {} & {} & {} & {} \\".format
TRIPLEROWFMT = r"{} & {} & {} & {} & {} & {} & {} \\".format
QUADROWFMT = r"{} & {} & {} & {} & {} & {} & {} & {} & {} \\".format




def sapt0_make_approx_frz(sapt_data, sapt_basis: str) -> float:
    return sapt_data[sapt_basis]["el"] + sapt_data[sapt_basis]["exch"]


def sapt0_make_approx_pol(sapt_data, sapt_basis: str) -> float:
    return (
        sapt_data[sapt_basis]["ind"]
        + sapt_data[sapt_basis]["exch-ind"]
        + sapt_data[sapt_basis]["ind_HO"]
    )


def sapt0_make_total_disp(sapt_data, sapt_basis: str) -> float:
    return sapt_data[sapt_basis]["disp"] + sapt_data[sapt_basis]["exch-disp"]


def sapt0_make_almo_like_quantities(sapt_data):
    ret = deepcopy(sapt_data)
    for basis in ret:
        if basis != "ct":
            ret[basis]['"frz"'] = sapt0_make_approx_frz(ret, basis)
            ret[basis]['"pol"'] = sapt0_make_approx_pol(ret, basis)
            ret[basis]['"disp"'] = sapt0_make_total_disp(ret, basis)
    return ret


def almo_calculate_percentages(almo_data):
    return {k: 100 * v / almo_data["int"] for k, v in almo_data.items()}


def sapt0_calculate_percentages(sapt_data):
    ret = dict()
    for k, v in sapt_data.items():
        if not isinstance(v, float):
            ret[k] = {k2: 100 * v2 / v["total"] for k2, v2 in v.items()}
    return ret


_GEOM_KEY = "geom"

MAP_LATEX_EDA = {
    "sol": r"\(\Delta E_{\textrm{solv}}\)",
    "elec": r"\(\Delta E_{\textrm{elec}}\)",
    "pauli": r"\(\Delta E_{\textrm{Pauli}}\)",
    "disp": r"\(\Delta E_{\textrm{disp}}\)",
    "pol": r"\(\Delta E_{\textrm{pol}}\)",
    # "pct": r"\(\Delta E_{\textrm{CT}}^{\textrm{pert}}\)",
    # "HO": r"\(\Delta E_{\textrm{CT}}^{\textrm{HO}}\)",
    "vct": r"\(\Delta E_{\textrm{CT}}\)",
    # "int": r"\(\Delta E_{\textrm{int}}\)",
}

MAP_LATEX_SAPT = {
    "el": r"\(\Delta E_{\textrm{elec}}\)",
    "exch": r"\(\Delta E_{\textrm{Pauli}}\)",
    '"disp"': r"\(\Delta E_{\textrm{disp}}\)",
    '"pol"': r"\(\Delta E_{\textrm{pol}}\)",
    # "total": r"\(\Delta E_{\textrm{int}}\)",
    # "ct": r"\(\Delta E_{\textrm{CT}}\)",
}


class SAPTCTEnergies(BaseModel):
    supersystem: float = Field(...)
    a_dimer: float = Field(...)
    b_dimer: float = Field(...)
    a_monomer: float = Field(...)
    b_monomer: float = Field(...)


def read_psi4_sapt0_obj(filename: Union[str, Path]) -> Tuple[Any, SAPTCTEnergies]:
    sapt_data, frgm_energies = read_psi4_sapt0(filename)
    for basis in {"monomer", "dimer"}:
        sapt_data[basis]["ct"] = sapt_data["ct"]
    return sapt_data, SAPTCTEnergies(
        supersystem=frgm_energies[0],
        a_dimer=frgm_energies[1],
        b_dimer=frgm_energies[2],
        a_monomer=frgm_energies[3],
        b_monomer=frgm_energies[4],
    )


def recalculate_eda_int(almo_eda_data):
    """Recalculate the total ALMO-EDA interaction energy based on existing
    components.

    This is probably necessary when the geometric distortion term has been added.
    """
    ret = deepcopy(almo_eda_data)
    ret["int"] = (
        ret.get(_GEOM_KEY, 0.0)
        + ret.get("sol", 0.0)
        # + ret["elec"]
        # + ret["pauli"]
        + ret["cls_elec"]
        + ret["mod_pauli"]
        + ret["disp"]
        + ret["pol"]
        + ret["vct"]
    )
    # The "original" interaction energy.  This will probably be slightly
    # different from the printed value (5th decimal place) due to at least one
    # unit conversion (kJ/mol to kcal/mol) and limited precision in the
    # original value parsed from the output file.
    #
    # ret["int"] = (
    #     ret.get("sol", 0.0)
    #     + ret["elec"]
    #     + ret["pauli"]
    #     + ret["disp"]
    #     + ret["pol"]
    #     + ret["vct"]
    # )
    return ret


def recalculate_sapt0_int(sapt_data):
    """Recalculate the total SAPT0 interaction energy based on existing
    components.

    This is probably necessary when the geometric distortion term has been
    added, even though it technically isn't a part of the true SAPT
    interaction.
    """
    ret = deepcopy(sapt_data)
    for basis in {"monomer", "dimer"}:
        # This is the "original" interaction energy plus the geometric
        # distortion term.
        ret[basis]["total"] = (
            ret[basis].get(_GEOM_KEY, 0.0)
            + ret[basis]["el"]
            + ret[basis]["exch"]
            + ret[basis]["ind"]
            + ret[basis]["exch-ind"]
            + ret[basis]["ind_HO"]
            + ret[basis]["disp"]
            + ret[basis]["exch-disp"]
        )
    return ret


def ratio_sapt0_data(left, right):
    """Take the ratio of the 'left' SAPT0 data to the 'right' and multiply by
    100.
    """
    ret = deepcopy(left)
    assert left.keys() == right.keys()
    for k, v in left.items():
        if not isinstance(v, float):
            for k2 in v:
                ret[k][k2] = 100 * left[k][k2] / right[k][k2]
    return ret


def ratio_eda_data(left, right):
    """Take the ratio of the 'left' ALMO-EDA data to the 'right' and multiply by
    100.
    """
    ret = deepcopy(left)
    assert left.keys() == right.keys()
    for k in left:
        ret[k] = 100 * left[k] / right[k]
    return ret


if __name__ == "__main__":

    mg_ct_sapt, mg_ct_sapt_frgm = read_psi4_sapt0_obj("mg_ct_sapt_gas.out")
    mg_ct_sapt = sapt0_make_almo_like_quantities(mg_ct_sapt)
    ca_ct_sapt, ca_ct_sapt_frgm = read_psi4_sapt0_obj("ca_ct_sapt_gas.out")
    ca_ct_sapt = sapt0_make_almo_like_quantities(ca_ct_sapt)
    mg_ct_eda_data, mg_ct_eda_frgm_energies = read_qchem_eda_v2("mg_ct_eda.out")
    ca_ct_eda_data, ca_ct_eda_frgm_energies = read_qchem_eda_v2("ca_ct_eda.out")
    pprint(mg_ct_eda_data)
    pprint(ca_ct_eda_data)

    # Geometric distortion term: the energy raising going from the isolated
    # fragment geometry (apo form with no metal bound) to the bound complex
    # geometry
    #
    # First put in hartrees, then kcal/mol
    data_edta_apo_freq = ccread("apo_freq.out")
    energy_eda_edta_apo_au = convertor(
        data_edta_apo_freq.scfenergies[0], "eV", "hartree"
    )
    energy_eda_edta_mg_au = mg_ct_eda_frgm_energies[1]
    energy_eda_edta_ca_au = ca_ct_eda_frgm_energies[1]

    energy_eda_edta_apo_kcal_mol = convertor(
        energy_eda_edta_apo_au, "hartree", "kcal/mol"
    )
    energy_eda_edta_mg_kcal_mol = convertor(
        energy_eda_edta_mg_au, "hartree", "kcal/mol"
    )
    energy_eda_edta_ca_kcal_mol = convertor(
        energy_eda_edta_ca_au, "hartree", "kcal/mol"
    )

    mg_ct_eda_data[_GEOM_KEY] = (
        energy_eda_edta_mg_kcal_mol - energy_eda_edta_apo_kcal_mol
    )
    ca_ct_eda_data[_GEOM_KEY] = (
        energy_eda_edta_ca_kcal_mol - energy_eda_edta_apo_kcal_mol
    )

    # Take into account the geometric distortion term.
    mg_ct_eda_data = recalculate_eda_int(mg_ct_eda_data)
    ca_ct_eda_data = recalculate_eda_int(ca_ct_eda_data)

    # need to do the same for SAPT, which is based on Hartree-Fock
    # wavefunctions rather than wB97X-D.
    data_edta_apo_sp_hf = ccread("apo_opt_sp_hf.out")
    energy_hf_edta_apo_au = convertor(
        data_edta_apo_sp_hf.scfenergies[0], "eV", "hartree"
    )

    energy_hf_edta_apo_kcal_mol = convertor(
        energy_hf_edta_apo_au, "hartree", "kcal/mol"
    )
    energy_hf_edta_mg_monomer_kcal_mol = convertor(
        mg_ct_sapt_frgm.b_monomer, "hartree", "kcal/mol"
    )
    energy_hf_edta_ca_monomer_kcal_mol = convertor(
        ca_ct_sapt_frgm.b_monomer, "hartree", "kcal/mol"
    )
    energy_hf_edta_mg_dimer_kcal_mol = convertor(
        mg_ct_sapt_frgm.b_dimer, "hartree", "kcal/mol"
    )
    energy_hf_edta_ca_dimer_kcal_mol = convertor(
        ca_ct_sapt_frgm.b_dimer, "hartree", "kcal/mol"
    )

    mg_ct_sapt["dimer"][_GEOM_KEY] = (
        energy_hf_edta_mg_dimer_kcal_mol - energy_hf_edta_apo_kcal_mol
    )
    ca_ct_sapt["dimer"][_GEOM_KEY] = (
        energy_hf_edta_ca_dimer_kcal_mol - energy_hf_edta_apo_kcal_mol
    )
    mg_ct_sapt["monomer"][_GEOM_KEY] = (
        energy_hf_edta_mg_monomer_kcal_mol - energy_hf_edta_apo_kcal_mol
    )
    ca_ct_sapt["monomer"][_GEOM_KEY] = (
        energy_hf_edta_ca_monomer_kcal_mol - energy_hf_edta_apo_kcal_mol
    )

    # Take into account the geometric distortion term.
    mg_ct_sapt = recalculate_sapt0_int(mg_ct_sapt)
    ca_ct_sapt = recalculate_sapt0_int(ca_ct_sapt)

    mg_ca_ratio_eda = ratio_eda_data(mg_ct_eda_data, ca_ct_eda_data)
    mg_ca_ratio_sapt = ratio_sapt0_data(mg_ct_sapt, ca_ct_sapt)

    mg_ct_eda_almo_pct = almo_calculate_percentages(mg_ct_eda_data)
    ca_ct_eda_almo_pct = almo_calculate_percentages(ca_ct_eda_data)
    mg_ct_sapt_pct = sapt0_calculate_percentages(mg_ct_sapt)
    ca_ct_sapt_pct = sapt0_calculate_percentages(ca_ct_sapt)

    ##########

    swap_metals_base = Path(".").resolve() / "swap_metals"
    cainmg_ct_eda_filename = swap_metals_base / "cainmg_ct_eda.out"
    mginca_ct_eda_filename = swap_metals_base / "mginca_ct_eda.out"
    cainmg_ct_sapt_gas_filename = swap_metals_base / "cainmg_ct_sapt_gas.out"
    mginca_ct_sapt_gas_filename = swap_metals_base / "mginca_ct_sapt_gas.out"

    mginca_ct_eda_data, mginca_ct_eda_frgm_energies = read_qchem_eda_v2(
        mginca_ct_eda_filename
    )
    cainmg_ct_eda_data, cainmg_ct_eda_frgm_energies = read_qchem_eda_v2(
        cainmg_ct_eda_filename
    )
    mginca_ct_sapt, mginca_ct_sapt_frgm = read_psi4_sapt0_obj(
        mginca_ct_sapt_gas_filename
    )
    mginca_ct_sapt = sapt0_make_almo_like_quantities(mginca_ct_sapt)
    cainmg_ct_sapt, cainmg_ct_sapt_frgm = read_psi4_sapt0_obj(
        cainmg_ct_sapt_gas_filename
    )
    cainmg_ct_sapt = sapt0_make_almo_like_quantities(cainmg_ct_sapt)

    mginca_ct_eda_data[_GEOM_KEY] = ca_ct_eda_data[_GEOM_KEY]
    cainmg_ct_eda_data[_GEOM_KEY] = mg_ct_eda_data[_GEOM_KEY]
    for basis in ("monomer", "dimer"):
        mginca_ct_sapt[basis][_GEOM_KEY] = ca_ct_sapt[basis][_GEOM_KEY]
        cainmg_ct_sapt[basis][_GEOM_KEY] = mg_ct_sapt[basis][_GEOM_KEY]

    # Take into account the geometric distortion term.
    mginca_ct_eda_data = recalculate_eda_int(mginca_ct_eda_data)
    cainmg_ct_eda_data = recalculate_eda_int(cainmg_ct_eda_data)
    mginca_ct_sapt = recalculate_sapt0_int(mginca_ct_sapt)
    cainmg_ct_sapt = recalculate_sapt0_int(cainmg_ct_sapt)

    ##########

    rows_combined = [
        r"\begin{tabular}{cSSSS}",
        r"\toprule",
        r"contribution & \multicolumn{2}{c}{ALMO-EDA} & \multicolumn{2}{c}{SAPT0 (monomer basis)} \\",
        r"\cmidrule{2-5}",
        r"(\si{\kcal\per\mol}) & \mg{} & \ca{} & \mg{} & \ca{} \\",
        r"\midrule",
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{geom}}\)",
            mg_ct_eda_data[_GEOM_KEY],
            ca_ct_eda_data[_GEOM_KEY],
            mg_ct_sapt["monomer"][_GEOM_KEY],
            ca_ct_sapt["monomer"][_GEOM_KEY],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{solv}}\)",
            mg_ct_eda_data["sol_elec"],
            ca_ct_eda_data["sol_elec"],
            0.0,
            0.0,
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{elec}}\)",
            mg_ct_eda_data["cls_elec"],
            ca_ct_eda_data["cls_elec"],
            mg_ct_sapt["monomer"]["el"],
            ca_ct_sapt["monomer"]["el"],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{Pauli}}\)",
            mg_ct_eda_data["mod_pauli(solv)"],
            ca_ct_eda_data["mod_pauli(solv)"],
            mg_ct_sapt["monomer"]["exch"],
            ca_ct_sapt["monomer"]["exch"],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{disp}}\)",
            mg_ct_eda_data["disp"],
            ca_ct_eda_data["disp"],
            mg_ct_sapt["monomer"]['"disp"'],
            ca_ct_sapt["monomer"]['"disp"'],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{pol}}\)",
            mg_ct_eda_data["pol"],
            ca_ct_eda_data["pol"],
            mg_ct_sapt["monomer"]['"pol"'],
            ca_ct_sapt["monomer"]['"pol"'],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{CT}}\)",
            mg_ct_eda_data["vct"],
            ca_ct_eda_data["vct"],
            mg_ct_sapt["monomer"]["ct"],
            ca_ct_sapt["monomer"]["ct"],
        ),
        r"\midrule",
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{int}}\)",
            mg_ct_eda_data["int"],
            ca_ct_eda_data["int"],
            mg_ct_sapt["monomer"]["total"],
            ca_ct_sapt["monomer"]["total"],
        ),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    table_combined = NoEscape("\n".join(rows_combined))

    with open("table_combined.tex", "w") as fh:
        fh.write(table_combined)
        fh.write("\n")

    rows_ratio = [
        r"\begin{tabular}{cSS}",
        r"\toprule",
        r"contribution & {ALMO-EDA} & {SAPT0 (monomer basis)} \\",
        r"\midrule",
        SINGLEROWFMT(
            r"\(\Delta E_{\textrm{geom}}\)",
            mg_ca_ratio_eda[_GEOM_KEY],
            mg_ca_ratio_sapt["monomer"][_GEOM_KEY],
        ),
        SINGLEROWFMT(
            r"\(\Delta E_{\textrm{solv}}\)",
            mg_ca_ratio_eda["sol_elec"],
            "N/A",
        ),
        SINGLEROWFMT(
            r"\(\Delta E_{\textrm{elec}}\)",
            mg_ca_ratio_eda["cls_elec"],
            mg_ca_ratio_sapt["monomer"]["el"],
        ),
        SINGLEROWFMT(
            r"\(\Delta E_{\textrm{Pauli}}\)",
            mg_ca_ratio_eda["mod_pauli(solv)"],
            mg_ca_ratio_sapt["monomer"]["exch"],
        ),
        SINGLEROWFMT(
            r"\(\Delta E_{\textrm{disp}}\)",
            mg_ca_ratio_eda["disp"],
            mg_ca_ratio_sapt["monomer"]['"disp"'],
        ),
        SINGLEROWFMT(
            r"\(\Delta E_{\textrm{pol}}\)",
            mg_ca_ratio_eda["pol"],
            mg_ca_ratio_sapt["monomer"]['"pol"'],
        ),
        SINGLEROWFMT(
            r"\(\Delta E_{\textrm{CT}}\)",
            mg_ca_ratio_eda["vct"],
            mg_ca_ratio_sapt["monomer"]["ct"],
        ),
        r"\midrule",
        SINGLEROWFMT(
            r"\(\Delta E_{\textrm{int}}\)",
            mg_ca_ratio_eda["int"],
            mg_ca_ratio_sapt["monomer"]["total"],
        ),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    table_ratio = NoEscape("\n".join(rows_ratio))

    with open("table_ratio.tex", "w") as fh:
        fh.write(table_ratio)
        fh.write("\n")

    rows_combined_nosolv = [
        r"\begin{tabular}{cSSSS}",
        r"\toprule",
        r"contribution & \multicolumn{2}{c}{ALMO-EDA} & \multicolumn{2}{c}{SAPT0 (monomer basis)} \\",
        r"\cmidrule{2-5}",
        r"(\si{\kcal\per\mol}) & \mg{} & \ca{} & \mg{} & \ca{} \\",
        r"\midrule",
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{geom}}\)",
            mg_ct_eda_data[_GEOM_KEY],
            ca_ct_eda_data[_GEOM_KEY],
            mg_ct_sapt["monomer"][_GEOM_KEY],
            ca_ct_sapt["monomer"][_GEOM_KEY],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{elec}}\)",
            mg_ct_eda_data["cls_elec"],
            ca_ct_eda_data["cls_elec"],
            mg_ct_sapt["monomer"]["el"],
            ca_ct_sapt["monomer"]["el"],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{Pauli}}\)",
            mg_ct_eda_data["mod_pauli"],
            ca_ct_eda_data["mod_pauli"],
            mg_ct_sapt["monomer"]["exch"],
            ca_ct_sapt["monomer"]["exch"],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{disp}}\)",
            mg_ct_eda_data["disp"],
            ca_ct_eda_data["disp"],
            mg_ct_sapt["monomer"]['"disp"'],
            ca_ct_sapt["monomer"]['"disp"'],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{pol}}\)",
            mg_ct_eda_data["pol"],
            ca_ct_eda_data["pol"],
            mg_ct_sapt["monomer"]['"pol"'],
            ca_ct_sapt["monomer"]['"pol"'],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{CT}}\)",
            mg_ct_eda_data["vct"],
            ca_ct_eda_data["vct"],
            mg_ct_sapt["monomer"]["ct"],
            ca_ct_sapt["monomer"]["ct"],
        ),
        r"\midrule",
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{int}}\)",
            mg_ct_eda_data["int"] - mg_ct_eda_data["sol"],
            ca_ct_eda_data["int"] - ca_ct_eda_data["sol"],
            mg_ct_sapt["monomer"]["total"],
            ca_ct_sapt["monomer"]["total"],
        ),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    table_combined_nosolv = NoEscape("\n".join(rows_combined_nosolv))

    with open("table_nosolv.tex", "w") as fh:
        fh.write(table_combined_nosolv)
        fh.write("\n")

    rows_percentages = [
        r"\begin{tabular}{cSSSS}",
        r"\toprule",
        r"contribution & \multicolumn{2}{c}{ALMO-EDA} & \multicolumn{2}{c}{SAPT0 (monomer basis)} \\",
        r"\cmidrule{2-5}",
        r"(\% of \(\Delta E_{\textrm{tot}}\)) & \mg{} & \ca{} & \mg{} & \ca{} \\",
        r"\midrule",
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{geom}}\)",
            mg_ct_eda_almo_pct[_GEOM_KEY],
            ca_ct_eda_almo_pct[_GEOM_KEY],
            mg_ct_sapt_pct["monomer"][_GEOM_KEY],
            ca_ct_sapt_pct["monomer"][_GEOM_KEY],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{solv}}\)",
            mg_ct_eda_almo_pct["sol_elec"],
            ca_ct_eda_almo_pct["sol_elec"],
            0.0,
            0.0,
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{elec}}\)",
            mg_ct_eda_almo_pct["cls_elec"],
            ca_ct_eda_almo_pct["cls_elec"],
            mg_ct_sapt_pct["monomer"]["el"],
            ca_ct_sapt_pct["monomer"]["el"],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{Pauli}}\)",
            mg_ct_eda_almo_pct["mod_pauli(solv)"],
            ca_ct_eda_almo_pct["mod_pauli(solv)"],
            mg_ct_sapt_pct["monomer"]["exch"],
            ca_ct_sapt_pct["monomer"]["exch"],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{disp}}\)",
            mg_ct_eda_almo_pct["disp"],
            ca_ct_eda_almo_pct["disp"],
            mg_ct_sapt_pct["monomer"]['"disp"'],
            ca_ct_sapt_pct["monomer"]['"disp"'],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{pol}}\)",
            mg_ct_eda_almo_pct["pol"],
            ca_ct_eda_almo_pct["pol"],
            mg_ct_sapt_pct["monomer"]['"pol"'],
            ca_ct_sapt_pct["monomer"]['"pol"'],
        ),
        DOUBLEROWFMT(
            r"\(\Delta E_{\textrm{CT}}\)",
            mg_ct_eda_almo_pct["vct"],
            ca_ct_eda_almo_pct["vct"],
            mg_ct_sapt_pct["monomer"]["ct"],
            ca_ct_sapt_pct["monomer"]["ct"],
        ),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    table_percentages = NoEscape("\n".join(rows_percentages))

    with open("table_percentages.tex", "w") as fh:
        fh.write(table_percentages)
        fh.write("\n")

    rows_dimer_combined = [
        r"\begin{tabular}{cSSSSSS}",
        r"\toprule",
        r"contribution & \multicolumn{2}{c}{ALMO-EDA} & \multicolumn{2}{c}{SAPT0 (monomer basis)} & \multicolumn{2}{c}{SAPT0 (dimer basis)} \\",
        r"\cmidrule{2-7}",
        r"(\si{\kcal\per\mol}) & \mg{} & \ca{} & \mg{} & \ca{} & \mg{} & \ca{} \\",
        r"\midrule",
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{geom}}\)",
            mg_ct_eda_data[_GEOM_KEY],
            ca_ct_eda_data[_GEOM_KEY],
            mg_ct_sapt["monomer"][_GEOM_KEY],
            ca_ct_sapt["monomer"][_GEOM_KEY],
            mg_ct_sapt["dimer"][_GEOM_KEY],
            ca_ct_sapt["dimer"][_GEOM_KEY],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{solv}}\)",
            mg_ct_eda_data["sol_elec"],
            ca_ct_eda_data["sol_elec"],
            0.0,
            0.0,
            0.0,
            0.0,
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{elec}}\)",
            mg_ct_eda_data["cls_elec"],
            ca_ct_eda_data["cls_elec"],
            mg_ct_sapt["monomer"]["el"],
            ca_ct_sapt["monomer"]["el"],
            mg_ct_sapt["dimer"]["el"],
            ca_ct_sapt["dimer"]["el"],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{Pauli}}\)",
            mg_ct_eda_data["mod_pauli(solv)"],
            ca_ct_eda_data["mod_pauli(solv)"],
            mg_ct_sapt["monomer"]["exch"],
            ca_ct_sapt["monomer"]["exch"],
            mg_ct_sapt["dimer"]["exch"],
            ca_ct_sapt["dimer"]["exch"],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{disp}}\)",
            mg_ct_eda_data["disp"],
            ca_ct_eda_data["disp"],
            mg_ct_sapt["monomer"]['"disp"'],
            ca_ct_sapt["monomer"]['"disp"'],
            mg_ct_sapt["dimer"]['"disp"'],
            ca_ct_sapt["dimer"]['"disp"'],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{pol}}\)",
            mg_ct_eda_data["pol"],
            ca_ct_eda_data["pol"],
            mg_ct_sapt["monomer"]['"pol"'],
            ca_ct_sapt["monomer"]['"pol"'],
            mg_ct_sapt["dimer"]['"pol"'],
            ca_ct_sapt["dimer"]['"pol"'],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{CT}}\)",
            mg_ct_eda_data["vct"],
            ca_ct_eda_data["vct"],
            mg_ct_sapt["monomer"]["ct"],
            ca_ct_sapt["monomer"]["ct"],
            mg_ct_sapt["dimer"]["ct"],
            ca_ct_sapt["dimer"]["ct"],
        ),
        r"\midrule",
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{int}}\)",
            mg_ct_eda_data["int"],
            ca_ct_eda_data["int"],
            mg_ct_sapt["monomer"]["total"],
            ca_ct_sapt["monomer"]["total"],
            mg_ct_sapt["dimer"]["total"],
            ca_ct_sapt["dimer"]["total"],
        ),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    table_dimer_combined = NoEscape("\n".join(rows_dimer_combined))

    with open("table_dimer_combined.tex", "w") as fh:
        fh.write(table_dimer_combined)
        fh.write("\n")

    rows_dimer_combined_nosolv = [
        r"\begin{tabular}{cSSSSSS}",
        r"\toprule",
        r"contribution & \multicolumn{2}{c}{ALMO-EDA} & \multicolumn{2}{c}{SAPT0 (monomer basis)} & \multicolumn{2}{c}{SAPT0 (dimer basis)} \\",
        r"\cmidrule{2-7}",
        r"(\si{\kcal\per\mol}) & \mg{} & \ca{} & \mg{} & \ca{} & \mg{} & \ca{} \\",
        r"\midrule",
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{geom}}\)",
            mg_ct_eda_data[_GEOM_KEY],
            ca_ct_eda_data[_GEOM_KEY],
            mg_ct_sapt["monomer"][_GEOM_KEY],
            ca_ct_sapt["monomer"][_GEOM_KEY],
            mg_ct_sapt["dimer"][_GEOM_KEY],
            ca_ct_sapt["dimer"][_GEOM_KEY],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{elec}}\)",
            mg_ct_eda_data["cls_elec"],
            ca_ct_eda_data["cls_elec"],
            mg_ct_sapt["monomer"]["el"],
            ca_ct_sapt["monomer"]["el"],
            mg_ct_sapt["dimer"]["el"],
            ca_ct_sapt["dimer"]["el"],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{Pauli}}\)",
            mg_ct_eda_data["mod_pauli"],
            ca_ct_eda_data["mod_pauli"],
            mg_ct_sapt["monomer"]["exch"],
            ca_ct_sapt["monomer"]["exch"],
            mg_ct_sapt["dimer"]["exch"],
            ca_ct_sapt["dimer"]["exch"],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{disp}}\)",
            mg_ct_eda_data["disp"],
            ca_ct_eda_data["disp"],
            mg_ct_sapt["monomer"]['"disp"'],
            ca_ct_sapt["monomer"]['"disp"'],
            mg_ct_sapt["dimer"]['"disp"'],
            ca_ct_sapt["dimer"]['"disp"'],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{pol}}\)",
            mg_ct_eda_data["pol"],
            ca_ct_eda_data["pol"],
            mg_ct_sapt["monomer"]['"pol"'],
            ca_ct_sapt["monomer"]['"pol"'],
            mg_ct_sapt["dimer"]['"pol"'],
            ca_ct_sapt["dimer"]['"pol"'],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{CT}}\)",
            mg_ct_eda_data["vct"],
            ca_ct_eda_data["vct"],
            mg_ct_sapt["monomer"]["ct"],
            ca_ct_sapt["monomer"]["ct"],
            mg_ct_sapt["dimer"]["ct"],
            ca_ct_sapt["dimer"]["ct"],
        ),
        r"\midrule",
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{int}}\)",
            mg_ct_eda_data["int"] - mg_ct_eda_data["sol"],
            ca_ct_eda_data["int"] - ca_ct_eda_data["sol"],
            mg_ct_sapt["monomer"]["total"],
            ca_ct_sapt["monomer"]["total"],
            mg_ct_sapt["dimer"]["total"],
            ca_ct_sapt["dimer"]["total"],
        ),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    table_dimer_combined_nosolv = NoEscape("\n".join(rows_dimer_combined_nosolv))

    with open("table_dimer_nosolv.tex", "w") as fh:
        fh.write(table_dimer_combined_nosolv)
        fh.write("\n")

    rows_dimer_percentages = [
        r"\begin{tabular}{cSSSSSS}",
        r"\toprule",
        r"contribution & \multicolumn{2}{c}{ALMO-EDA} & \multicolumn{2}{c}{SAPT0 (monomer basis)} & \multicolumn{2}{c}{SAPT0 (dimer basis)} \\",
        r"\cmidrule{2-7}",
        r"(\% of \(\Delta E_{\textrm{tot}}\)) & \mg{} & \ca{} & \mg{} & \ca{} & \mg{} & \ca{} \\",
        r"\midrule",
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{geom}}\)",
            mg_ct_eda_almo_pct[_GEOM_KEY],
            ca_ct_eda_almo_pct[_GEOM_KEY],
            mg_ct_sapt_pct["monomer"][_GEOM_KEY],
            ca_ct_sapt_pct["monomer"][_GEOM_KEY],
            mg_ct_sapt_pct["dimer"][_GEOM_KEY],
            ca_ct_sapt_pct["dimer"][_GEOM_KEY],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{solv}}\)",
            mg_ct_eda_almo_pct["sol_elec"],
            ca_ct_eda_almo_pct["sol_elec"],
            0.0,
            0.0,
            0.0,
            0.0,
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{elec}}\)",
            mg_ct_eda_almo_pct["cls_elec"],
            ca_ct_eda_almo_pct["cls_elec"],
            mg_ct_sapt_pct["monomer"]["el"],
            ca_ct_sapt_pct["monomer"]["el"],
            mg_ct_sapt_pct["dimer"]["el"],
            ca_ct_sapt_pct["dimer"]["el"],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{Pauli}}\)",
            mg_ct_eda_almo_pct["mod_pauli(solv)"],
            ca_ct_eda_almo_pct["mod_pauli(solv)"],
            mg_ct_sapt_pct["monomer"]["exch"],
            ca_ct_sapt_pct["monomer"]["exch"],
            mg_ct_sapt_pct["dimer"]["exch"],
            ca_ct_sapt_pct["dimer"]["exch"],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{disp}}\)",
            mg_ct_eda_almo_pct["disp"],
            ca_ct_eda_almo_pct["disp"],
            mg_ct_sapt_pct["monomer"]['"disp"'],
            ca_ct_sapt_pct["monomer"]['"disp"'],
            mg_ct_sapt_pct["dimer"]['"disp"'],
            ca_ct_sapt_pct["dimer"]['"disp"'],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{pol}}\)",
            mg_ct_eda_almo_pct["pol"],
            ca_ct_eda_almo_pct["pol"],
            mg_ct_sapt_pct["monomer"]['"pol"'],
            ca_ct_sapt_pct["monomer"]['"pol"'],
            mg_ct_sapt_pct["dimer"]['"pol"'],
            ca_ct_sapt_pct["dimer"]['"pol"'],
        ),
        TRIPLEROWFMT(
            r"\(\Delta E_{\textrm{CT}}\)",
            mg_ct_eda_almo_pct["vct"],
            ca_ct_eda_almo_pct["vct"],
            mg_ct_sapt_pct["monomer"]["ct"],
            ca_ct_sapt_pct["monomer"]["ct"],
            mg_ct_sapt_pct["dimer"]["ct"],
            ca_ct_sapt_pct["dimer"]["ct"],
        ),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    table_dimer_percentages = NoEscape("\n".join(rows_dimer_percentages))

    with open("table_dimer_percentages.tex", "w") as fh:
        fh.write(table_dimer_percentages)
        fh.write("\n")

    rows_combined_swap = [
        r"\begin{tabular}{cSSSSSSSS}",
        r"\toprule",
        r"contribution & \multicolumn{4}{c}{ALMO-EDA} & \multicolumn{4}{c}{SAPT0 (monomer basis)} \\",
        r"\cmidrule{2-9}",
        r"(\si{\kcal\per\mol}) & \mg{} & \mg{}* & \ca{} & \ca{}* & \mg{} & \mg{}* & \ca{} & \ca{}* \\",
        r"\midrule",
        QUADROWFMT(
            r"\(\Delta E_{\textrm{geom}}\)",
            mg_ct_eda_data[_GEOM_KEY],
            mginca_ct_eda_data[_GEOM_KEY],
            ca_ct_eda_data[_GEOM_KEY],
            cainmg_ct_eda_data[_GEOM_KEY],
            mg_ct_sapt["monomer"][_GEOM_KEY],
            mginca_ct_sapt["monomer"][_GEOM_KEY],
            ca_ct_sapt["monomer"][_GEOM_KEY],
            cainmg_ct_sapt["monomer"][_GEOM_KEY],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{solv}}\)",
            mg_ct_eda_data["sol_elec"],
            mginca_ct_eda_data["sol_elec"],
            ca_ct_eda_data["sol_elec"],
            cainmg_ct_eda_data["sol_elec"],
            0.0,
            0.0,
            0.0,
            0.0,
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{elec}}\)",
            mg_ct_eda_data["cls_elec"],
            mginca_ct_eda_data["cls_elec"],
            ca_ct_eda_data["cls_elec"],
            cainmg_ct_eda_data["cls_elec"],
            mg_ct_sapt["monomer"]["el"],
            mginca_ct_sapt["monomer"]["el"],
            ca_ct_sapt["monomer"]["el"],
            cainmg_ct_sapt["monomer"]["el"],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{Pauli}}\)",
            mg_ct_eda_data["mod_pauli(solv)"],
            mginca_ct_eda_data["mod_pauli(solv)"],
            ca_ct_eda_data["mod_pauli(solv)"],
            cainmg_ct_eda_data["mod_pauli(solv)"],
            mg_ct_sapt["monomer"]["exch"],
            mginca_ct_sapt["monomer"]["exch"],
            ca_ct_sapt["monomer"]["exch"],
            cainmg_ct_sapt["monomer"]["exch"],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{disp}}\)",
            mg_ct_eda_data["disp"],
            mginca_ct_eda_data["disp"],
            ca_ct_eda_data["disp"],
            cainmg_ct_eda_data["disp"],
            mg_ct_sapt["monomer"]['"disp"'],
            mginca_ct_sapt["monomer"]['"disp"'],
            ca_ct_sapt["monomer"]['"disp"'],
            cainmg_ct_sapt["monomer"]['"disp"'],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{pol}}\)",
            mg_ct_eda_data["pol"],
            mginca_ct_eda_data["pol"],
            ca_ct_eda_data["pol"],
            cainmg_ct_eda_data["pol"],
            mg_ct_sapt["monomer"]['"pol"'],
            mginca_ct_sapt["monomer"]['"pol"'],
            ca_ct_sapt["monomer"]['"pol"'],
            cainmg_ct_sapt["monomer"]['"pol"'],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{CT}}\)",
            mg_ct_eda_data["vct"],
            mginca_ct_eda_data["vct"],
            ca_ct_eda_data["vct"],
            cainmg_ct_eda_data["vct"],
            mg_ct_sapt["monomer"]["ct"],
            mginca_ct_sapt["monomer"]["ct"],
            ca_ct_sapt["monomer"]["ct"],
            cainmg_ct_sapt["monomer"]["ct"],
        ),
        r"\midrule",
        QUADROWFMT(
            r"\(\Delta E_{\textrm{int}}\)",
            mg_ct_eda_data["int"],
            mginca_ct_eda_data["int"],
            ca_ct_eda_data["int"],
            cainmg_ct_eda_data["int"],
            mg_ct_sapt["monomer"]["total"],
            mginca_ct_sapt["monomer"]["total"],
            ca_ct_sapt["monomer"]["total"],
            cainmg_ct_sapt["monomer"]["total"],
        ),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    table_combined_swap = NoEscape("\n".join(rows_combined_swap))

    with open("table_swap.tex", "w") as fh:
        fh.write(table_combined_swap)
        fh.write("\n")

    rows_combined_swap_nosolv = [
        r"\begin{tabular}{cSSSSSSSS}",
        r"\toprule",
        r"contribution & \multicolumn{4}{c}{ALMO-EDA} & \multicolumn{4}{c}{SAPT0 (monomer basis)} \\",
        r"\cmidrule{2-9}",
        r"(\si{\kcal\per\mol}) & \mg{} & \mg{}* & \ca{} & \ca{}* & \mg{} & \mg{}* & \ca{} & \ca{}* \\",
        r"\midrule",
        QUADROWFMT(
            r"\(\Delta E_{\textrm{geom}}\)",
            mg_ct_eda_data[_GEOM_KEY],
            mginca_ct_eda_data[_GEOM_KEY],
            ca_ct_eda_data[_GEOM_KEY],
            cainmg_ct_eda_data[_GEOM_KEY],
            mg_ct_sapt["monomer"][_GEOM_KEY],
            mginca_ct_sapt["monomer"][_GEOM_KEY],
            ca_ct_sapt["monomer"][_GEOM_KEY],
            cainmg_ct_sapt["monomer"][_GEOM_KEY],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{elec}}\)",
            mg_ct_eda_data["cls_elec"],
            mginca_ct_eda_data["cls_elec"],
            ca_ct_eda_data["cls_elec"],
            cainmg_ct_eda_data["cls_elec"],
            mg_ct_sapt["monomer"]["el"],
            mginca_ct_sapt["monomer"]["el"],
            ca_ct_sapt["monomer"]["el"],
            cainmg_ct_sapt["monomer"]["el"],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{Pauli}}\)",
            mg_ct_eda_data["mod_pauli"],
            mginca_ct_eda_data["mod_pauli"],
            ca_ct_eda_data["mod_pauli"],
            cainmg_ct_eda_data["mod_pauli"],
            mg_ct_sapt["monomer"]["exch"],
            mginca_ct_sapt["monomer"]["exch"],
            ca_ct_sapt["monomer"]["exch"],
            cainmg_ct_sapt["monomer"]["exch"],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{disp}}\)",
            mg_ct_eda_data["disp"],
            mginca_ct_eda_data["disp"],
            ca_ct_eda_data["disp"],
            cainmg_ct_eda_data["disp"],
            mg_ct_sapt["monomer"]['"disp"'],
            mginca_ct_sapt["monomer"]['"disp"'],
            ca_ct_sapt["monomer"]['"disp"'],
            cainmg_ct_sapt["monomer"]['"disp"'],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{pol}}\)",
            mg_ct_eda_data["pol"],
            mginca_ct_eda_data["pol"],
            ca_ct_eda_data["pol"],
            cainmg_ct_eda_data["pol"],
            mg_ct_sapt["monomer"]['"pol"'],
            mginca_ct_sapt["monomer"]['"pol"'],
            ca_ct_sapt["monomer"]['"pol"'],
            cainmg_ct_sapt["monomer"]['"pol"'],
        ),
        QUADROWFMT(
            r"\(\Delta E_{\textrm{CT}}\)",
            mg_ct_eda_data["vct"],
            mginca_ct_eda_data["vct"],
            ca_ct_eda_data["vct"],
            cainmg_ct_eda_data["vct"],
            mg_ct_sapt["monomer"]["ct"],
            mginca_ct_sapt["monomer"]["ct"],
            ca_ct_sapt["monomer"]["ct"],
            cainmg_ct_sapt["monomer"]["ct"],
        ),
        r"\midrule",
        QUADROWFMT(
            r"\(\Delta E_{\textrm{int}}\)",
            mg_ct_eda_data["int"] - mg_ct_eda_data["sol"],
            mginca_ct_eda_data["int"] - mginca_ct_eda_data["sol"],
            ca_ct_eda_data["int"] - ca_ct_eda_data["sol"],
            cainmg_ct_eda_data["int"] - cainmg_ct_eda_data["sol"],
            mg_ct_sapt["monomer"]["total"],
            mginca_ct_sapt["monomer"]["total"],
            ca_ct_sapt["monomer"]["total"],
            cainmg_ct_sapt["monomer"]["total"],
        ),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    table_combined_swap_nosolv = NoEscape("\n".join(rows_combined_swap_nosolv))

    with open("table_swap_nosolv.tex", "w") as fh:
        fh.write(table_combined_swap_nosolv)
        fh.write("\n")
