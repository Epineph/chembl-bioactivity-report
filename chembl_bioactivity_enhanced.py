"""Bioactivity cleanup and transparent structure-based PK estimates.

The pharmacokinetic estimates in this module are screening heuristics, not
clinical predictions. They intentionally expose the descriptors, rules, and
confidence level so a result can be interpreted for compounds without human
trial data while still making the uncertainty visible.
"""

from __future__ import annotations

import math
from typing import Any

import pandas as pd
import requests

try:  # Optional: available in the Binder environment, but tests should not require it.
    import pubchempy as pcp
except Exception:  # pragma: no cover - optional dependency
    pcp = None

try:  # Optional: used for SMILES-derived descriptors when available.
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
except Exception:  # pragma: no cover - optional dependency
    Chem = None
    Crippen = None
    Descriptors = None
    Lipinski = None
    rdMolDescriptors = None


PK_ESTIMATE_NOTE = (
    "Structure-based PK estimates are approximate screening values, not dosing "
    "or clinical safety guidance. Absorption and BBB predictions have stronger "
    "descriptor evidence than clearance and half-life, which are lower-confidence "
    "heuristics without in vitro or human PK data."
)

_PUG_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
    "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,"
    "RotatableBondCount,Charge,CanonicalSMILES,IsomericSMILES/JSON"
)


def _as_float(value: Any, default: float | None = None) -> float | None:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _clip(value: float, low: float, high: float) -> float:
    return max(low, min(high, value))


def _sigmoid(value: float) -> float:
    if value >= 0:
        z = math.exp(-value)
        return 1.0 / (1.0 + z)
    z = math.exp(value)
    return z / (1.0 + z)


def _first_present(data: dict[str, Any], *keys: str) -> Any:
    lower_map = {str(k).lower(): v for k, v in data.items()}
    for key in keys:
        if key in data and data[key] not in (None, ""):
            return data[key]
        value = lower_map.get(key.lower())
        if value not in (None, ""):
            return value
    return None


def _format_number(value: float | int | None, digits: int = 2) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return ""
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def _activity_value_to_nm(activity_type: str, value: Any, unit: str) -> float | None:
    numeric = _as_float(value)
    if numeric is None:
        return None

    typ = (activity_type or "").strip().upper()
    normalized_unit = (unit or "").strip().replace("μ", "u").replace("µ", "u")
    unit_lower = normalized_unit.lower()

    if typ == "KA" and unit_lower in {"m^-1", "m-1", "1/m"} and numeric != 0:
        return (1.0 / numeric) * 1e9

    multipliers = {
        "pm": 1e-3,
        "nm": 1.0,
        "um": 1e3,
        "µm": 1e3,
        "mm": 1e6,
        "m": 1e9,
    }
    return numeric * multipliers.get(unit_lower) if unit_lower in multipliers else None


def clean_activity_df(activities: list[dict[str, Any]]) -> pd.DataFrame:
    """Normalize ChEMBL activity rows and add a numeric nM column when possible."""
    rows: list[dict[str, Any]] = []
    for activity in activities:
        target = activity.get("target_pref_name") or activity.get("target_chembl_id") or "Unknown"
        activity_type = activity.get("standard_type") or ""
        value = activity.get("standard_value") or ""
        units = activity.get("standard_units") or ""
        value_nm = _activity_value_to_nm(activity_type, value, units)
        kd_nm = value_nm if (activity_type or "").strip().upper() == "KA" else None

        rows.append(
            {
                "Target": target,
                "Activity": activity_type,
                "Value": value,
                "Units": units,
                "Value_nM": round(value_nm, 3) if value_nm is not None else None,
                "Kd (nM) (from KA)": round(kd_nm, 3) if kd_nm is not None else "",
            }
        )

    columns = ["Target", "Activity", "Value", "Units", "Value_nM", "Kd (nM) (from KA)"]
    df = pd.DataFrame(rows, columns=columns)
    if df.empty:
        return df
    return df[df["Value"].astype(str).str.len() > 0].reset_index(drop=True)


def filter_by_threshold(df: pd.DataFrame, column: str = "Value_nM", max_nM: float = 10_000) -> pd.DataFrame:
    """Keep rows at or below an nM threshold while preserving rows without numeric values."""
    if df is None or df.empty or column not in df.columns:
        return df
    numeric = pd.to_numeric(df[column], errors="coerce")
    return df[numeric.isna() | (numeric <= max_nM)].reset_index(drop=True)


def descriptors_from_smiles(smiles: str) -> dict[str, Any]:
    """Calculate core PK descriptors from a SMILES string using RDKit."""
    if Chem is None:
        return {}
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return {}
    return {
        "MolecularWeight": Descriptors.MolWt(mol),
        "XLogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBondDonorCount": Lipinski.NumHDonors(mol),
        "HBondAcceptorCount": Lipinski.NumHAcceptors(mol),
        "RotatableBondCount": Lipinski.NumRotatableBonds(mol),
        "Charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        "CanonicalSMILES": Chem.MolToSmiles(mol),
        "DescriptorSource": "RDKit from SMILES",
    }


def descriptors_from_pubchem_cid(cid: int, smiles: str | None = None) -> dict[str, Any]:
    """Fetch PubChem descriptors, with RDKit used to fill gaps when possible."""
    descriptors: dict[str, Any] = {}

    if pcp is not None:
        try:
            props = [
                "MolecularWeight",
                "XLogP",
                "TPSA",
                "HBondDonorCount",
                "HBondAcceptorCount",
                "RotatableBondCount",
                "Charge",
                "CanonicalSMILES",
                "IsomericSMILES",
            ]
            df = pcp.get_properties(props, int(cid), as_dataframe=True)
            if not df.empty:
                descriptors.update(df.iloc[0].dropna().to_dict())
                descriptors["DescriptorSource"] = "PubChem computed properties"
        except Exception:
            pass

    if not descriptors:
        try:
            response = requests.get(_PUG_PROPERTY_URL.format(cid=int(cid)), timeout=20)
            if response.ok:
                properties = response.json().get("PropertyTable", {}).get("Properties", [])
                if properties:
                    descriptors.update(properties[0])
                    descriptors["DescriptorSource"] = "PubChem PUG REST computed properties"
        except Exception:
            pass

    smiles = smiles or descriptors.get("IsomericSMILES") or descriptors.get("CanonicalSMILES")
    if smiles:
        rdkit_descriptors = descriptors_from_smiles(str(smiles))
        filled_with_rdkit = False
        source_before_rdkit = descriptors.get("DescriptorSource")
        for key, value in rdkit_descriptors.items():
            if descriptors.get(key) in (None, ""):
                descriptors[key] = value
                filled_with_rdkit = True
        if filled_with_rdkit and source_before_rdkit:
            descriptors["DescriptorSource"] = f"{descriptors['DescriptorSource']} + RDKit fallback"

    return descriptors


def normalize_pk_descriptors(descriptors: dict[str, Any]) -> dict[str, float | str | None]:
    """Map PubChem/RDKit descriptor names to a stable internal schema."""
    normalized: dict[str, float | str | None] = {
        "mw": _as_float(_first_present(descriptors, "MolecularWeight", "MolWt", "MW", "molecular_weight"), 0.0),
        "logp": _as_float(_first_present(descriptors, "XLogP", "MolLogP", "LogP", "cLogP", "logp"), 0.0),
        "tpsa": _as_float(_first_present(descriptors, "TPSA", "tpsa", "TopologicalPolarSurfaceArea"), 0.0),
        "hbd": _as_float(_first_present(descriptors, "HBondDonorCount", "HBD", "NumHDonors"), 0.0),
        "hba": _as_float(_first_present(descriptors, "HBondAcceptorCount", "HBA", "NumHAcceptors"), 0.0),
        "rotb": _as_float(_first_present(descriptors, "RotatableBondCount", "RotB", "NumRotatableBonds"), 0.0),
        "charge": _as_float(_first_present(descriptors, "Charge", "FormalCharge", "formal_charge"), 0.0),
        "smiles": _first_present(descriptors, "IsomericSMILES", "CanonicalSMILES", "SMILES"),
        "source": _first_present(descriptors, "DescriptorSource", "Source") or "Provided descriptors",
    }

    if normalized["logp"] == 0.0 and normalized["smiles"]:
        rdkit_descriptors = descriptors_from_smiles(str(normalized["smiles"]))
        if rdkit_descriptors:
            normalized["logp"] = _as_float(rdkit_descriptors.get("XLogP"), normalized["logp"])
            normalized["source"] = f"{normalized['source']} + RDKit"

    return normalized


def predict_pk_from_descriptors(descriptors: dict[str, Any]) -> dict[str, Any]:
    """Return numeric and categorical PK estimates from molecular descriptors.

    The formulas combine widely used small-molecule ADME cutoffs with empirical
    descriptor relationships. They are intended for ranking and hypothesis
    generation, not individual-patient prediction.
    """
    d = normalize_pk_descriptors(descriptors)
    mw = float(d["mw"] or 0.0)
    logp = float(d["logp"] or 0.0)
    tpsa = float(d["tpsa"] or 0.0)
    hbd = float(d["hbd"] or 0.0)
    hba = float(d["hba"] or 0.0)
    rotb = float(d["rotb"] or 0.0)
    charge = float(d["charge"] or 0.0)

    lipinski_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    veber_pass = tpsa <= 140 and rotb <= 10
    egan_pass = tpsa <= 132 and logp <= 5.88

    absorption_penalty = (
        max(0.0, (tpsa - 90.0) / 55.0)
        + max(0.0, (rotb - 7.0) / 5.0)
        + max(0.0, (mw - 450.0) / 180.0)
        + max(0.0, (logp - 4.5) / 1.8)
        + max(0.0, (-1.0 - logp) / 2.0)
        + max(0.0, (hbd - 3.0) / 3.0)
        + 0.35 * lipinski_violations
    )
    fraction_absorbed_pct = _clip(100.0 * _sigmoid(2.2 - absorption_penalty), 5.0, 98.0)

    first_pass_factor = _clip(
        0.82
        - 0.05 * lipinski_violations
        - 0.06 * max(0.0, logp - 3.0)
        - 0.035 * max(0.0, rotb - 7.0)
        - 0.08 * max(0.0, (tpsa - 140.0) / 50.0)
        - 0.05 * (abs(charge) >= 2),
        0.18,
        0.9,
    )
    oral_bioavailability_pct = _clip(fraction_absorbed_pct * first_pass_factor, 1.0, 95.0)

    # Clark-style descriptor model: logBB = 0.152*cLogP - 0.0148*PSA + 0.139.
    logbb = 0.152 * logp - 0.0148 * tpsa + 0.139
    if logbb > 0.3:
        bbb_class = "High CNS penetration tendency"
    elif logbb >= -1.0:
        bbb_class = "Possible/moderate CNS penetration"
    else:
        bbb_class = "Low CNS penetration tendency"

    charge_adjustment = 0.18 if charge > 0 else (-0.18 if charge < 0 else 0.0)
    log10_vd = (
        0.16 * (logp - 1.5)
        - 0.0035 * (tpsa - 75.0)
        - 0.0008 * (mw - 350.0)
        + charge_adjustment
    )
    vd_l_kg = 10 ** _clip(log10_vd, -1.3, 1.2)

    hepatic_score = _sigmoid(-2.0 + 0.45 * logp + 0.006 * (mw - 300.0) + 0.08 * rotb - 0.008 * tpsa)
    hepatic_extraction = 0.02 + 0.55 * hepatic_score
    hepatic_cl_ml_min_kg = 20.7 * hepatic_extraction
    renal_score = _sigmoid(0.8 - 0.75 * logp + 0.012 * (tpsa - 70.0) - 0.006 * (mw - 300.0))
    renal_cl_ml_min_kg = 1.0 * renal_score
    clearance_ml_min_kg = _clip(hepatic_cl_ml_min_kg + renal_cl_ml_min_kg, 0.2, 25.0)

    half_life_h = 0.693 * vd_l_kg / (clearance_ml_min_kg * 0.06)

    if clearance_ml_min_kg < 3:
        clearance_class = "Low clearance"
    elif clearance_ml_min_kg <= 12:
        clearance_class = "Moderate clearance"
    else:
        clearance_class = "High clearance"

    if half_life_h < 3:
        half_life_class = "Short predicted half-life"
    elif half_life_h <= 12:
        half_life_class = "Moderate predicted half-life"
    else:
        half_life_class = "Long predicted half-life"

    applicability_flags = []
    if not (100 <= mw <= 900):
        applicability_flags.append("MW outside common small-molecule range")
    if not (-3 <= logp <= 8):
        applicability_flags.append("logP outside calibrated screening range")
    if tpsa > 250:
        applicability_flags.append("TPSA very high")
    if lipinski_violations >= 2:
        applicability_flags.append("multiple Lipinski violations")
    applicability = "; ".join(applicability_flags) if applicability_flags else "Within broad small-molecule descriptor range"

    return {
        "mw": mw,
        "logp": logp,
        "tpsa": tpsa,
        "hbd": hbd,
        "hba": hba,
        "rotb": rotb,
        "charge": charge,
        "source": d["source"],
        "lipinski_violations": lipinski_violations,
        "veber_pass": veber_pass,
        "egan_pass": egan_pass,
        "fraction_absorbed_pct": fraction_absorbed_pct,
        "oral_bioavailability_pct": oral_bioavailability_pct,
        "logbb": logbb,
        "bbb_class": bbb_class,
        "vd_l_kg": vd_l_kg,
        "clearance_ml_min_kg": clearance_ml_min_kg,
        "clearance_class": clearance_class,
        "half_life_h": half_life_h,
        "half_life_class": half_life_class,
        "applicability": applicability,
    }


def pk_profile_from_descriptors(descriptors: dict[str, Any], label: str = "Compound") -> pd.DataFrame:
    """Build a display-ready pharmacokinetic estimate table."""
    if not descriptors:
        return pd.DataFrame(
            [
                {
                    "Parameter": "Predicted pharmacokinetics unavailable",
                    "Estimate": "",
                    "Units": "",
                    "Interpretation": "No usable structure descriptors were available.",
                    "Evidence / formula": "Requires MW, logP, TPSA, HBD, HBA, rotatable bonds, and charge.",
                    "Confidence": "Unavailable",
                }
            ]
        )

    p = predict_pk_from_descriptors(descriptors)
    rows = [
        (
            "Descriptor source",
            p["source"],
            "",
            label,
            "Descriptors are read from PubChem or calculated from SMILES with RDKit.",
            "Input",
        ),
        ("Molecular weight", _format_number(p["mw"], 2), "g/mol", "", "Lipinski RO5 uses MW <= 500.", "Input"),
        ("XLogP / MolLogP", _format_number(p["logp"], 2), "log10", "", "Lipinski RO5 uses logP <= 5.", "Input"),
        ("Topological polar surface area", _format_number(p["tpsa"], 1), "A^2", "", "Veber oral exposure threshold commonly uses TPSA <= 140 A^2.", "Input"),
        ("H-bond donors", _format_number(p["hbd"], 0), "count", "", "Lipinski RO5 uses HBD <= 5.", "Input"),
        ("H-bond acceptors", _format_number(p["hba"], 0), "count", "", "Lipinski RO5 uses HBA <= 10.", "Input"),
        ("Rotatable bonds", _format_number(p["rotb"], 0), "count", "", "Veber oral exposure threshold commonly uses rotatable bonds <= 10.", "Input"),
        (
            "Lipinski RO5 violations",
            str(p["lipinski_violations"]),
            "count",
            "Lower is generally better for oral drug-likeness.",
            "MW > 500, logP > 5, HBD > 5, HBA > 10.",
            "Medium",
        ),
        (
            "Veber oral permeability rule",
            "Pass" if p["veber_pass"] else "Fail",
            "",
            "Pass supports oral absorption potential.",
            "TPSA <= 140 A^2 and rotatable bonds <= 10.",
            "Medium",
        ),
        (
            "Predicted fraction absorbed orally",
            _format_number(p["fraction_absorbed_pct"], 0),
            "%",
            "Descriptor-limited estimate of intestinal absorption.",
            "Logistic penalty from TPSA, rotatable bonds, MW, logP, HBD, and RO5 violations.",
            "Medium-low",
        ),
        (
            "Predicted oral bioavailability potential",
            _format_number(p["oral_bioavailability_pct"], 0),
            "%",
            "Relative oral bioavailability potential before compound-specific metabolism data.",
            "Fraction absorbed scaled by first-pass/metabolic-liability descriptor penalties.",
            "Low",
        ),
        (
            "Predicted logBB",
            _format_number(p["logbb"], 2),
            "log10 brain:blood",
            p["bbb_class"],
            "Clark-style model: 0.152 * logP - 0.0148 * TPSA + 0.139.",
            "Medium-low",
        ),
        (
            "Predicted steady-state volume of distribution",
            _format_number(p["vd_l_kg"], 2),
            "L/kg",
            "Higher values suggest more tissue distribution.",
            "Heuristic log10(Vd) from logP, TPSA, MW, and formal charge.",
            "Low",
        ),
        (
            "Predicted total clearance",
            _format_number(p["clearance_ml_min_kg"], 2),
            "mL/min/kg",
            p["clearance_class"],
            "Hepatic extraction proxy plus renal polarity proxy, capped below human hepatic blood flow.",
            "Low",
        ),
        (
            "Predicted terminal half-life",
            _format_number(p["half_life_h"], 2),
            "h",
            p["half_life_class"],
            "t1/2 = 0.693 * Vd / CL using predicted Vd and clearance.",
            "Low",
        ),
        (
            "Applicability domain",
            p["applicability"],
            "",
            PK_ESTIMATE_NOTE,
            "Best used for neutral or simply ionized small molecules with complete descriptors.",
            "Context",
        ),
    ]

    return pd.DataFrame(
        rows,
        columns=["Parameter", "Estimate", "Units", "Interpretation", "Evidence / formula", "Confidence"],
    )


def pk_profile_for_compound(compound: int | str | None, smiles: str | None = None) -> pd.DataFrame:
    """Return a PK profile for a PubChem CID or SMILES-like structure input."""
    descriptors: dict[str, Any] = {}
    label = "Compound"

    if smiles:
        descriptors.update(descriptors_from_smiles(smiles))

    if compound is not None:
        compound_text = str(compound).strip()
        if compound_text.isdigit():
            label = f"PubChem CID {compound_text}"
            fetched = descriptors_from_pubchem_cid(int(compound_text), smiles=smiles)
            for key, value in descriptors.items():
                if value not in (None, "") and fetched.get(key) in (None, ""):
                    fetched[key] = value
            descriptors = fetched
        elif not descriptors:
            label = "SMILES input"
            descriptors = descriptors_from_smiles(compound_text)

    return pk_profile_from_descriptors(descriptors, label=label)
