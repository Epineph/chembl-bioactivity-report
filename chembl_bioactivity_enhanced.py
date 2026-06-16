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

PK_SIMULATION_NOTE = (
    "Dose-curve simulations assume a healthy 70 kg young adult by default, a "
    "one-compartment model, first-order absorption for non-IV routes, normal renal "
    "and CYP activity, and no inhibitors or inducers. They are for exploratory "
    "comparison only and must not be used for prescribing or self-dosing."
)

_PUG_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
    "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,"
    "RotatableBondCount,Charge,CanonicalSMILES,IsomericSMILES/JSON"
)

_DOSE_UNIT_TO_MG = {
    "microgram": 0.001,
    "micrograms": 0.001,
    "ug": 0.001,
    "mcg": 0.001,
    "milligram": 1.0,
    "milligrams": 1.0,
    "mg": 1.0,
    "gram": 1000.0,
    "grams": 1000.0,
    "g": 1000.0,
}

ROUTE_DISPLAY_NAMES = {
    "oral": "Oral",
    "sublingual": "Sublingual",
    "intranasal": "Intranasal",
    "subcutaneous": "Subcutaneous",
    "intravenous": "Intravenous",
}

ACTIVE_METABOLITE_NOTES = {
    "codeine": {
        "Active metabolite(s)": "Morphine; codeine-6-glucuronide may also contribute",
        "Main pathway": "CYP2D6 O-demethylation to morphine; CYP3A4 forms norcodeine",
        "PK implication": "Parent-codeine concentration is a poor analgesic-effect proxy. IV codeine bypasses gut absorption but does not remove the need for CYP2D6 activation.",
    },
    "tramadol": {
        "Active metabolite(s)": "O-desmethyltramadol (M1)",
        "Main pathway": "CYP2D6 forms M1; CYP3A4/CYP2B6 form less active N-desmethyltramadol",
        "PK implication": "Analgesic effect depends partly on active metabolite formation; parent-only curves can misrepresent onset and potency.",
    },
    "diazepam": {
        "Active metabolite(s)": "Nordazepam/desmethyldiazepam, temazepam, oxazepam",
        "Main pathway": "CYP3A4/CYP2C19 N-demethylation and hydroxylation",
        "PK implication": "Parent diazepam curve can underestimate total active-moiety duration because active metabolites persist.",
    },
    "clopidogrel": {
        "Active metabolite(s)": "Thiol active metabolite",
        "Main pathway": "Two-step CYP oxidation, including CYP2C19",
        "PK implication": "Parent concentration is not a useful antiplatelet-effect proxy; activation is required.",
    },
    "prednisone": {
        "Active metabolite(s)": "Prednisolone",
        "Main pathway": "Hepatic 11-beta-hydroxysteroid dehydrogenase activation",
        "PK implication": "Therapeutic effect mainly follows prednisolone exposure rather than parent prednisone exposure.",
    },
    "enalapril": {
        "Active metabolite(s)": "Enalaprilat",
        "Main pathway": "Esterase hydrolysis",
        "PK implication": "Prodrug activation is required; route changes can alter the parent-to-active exposure relationship.",
    },
    "lisdexamfetamine": {
        "Active metabolite(s)": "Dextroamphetamine",
        "Main pathway": "Enzymatic cleavage after absorption, mainly in blood",
        "PK implication": "Parent concentration is not the desired stimulant-effect curve; active amphetamine formation controls onset.",
    },
    "risperidone": {
        "Active metabolite(s)": "9-hydroxyrisperidone/paliperidone",
        "Main pathway": "CYP2D6 hydroxylation",
        "PK implication": "Parent plus active metabolite better represent active moiety than parent alone.",
    },
    "fluoxetine": {
        "Active metabolite(s)": "Norfluoxetine",
        "Main pathway": "CYP-mediated N-demethylation",
        "PK implication": "Long-lived active metabolite can dominate duration after repeated dosing; single-dose parent curve is incomplete.",
    },
}


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
    text = f"{value:.{digits}f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


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


def dose_to_mg(amount: Any, unit: str) -> float:
    """Convert a user-entered dose amount to mg."""
    dose = _as_float(amount, 0.0) or 0.0
    multiplier = _DOSE_UNIT_TO_MG.get((unit or "mg").strip().lower())
    if multiplier is None:
        raise ValueError(f"Unsupported dose unit: {unit}")
    return max(0.0, dose * multiplier)


def _normalize_route(route: str) -> str:
    text = (route or "oral").strip().lower().replace("-", " ").replace("_", " ")
    aliases = {
        "po": "oral",
        "by mouth": "oral",
        "sl": "sublingual",
        "sublingual": "sublingual",
        "in": "intranasal",
        "nasal": "intranasal",
        "intranasal": "intranasal",
        "sc": "subcutaneous",
        "subcutaneous": "subcutaneous",
        "sub cutaneous": "subcutaneous",
        "iv": "intravenous",
        "intravenous": "intravenous",
    }
    return aliases.get(text, text)


def route_parameters(route: str, prediction: dict[str, Any]) -> dict[str, Any]:
    """Estimate route-specific bioavailability, absorption rate, and lag."""
    route_key = _normalize_route(route)
    if route_key not in ROUTE_DISPLAY_NAMES:
        raise ValueError(f"Unsupported route: {route}")

    fa = _clip(float(prediction.get("fraction_absorbed_pct", 0.0)) / 100.0, 0.01, 0.98)
    oral_f = _clip(float(prediction.get("oral_bioavailability_pct", 0.0)) / 100.0, 0.01, 0.95)
    tpsa = float(prediction.get("tpsa", 0.0) or 0.0)
    rotb = float(prediction.get("rotb", 0.0) or 0.0)
    mw = float(prediction.get("mw", 0.0) or 0.0)
    charge = float(prediction.get("charge", 0.0) or 0.0)
    permeability = _clip(1.15 - max(0.0, tpsa - 75.0) / 220.0 - max(0.0, rotb - 6.0) / 18.0, 0.35, 1.35)

    if route_key == "oral":
        return {
            "route": "Oral",
            "F": oral_f,
            "ka_h": _clip(1.05 * permeability, 0.25, 2.5),
            "lag_h": 0.25,
            "assumption": "GI absorption with first-pass extraction reflected in predicted oral bioavailability.",
        }
    if route_key == "sublingual":
        return {
            "route": "Sublingual",
            "F": _clip(0.15 + 0.68 * fa, 0.05, 0.85),
            "ka_h": _clip(2.4 * permeability, 0.6, 4.0),
            "lag_h": 0.05,
            "assumption": "Partial mucosal absorption with partial first-pass bypass; swallowed fraction is not separately modeled.",
        }
    if route_key == "intranasal":
        return {
            "route": "Intranasal",
            "F": _clip(0.10 + 0.62 * fa, 0.05, 0.80),
            "ka_h": _clip(2.0 * permeability, 0.5, 3.5),
            "lag_h": 0.05,
            "assumption": "Nasal mucosal absorption proxy; formulation, pH, and local tolerability are not modeled.",
        }
    if route_key == "subcutaneous":
        bioavailability = 0.92 - 0.18 * max(0.0, mw - 800.0) / 800.0 - 0.06 * (abs(charge) >= 2)
        return {
            "route": "Subcutaneous",
            "F": _clip(bioavailability, 0.45, 0.98),
            "ka_h": _clip(0.45 * permeability, 0.08, 1.2),
            "lag_h": 0.10,
            "assumption": "Depot-like first-order absorption; injection-site blood flow and formulation are not modeled.",
        }
    return {
        "route": "Intravenous",
        "F": 1.0,
        "ka_h": None,
        "lag_h": 0.0,
        "assumption": "Instantaneous IV bolus into central compartment; infusion kinetics are not modeled.",
    }


def _linear_crossing_time(t0: float, c0: float, t1: float, c1: float, threshold: float) -> float:
    if c1 == c0:
        return t1
    return t0 + (threshold - c0) * (t1 - t0) / (c1 - c0)


def _threshold_times(times: list[float], concentrations: list[float], threshold: float | None) -> dict[str, Any]:
    if threshold is None or threshold <= 0:
        return {"first_h": None, "last_h": None, "duration_h": None}

    first_h = None
    fall_h = None
    for idx, concentration in enumerate(concentrations):
        if concentration >= threshold:
            if first_h is None:
                if idx == 0:
                    first_h = times[idx]
                else:
                    first_h = _linear_crossing_time(
                        times[idx - 1], concentrations[idx - 1], times[idx], concentration, threshold
                    )
        elif first_h is not None and idx > 0 and concentrations[idx - 1] >= threshold:
            fall_h = _linear_crossing_time(times[idx - 1], concentrations[idx - 1], times[idx], concentration, threshold)
            break

    duration_h = None if first_h is None or fall_h is None else max(0.0, fall_h - first_h)
    return {"first_h": first_h, "last_h": fall_h, "duration_h": duration_h}


def _therapeutic_window_duration(
    times: list[float], concentrations: list[float], mec_mg_l: float | None, mtc_mg_l: float | None
) -> float | None:
    if mec_mg_l is None or mec_mg_l <= 0 or mtc_mg_l is None or mtc_mg_l <= mec_mg_l:
        return None
    total = 0.0
    for idx in range(1, len(times)):
        c_mid = (concentrations[idx - 1] + concentrations[idx]) / 2.0
        if mec_mg_l <= c_mid < mtc_mg_l:
            total += times[idx] - times[idx - 1]
    return total


def simulate_pk_for_route(
    prediction: dict[str, Any],
    dose_mg: float,
    route: str,
    body_weight_kg: float = 70.0,
    mec_mg_l: float | None = None,
    mtc_mg_l: float | None = None,
    duration_h: float = 24.0,
    time_step_h: float = 0.05,
    injection_concentration_mg_ml: float | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Simulate a one-compartment concentration-time curve for one route."""
    route_params = route_parameters(route, prediction)
    body_weight_kg = max(1.0, float(body_weight_kg or 70.0))
    dose_mg = max(0.0, float(dose_mg or 0.0))
    duration_h = max(0.25, float(duration_h or 24.0))
    time_step_h = _clip(float(time_step_h or 0.05), 0.01, 1.0)

    vd_l_kg = max(0.01, float(prediction.get("vd_l_kg", 0.7) or 0.7))
    clearance_ml_min_kg = max(0.01, float(prediction.get("clearance_ml_min_kg", 5.0) or 5.0))
    vd_l = vd_l_kg * body_weight_kg
    clearance_l_h = clearance_ml_min_kg * 0.06 * body_weight_kg
    ke_h = clearance_l_h / vd_l
    half_life_h = math.log(2) / ke_h if ke_h > 0 else None
    absorbed_dose_mg = route_params["F"] * dose_mg

    steps = int(math.ceil(duration_h / time_step_h))
    times = [round(idx * time_step_h, 6) for idx in range(steps + 1)]
    concentrations: list[float] = []
    ka_h = route_params["ka_h"]
    lag_h = float(route_params["lag_h"] or 0.0)

    for time_h in times:
        if dose_mg <= 0:
            concentration = 0.0
        elif ka_h is None:
            concentration = absorbed_dose_mg / vd_l * math.exp(-ke_h * time_h)
        else:
            t_eff = time_h - lag_h
            if t_eff <= 0:
                concentration = 0.0
            elif abs(ka_h - ke_h) < 1e-6:
                concentration = absorbed_dose_mg / vd_l * ka_h * t_eff * math.exp(-ke_h * t_eff)
            else:
                concentration = (
                    absorbed_dose_mg
                    * ka_h
                    / (vd_l * (ka_h - ke_h))
                    * (math.exp(-ke_h * t_eff) - math.exp(-ka_h * t_eff))
                )
        concentrations.append(max(0.0, concentration))

    cmax_mg_l = max(concentrations)
    tmax_h = times[concentrations.index(cmax_mg_l)]
    auc_0_inf_mg_h_l = absorbed_dose_mg / clearance_l_h if clearance_l_h > 0 else None
    auc_last = sum(
        (concentrations[idx - 1] + concentrations[idx]) / 2.0 * (times[idx] - times[idx - 1])
        for idx in range(1, len(times))
    )
    mec_times = _threshold_times(times, concentrations, mec_mg_l)
    mtc_times = _threshold_times(times, concentrations, mtc_mg_l)
    therapeutic_duration_h = _therapeutic_window_duration(times, concentrations, mec_mg_l, mtc_mg_l)
    volume_to_administer_ml = None
    concentration = _as_float(injection_concentration_mg_ml)
    if concentration and concentration > 0 and route_params["route"] in {"Intravenous", "Subcutaneous"}:
        volume_to_administer_ml = dose_mg / concentration

    curve = pd.DataFrame(
        {
            "Route": route_params["route"],
            "Time (h)": times,
            "Concentration (mg/L)": concentrations,
        }
    )
    summary = {
        "Route": route_params["route"],
        "Dose (mg)": dose_mg,
        "Bioavailability F": route_params["F"],
        "Absorption ka (1/h)": ka_h,
        "Lag time (h)": lag_h,
        "Vd (L/kg)": vd_l_kg,
        "CL (mL/min/kg)": clearance_ml_min_kg,
        "ke (1/h)": ke_h,
        "Half-life (h)": half_life_h,
        "Cmax (mg/L)": cmax_mg_l,
        "Tmax (h)": tmax_h,
        "AUC 0-inf (mg*h/L)": auc_0_inf_mg_h_l,
        "AUC 0-last (mg*h/L)": auc_last,
        "Onset to MEC (h)": mec_times["first_h"],
        "Falls below MEC (h)": mec_times["last_h"],
        "Duration above MEC (h)": mec_times["duration_h"],
        "Therapeutic-window duration (h)": therapeutic_duration_h,
        "First exceeds MTC (h)": mtc_times["first_h"],
        "Time above MTC (h)": mtc_times["duration_h"],
        "Injection concentration (mg/mL)": concentration,
        "Volume to administer (mL)": volume_to_administer_ml,
        "Assumption": route_params["assumption"],
    }
    return curve, summary


def simulate_pk_curves(
    prediction: dict[str, Any],
    dose_amount: Any,
    dose_unit: str,
    routes: list[str] | tuple[str, ...],
    body_weight_kg: float = 70.0,
    mec_mg_l: float | None = None,
    mtc_mg_l: float | None = None,
    duration_h: float = 24.0,
    injection_concentration_mg_ml: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Simulate concentration-time curves for one or more routes."""
    dose_mg = dose_to_mg(dose_amount, dose_unit)
    selected_routes = list(routes) or ["Oral"]
    curves = []
    summaries = []
    for route in selected_routes:
        curve, summary = simulate_pk_for_route(
            prediction,
            dose_mg=dose_mg,
            route=route,
            body_weight_kg=body_weight_kg,
            mec_mg_l=mec_mg_l if mec_mg_l and mec_mg_l > 0 else None,
            mtc_mg_l=mtc_mg_l if mtc_mg_l and mtc_mg_l > 0 else None,
            duration_h=duration_h,
            injection_concentration_mg_ml=injection_concentration_mg_ml,
        )
        curves.append(curve)
        summaries.append(summary)
    curve_df = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()
    summary_df = pd.DataFrame(summaries)
    return curve_df, summary_df


def active_metabolite_notes_for_compound(compound_name: str) -> pd.DataFrame:
    """Return curated active-metabolite caveats for common examples."""
    name = (compound_name or "").strip().casefold()
    note = ACTIVE_METABOLITE_NOTES.get(name)
    if note is None:
        note = {
            "Active metabolite(s)": "No curated active-metabolite model available",
            "Main pathway": "Unknown from the current lightweight rules",
            "PK implication": "The concentration-time curve is parent-compound only. If efficacy depends on an active metabolite or prodrug activation, the therapeutic-effect curve may differ substantially.",
        }
    return pd.DataFrame([{ "Compound": compound_name or "Compound", **note }])


def pk_formula_markdown() -> str:
    """LaTeX formulas shown in the notebook for the simulated PK model."""
    return r"""
**One-compartment simulation formulas**

Elimination and half-life:

$$k_e = \frac{CL}{V_d}$$

$$t_{1/2} = \frac{\ln(2)}{k_e}$$

IV bolus concentration:

$$C(t) = \frac{F \cdot Dose}{V_d} e^{-k_e t},\quad F=1$$

First-order absorption for oral, sublingual, intranasal, and subcutaneous routes:

$$C(t) = \frac{F \cdot Dose \cdot k_a}{V_d(k_a-k_e)}\left(e^{-k_e(t-t_{lag})} - e^{-k_a(t-t_{lag})}\right)$$

Exposure:

$$AUC_{0-\infty} = \frac{F \cdot Dose}{CL}$$

Injection volume, when a concentration is supplied:

$$Volume\ to\ administer = \frac{Prescribed\ dose}{Concentration}$$

Onset is the first time the predicted concentration reaches the entered MEC. Duration is the time above MEC; if MTC is supplied, the table also reports time above MTC and time inside the approximate therapeutic window.
""".strip()


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
