"""Bioactivity cleanup and transparent structure-based PK estimates.

The pharmacokinetic estimates in this module are screening heuristics, not
clinical predictions. They intentionally expose the descriptors, rules, and
confidence level so a result can be interpreted for compounds without human
trial data while still making the uncertainty visible.
"""

from __future__ import annotations

import base64
import html
import io
import json
import math
import os
import time
from pathlib import Path
from typing import Any, Sequence
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests

try:  # Optional: available in the Binder environment, but tests should not require it.
    import pubchempy as pcp
except Exception:  # pragma: no cover - optional dependency
    pcp = None

try:  # Optional: used for ChEMBL multitask target prediction.
    import onnxruntime as ort
except Exception:  # pragma: no cover - optional dependency
    ort = None

try:  # Optional: used for SMILES-derived descriptors when available.
    from rdkit import Chem, DataStructs
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
except Exception:  # pragma: no cover - optional dependency
    Chem = None
    DataStructs = None
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
    "comparison only and must not be used for prescribing or self-dosing. "
    "Elimination half-life is route-independent in this linear model; apparent "
    "terminal half-life becomes route-dependent only when absorption is slower "
    "than elimination."
)

PD_TARGET_PREDICTION_NOTE = (
    "Predicted target hypotheses use ChEMBL's ligand-based multitask model. "
    "Scores rank possible molecular targets and off-targets; they do not infer "
    "agonism, antagonism, efficacy, tissue response, clinical effect, or safety."
)

_PUG_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
    "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,"
    "RotatableBondCount,Charge,CanonicalSMILES,IsomericSMILES/JSON"
)
_CHEMBL_TARGET_URL = "https://www.ebi.ac.uk/chembl/api/data/target.json"
_CHEMBL_MULTITASK_VERSION = "chembl_36_model"
_CHEMBL_MULTITASK_BASE_URL = (
    "https://raw.githubusercontent.com/chembl/chembl_multitask_model/main/"
    f"trained_models/{_CHEMBL_MULTITASK_VERSION}"
)
_CHEMBL_MULTITASK_MODEL_FILE = "chembl_36_multitask_q8.onnx"
_CHEMBL_MULTITASK_TARGETS_FILE = "targets_36_all.json"
_CHEMBL_MULTITASK_MODEL_URL = (
    f"{_CHEMBL_MULTITASK_BASE_URL}/{_CHEMBL_MULTITASK_MODEL_FILE}"
)
_CHEMBL_MULTITASK_TARGETS_URL = (
    f"{_CHEMBL_MULTITASK_BASE_URL}/{_CHEMBL_MULTITASK_TARGETS_FILE}"
)
_CHEMBL_MULTITASK_MODEL_LABEL = "ChEMBL 36 multitask q8 ONNX"
_CHEMBL_MULTITASK_FP_SIZE = 1024
_CHEMBL_MULTITASK_FP_RADIUS = 2

_PD_TARGET_COLUMNS = [
    "Rank",
    "Target ChEMBL ID",
    "Target",
    "Organism",
    "Target Type",
    "Prediction Score",
    "Model",
    "Evidence / limitation",
]
_CHEMBL_TARGET_METADATA_CACHE: dict[str, dict[str, Any]] = {}
_CHEMBL_MULTITASK_SESSION_CACHE: dict[str, Any] = {}

_HTTP_SESSION = requests.Session()
_HTTP_SESSION.headers.update(
    {
        "User-Agent": "chembl-bioactivity-report/0.2 (+https://github.com/Epineph/chembl-bioactivity-report)",
        "Accept": "application/json, text/plain;q=0.5",
    }
)

_PK_DESCRIPTOR_ALIASES = {
    "mw": ("MolecularWeight", "MolWt", "MW", "molecular_weight"),
    "logp": ("XLogP", "MolLogP", "LogP", "cLogP", "logp"),
    "tpsa": ("TPSA", "tpsa", "TopologicalPolarSurfaceArea"),
    "hbd": ("HBondDonorCount", "HBD", "NumHDonors"),
    "hba": ("HBondAcceptorCount", "HBA", "NumHAcceptors"),
    "rotb": ("RotatableBondCount", "RotB", "NumRotatableBonds"),
    "charge": ("Charge", "FormalCharge", "formal_charge"),
}

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
        "Metabolite model": [
            {
                "name": "Morphine",
                "formation_fraction": 0.10,
                "half_life_h": 2.5,
                "vd_multiplier": 0.9,
                "potency_relative_to_parent": 15.0,
                "basis": "Approximate CYP2D6 normal-metabolizer activation proxy; parent-codeine curve is not an analgesic-effect proxy.",
            }
        ],
    },
    "tramadol": {
        "Active metabolite(s)": "O-desmethyltramadol (M1)",
        "Main pathway": "CYP2D6 forms M1; CYP3A4/CYP2B6 form less active N-desmethyltramadol",
        "PK implication": "Analgesic effect depends on both parent tramadol and M1. M1 dominates the opioid-receptor component because its MOR Ki is much lower than parent tramadol.",
        "Metabolite model": [
            {
                "name": "O-desmethyltramadol (M1)",
                "formation_fraction": 0.20,
                "half_life_h": 7.0,
                "vd_multiplier": 1.0,
                "potency_relative_to_parent": 2.4 / 0.0034,
                "basis": "MOR affinity proxy from Ki parent 2.4 uM vs M1 0.0034 uM; CYP2D6 normal-metabolizer assumption.",
            }
        ],
    },
    "diazepam": {
        "Active metabolite(s)": "Nordazepam/desmethyldiazepam, temazepam, oxazepam",
        "Main pathway": "CYP3A4/CYP2C19 N-demethylation and hydroxylation",
        "PK implication": "Parent diazepam curve can underestimate total active-moiety duration because active metabolites persist.",
        "Metabolite model": [
            {
                "name": "Nordazepam/desmethyldiazepam",
                "formation_fraction": 0.35,
                "half_life_h": 60.0,
                "vd_multiplier": 1.1,
                "potency_relative_to_parent": 0.7,
                "basis": "Approximate long-lived active-metabolite duration proxy; receptor potency differences are smaller than tramadol/M1.",
            }
        ],
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


def _http_get(url: str, timeout: float = 30.0, retries: int = 3, backoff: float = 1.6):
    """GET with a small retry loop for PubChem/PUG calls."""
    for attempt in range(retries):
        try:
            response = _HTTP_SESSION.get(url, timeout=timeout, allow_redirects=True)
            if response.status_code == 429 or 500 <= response.status_code < 600:
                time.sleep(backoff**attempt)
                continue
            return response
        except requests.RequestException:
            time.sleep(backoff**attempt)
    return None


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


def has_complete_pk_descriptors(descriptors: dict[str, Any] | None) -> bool:
    """Return whether all required PK screening descriptors are present."""
    if not descriptors:
        return False
    return all(
        _first_present(descriptors, *aliases) is not None
        for aliases in _PK_DESCRIPTOR_ALIASES.values()
    )


def missing_pk_descriptor_names(descriptors: dict[str, Any] | None) -> list[str]:
    """Return display names for required PK descriptors missing from a mapping."""
    labels = {
        "mw": "molecular weight",
        "logp": "logP",
        "tpsa": "TPSA",
        "hbd": "H-bond donors",
        "hba": "H-bond acceptors",
        "rotb": "rotatable bonds",
        "charge": "formal charge",
    }
    descriptors = descriptors or {}
    return [
        labels[key]
        for key, aliases in _PK_DESCRIPTOR_ALIASES.items()
        if _first_present(descriptors, *aliases) is None
    ]


def merge_pk_descriptor_fallback(
    primary: dict[str, Any] | None, fallback: dict[str, Any] | None
) -> dict[str, Any]:
    """Fill missing descriptor keys from a fallback source without overwriting values."""
    merged = dict(primary or {})
    if not fallback:
        return merged

    original_source = merged.get("DescriptorSource")
    fallback_source = fallback.get("DescriptorSource")
    changed = False
    for key, value in fallback.items():
        if key == "DescriptorSource" or value in (None, ""):
            continue
        if merged.get(key) in (None, ""):
            merged[key] = value
            changed = True

    if fallback_source and changed:
        if original_source:
            merged["DescriptorSource"] = (
                f"{original_source} + {fallback_source} fallback"
            )
        else:
            merged["DescriptorSource"] = fallback_source
    elif fallback_source and not original_source:
        merged["DescriptorSource"] = fallback_source
    return merged


def clean_activity_df(activities: list[dict[str, Any]]) -> pd.DataFrame:
    """Normalize ChEMBL activity rows and add a numeric nM column when possible."""
    rows: list[dict[str, Any]] = []
    for activity in activities:
        target = (
            activity.get("target_pref_name")
            or activity.get("target_chembl_id")
            or "Unknown"
        )
        activity_type = activity.get("standard_type") or ""
        value = activity.get("standard_value") or ""
        units = activity.get("standard_units") or ""
        value_nm = _activity_value_to_nm(activity_type, value, units)
        kd_nm = value_nm if (activity_type or "").strip().upper() == "KA" else None
        relation = _first_present(activity, "standard_relation", "relation") or ""
        pchembl = _first_present(activity, "pchembl_value", "pChEMBL") or ""
        assay_id = _first_present(activity, "assay_chembl_id", "assay_id") or ""
        assay_type = _first_present(activity, "assay_type") or ""
        assay_description = _first_present(activity, "assay_description") or ""
        validity = _first_present(activity, "data_validity_comment") or ""
        duplicate = _first_present(activity, "potential_duplicate")
        duplicate = "" if duplicate is None else duplicate

        rows.append(
            {
                "Target": target,
                "Activity": activity_type,
                "Relation": relation,
                "Value": value,
                "Units": units,
                "Value_nM": round(value_nm, 3) if value_nm is not None else None,
                "pChEMBL": pchembl,
                "Kd (nM) (from KA)": round(kd_nm, 3) if kd_nm is not None else "",
                "Assay ChEMBL ID": assay_id,
                "Assay Type": assay_type,
                "Assay Description": assay_description,
                "Data Validity": validity,
                "Potential Duplicate": duplicate,
            }
        )

    columns = [
        "Target",
        "Activity",
        "Relation",
        "Value",
        "Units",
        "Value_nM",
        "pChEMBL",
        "Kd (nM) (from KA)",
        "Assay ChEMBL ID",
        "Assay Type",
        "Assay Description",
        "Data Validity",
        "Potential Duplicate",
    ]
    df = pd.DataFrame(rows, columns=columns)
    if df.empty:
        return df
    return df[df["Value"].astype(str).str.len() > 0].reset_index(drop=True)


def aggregate_activity_replicates(
    df: pd.DataFrame,
    group_cols: tuple[str, ...] = (
        "Target",
        "Activity",
        "Assay ChEMBL ID",
        "Assay Type",
        "Assay Description",
        "Relation",
        "Data Validity",
        "Potential Duplicate",
    ),
    value_col: str = "Value_nM",
    sigma_cutoff: float | None = None,
) -> pd.DataFrame:
    """Average repeated pharmacodynamic values, optionally with SD clipping.

    Numeric rows are summarized within each target/activity group after
    conversion to nM. Assay and relation fields are included in the default
    grouping when present so heterogeneous rows are not averaged together.
    Non-numeric rows are retained because their units cannot be compared safely.
    """
    if df is None or df.empty or value_col not in df.columns:
        return df

    effective_group_cols = tuple(col for col in group_cols if col in df.columns)
    if not effective_group_cols:
        return df.copy()

    out_rows: list[dict[str, Any]] = []
    numeric_values = pd.to_numeric(df[value_col], errors="coerce")
    cutoff = _as_float(sigma_cutoff)
    cutoff = cutoff if cutoff is not None and cutoff > 0 else None

    grouped = df.assign(__numeric_value__=numeric_values).groupby(
        list(effective_group_cols), dropna=False, sort=False
    )
    for keys, group in grouped:
        numeric_group = group[group["__numeric_value__"].notna()]
        non_numeric_group = group[group["__numeric_value__"].isna()]

        if not numeric_group.empty:
            values = numeric_group["__numeric_value__"].astype(float)
            keep_mask = pd.Series(True, index=values.index)
            if cutoff is not None and len(values) >= 3:
                center = values.mean()
                sd = values.std(ddof=1)
                if sd > 0 and not math.isnan(sd):
                    keep_mask = (values - center).abs() <= cutoff * sd
                    if not keep_mask.any():
                        keep_mask = pd.Series(True, index=values.index)

            retained = values[keep_mask]
            mean_nm = retained.mean()
            sd_nm = retained.std(ddof=1) if len(retained) >= 2 else None
            excluded = int(len(values) - len(retained))

            first = numeric_group.iloc[0].drop(labels="__numeric_value__").to_dict()
            if not isinstance(keys, tuple):
                keys = (keys,)
            for col, key in zip(effective_group_cols, keys):
                first[col] = key

            activity_type = str(first.get("Activity", "")).strip().upper()
            first.update(
                {
                    "Value": round(mean_nm, 3),
                    "Units": "nM",
                    value_col: round(mean_nm, 3),
                    "Kd (nM) (from KA)": (
                        round(mean_nm, 3) if activity_type == "KA" else ""
                    ),
                    "Replicates": int(len(values)),
                    "Retained": int(len(retained)),
                    "Excluded": excluded,
                    "SD_nM": round(sd_nm, 3) if sd_nm is not None else "",
                }
            )
            out_rows.append(first)

        for _, row in non_numeric_group.drop(columns="__numeric_value__").iterrows():
            retained_row = row.to_dict()
            retained_row.update(
                {
                    "Replicates": "",
                    "Retained": "",
                    "Excluded": "",
                    "SD_nM": "",
                }
            )
            out_rows.append(retained_row)

    columns = list(df.columns)
    for col in ["Replicates", "Retained", "Excluded", "SD_nM"]:
        if col not in columns:
            columns.append(col)
    return pd.DataFrame(out_rows, columns=columns).reset_index(drop=True)


def filter_by_threshold(
    df: pd.DataFrame, column: str = "Value_nM", max_nM: float = 10_000
) -> pd.DataFrame:
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
            response = _http_get(_PUG_PROPERTY_URL.format(cid=int(cid)), timeout=20)
            if response and response.ok:
                properties = (
                    response.json().get("PropertyTable", {}).get("Properties", [])
                )
                if properties:
                    descriptors.update(properties[0])
                    descriptors["DescriptorSource"] = (
                        "PubChem PUG REST computed properties"
                    )
        except Exception:
            pass

    smiles = (
        smiles
        or descriptors.get("IsomericSMILES")
        or descriptors.get("CanonicalSMILES")
    )
    if smiles:
        rdkit_descriptors = descriptors_from_smiles(str(smiles))
        filled_with_rdkit = False
        source_before_rdkit = descriptors.get("DescriptorSource")
        for key, value in rdkit_descriptors.items():
            if descriptors.get(key) in (None, ""):
                descriptors[key] = value
                filled_with_rdkit = True
        if filled_with_rdkit and source_before_rdkit:
            descriptors["DescriptorSource"] = (
                f"{descriptors['DescriptorSource']} + RDKit fallback"
            )

    return descriptors


def complete_pk_descriptors(
    descriptors: dict[str, Any] | None,
    cid: int | None = None,
    smiles: str | None = None,
) -> dict[str, Any]:
    """Fill incomplete PK descriptors from PubChem and RDKit when possible."""
    completed = dict(descriptors or {})
    if has_complete_pk_descriptors(completed):
        return completed

    if cid is not None:
        completed = merge_pk_descriptor_fallback(
            completed, descriptors_from_pubchem_cid(int(cid), smiles=smiles)
        )
        smiles = (
            smiles
            or completed.get("IsomericSMILES")
            or completed.get("CanonicalSMILES")
        )
        if has_complete_pk_descriptors(completed):
            return completed

    if smiles:
        completed = merge_pk_descriptor_fallback(
            completed, descriptors_from_smiles(str(smiles))
        )
    return completed


def _chembl_multitask_cache_dir(cache_dir: str | os.PathLike[str] | None = None) -> Path:
    if cache_dir is not None:
        return Path(cache_dir).expanduser()
    env_dir = os.environ.get("CHEMBL_MULTITASK_MODEL_DIR")
    if env_dir:
        return Path(env_dir).expanduser()
    return (
        Path.home()
        / ".cache"
        / "chembl_bioactivity_report"
        / _CHEMBL_MULTITASK_VERSION
    )


def chembl_multitask_asset_paths(
    cache_dir: str | os.PathLike[str] | None = None,
) -> dict[str, Path]:
    """Return local paths for the ChEMBL multitask ONNX model assets."""
    base_dir = _chembl_multitask_cache_dir(cache_dir)
    return {
        "model": base_dir / _CHEMBL_MULTITASK_MODEL_FILE,
        "targets": base_dir / _CHEMBL_MULTITASK_TARGETS_FILE,
    }


def _download_file(url: str, destination: Path, timeout: float = 90.0) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    response = _http_get(url, timeout=timeout)
    if response is None or not response.ok:
        status = "no response" if response is None else f"HTTP {response.status_code}"
        raise RuntimeError(f"Could not download {url}: {status}.")

    temporary_path = destination.with_suffix(f"{destination.suffix}.tmp")
    temporary_path.write_bytes(response.content)
    temporary_path.replace(destination)


def ensure_chembl_multitask_assets(
    cache_dir: str | os.PathLike[str] | None = None,
    download: bool = True,
) -> dict[str, Path]:
    """Ensure the ChEMBL multitask model and ordered target list are available."""
    paths = chembl_multitask_asset_paths(cache_dir)
    missing = [name for name, path in paths.items() if not path.exists()]
    if missing and not download:
        missing_files = ", ".join(str(paths[name]) for name in missing)
        raise FileNotFoundError(
            "ChEMBL multitask model asset(s) are missing and download=False: "
            f"{missing_files}."
        )

    if "model" in missing:
        _download_file(_CHEMBL_MULTITASK_MODEL_URL, paths["model"])
    if "targets" in missing:
        _download_file(_CHEMBL_MULTITASK_TARGETS_URL, paths["targets"])
    return paths


def load_chembl_multitask_targets(
    cache_dir: str | os.PathLike[str] | None = None,
    download: bool = True,
) -> list[str]:
    """Load the ordered ChEMBL target IDs matching the multitask ONNX outputs."""
    paths = ensure_chembl_multitask_assets(cache_dir=cache_dir, download=download)
    try:
        payload = json.loads(paths["targets"].read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Could not parse ChEMBL target list: {paths['targets']}."
        ) from exc

    if isinstance(payload, list):
        target_ids = [str(item) for item in payload]
    elif isinstance(payload, dict):
        target_ids = []
        for key in ("targets", "target_ids", "target_chembl_ids"):
            value = payload.get(key)
            if isinstance(value, list):
                target_ids = [str(item) for item in value]
                break
    else:
        target_ids = []

    if not target_ids:
        raise ValueError("ChEMBL multitask target list is empty or unsupported.")
    return target_ids


def _chembl_multitask_fingerprint(smiles: str) -> np.ndarray:
    if Chem is None or DataStructs is None or rdMolDescriptors is None:
        raise RuntimeError(
            "RDKit is required for ChEMBL multitask target prediction."
        )
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        raise ValueError("Cannot predict targets from an invalid SMILES string.")

    fingerprint = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol,
        _CHEMBL_MULTITASK_FP_RADIUS,
        nBits=_CHEMBL_MULTITASK_FP_SIZE,
    )
    values = np.zeros((_CHEMBL_MULTITASK_FP_SIZE,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fingerprint, values)
    return values.astype(np.float32)


def _chembl_multitask_session(model_path: Path):
    if ort is None:
        raise RuntimeError(
            "onnxruntime is required for ChEMBL multitask target prediction."
        )

    key = str(model_path)
    if key not in _CHEMBL_MULTITASK_SESSION_CACHE:
        _CHEMBL_MULTITASK_SESSION_CACHE[key] = ort.InferenceSession(
            key,
            providers=["CPUExecutionProvider"],
        )
    return _CHEMBL_MULTITASK_SESSION_CACHE[key]


def _default_target_metadata(target_id: str) -> dict[str, Any]:
    return {
        "target_chembl_id": target_id,
        "pref_name": target_id,
        "organism": "Unknown",
        "target_type": "Unknown",
    }


def fetch_chembl_target_metadata(
    target_ids: Sequence[str],
    chunk_size: int = 100,
    timeout: float = 8.0,
) -> dict[str, dict[str, Any]]:
    """Fetch compact target metadata for ChEMBL target IDs.

    Metadata lookup is deliberately separate from model inference. A temporary
    ChEMBL metadata failure should not corrupt model scores or target ordering.
    """
    ordered_ids = [str(target_id) for target_id in target_ids if target_id]
    missing_ids = [
        target_id
        for target_id in dict.fromkeys(ordered_ids)
        if target_id not in _CHEMBL_TARGET_METADATA_CACHE
    ]

    chunk_size = max(1, int(chunk_size or 100))
    for start in range(0, len(missing_ids), chunk_size):
        chunk = missing_ids[start : start + chunk_size]
        query = urlencode(
            {"target_chembl_id__in": ",".join(chunk), "limit": len(chunk)},
            safe=",",
        )
        response = _http_get(
            f"{_CHEMBL_TARGET_URL}?{query}",
            timeout=timeout,
            retries=1,
        )
        if response is None or not response.ok:
            continue
        try:
            payload = response.json()
        except ValueError:
            continue

        for record in payload.get("targets", []):
            target_id = record.get("target_chembl_id")
            if not target_id:
                continue
            _CHEMBL_TARGET_METADATA_CACHE[str(target_id)] = {
                "target_chembl_id": str(target_id),
                "pref_name": record.get("pref_name") or str(target_id),
                "organism": record.get("organism") or "Unknown",
                "target_type": record.get("target_type") or "Unknown",
            }

    return {
        target_id: _CHEMBL_TARGET_METADATA_CACHE.get(
            target_id, _default_target_metadata(target_id)
        )
        for target_id in ordered_ids
    }


def _flatten_chembl_multitask_outputs(outputs: Sequence[Any]) -> np.ndarray:
    if not outputs:
        return np.array([], dtype=np.float32)

    if len(outputs) == 1:
        return np.asarray(outputs[0], dtype=np.float32).reshape(-1)

    values = []
    for output in outputs:
        arr = np.asarray(output, dtype=np.float32).reshape(-1)
        if arr.size:
            values.append(float(arr[0]))
    return np.asarray(values, dtype=np.float32)


def build_pd_target_prediction_df(
    scores: Sequence[float],
    target_ids: Sequence[str],
    target_metadata: dict[str, dict[str, Any]] | None = None,
    top_n: int = 25,
    human_only: bool = True,
    min_score: float | None = None,
    model_label: str = _CHEMBL_MULTITASK_MODEL_LABEL,
) -> pd.DataFrame:
    """Format multitask target scores as a pharmacodynamic hypothesis table."""
    scores_array = np.asarray(scores, dtype=float).reshape(-1)
    targets = [str(target_id) for target_id in target_ids]
    if scores_array.size != len(targets):
        raise ValueError(
            "Target prediction score count does not match target ID count "
            f"({scores_array.size} scores vs {len(targets)} target IDs)."
        )
    if scores_array.size == 0 or top_n <= 0:
        return pd.DataFrame(columns=_PD_TARGET_COLUMNS)

    metadata = target_metadata or {}
    records = []
    for index in np.argsort(-scores_array):
        target_id = targets[int(index)]
        score = float(scores_array[int(index)])
        if min_score is not None and score < min_score:
            continue

        target_info = metadata.get(target_id, _default_target_metadata(target_id))
        organism = target_info.get("organism") or "Unknown"
        if human_only and organism not in {"Homo sapiens", "Unknown"}:
            continue

        limitation = (
            "Ligand-based target/off-target hypothesis; not activity direction, "
            "efficacy, tissue response, clinical effect, or safety."
        )
        if organism == "Unknown":
            limitation = (
                f"{limitation} Target metadata unavailable; organism is not "
                "confirmed."
            )

        records.append(
            {
                "Rank": len(records) + 1,
                "Target ChEMBL ID": target_id,
                "Target": target_info.get("pref_name") or target_id,
                "Organism": organism,
                "Target Type": target_info.get("target_type") or "Unknown",
                "Prediction Score": round(score, 4),
                "Model": model_label,
                "Evidence / limitation": limitation,
            }
        )
        if len(records) >= int(top_n):
            break

    return pd.DataFrame.from_records(records, columns=_PD_TARGET_COLUMNS)


def predict_pd_targets_from_smiles(
    smiles: str,
    top_n: int = 25,
    human_only: bool = True,
    min_score: float | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
    download: bool = True,
) -> pd.DataFrame:
    """Predict likely ChEMBL targets/off-targets for a molecule SMILES string."""
    if not smiles or not str(smiles).strip():
        raise ValueError("A SMILES string is required for target prediction.")

    paths = ensure_chembl_multitask_assets(cache_dir=cache_dir, download=download)
    target_ids = load_chembl_multitask_targets(cache_dir=cache_dir, download=False)
    session = _chembl_multitask_session(paths["model"])

    model_inputs = session.get_inputs()
    if not model_inputs:
        raise RuntimeError("ChEMBL multitask ONNX model has no input nodes.")
    fingerprint = _chembl_multitask_fingerprint(smiles)
    outputs = session.run(None, {model_inputs[0].name: fingerprint})
    scores = _flatten_chembl_multitask_outputs(outputs)
    output_names = [output.name for output in session.get_outputs()]

    if scores.size == len(output_names) and all(
        name.startswith("CHEMBL") for name in output_names
    ):
        target_ids = output_names
    elif scores.size != len(target_ids):
        raise RuntimeError(
            "ChEMBL multitask ONNX output count does not match the target "
            f"list ({scores.size} scores vs {len(target_ids)} target IDs)."
        )

    candidate_count = min(len(target_ids), max(int(top_n) * 4, int(top_n) + 25))
    candidate_indices = np.argsort(-scores)[:candidate_count]
    candidate_targets = [target_ids[int(index)] for index in candidate_indices]
    target_metadata = fetch_chembl_target_metadata(candidate_targets)
    return build_pd_target_prediction_df(
        scores,
        target_ids,
        target_metadata=target_metadata,
        top_n=top_n,
        human_only=human_only,
        min_score=min_score,
    )


def normalize_pk_descriptors(
    descriptors: dict[str, Any],
) -> dict[str, float | str | None]:
    """Map PubChem/RDKit descriptor names to a stable internal schema."""
    normalized: dict[str, float | str | None] = {
        "mw": _as_float(
            _first_present(descriptors, *_PK_DESCRIPTOR_ALIASES["mw"]), 0.0
        ),
        "logp": _as_float(
            _first_present(descriptors, *_PK_DESCRIPTOR_ALIASES["logp"]), 0.0
        ),
        "tpsa": _as_float(
            _first_present(descriptors, *_PK_DESCRIPTOR_ALIASES["tpsa"]), 0.0
        ),
        "hbd": _as_float(
            _first_present(descriptors, *_PK_DESCRIPTOR_ALIASES["hbd"]), 0.0
        ),
        "hba": _as_float(
            _first_present(descriptors, *_PK_DESCRIPTOR_ALIASES["hba"]), 0.0
        ),
        "rotb": _as_float(
            _first_present(descriptors, *_PK_DESCRIPTOR_ALIASES["rotb"]), 0.0
        ),
        "charge": _as_float(
            _first_present(descriptors, *_PK_DESCRIPTOR_ALIASES["charge"]), 0.0
        ),
        "smiles": _first_present(
            descriptors, "IsomericSMILES", "CanonicalSMILES", "SMILES"
        ),
        "source": _first_present(descriptors, "DescriptorSource", "Source")
        or "Provided descriptors",
    }

    if not has_complete_pk_descriptors(descriptors) and normalized["smiles"]:
        rdkit_descriptors = descriptors_from_smiles(str(normalized["smiles"]))
        if rdkit_descriptors:
            rdkit_values = {
                "mw": rdkit_descriptors.get("MolecularWeight"),
                "logp": rdkit_descriptors.get("XLogP"),
                "tpsa": rdkit_descriptors.get("TPSA"),
                "hbd": rdkit_descriptors.get("HBondDonorCount"),
                "hba": rdkit_descriptors.get("HBondAcceptorCount"),
                "rotb": rdkit_descriptors.get("RotatableBondCount"),
                "charge": rdkit_descriptors.get("Charge"),
            }
            for key, value in rdkit_values.items():
                if _first_present(descriptors, *_PK_DESCRIPTOR_ALIASES[key]) is None:
                    normalized[key] = _as_float(value, normalized[key])
            normalized["source"] = f"{normalized['source']} + RDKit"

    return normalized


def predict_pk_from_descriptors(descriptors: dict[str, Any]) -> dict[str, Any]:
    """Return numeric and categorical PK estimates from molecular descriptors.

    The formulas combine widely used small-molecule ADME cutoffs with empirical
    descriptor relationships. They are intended for ranking and hypothesis
    generation, not individual-patient prediction.
    """
    completed_descriptors = complete_pk_descriptors(descriptors)
    if not has_complete_pk_descriptors(completed_descriptors):
        missing = ", ".join(missing_pk_descriptor_names(completed_descriptors))
        raise ValueError(
            "Cannot predict PK from incomplete descriptors; missing required "
            f"descriptor(s): {missing}."
        )

    d = normalize_pk_descriptors(completed_descriptors)
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
    oral_bioavailability_pct = _clip(
        fraction_absorbed_pct * first_pass_factor, 1.0, 95.0
    )

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

    hepatic_score = _sigmoid(
        -2.0 + 0.45 * logp + 0.006 * (mw - 300.0) + 0.08 * rotb - 0.008 * tpsa
    )
    hepatic_extraction = 0.02 + 0.55 * hepatic_score
    hepatic_cl_ml_min_kg = 20.7 * hepatic_extraction
    renal_score = _sigmoid(
        0.8 - 0.75 * logp + 0.012 * (tpsa - 70.0) - 0.006 * (mw - 300.0)
    )
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
    applicability = (
        "; ".join(applicability_flags)
        if applicability_flags
        else "Within broad small-molecule descriptor range"
    )

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
    oral_f = _clip(
        float(prediction.get("oral_bioavailability_pct", 0.0)) / 100.0, 0.01, 0.95
    )
    tpsa = float(prediction.get("tpsa", 0.0) or 0.0)
    rotb = float(prediction.get("rotb", 0.0) or 0.0)
    mw = float(prediction.get("mw", 0.0) or 0.0)
    charge = float(prediction.get("charge", 0.0) or 0.0)
    permeability = _clip(
        1.15 - max(0.0, tpsa - 75.0) / 220.0 - max(0.0, rotb - 6.0) / 18.0, 0.35, 1.35
    )

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
        bioavailability = (
            0.92 - 0.18 * max(0.0, mw - 800.0) / 800.0 - 0.06 * (abs(charge) >= 2)
        )
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


def _linear_crossing_time(
    t0: float, c0: float, t1: float, c1: float, threshold: float
) -> float:
    if c1 == c0:
        return t1
    return t0 + (threshold - c0) * (t1 - t0) / (c1 - c0)


def _threshold_times(
    times: list[float], concentrations: list[float], threshold: float | None
) -> dict[str, Any]:
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
                        times[idx - 1],
                        concentrations[idx - 1],
                        times[idx],
                        concentration,
                        threshold,
                    )
        elif first_h is not None and idx > 0 and concentrations[idx - 1] >= threshold:
            fall_h = _linear_crossing_time(
                times[idx - 1],
                concentrations[idx - 1],
                times[idx],
                concentration,
                threshold,
            )
            break

    duration_h = (
        None if first_h is None or fall_h is None else max(0.0, fall_h - first_h)
    )
    return {"first_h": first_h, "last_h": fall_h, "duration_h": duration_h}


def _therapeutic_window_duration(
    times: list[float],
    concentrations: list[float],
    mec_mg_l: float | None,
    mtc_mg_l: float | None,
) -> float | None:
    if mec_mg_l is None or mec_mg_l <= 0 or mtc_mg_l is None or mtc_mg_l <= mec_mg_l:
        return None
    total = 0.0
    for idx in range(1, len(times)):
        c_mid = (concentrations[idx - 1] + concentrations[idx]) / 2.0
        if mec_mg_l <= c_mid < mtc_mg_l:
            total += times[idx] - times[idx - 1]
    return total


def _post_peak_half_time(
    times: list[float], concentrations: list[float]
) -> float | None:
    """Time after Cmax for concentration to fall to half of Cmax."""
    if not times or not concentrations:
        return None
    cmax = max(concentrations)
    if cmax <= 0:
        return None
    peak_idx = concentrations.index(cmax)
    threshold = cmax / 2.0
    for idx in range(peak_idx + 1, len(times)):
        if concentrations[idx] <= threshold <= concentrations[idx - 1]:
            crossing = _linear_crossing_time(
                times[idx - 1],
                concentrations[idx - 1],
                times[idx],
                concentrations[idx],
                threshold,
            )
            return max(0.0, crossing - times[peak_idx])
    return None


def estimate_threshold_from_reference_dose(
    prediction: dict[str, Any],
    reference_dose_amount: Any,
    reference_dose_unit: str,
    reference_route: str,
    body_weight_kg: float = 70.0,
    cmax_fraction: float = 0.5,
    duration_h: float = 24.0,
) -> tuple[float | None, dict[str, Any] | None]:
    """Estimate an exposure threshold from a reference active or toxic dose.

    If a dose is known to be minimally active or minimally toxic, the model uses
    a fraction of that route's predicted Cmax as the concentration threshold.
    The default 50% of Cmax is deliberately conservative and user-adjustable in
    the UI.
    """
    dose_mg = dose_to_mg(reference_dose_amount, reference_dose_unit)
    if dose_mg <= 0:
        return None, None
    _, summary = simulate_pk_for_route(
        prediction,
        dose_mg=dose_mg,
        route=reference_route,
        body_weight_kg=body_weight_kg,
        duration_h=duration_h,
    )
    cmax = summary.get("Cmax (mg/L)")
    if cmax is None or cmax <= 0:
        return None, summary
    fraction = _clip(float(cmax_fraction or 0.5), 0.05, 1.0)
    return cmax * fraction, summary


def estimate_exposure_thresholds(
    prediction: dict[str, Any],
    active_dose_amount: Any = 0,
    toxic_dose_amount: Any = 0,
    reference_dose_unit: str = "milligram",
    reference_route: str = "Oral",
    body_weight_kg: float = 70.0,
    cmax_fraction: float = 0.5,
    duration_h: float = 24.0,
) -> tuple[float | None, float | None, pd.DataFrame]:
    """Estimate MEC/MTC from reference active/toxic doses when supplied."""
    rows = []
    mec, mec_summary = estimate_threshold_from_reference_dose(
        prediction,
        active_dose_amount,
        reference_dose_unit,
        reference_route,
        body_weight_kg=body_weight_kg,
        cmax_fraction=cmax_fraction,
        duration_h=duration_h,
    )
    if mec is not None and mec_summary is not None:
        rows.append(
            {
                "Threshold": "Estimated MEC",
                "Estimate (mg/L)": mec,
                "Reference route": mec_summary["Route"],
                "Reference dose (mg)": mec_summary["Dose (mg)"],
                "Reference Cmax (mg/L)": mec_summary["Cmax (mg/L)"],
                "Cmax fraction": _clip(float(cmax_fraction or 0.5), 0.05, 1.0),
                "Assumption": "Minimum active dose reaches a concentration near the effect threshold.",
            }
        )

    mtc, mtc_summary = estimate_threshold_from_reference_dose(
        prediction,
        toxic_dose_amount,
        reference_dose_unit,
        reference_route,
        body_weight_kg=body_weight_kg,
        cmax_fraction=cmax_fraction,
        duration_h=duration_h,
    )
    if mtc is not None and mtc_summary is not None:
        rows.append(
            {
                "Threshold": "Estimated MTC",
                "Estimate (mg/L)": mtc,
                "Reference route": mtc_summary["Route"],
                "Reference dose (mg)": mtc_summary["Dose (mg)"],
                "Reference Cmax (mg/L)": mtc_summary["Cmax (mg/L)"],
                "Cmax fraction": _clip(float(cmax_fraction or 0.5), 0.05, 1.0),
                "Assumption": "Minimum toxic dose reaches a concentration near the toxicity threshold.",
            }
        )

    return mec, mtc, pd.DataFrame(rows)


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
    clearance_ml_min_kg = max(
        0.01, float(prediction.get("clearance_ml_min_kg", 5.0) or 5.0)
    )
    vd_l = vd_l_kg * body_weight_kg
    clearance_l_h = clearance_ml_min_kg * 0.06 * body_weight_kg
    ke_h = clearance_l_h / vd_l
    ka_h = route_params["ka_h"]
    half_life_h = math.log(2) / ke_h if ke_h > 0 else None
    absorption_half_life_h = math.log(2) / ka_h if ka_h else None
    if ka_h is None:
        apparent_terminal_half_life_h = half_life_h
        terminal_phase_driver = "Elimination"
    elif ka_h < ke_h:
        apparent_terminal_half_life_h = math.log(2) / ka_h
        terminal_phase_driver = "Absorption-limited (flip-flop)"
    else:
        apparent_terminal_half_life_h = half_life_h
        terminal_phase_driver = "Elimination"
    absorbed_dose_mg = route_params["F"] * dose_mg

    steps = int(math.ceil(duration_h / time_step_h))
    times = [round(idx * time_step_h, 6) for idx in range(steps + 1)]
    concentrations: list[float] = []
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
                concentration = (
                    absorbed_dose_mg / vd_l * ka_h * t_eff * math.exp(-ke_h * t_eff)
                )
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
    post_peak_half_time_h = _post_peak_half_time(times, concentrations)
    auc_0_inf_mg_h_l = absorbed_dose_mg / clearance_l_h if clearance_l_h > 0 else None
    auc_last = sum(
        (concentrations[idx - 1] + concentrations[idx])
        / 2.0
        * (times[idx] - times[idx - 1])
        for idx in range(1, len(times))
    )
    mec_times = _threshold_times(times, concentrations, mec_mg_l)
    mtc_times = _threshold_times(times, concentrations, mtc_mg_l)
    therapeutic_duration_h = _therapeutic_window_duration(
        times, concentrations, mec_mg_l, mtc_mg_l
    )
    volume_to_administer_ml = None
    concentration = _as_float(injection_concentration_mg_ml)
    if (
        concentration
        and concentration > 0
        and route_params["route"] in {"Intravenous", "Subcutaneous"}
    ):
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
        "Elimination t1/2 (h)": half_life_h,
        "Absorption t1/2 (h)": absorption_half_life_h,
        "Apparent terminal t1/2 (h)": apparent_terminal_half_life_h,
        "Terminal phase driver": terminal_phase_driver,
        "Post-peak 50% decline (h)": post_peak_half_time_h,
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


def active_metabolite_model_for_compound(compound_name: str) -> list[dict[str, Any]]:
    """Return curated metabolite model parameters when available."""
    note = ACTIVE_METABOLITE_NOTES.get((compound_name or "").strip().casefold(), {})
    return list(note.get("Metabolite model", []))


def simulate_active_metabolite_curves(
    compound_name: str,
    parent_curve_df: pd.DataFrame,
    parent_summary_df: pd.DataFrame,
    parent_prediction: dict[str, Any],
    body_weight_kg: float = 70.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Approximate active-metabolite and potency-weighted active-moiety curves.

    This is intentionally restricted to curated examples. It is a rough
    formation-limited model, not a PBPK or enzyme-genotype model.
    """
    models = active_metabolite_model_for_compound(compound_name)
    if (
        not models
        or parent_curve_df is None
        or parent_curve_df.empty
        or parent_summary_df is None
        or parent_summary_df.empty
    ):
        return pd.DataFrame(), pd.DataFrame()

    body_weight_kg = max(1.0, float(body_weight_kg or 70.0))
    parent_vd_l_kg = max(0.01, float(parent_prediction.get("vd_l_kg", 0.7) or 0.7))
    rows = []
    summary_rows = []

    for route, route_curve in parent_curve_df.groupby("Route"):
        route_curve = route_curve.sort_values("Time (h)")
        times = list(route_curve["Time (h)"].astype(float))
        parent_concs = list(route_curve["Concentration (mg/L)"].astype(float))
        summary_match = parent_summary_df[parent_summary_df["Route"] == route]
        if summary_match.empty:
            continue
        route_summary = summary_match.iloc[0]
        parent_cl_l_h = max(
            0.001, float(route_summary["CL (mL/min/kg)"]) * 0.06 * body_weight_kg
        )

        for idx, time_h in enumerate(times):
            rows.append(
                {
                    "Route": route,
                    "Time (h)": time_h,
                    "Curve": "Parent active-moiety index",
                    "Relative active-moiety index": parent_concs[idx],
                    "Basis": "Parent concentration, weight 1.0",
                }
            )

        total_by_time = {time_h: parent_concs[idx] for idx, time_h in enumerate(times)}
        for model in models:
            amount_mg = 0.0
            metabolite_concs = []
            half_life_h = max(0.05, float(model.get("half_life_h", 6.0) or 6.0))
            ke_m = math.log(2) / half_life_h
            vd_m_l = (
                parent_vd_l_kg
                * float(model.get("vd_multiplier", 1.0) or 1.0)
                * body_weight_kg
            )
            formation_fraction = _clip(
                float(model.get("formation_fraction", 0.1) or 0.1), 0.0, 1.0
            )
            potency = max(
                0.0, float(model.get("potency_relative_to_parent", 1.0) or 1.0)
            )
            previous_time = times[0] if times else 0.0

            for idx, time_h in enumerate(times):
                dt = 0.0 if idx == 0 else max(0.0, time_h - previous_time)
                input_rate_mg_h = formation_fraction * parent_cl_l_h * parent_concs[idx]
                amount_mg = max(
                    0.0, amount_mg + (input_rate_mg_h - ke_m * amount_mg) * dt
                )
                concentration = amount_mg / vd_m_l if vd_m_l > 0 else 0.0
                metabolite_concs.append(concentration)
                effect_index = concentration * potency
                total_by_time[time_h] = total_by_time.get(time_h, 0.0) + effect_index
                rows.append(
                    {
                        "Route": route,
                        "Time (h)": time_h,
                        "Curve": f"{model['name']} active-moiety index",
                        "Relative active-moiety index": effect_index,
                        "Basis": model.get(
                            "basis", "Curated active-metabolite approximation"
                        ),
                    }
                )
                previous_time = time_h

            peak_index = max((c * potency for c in metabolite_concs), default=0.0)
            peak_time = (
                times[[c * potency for c in metabolite_concs].index(peak_index)]
                if metabolite_concs
                else None
            )
            summary_rows.append(
                {
                    "Route": route,
                    "Metabolite": model["name"],
                    "Formation fraction": formation_fraction,
                    "Metabolite half-life (h)": half_life_h,
                    "Relative potency vs parent": potency,
                    "Peak active-moiety index": peak_index,
                    "Peak time (h)": peak_time,
                    "Basis": model.get(
                        "basis", "Curated active-metabolite approximation"
                    ),
                }
            )

        for time_h, total_index in total_by_time.items():
            rows.append(
                {
                    "Route": route,
                    "Time (h)": time_h,
                    "Curve": "Total modeled active-moiety index",
                    "Relative active-moiety index": total_index,
                    "Basis": "Parent plus curated active-metabolite potency-weighted indices",
                }
            )

    return pd.DataFrame(rows), pd.DataFrame(summary_rows)


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
    display_note = {
        key: value for key, value in note.items() if key != "Metabolite model"
    }
    if "Metabolite model" in note:
        display_note["Modeled metabolite curve"] = ", ".join(
            model["name"] for model in note["Metabolite model"]
        )
    return pd.DataFrame([{"Compound": compound_name or "Compound", **display_note}])


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


def pk_formula_items() -> list[tuple[str, str]]:
    """Formula labels and LaTeX equations for rendered notebook display."""
    return [
        ("Elimination rate", r"k_e = \frac{CL}{V_d}"),
        ("Half-life", r"t_{1/2} = \frac{\ln(2)}{k_e}"),
        (
            "IV bolus concentration",
            r"C(t) = \frac{F \cdot D}{V_d} e^{-k_e t},\quad F=1",
        ),
        (
            "First-order absorption concentration",
            r"C(t) = \frac{F \cdot D \cdot k_a}{V_d(k_a-k_e)}\left(e^{-k_e(t-t_{lag})} - e^{-k_a(t-t_{lag})}\right)",
        ),
        ("Exposure", r"AUC_{0-\infty} = \frac{F \cdot D}{CL}"),
        ("Injection volume", r"V_{admin} = \frac{D_{prescribed}}{C_{solution}}"),
        (
            "Reference-dose threshold calibration",
            r"MEC \approx f_{Cmax} \cdot C_{max,ref\ active}",
        ),
        (
            "Reference toxic threshold calibration",
            r"MTC \approx f_{Cmax} \cdot C_{max,ref\ toxic}",
        ),
    ]


def pk_estimate_formula_items() -> list[tuple[str, str]]:
    """Formula labels for descriptor-derived estimates."""
    return [
        (
            "Absorption descriptor penalty",
            r"P_{abs}=P_{TPSA}+P_{rotB}+P_{MW}+P_{logP}+P_{HBD}+0.35\cdot RO5_{viol}",
        ),
        ("Predicted fraction absorbed", r"F_a = 100\cdot \sigma(2.2-P_{abs})"),
        ("Oral bioavailability potential", r"F_{oral}=F_a\cdot FP_{factor}"),
        (
            "Clark-style brain:blood estimate",
            r"logBB=0.152\cdot logP-0.0148\cdot TPSA+0.139",
        ),
        (
            "Volume of distribution heuristic",
            r"log_{10}(V_d)=0.16(logP-1.5)-0.0035(TPSA-75)-0.0008(MW-350)+q_{adj}",
        ),
        ("Hepatic clearance proxy", r"CL_H=20.7\cdot(0.02+0.55\cdot \sigma(S_H))"),
        ("Renal clearance proxy", r"CL_R=1.0\cdot \sigma(S_R)"),
        ("Total clearance proxy", r"CL=CL_H+CL_R"),
        ("Elimination half-life", r"t_{1/2,elim}=\frac{\ln(2)\cdot V_d}{CL}"),
        ("Absorption half-life", r"t_{1/2,abs}=\frac{\ln(2)}{k_a}"),
        ("Apparent terminal half-life", r"t_{1/2,term}=\frac{\ln(2)}{\min(k_e,k_a)}"),
        (
            "Post-peak half-time",
            r"t_{50,postpeak}=t(C=0.5\cdot C_{max},\ t>T_{max})-T_{max}",
        ),
    ]


def pk_profile_from_descriptors(
    descriptors: dict[str, Any], label: str = "Compound"
) -> pd.DataFrame:
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

    completed_descriptors = complete_pk_descriptors(
        descriptors,
        smiles=_first_present(
            descriptors, "IsomericSMILES", "CanonicalSMILES", "SMILES"
        ),
    )
    if not has_complete_pk_descriptors(completed_descriptors):
        missing = ", ".join(missing_pk_descriptor_names(completed_descriptors))
        return pd.DataFrame(
            [
                {
                    "Parameter": "Predicted pharmacokinetics unavailable",
                    "Estimate": "",
                    "Units": "",
                    "Interpretation": f"Incomplete structure descriptors for {label}.",
                    "Evidence / formula": f"Missing required descriptor(s): {missing}.",
                    "Confidence": "Unavailable",
                }
            ]
        )

    p = predict_pk_from_descriptors(completed_descriptors)
    rows = [
        (
            "Descriptor source",
            p["source"],
            "",
            label,
            "Descriptors are read from PubChem or calculated from SMILES with RDKit.",
            "Input",
        ),
        (
            "Molecular weight",
            _format_number(p["mw"], 2),
            "g/mol",
            "",
            "Lipinski RO5 uses MW <= 500.",
            "Input",
        ),
        (
            "XLogP / MolLogP",
            _format_number(p["logp"], 2),
            "log10",
            "",
            "Lipinski RO5 uses logP <= 5.",
            "Input",
        ),
        (
            "Topological polar surface area",
            _format_number(p["tpsa"], 1),
            "A^2",
            "",
            "Veber oral exposure threshold commonly uses TPSA <= 140 A^2.",
            "Input",
        ),
        (
            "H-bond donors",
            _format_number(p["hbd"], 0),
            "count",
            "",
            "Lipinski RO5 uses HBD <= 5.",
            "Input",
        ),
        (
            "H-bond acceptors",
            _format_number(p["hba"], 0),
            "count",
            "",
            "Lipinski RO5 uses HBA <= 10.",
            "Input",
        ),
        (
            "Rotatable bonds",
            _format_number(p["rotb"], 0),
            "count",
            "",
            "Veber oral exposure threshold commonly uses rotatable bonds <= 10.",
            "Input",
        ),
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
        columns=[
            "Parameter",
            "Estimate",
            "Units",
            "Interpretation",
            "Evidence / formula",
            "Confidence",
        ],
    )


def pk_profile_for_compound(
    compound: int | str | None, smiles: str | None = None
) -> pd.DataFrame:
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


def report_row_limit(value: Any) -> int | None:
    """Normalize report row-limit widget values."""
    if value is None:
        return 50
    text = str(value).strip().casefold()
    if not text:
        return 50
    if text in {"all", "none"}:
        return None
    try:
        limit = int(float(text))
    except (TypeError, ValueError):
        return 50
    return None if limit <= 0 else limit


def limited_report_dataframe(
    df: pd.DataFrame, row_limit: int | None = 50
) -> pd.DataFrame:
    """Return a display-safe dataframe limited to the requested number of rows."""
    if df is None or df.empty:
        return pd.DataFrame()
    out = df.copy()
    if row_limit is not None:
        out = out.head(row_limit)
    return out.fillna("")


def _section_dataframe(section: dict[str, Any], row_limit: int | None) -> pd.DataFrame:
    return limited_report_dataframe(section.get("df"), row_limit=row_limit)


def _section_images(section: dict[str, Any]) -> list[dict[str, Any]]:
    images = section.get("images") or []
    if section.get("image"):
        images = [section["image"], *images]
    normalized = []
    for image in images:
        data = image.get("data") if isinstance(image, dict) else None
        if isinstance(data, str):
            data = data.encode("utf-8")
        if not data:
            continue
        normalized.append(
            {
                "title": str(image.get("title") or "Figure")
                if isinstance(image, dict)
                else "Figure",
                "data": data,
                "mime": str(image.get("mime") or "image/png")
                if isinstance(image, dict)
                else "image/png",
            }
        )
    return normalized


def build_interactive_html_report(
    title: str, sections: list[dict[str, Any]], row_limit: Any = 50
) -> str:
    """Build a self-contained HTML report with collapsible table sections."""
    limit = report_row_limit(row_limit)
    safe_title = html.escape(title or "ChEMBL Bioactivity Report")
    body = [
        "<!doctype html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        f"<title>{safe_title}</title>",
        "<style>",
        "body{font-family:system-ui,-apple-system,Segoe UI,sans-serif;margin:2rem;color:#17202a;line-height:1.45}",
        "h1{margin-bottom:.2rem} .meta{color:#5d6d7e;margin-bottom:1.5rem}",
        "details{border:1px solid #d5d8dc;border-radius:8px;margin:1rem 0;padding:.75rem;background:#fbfcfc}",
        "summary{font-weight:700;cursor:pointer;font-size:1.05rem}",
        "table{border-collapse:collapse;width:100%;font-size:.88rem;margin-top:.75rem}",
        "th,td{border:1px solid #d6dbdf;padding:.35rem .45rem;vertical-align:top}",
        "th{background:#eaf2f8;text-align:left;position:sticky;top:0}",
        "tr:nth-child(even){background:#f8f9f9}",
        ".report-image{margin:.8rem 0 1rem}",
        ".report-image img{max-width:100%;height:auto;border:1px solid #d5d8dc;border-radius:6px;background:white}",
        ".note{color:#566573;margin:.4rem 0}",
        ".small{font-size:.82rem;color:#5d6d7e}",
        "@media print{details{break-inside:avoid} summary{cursor:default} body{margin:1cm}}",
        "</style>",
        "</head><body>",
        f"<h1>{safe_title}</h1>",
        '<div class="meta">Generated from the ChEMBL Bioactivity and Predicted Pharmacokinetics app. Use browser print/save-as-PDF for a printable copy.</div>',
    ]

    if not sections:
        body.append("<p>No report sections were available.</p>")

    for section in sections:
        df = _section_dataframe(section, limit)
        images = _section_images(section)
        if df.empty and not images:
            continue
        section_title = html.escape(str(section.get("title") or "Table"))
        original_rows = (
            len(section.get("df")) if section.get("df") is not None else len(df)
        )
        shown_rows = len(df)
        note = section.get("note")
        open_attr = " open" if section.get("open", True) else ""
        table_meta = (
            f' <span class="small">({shown_rows} of {original_rows} rows)</span>'
            if not df.empty
            else ""
        )
        body.append(f"<details{open_attr}>")
        body.append(f"<summary>{section_title}{table_meta}</summary>")
        if note:
            body.append(f'<div class="note">{html.escape(str(note))}</div>')
        for image in images:
            b64 = base64.b64encode(image["data"]).decode("ascii")
            image_title = html.escape(image["title"])
            mime = html.escape(image["mime"])
            body.append(
                f'<figure class="report-image"><img src="data:{mime};base64,{b64}" alt="{image_title}">'
                f'<figcaption class="small">{image_title}</figcaption></figure>'
            )
        if not df.empty:
            body.append(df.to_html(index=False, escape=True, border=0))
        body.append("</details>")

    body.append(
        '<p class="small"><strong>Disclaimer:</strong> Pharmacokinetic outputs are approximate screening estimates, not dosing or clinical safety guidance.</p>'
    )
    body.append("</body></html>")
    return "\n".join(body)


def _pdf_cell(value: Any):
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import Paragraph

    styles = getSampleStyleSheet()
    style = styles["BodyText"]
    style.fontSize = 6
    style.leading = 7
    if value is None or (isinstance(value, float) and math.isnan(value)):
        text = ""
    else:
        text = html.escape(str(value))
    return Paragraph(text, style)


def build_pdf_report(
    title: str, sections: list[dict[str, Any]], row_limit: Any = 50
) -> bytes:
    """Build a PDF report of current tables. Requires reportlab."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.lib.utils import ImageReader
    from reportlab.platypus import Image as PdfImage
    from reportlab.platypus import (
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    limit = report_row_limit(row_limit)
    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf,
        pagesize=landscape(A4),
        leftMargin=0.8 * cm,
        rightMargin=0.8 * cm,
        topMargin=0.8 * cm,
        bottomMargin=0.8 * cm,
    )
    styles = getSampleStyleSheet()
    story = [
        Paragraph(html.escape(title or "ChEMBL Bioactivity Report"), styles["Title"])
    ]
    story.append(
        Paragraph(
            "Approximate screening report generated by the ChEMBL Bioactivity and Predicted Pharmacokinetics app.",
            styles["BodyText"],
        )
    )
    story.append(Spacer(1, 0.25 * cm))

    for section in sections:
        df = _section_dataframe(section, limit)
        images = _section_images(section)
        if df.empty and not images:
            continue
        section_title = str(section.get("title") or "Table")
        original_rows = (
            len(section.get("df")) if section.get("df") is not None else len(df)
        )
        table_meta = f" ({len(df)} of {original_rows} rows)" if not df.empty else ""
        story.append(
            Paragraph(html.escape(f"{section_title}{table_meta}"), styles["Heading2"])
        )
        if section.get("note"):
            story.append(
                Paragraph(html.escape(str(section["note"])), styles["BodyText"])
            )

        for image in images:
            try:
                image_buf = io.BytesIO(image["data"])
                width_px, height_px = ImageReader(image_buf).getSize()
                scale = min(doc.width / width_px, (10.5 * cm) / height_px, 1.0)
                image_buf.seek(0)
                story.append(
                    PdfImage(
                        image_buf, width=width_px * scale, height=height_px * scale
                    )
                )
                story.append(Paragraph(html.escape(image["title"]), styles["BodyText"]))
                story.append(Spacer(1, 0.2 * cm))
            except Exception:
                continue

        if not df.empty:
            columns = list(df.columns)
            data = [[_pdf_cell(col) for col in columns]]
            for _, row in df.iterrows():
                data.append([_pdf_cell(row.get(col, "")) for col in columns])
            col_count = max(1, len(columns))
            table = Table(
                data, repeatRows=1, colWidths=[doc.width / col_count] * col_count
            )
            table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D6EAF8")),
                        ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#AAB7B8")),
                        ("VALIGN", (0, 0), (-1, -1), "TOP"),
                        (
                            "ROWBACKGROUNDS",
                            (0, 1),
                            (-1, -1),
                            [colors.white, colors.HexColor("#F8F9F9")],
                        ),
                        ("LEFTPADDING", (0, 0), (-1, -1), 2),
                        ("RIGHTPADDING", (0, 0), (-1, -1), 2),
                        ("TOPPADDING", (0, 0), (-1, -1), 2),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
                    ]
                )
            )
            story.append(table)
        story.append(Spacer(1, 0.35 * cm))

    story.append(
        Paragraph(
            "Disclaimer: PK estimates are approximate screening values and must not be used for prescribing, self-dosing, or clinical safety decisions.",
            styles["BodyText"],
        )
    )
    doc.build(story)
    return buf.getvalue()
