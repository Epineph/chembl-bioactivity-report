# ===========================
# Neat ChEMBL Bioactivity Report (Integrated, Binder-hardened)
# ===========================
# Features
# - Pharmacodynamics from ChEMBL (human targets)
# - Search/sort/filter table (itables/DataTables)
# - Download CSV/Excel (separator selectable)
# - 2D structure (RDKit if available; else PubChem PNG)
# - 3D interactive structure (py3Dmol; requires PubChem 3D)
# - PubChem properties (robust walker for Experimental/Computed)
# - Predicted pharmacokinetics from transparent structure-based heuristics
# - Binder/Voila-hardened networking (UA + retries + TXT fallback)
#
# Notes
# * RDKit via conda-forge is recommended (see environment.yml).
# * 3D model depends on PubChem having a 3D conformer.
# * Experimental properties are heterogeneous; availability varies.
# * This cell only builds the UI; no network calls occur until you press "Search".

from chembl_bioactivity_enhanced import (
    PK_ESTIMATE_NOTE,
    PK_SIMULATION_NOTE,
    active_metabolite_notes_for_compound,
    build_interactive_html_report,
    build_pdf_report,
    descriptors_from_pubchem_cid,
    estimate_exposure_thresholds,
    pk_estimate_formula_items,
    pk_formula_items,
    pk_profile_from_descriptors,
    predict_pk_from_descriptors,
    report_row_limit,
    simulate_active_metabolite_curves,
    simulate_pk_curves,
)
import html, io, base64, json, re, requests, contextlib, warnings, time, unicodedata
from urllib.parse import quote
import pandas as pd
from IPython.display import display, clear_output, Markdown, HTML, Image
import ipywidgets as widgets

# Silence known deprecation warning from chembl_webresource_client
warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated as an API",
    category=UserWarning,
    module=r"chembl_webresource_client.*"
)

from chembl_webresource_client.new_client import new_client
from itables import init_notebook_mode, show, options as itbl_options
itbl_options.warn_on_undocumented_option = False
# Default: 50 rows per page + common sizes incl. “All”
itbl_options.pageLength = 50
itbl_options.lengthMenu = [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]]
init_notebook_mode(all_interactive=True)

# ----------------------------
# Optional imports (quiet, graceful)
# ----------------------------
def _silent_import(module_name, attrs=()):
    try:
        buf = io.StringIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(buf):
                mod = __import__(module_name, fromlist=list(attrs))
                for a in attrs:
                    getattr(mod, a)
        return True, ""
    except Exception as e:
        return False, f"{module_name} disabled ({e.__class__.__name__}: {e})"

def _silent_import_rdkit():
    try:
        buf = io.StringIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(buf):
                from rdkit import Chem  # noqa: F401
                from rdkit.Chem import Draw, AllChem  # noqa: F401
        return True, ""
    except Exception as e:
        return False, f"RDKit disabled ({e.__class__.__name__}: {e})"

_PCP, _pcp_msg = _silent_import("pubchempy")
_P3D, _p3d_msg = _silent_import("py3Dmol")
_MPL, _mpl_msg = _silent_import("matplotlib.pyplot")
_RDKIT, _rdkit_msg = _silent_import_rdkit()

if _RDKIT:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
if _PCP:
    import pubchempy as pcp
if _P3D:
    import py3Dmol
if _MPL:
    import matplotlib.pyplot as plt

# ----------------------------
# Robust PubChem networking for Binder/Voila
# ----------------------------
_UA = "chembl-bioactivity-report/0.2 (+https://github.com/your-org/chembl-bioactivity-report)"
_SESS = requests.Session()
_SESS.headers.update({
    "User-Agent": _UA,
    "Accept": "application/json, text/plain;q=0.5",
})

def _http_get(url: str, timeout: float = 30.0, retries: int = 3, backoff: float = 1.6):
    """GET with tiny exponential backoff (handles 429/5xx/DNS hiccups)."""
    for attempt in range(retries):
        try:
            r = _SESS.get(url, timeout=timeout, allow_redirects=True)
            if r.status_code in (429,) or (500 <= r.status_code < 600):
                time.sleep(backoff ** attempt)
                continue
            return r
        except requests.RequestException:
            time.sleep(backoff ** attempt)
    return None

# ----------------------------
# ChEMBL helpers
# ----------------------------
def get_chembl_id(compound: str) -> str:
    mol_client = new_client.molecule
    res = mol_client.filter(pref_name__iexact=compound)
    if not res:
        raise ValueError(f"No ChEMBL entry for '{compound}'")
    return res[0]['molecule_chembl_id']

def fetch_activities(chembl_id: str) -> list[dict]:
    act_client = new_client.activity
    acts = act_client.filter(
        molecule_chembl_id=chembl_id,
        target_organism__iexact='Homo sapiens'
    ).only([
        'target_chembl_id',
        'standard_type',
        'standard_value',
        'standard_units'
    ])
    return list(acts)

def fetch_target_names(target_ids: set[str]) -> dict[str, str]:
    tgt_client = new_client.target
    names = {}
    for tid in target_ids:
        rec = tgt_client.filter(target_chembl_id=tid).only(['pref_name'])
        names[tid] = rec[0]['pref_name'] if rec else tid
    return names

def build_activity_df(acts: list[dict]) -> pd.DataFrame:
    rows = []
    for a in acts:
        tid   = a.get('target_chembl_id') or 'Unknown'
        typ   = a.get('standard_type')    or ''
        val   = a.get('standard_value')   or ''
        unit  = a.get('standard_units')   or ''
        kd_nm = ''
        if typ and isinstance(typ, str) and typ.upper() == 'KA' and val:
            u = (unit or '').strip()
            if u in {'M^-1', 'M-1', '1/M'}:
                try:
                    kd_m  = 1.0 / float(val)
                    kd_nm = round(kd_m * 1e9, 3)
                except Exception:
                    kd_nm = ''
        rows.append({
            'Target': tid,
            'Activity': typ,
            'Value': val,
            'Units': unit,
            'Kd (nM) (from KA)': kd_nm
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    unique_tids = set(df['Target'])
    name_map = fetch_target_names(unique_tids)
    df['Target'] = df['Target'].map(name_map)

    df = df[['Target', 'Activity', 'Value', 'Units', 'Kd (nM) (from KA)']]
    df = df[df['Value'].astype(str).str.len() > 0].reset_index(drop=True)
    return df

# ----------------------------
# PubChem helpers (ID + basic properties)
# ----------------------------
def _normalize_text(s: str) -> str:
    return unicodedata.normalize("NFKC", (s or "")).strip()

def pubchem_cid_from_name(name: str) -> int | None:
    """Resolve name -> PubChem CID (PubChemPy → PUG JSON → TXT fallback)."""
    nm = _normalize_text(name)
    if not nm:
        return None

    if _PCP:
        try:
            cids = pcp.get_cids(nm, namespace='name')
            if cids:
                return int(cids[0])
        except Exception:
            pass

    # JSON first
    url_json = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(nm)}/cids/JSON"
    r = _http_get(url_json, timeout=35.0)
    if r and r.ok:
        try:
            js = r.json()
            ids = js.get('IdentifierList', {}).get('CID', [])
            if ids:
                return int(ids[0])
        except Exception:
            pass

    # TXT fallback
    url_txt = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(nm)}/cids/TXT"
    r2 = _http_get(url_txt, timeout=35.0)
    if r2 and r2.ok and r2.text:
        try:
            for line in r2.text.splitlines():
                line = line.strip()
                if not line:
                    continue
                m = re.search(r"\d+", line)
                if m:
                    return int(m.group(0))
        except Exception:
            pass

    return None

def pubchem_basic_props_df(cid: int) -> pd.DataFrame:
    """Compact identifiers/descriptors table. Prefers PubChemPy; falls back to PUG-View."""
    if _PCP:
        try:
            props = [
                'IUPACName', 'MolecularFormula', 'MolecularWeight',
                'CanonicalSMILES', 'IsomericSMILES', 'InChIKey',
                'XLogP', 'ExactMass', 'TPSA',
                'HBondDonorCount', 'HBondAcceptorCount',
                'RotatableBondCount', 'FormalCharge'
            ]
            df = pcp.get_properties(props, cid, as_dataframe=True)
            df.insert(0, 'Source', 'PubChem (computed)')
            tidy = df.T.reset_index().rename(columns={'index': 'Property'})
            value_cols = [col for col in tidy.columns if col != 'Property']
            if value_cols:
                tidy = tidy.rename(columns={value_cols[0]: 'Value'})
            return tidy
        except Exception:
            pass
    # PUG-View fallback
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
        r = _http_get(url, timeout=35.0)
        if not r or not r.ok:
            return pd.DataFrame(columns=['Source','Property','Value'])
        data = r.json()
        props = []
        record = data.get('Record', {})
        def pick_strings(sections, wanted):
            for s in sections:
                if 'Information' in s:
                    for info in s['Information']:
                        name = info.get('Name', '')
                        if name in wanted and 'Value' in info and 'StringWithMarkup' in info['Value']:
                            val = ' '.join(v.get('String', '') for v in info['Value']['StringWithMarkup'])
                            props.append((name, val))
                if 'Section' in s:
                    pick_strings(s['Section'], wanted)
        wanted = {'IUPAC Name', 'Molecular Formula', 'InChIKey', 'Canonical SMILES', 'Isomeric SMILES'}
        pick_strings(record.get('Section', []), wanted)
        df = pd.DataFrame(props, columns=['Property', 'Value'])
        if not df.empty:
            df.insert(0, 'Source', 'PubChem')
        return df
    except Exception:
        return pd.DataFrame(columns=['Source','Property','Value'])

# ----------------------------
# PubChem helpers (robust property walker)
# ----------------------------
def _flatten_string_with_markup(obj) -> str:
    out = []
    if isinstance(obj, dict):
        if 'String' in obj:
            out.append(str(obj.get('String', '')))
        if 'StringWithMarkup' in obj:
            for item in obj['StringWithMarkup']:
                out.append(_flatten_string_with_markup(item))
        if 'List' in obj:
            for item in obj['List']:
                out.append(_flatten_string_with_markup(item))
    elif isinstance(obj, list):
        for item in obj:
            out.append(_flatten_string_with_markup(item))
    elif obj is not None:
        out.append(str(obj))
    return " ".join(s for s in out if s).strip()

def _flatten_table(tbl: dict) -> str:
    rows_text = []
    headers = []
    if 'Columns' in tbl and isinstance(tbl['Columns'], list):
        headers = [c.get('Name', '').strip() for c in tbl['Columns']]
    for row in tbl.get('Row', []):
        cells = row.get('Cell', [])
        cell_strs = [_flatten_string_with_markup(c) for c in cells]
        if headers and len(headers) == len(cell_strs):
            row_str = " | ".join(f"{h}: {v}" for h, v in zip(headers, cell_strs))
        else:
            row_str = " | ".join(cell_strs)
        rows_text.append(row_str.strip())
    title = tbl.get('Title') or ""
    body = "; ".join(r for r in rows_text if r)
    return f"{title}: {body}".strip() if title and body else (body or title)

def _extract_value_from_information(info: dict) -> str | None:
    v = info.get('Value')
    if v:
        if 'StringWithMarkup' in v or 'String' in v or 'List' in v:
            text = _flatten_string_with_markup(v)
            if text:
                return text
        if 'Number' in v:
            num = v.get('Number')
            unit = v.get('Unit') or v.get('Units') or ""
            return f"{num} {unit}".strip() if num is not None else None
    tbl = info.get('Table')
    if tbl:
        return _flatten_table(tbl)
    return None

def pubchem_properties_all(cid: int) -> pd.DataFrame:
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
        r = _http_get(url, timeout=40.0)
        if not r or not r.ok:
            return pd.DataFrame(columns=['Property','Value','Source'])
        js = r.json()
    except Exception:
        return pd.DataFrame(columns=['Property','Value','Source'])

    out = []
    def walk(sections, path):
        for s in sections:
            heading = s.get('TOCHeading') or s.get('Name') or ""
            new_path = path + ([heading] if heading else [])
            for info in s.get('Information', []):
                name = (info.get('Name') or heading or "Property").strip()
                val = _extract_value_from_information(info)
                if val:
                    src = " > ".join(p for p in new_path if p)
                    out.append((name, val.strip(), src))
            if 'Table' in s and isinstance(s['Table'], dict):
                t_str = _flatten_table(s['Table'])
                if t_str:
                    src = " > ".join(p for p in new_path if p)
                    title = s['Table'].get('Title') or heading or "Table"
                    out.append((title.strip(), t_str.strip(), src))
            if 'Section' in s:
                walk(s['Section'], new_path)
    record = js.get('Record', {})
    walk(record.get('Section', []), path=[])
    df = pd.DataFrame(out, columns=['Property','Value','Source']).drop_duplicates()
    return df

def pubchem_experimental_props_df(cid: int) -> pd.DataFrame:
    """
    Return prioritized subset of properties; if none found, return the full set.
    Uses non-capturing groups to avoid pandas warnings and sets na=False.
    """
    df_all = pubchem_properties_all(cid)
    if df_all.empty:
        return df_all

    patterns = [
        r'\bmelting point\b',
        r'\bboiling point\b',
        r'\bsolubilit(?:y|ies)\b',
        r'\bpK(?:a|A)\b',
        r'\bpH\b',
        r'\blog\s*P\b',
        r'\bX?logP.*',
        r'\bdensity\b',
        r'\bvapou?r pressure\b',
        r'\bflash point\b',
        r'\bappearance\b',
        r'\bcolor/?form\b',
    ]
    rx = re.compile("|".join(patterns), re.IGNORECASE)

    prop_hit = df_all["Property"].str.contains(rx, na=False)
    src_hit  = df_all    ["Source"].str.contains(rx, na=False)
    df_sub = df_all[prop_hit | src_hit]

    if not df_sub.empty:
        return df_sub

    fallback_patterns = [
        r'\bdescriptor\b', r'\bphysical\b', r'\bchemical\b',
        r'\bpartition\b', r'\bacid dissociation\b'
    ]
    rx2 = re.compile("|".join(fallback_patterns), re.IGNORECASE)
    fb_hit = df_all["Source"].str.contains(rx2, na=False)
    df_fb = df_all[fb_hit]

    return df_fb if not df_fb.empty else df_all

# ----------------------------
# Rendering helpers
# ----------------------------
def rdkit_image_from_smiles(smiles: str, size=(320, 320)):
    if not _RDKIT:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        AllChem.Compute2DCoords(mol)
        return Draw.MolToImage(mol, size=size)
    except Exception:
        return None

def pubchem_png_image(cid: int, size=300) -> bytes | None:
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?image_size={size}x{size}"
        r = _http_get(url, timeout=35.0)
        return r.content if r and r.ok else None
    except Exception:
        return None

def py3dmol_view_from_pubchem_cid(cid: int, width=420, height=320):
    if not _P3D:
        return None
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        r = _http_get(url, timeout=40.0)
        if not r or not r.ok or not r.text or ("V2000" not in r.text and "V3000" not in r.text):
            return None
        v = py3Dmol.view(width=width, height=height)
        v.addModel(r.text, 'sdf')
        v.setStyle({'stick': {}})
        v.zoomTo()
        return v
    except Exception:
        return None

# ----------------------------
# Sorting helper (non-deprecated)
# ----------------------------
def sort_dataframe(df: pd.DataFrame, by: str, ascending: bool) -> pd.DataFrame:
    s = df[by]
    sn = pd.to_numeric(s, errors='coerce')
    if sn.notna().sum() >= max(1, len(s) // 2):
        df = df.assign(__sort_key__=sn)
        df = df.sort_values(by="__sort_key__", ascending=ascending, na_position="last").drop(columns="__sort_key__")
    else:
        df = df.sort_values(by=by, ascending=ascending, key=lambda x: x.astype(str).str.casefold())
    return df.reset_index(drop=True)

# ----------------------------
# Download links
# ----------------------------
def make_download_link(df: pd.DataFrame, filename: str, filetype: str = "csv", sep: str = ",") -> str:
    if df is None or df.empty:
        return "<em>No data to download</em>"
    if filetype == "csv":
        buf = io.StringIO(); df.to_csv(buf, index=False, sep=sep)
        data, mime = buf.getvalue(), "text/csv"
    elif filetype == "xlsx":
        buf = io.BytesIO(); df.to_excel(buf, index=False)
        data, mime = buf.getvalue(), "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    else:
        raise ValueError("Unsupported type")
    b64 = base64.b64encode(data.encode() if isinstance(data, str) else data).decode()
    return f'<a download="{filename}" href="data:{mime};base64,{b64}">⬇️ Download {filename}</a>'

def make_blob_download_link(data: str | bytes, filename: str, mime: str, label: str | None = None) -> str:
    payload = data.encode("utf-8") if isinstance(data, str) else data
    if not payload:
        return "<em>No data to download</em>"
    b64 = base64.b64encode(payload).decode("ascii")
    text = html.escape(label or f"Download {filename}")
    safe_filename = html.escape(filename)
    return f'<a download="{safe_filename}" href="data:{mime};base64,{b64}">{text}</a>'

def safe_report_filename(compound: str, extension: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", compound or "compound").strip("_") or "compound"
    return f"{slug}_chembl_pk_report.{extension}"

def display_report_downloads(compound: str, sections: list[dict], row_limit):
    if not sections:
        return
    title = f"ChEMBL Bioactivity and Predicted Pharmacokinetics: {compound}"
    display(Markdown("#### Report Downloads"))
    html_report = build_interactive_html_report(title, sections, row_limit=row_limit)
    display(HTML(make_blob_download_link(html_report, safe_report_filename(compound, "html"), "text/html", "Download interactive HTML report")))
    try:
        pdf_report = build_pdf_report(title, sections, row_limit=row_limit)
        display(HTML(make_blob_download_link(pdf_report, safe_report_filename(compound, "pdf"), "application/pdf", "Download PDF report")))
    except Exception as e:
        display(Markdown(f"> PDF report unavailable: {e}. The HTML report can still be printed or saved as PDF from the browser."))

def value_from_property_df(df: pd.DataFrame, *property_names: str):
    if df is None or df.empty or "Property" not in df.columns or "Value" not in df.columns:
        return None
    wanted = {name.casefold() for name in property_names}
    hits = df[df["Property"].astype(str).str.casefold().isin(wanted)]
    if hits.empty:
        return None
    value = hits.iloc[0]["Value"]
    return None if pd.isna(value) else str(value)

def descriptors_from_property_df(df: pd.DataFrame) -> dict:
    if df is None or df.empty or "Property" not in df.columns or "Value" not in df.columns:
        return {}
    out = {}
    for _, row in df.iterrows():
        prop = str(row.get("Property", "")).strip()
        value = row.get("Value")
        if prop and not pd.isna(value):
            out[prop] = value
    if out:
        out["DescriptorSource"] = "PubChem basic properties"
    return out

def format_pk_simulation_summary(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    columns = [
        "Route",
        "Dose (mg)",
        "Bioavailability F",
        "Absorption ka (1/h)",
        "Lag time (h)",
        "Cmax (mg/L)",
        "Tmax (h)",
        "AUC 0-inf (mg*h/L)",
        "Elimination t1/2 (h)",
        "Absorption t1/2 (h)",
        "Post-peak 50% decline (h)",
        "Onset to MEC (h)",
        "Falls below MEC (h)",
        "Duration above MEC (h)",
        "Therapeutic-window duration (h)",
        "First exceeds MTC (h)",
        "Time above MTC (h)",
        "Volume to administer (mL)",
        "Assumption",
    ]
    out = df[[col for col in columns if col in df.columns]].copy()
    numeric_cols = [col for col in out.columns if col != "Assumption" and col != "Route"]
    for col in numeric_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce").round(3)
    return out.fillna("")

def plot_pk_curves(curve_df: pd.DataFrame, mec_mg_l: float | None = None, mtc_mg_l: float | None = None):
    if not _MPL or curve_df is None or curve_df.empty:
        return None
    fig, ax = plt.subplots(figsize=(9, 4.8))
    for route, route_df in curve_df.groupby("Route"):
        ax.plot(route_df["Time (h)"], route_df["Concentration (mg/L)"], label=route, linewidth=2)
    if mec_mg_l and mec_mg_l > 0:
        ax.axhline(mec_mg_l, color="green", linestyle="--", linewidth=1.4, label="MEC")
    if mtc_mg_l and mtc_mg_l > 0:
        ax.axhline(mtc_mg_l, color="red", linestyle="--", linewidth=1.4, label="MTC")
    ax.set_title("Predicted concentration-time curve")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Concentration (mg/L)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    return fig

def plot_active_metabolite_curves(curve_df: pd.DataFrame):
    if not _MPL or curve_df is None or curve_df.empty:
        return None
    fig, ax = plt.subplots(figsize=(9, 4.8))
    for (route, curve_name), curve_part in curve_df.groupby(["Route", "Curve"]):
        linestyle = "--" if "metabolite" in curve_name.casefold() else "-"
        linewidth = 2.6 if "Total" in curve_name else 1.8
        ax.plot(
            curve_part["Time (h)"],
            curve_part["Relative active-moiety index"],
            label=f"{route}: {curve_name}",
            linestyle=linestyle,
            linewidth=linewidth,
        )
    ax.set_title("Curated active metabolite / active-moiety approximation")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Relative active-moiety index")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize="small")
    fig.tight_layout()
    return fig

def formula_svg_html(label: str, formula: str) -> str:
    label_html = html.escape(label)
    formula_html = html.escape(formula)
    if not _MPL:
        return f'<div class="pk-formula"><strong>{label_html}</strong><pre>{formula_html}</pre></div>'
    try:
        fig = plt.figure(figsize=(0.01, 0.01))
        fig.text(0, 0, f"${formula}$", fontsize=15, color="black")
        buf = io.StringIO()
        fig.savefig(buf, format="svg", bbox_inches="tight", transparent=True, pad_inches=0.08)
        plt.close(fig)
        svg = buf.getvalue()
        b64 = base64.b64encode(svg.encode("utf-8")).decode("ascii")
        return (
            '<div class="pk-formula" style="margin:0.7rem 0;">'
            f'<div style="font-weight:600;margin-bottom:0.2rem;">{label_html}</div>'
            f'<img src="data:image/svg+xml;base64,{b64}" alt="{formula_html}" '
            'style="max-width:100%;height:auto;background:white;" />'
            '</div>'
        )
    except Exception:
        return f'<div class="pk-formula"><strong>{label_html}</strong><pre>{formula_html}</pre></div>'

def display_pk_formulas():
    blocks = ["<h4>Descriptor-derived estimate formulas</h4>"]
    blocks.extend(formula_svg_html(label, formula) for label, formula in pk_estimate_formula_items())
    blocks.append(
        "<p>These descriptor formulas are transparent heuristics. The penalty terms use clipped descriptor excesses; "
        "<code>sigma</code> is the logistic function, and <code>q_adj</code> is a formal-charge adjustment.</p>"
    )
    blocks.append("<h4>One-compartment simulation formulas</h4>")
    blocks.extend(formula_svg_html(label, formula) for label, formula in pk_formula_items())
    blocks.append(
        "<p>Onset is the first time the predicted concentration reaches MEC. "
        "Duration is the time above MEC; if MTC is supplied or estimated, the table also reports "
        "time above MTC and time inside the approximate therapeutic window.</p>"
    )
    display(HTML("\n".join(blocks)))

# ----------------------------
# Interactive UI
# ----------------------------
def interactive_mode():
    text = widgets.Text(
        value='scopolamine',
        description='Compound:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='420px')
    )
    act_filter = widgets.SelectMultiple(
        options=['IC50', 'Ki', 'KA', 'Kd'],
        value=['IC50', 'Ki'],
        description='Activity filter',
        layout=widgets.Layout(width='220px', height='110px')
    )
    sort_col = widgets.Dropdown(
        options=['Target', 'Activity', 'Value', 'Units', 'Kd (nM) (from KA)'],
        value='Target',
        description='Sort by:',
        layout=widgets.Layout(width='220px')
    )
    sort_asc = widgets.ToggleButtons(
        options=[('Asc', True), ('Desc', False)],
        value=True,
        description='Order:',
        layout=widgets.Layout(width='200px')
    )
    _sep_options = [('Comma ,', ','), ('Semicolon ;', ';'), ('Tab \\t', '\t')]
    sep_choice = widgets.Dropdown(
        options=_sep_options,
        value=_sep_options[0][1],
        description='CSV separator:',
        layout=widgets.Layout(width='240px')
    )
    report_rows = widgets.Dropdown(
        options=[('10 rows', 10), ('25 rows', 25), ('50 rows', 50), ('100 rows', 100), ('All rows', 'all')],
        value=50,
        description='Report rows/table:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='240px')
    )
    dose_amount = widgets.FloatText(
        value=10.0,
        description='Dose:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='180px')
    )
    dose_unit = widgets.Dropdown(
        options=[('micrograms', 'microgram'), ('milligrams', 'milligram'), ('grams', 'gram')],
        value='milligram',
        description='Unit:',
        layout=widgets.Layout(width='190px')
    )
    body_weight = widgets.FloatText(
        value=70.0,
        description='Body weight kg:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='220px')
    )
    route_select = widgets.SelectMultiple(
        options=['Oral', 'Sublingual', 'Intranasal', 'Subcutaneous', 'Intravenous'],
        value=['Oral'],
        description='PK routes',
        layout=widgets.Layout(width='260px', height='120px')
    )
    mec_input = widgets.FloatText(
        value=0.0,
        description='MEC mg/L:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='190px')
    )
    mtc_input = widgets.FloatText(
        value=0.0,
        description='MTC mg/L:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='190px')
    )
    duration_input = widgets.FloatText(
        value=24.0,
        description='Curve hours:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='200px')
    )
    injection_conc = widgets.FloatText(
        value=0.0,
        description='Injection mg/mL:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='230px')
    )
    ref_active_dose = widgets.FloatText(
        value=0.0,
        description='Min active dose:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='220px')
    )
    ref_toxic_dose = widgets.FloatText(
        value=0.0,
        description='Min toxic dose:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='220px')
    )
    ref_dose_unit = widgets.Dropdown(
        options=[('micrograms', 'microgram'), ('milligrams', 'milligram'), ('grams', 'gram')],
        value='milligram',
        description='Ref unit:',
        layout=widgets.Layout(width='200px')
    )
    ref_route = widgets.Dropdown(
        options=['Oral', 'Sublingual', 'Intranasal', 'Subcutaneous', 'Intravenous'],
        value='Oral',
        description='Ref route:',
        layout=widgets.Layout(width='220px')
    )
    ref_cmax_fraction = widgets.FloatSlider(
        value=0.5,
        min=0.05,
        max=1.0,
        step=0.05,
        description='Ref Cmax fraction:',
        readout_format='.2f',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='380px')
    )
    show_3d = widgets.Checkbox(
        value=False,
        description='Fetch 3D structure (slower)',
        indent=False,
        layout=widgets.Layout(width='250px')
    )
    show_pubchem_props = widgets.Checkbox(
        value=False,
        description='Fetch full PubChem properties (slower)',
        indent=False,
        layout=widgets.Layout(width='310px')
    )
    run_btn = widgets.Button(description="Search", button_style='primary')
    out = widgets.Output()

    controls = widgets.HBox([text, run_btn])
    filters  = widgets.HBox([act_filter, sort_col, sort_asc, sep_choice, report_rows])
    pk_controls = widgets.VBox([
        widgets.HBox([dose_amount, dose_unit, body_weight, duration_input]),
        widgets.HBox([route_select, widgets.VBox([mec_input, mtc_input, injection_conc])]),
        widgets.HBox([ref_active_dose, ref_toxic_dose, ref_dose_unit, ref_route]),
        ref_cmax_fraction,
        widgets.HBox([show_3d, show_pubchem_props]),
        widgets.HTML("<em>MEC/MTC are optional thresholds in mg/L. If left as 0, a minimum active/toxic reference dose estimates them as: threshold = Ref Cmax fraction × predicted Cmax for that reference dose and route. Example: 0.50 means half of the predicted reference-dose peak.</em>"),
    ])

    def on_click(_):
        with out:
            clear_output()

            compound = text.value
            compound = _normalize_text(compound)
            if not compound:
                print("❌ Please enter a compound name.")
                return

            display(Markdown(f"## Results for **{compound}**"))
            display(Markdown("Data sources: **ChEMBL** (bioactivity), **PubChem** (structure & properties), and transparent structure-based PK heuristics."))
            report_sections = []
            report_limit = report_row_limit(report_rows.value)
            datatable_page_length = -1 if report_limit is None else report_limit

            # --- ChEMBL PD table
            df_pd = pd.DataFrame()
            try:
                chembl_id = get_chembl_id(compound)
                acts = fetch_activities(chembl_id)
                df_pd = build_activity_df(acts)

                if not df_pd.empty:
                    selected = list(act_filter.value)
                    if selected:
                        df_pd = df_pd[df_pd['Activity'].isin(selected)]
                    df_pd = sort_dataframe(df_pd, by=sort_col.value, ascending=sort_asc.value)

                    display(Markdown("### Pharmacodynamic Bioactivities (Homo sapiens)"))
                    show(df_pd, classes="display compact cell-border", maxBytes=0, pageLength=datatable_page_length)
                    report_sections.append({"title": "Pharmacodynamic Bioactivities (Homo sapiens)", "df": df_pd})

                    display(HTML(make_download_link(df_pd, "bioactivity.csv", "csv", sep=sep_choice.value)))
                    display(HTML(make_download_link(df_pd, "bioactivity.xlsx", "xlsx")))
                else:
                    display(Markdown("> No human bioactivity rows returned by ChEMBL."))
            except Exception as e:
                display(Markdown(f"**ChEMBL error:** {e}"))

            # --- PubChem: CID, 2D image, 3D viewer, properties
            cid = None
            try:
                cid = pubchem_cid_from_name(compound)
            except Exception:
                cid = None

            if cid is None:
                display(Markdown("> PubChem lookup failed; structure/properties unavailable."))
                display_report_downloads(compound, report_sections, report_rows.value)
                return

            display(Markdown(f"### PubChem\nCID: **{cid}**"))

            df_basic = pubchem_basic_props_df(cid)
            if not df_basic.empty:
                show(df_basic, classes="display compact cell-border", maxBytes=0, pageLength=datatable_page_length)
                report_sections.append({"title": "PubChem Basic Properties", "df": df_basic})
                display(HTML(make_download_link(df_basic, "pubchem_basic.csv", "csv", sep=sep_choice.value)))

            # SMILES for RDKit rendering and predicted PK descriptors. Avoid deprecated PubChemPy accessors.
            smiles = value_from_property_df(
                df_basic,
                "IsomericSMILES",
                "CanonicalSMILES",
                "Isomeric SMILES",
                "Canonical SMILES",
            )
            descriptors = descriptors_from_property_df(df_basic)
            has_core_descriptors = any(key in descriptors for key in ["MolecularWeight", "XLogP", "TPSA"])
            if not has_core_descriptors:
                descriptors = descriptors_from_pubchem_cid(cid, smiles=smiles)
            smiles = smiles or descriptors.get("IsomericSMILES") or descriptors.get("CanonicalSMILES")

            display(Markdown("#### Predicted Pharmacokinetics (structure-based)"))
            df_pk = pk_profile_from_descriptors(descriptors, label=f"PubChem CID {cid}")
            if not df_pk.empty:
                show(df_pk, classes="display compact cell-border", maxBytes=0, pageLength=datatable_page_length)
                report_sections.append({"title": "Predicted Pharmacokinetics (structure-based)", "df": df_pk, "note": PK_ESTIMATE_NOTE})
                display(HTML(make_download_link(df_pk, "predicted_pharmacokinetics.csv", "csv", sep=sep_choice.value)))
                display(HTML(make_download_link(df_pk, "predicted_pharmacokinetics.xlsx", "xlsx")))
            display(Markdown(f"> {PK_ESTIMATE_NOTE}"))

            display(Markdown("#### Dose and Route Simulation"))
            if not descriptors:
                display(Markdown("> PK simulation unavailable because PubChem/RDKit descriptors were not available."))
            else:
                try:
                    pk_prediction = predict_pk_from_descriptors(descriptors)
                    estimated_mec, estimated_mtc, threshold_df = estimate_exposure_thresholds(
                        pk_prediction,
                        active_dose_amount=ref_active_dose.value,
                        toxic_dose_amount=ref_toxic_dose.value,
                        reference_dose_unit=ref_dose_unit.value,
                        reference_route=ref_route.value,
                        body_weight_kg=body_weight.value,
                        cmax_fraction=ref_cmax_fraction.value,
                        duration_h=duration_input.value,
                    )
                    mec = mec_input.value if mec_input.value and mec_input.value > 0 else estimated_mec
                    mtc = mtc_input.value if mtc_input.value and mtc_input.value > 0 else estimated_mtc
                    if not threshold_df.empty:
                        threshold_display = threshold_df.copy()
                        for col in ["Estimate (mg/L)", "Reference dose (mg)", "Reference Cmax (mg/L)", "Cmax fraction"]:
                            if col in threshold_display.columns:
                                threshold_display[col] = pd.to_numeric(threshold_display[col], errors="coerce").round(4)
                        display(Markdown("##### Estimated MEC/MTC From Reference Dose"))
                        show(threshold_display, classes="display compact cell-border", maxBytes=0, pageLength=datatable_page_length)
                        report_sections.append({"title": "Estimated MEC/MTC From Reference Dose", "df": threshold_display})
                    if mec_input.value and mec_input.value > 0:
                        display(Markdown(f"> Using manually entered MEC: **{mec:.4g} mg/L**."))
                    elif estimated_mec is not None:
                        display(Markdown(f"> Estimated MEC from reference active dose: **{estimated_mec:.4g} mg/L**."))
                    if mtc_input.value and mtc_input.value > 0:
                        display(Markdown(f"> Using manually entered MTC: **{mtc:.4g} mg/L**."))
                    elif estimated_mtc is not None:
                        display(Markdown(f"> Estimated MTC from reference toxic dose: **{estimated_mtc:.4g} mg/L**."))
                    curve_df, sim_summary = simulate_pk_curves(
                        pk_prediction,
                        dose_amount=dose_amount.value,
                        dose_unit=dose_unit.value,
                        routes=list(route_select.value),
                        body_weight_kg=body_weight.value,
                        mec_mg_l=mec,
                        mtc_mg_l=mtc,
                        duration_h=duration_input.value,
                        injection_concentration_mg_ml=injection_conc.value,
                    )
                    fig = plot_pk_curves(curve_df, mec_mg_l=mec, mtc_mg_l=mtc)
                    if fig is not None:
                        display(fig)
                        plt.close(fig)
                    else:
                        display(Markdown("> Matplotlib is unavailable; showing simulation tables only."))

                    sim_display = format_pk_simulation_summary(sim_summary)
                    if not sim_display.empty:
                        show(sim_display, classes="display compact cell-border", maxBytes=0, pageLength=datatable_page_length)
                        report_sections.append({"title": "Dose and Route Simulation Summary", "df": sim_display, "note": PK_SIMULATION_NOTE})
                        display(HTML(make_download_link(sim_display, "pk_simulation_summary.csv", "csv", sep=sep_choice.value)))
                        display(HTML(make_download_link(curve_df, "pk_concentration_curve.csv", "csv", sep=sep_choice.value)))
                    if mec is None:
                        display(Markdown("> Enter an MEC above 0 mg/L or a minimum active reference dose to calculate onset and duration above MEC."))
                    display(Markdown("> **Half-life note:** elimination half-life is compound/system dependent and usually route-independent in this one-compartment model. Route differences are shown as absorption half-life and post-peak 50% decline time."))
                    display(Markdown(f"> {PK_SIMULATION_NOTE}"))
                    display_pk_formulas()

                    metabolite_curve_df, metabolite_summary_df = simulate_active_metabolite_curves(
                        compound,
                        curve_df,
                        sim_summary,
                        pk_prediction,
                        body_weight_kg=body_weight.value,
                    )
                    if not metabolite_curve_df.empty:
                        display(Markdown("##### Active Metabolite / Active-Moiety Curve"))
                        metabolite_fig = plot_active_metabolite_curves(metabolite_curve_df)
                        if metabolite_fig is not None:
                            display(metabolite_fig)
                            plt.close(metabolite_fig)
                        if not metabolite_summary_df.empty:
                            metabolite_display = metabolite_summary_df.copy()
                            for col in ["Formation fraction", "Metabolite half-life (h)", "Relative potency vs parent", "Peak active-moiety index", "Peak time (h)"]:
                                metabolite_display[col] = pd.to_numeric(metabolite_display[col], errors="coerce").round(4)
                            show(metabolite_display, classes="display compact cell-border", maxBytes=0, pageLength=datatable_page_length)
                            report_sections.append({"title": "Active Metabolite / Active-Moiety Summary", "df": metabolite_display})
                            display(HTML(make_download_link(metabolite_display, "active_metabolite_summary.csv", "csv", sep=sep_choice.value)))
                            display(HTML(make_download_link(metabolite_curve_df, "active_metabolite_curve.csv", "csv", sep=sep_choice.value)))
                except Exception as e:
                    display(Markdown(f"**PK simulation error:** {e}"))

            display(Markdown("#### Active Metabolite / Prodrug Caveat"))
            df_metabolite = active_metabolite_notes_for_compound(compound)
            show(df_metabolite, classes="display compact cell-border", maxBytes=0, pageLength=datatable_page_length)
            report_sections.append({"title": "Active Metabolite / Prodrug Caveat", "df": df_metabolite})

            # 2D structure
            display(Markdown("#### 2D Structure"))
            img = rdkit_image_from_smiles(smiles) if smiles else None
            if img is not None:
                display(img)
            else:
                png = pubchem_png_image(cid, size=320)
                if png:
                    display(Image(png))
                else:
                    display(Markdown("> Unable to render 2D structure (RDKit and PNG fallback both failed)."))

            # 3D structure
            if show_3d.value:
                display(Markdown("#### 3D Structure (interactive)"))
                viewer = py3dmol_view_from_pubchem_cid(cid)
                if viewer is not None:
                    viewer.show()
                else:
                    msg = "py3Dmol not available" if not _P3D else "No PubChem 3D conformer found or retrieval failed"
                    display(Markdown(f"> 3D viewer unavailable: {msg}."))
            else:
                display(Markdown("> 3D structure skipped. Enable **Fetch 3D structure (slower)** to retrieve it."))

            # Experimental/Computed properties (robust)
            if show_pubchem_props.value:
                display(Markdown("#### Experimental / Computed Properties (PubChem)"))
                df_exp = pubchem_experimental_props_df(cid)
                if not df_exp.empty:
                    show(df_exp, classes="display compact cell-border", maxBytes=0, pageLength=datatable_page_length)
                    report_sections.append({"title": "Experimental / Computed Properties (PubChem)", "df": df_exp})
                    display(HTML(make_download_link(df_exp, "pubchem_properties.csv", "csv", sep=sep_choice.value)))
                    display(HTML(make_download_link(df_exp, "pubchem_properties.xlsx", "xlsx")))
                else:
                    display(Markdown("> No experimental/computed properties found (or parse failed)."))
            else:
                display(Markdown("> Full PubChem property scrape skipped. Enable **Fetch full PubChem properties (slower)** to retrieve it."))

            display_report_downloads(compound, report_sections, report_rows.value)

    run_btn.on_click(on_click)
    display(widgets.VBox([controls, filters, pk_controls, out]))

if __name__ == "__main__":
    interactive_mode()
