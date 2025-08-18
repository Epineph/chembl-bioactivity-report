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
# - Binder/Voila-hardened networking (UA + retries + TXT fallback)
#
# Notes
# * RDKit via conda-forge is recommended (see environment.yml).
# * 3D model depends on PubChem having a 3D conformer.
# * Experimental properties are heterogeneous; availability varies.
# * This cell only builds the UI; no network calls occur until you press "Search".

import io, base64, json, re, requests, contextlib, warnings, time, unicodedata
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
_RDKIT, _rdkit_msg = _silent_import_rdkit()

if _RDKIT:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
if _PCP:
    import pubchempy as pcp
if _P3D:
    import py3Dmol

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
            tidy = df.T.reset_index().rename(columns={'index': 'Property', 0: 'Value'})
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
    run_btn = widgets.Button(description="Search", button_style='primary')
    out = widgets.Output()

    controls = widgets.HBox([text, run_btn])
    filters  = widgets.HBox([act_filter, sort_col, sort_asc, sep_choice])

    def on_click(_):
        with out:
            clear_output()

            compound = text.value
            compound = _normalize_text(compound)
            if not compound:
                print("❌ Please enter a compound name.")
                return

            display(Markdown(f"## Results for **{compound}**"))
            display(Markdown("Data sources: **ChEMBL** (bioactivity), **PubChem** (structure & properties)."))

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
                    show(df_pd, classes="display compact cell-border", maxBytes=0, pageLength=50)

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
                return

            display(Markdown(f"### PubChem\nCID: **{cid}**"))

            df_basic = pubchem_basic_props_df(cid)
            if not df_basic.empty:
                show(df_basic, classes="display compact cell-border", maxBytes=0, pageLength=50)
                display(HTML(make_download_link(df_basic, "pubchem_basic.csv", "csv", sep=sep_choice.value)))

            # 2D structure
            smiles = None
            if _PCP:
                try:
                    comps = pcp.get_compounds([cid], 'cid')
                    if comps:
                        smiles = comps[0].isomeric_smiles or comps[0].canonical_smiles
                except Exception:
                    smiles = None

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
            display(Markdown("#### 3D Structure (interactive)"))
            viewer = py3dmol_view_from_pubchem_cid(cid)
            if viewer is not None:
                viewer.show()
            else:
                msg = "py3Dmol not available" if not _P3D else "No PubChem 3D conformer found or retrieval failed"
                display(Markdown(f"> 3D viewer unavailable: {msg}."))

            # Experimental/Computed properties (robust)
            display(Markdown("#### Experimental / Computed Properties (PubChem)"))
            df_exp = pubchem_experimental_props_df(cid)
            if not df_exp.empty:
                show(df_exp, classes="display compact cell-border", maxBytes=0, pageLength=50)
                display(HTML(make_download_link(df_exp, "pubchem_properties.csv", "csv", sep=sep_choice.value)))
                display(HTML(make_download_link(df_exp, "pubchem_properties.xlsx", "xlsx")))
            else:
                display(Markdown("> No experimental/computed properties found (or parse failed)."))

    run_btn.on_click(on_click)
    display(widgets.VBox([controls, filters, out]))

# Launch UI (Voila will execute the notebook top-down)
interactive_mode()
