# %% [markdown]

# notebook to fetch pharmacodynamic bioactivity of compounds

Neat ChEMBL Bioactivity Report

- **Fetch pharmacodynamic bioactivity**
- Pulls all human bioactivities (*IC50*, *Ki*, *KA*, etc.) for a given compound.
- Looks up each target’s *preferred name*.
- **Builds a DataFrame and prints it as a pretty GitHub‑style Markdown table** when run
  in a terminal, or as an HTML table in Jupyter.
- When run from terminal, if using linux, make sure the script is executable by: sudo
  chmod +x chembl_bioactivity.py anr run from terminal:
- python chembl_bioactivity.py <compound>

______________________________________________________________________

# %% [markdown]

## 2) Install the required libraries

The notebook **requires the following libraries:**

- chembl-webresource-client
- pandas
- tabulate

Install with pip, e.g.,:

python -m pip install chembl-webresource-client pandas tabulate

Or use your favorite packpage manager, e.g.,:

conda install -c conda-forge chembl-webresource-client pandas tabulate

or run clone the repository and run:

python -m pip install -r requirements.txt

# %% [markdown]

## 3) Imports & Setup

```python
import sys
import pandas as pd
from chembl_webresource_client.new_client import new_client
from tabulate import tabulate
```

______________________________________________________________________

# %% [markdown]

## 4) Pass compount name to function

Pass a string to to the function get_chembl_id, e.g.,

def get_chembl_id(compound: Literal['scopolamine']) -> str:

```python
from typing import Literal

def get_chembl_id(compound: Literal['aspirin']) -> str:
    """Lookup compound by preferred name in ChEMBL."""
    mol_client = new_client.molecule
    # Use the passed-in compound string, not a hard-coded value
    res = mol_client.filter(pref_name__iexact=compound)
    if not res:
        raise ValueError(f"No ChEMBL entry for '{compound}'")
    return res[0]['molecule_chembl_id']
```

______________________________________________________________________

```python
def fetch_activities(chembl_id: str) -> list[dict]:
    """Fetch all Homo sapiens bioactivities for the given ChEMBL ID."""
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
```

______________________________________________________________________

```python
def fetch_target_names(target_ids: set[str]) -> dict[str, str]:
    """Map each ChEMBL target ID to its preferred name."""
    tgt_client = new_client.target
    names = {}
    for tid in target_ids:
        rec = tgt_client.filter(target_chembl_id=tid).only(['pref_name'])
        names[tid] = rec[0]['pref_name'] if rec else tid
    return names
```

______________________________________________________________________

```python
def build_activity_df(acts: list[dict]) -> pd.DataFrame:
    """Build a tidy DataFrame and compute Kd (nM) for KA entries."""
    rows = []
    for a in acts:
        tid   = a.get('target_chembl_id') or 'Unknown'
        typ   = a.get('standard_type')    or ''
        val   = a.get('standard_value')   or ''
        unit  = a.get('standard_units')   or ''
        kd_nm = ''
        # Compute Kd in nanomolar for association constants
        if typ.upper() == 'KA' and val and unit.strip() in ['M^-1','M-1','1/M']:
            try:
                kd_m  = 1.0 / float(val)
                kd_nm = round(kd_m * 1e9, 2)
            except Exception:
                kd_nm = ''
        rows.append({
            'Target (ChEMBL)': tid,
            'Activity':        typ,
            'Value':           val,
            'Units':           unit,
            'Kd (nM)':         kd_nm
        })
    df = pd.DataFrame(rows)
    # Replace ChEMBL IDs with human-readable names
    unique_tids = set(df['Target (ChEMBL)'])
    name_map = fetch_target_names(unique_tids)
    df['Target (ChEMBL)'] = df['Target (ChEMBL)'].map(name_map)
    # Keep only non-empty bioactivity entries
    df = df[['Target (ChEMBL)', 'Activity', 'Value', 'Units', 'Kd (nM)']]
    df.dropna(how='all', subset=['Value'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df
```

______________________________________________________________________

```python
def display_df(df: pd.DataFrame):
    """Display as HTML table in Jupyter, or Markdown table on CLI."""
    try:
        from IPython.display import display
        styler = df.style
        # Use hide_index if available, otherwise fallback to hide
        if hasattr(styler, 'hide_index'):
            styler = styler.hide_index()
        else:
            styler = styler.hide(axis="index")
        display(styler)
    except ImportError:
        print(tabulate(df, headers='keys', tablefmt='github', showindex=False))
```

______________________________________________________________________

```python
def display_df(df: pd.DataFrame):
    """Display as HTML table in Jupyter, or Markdown table on CLI."""
    try:
        from IPython.display import display
        styler = df.style
        # Use hide_index if available, otherwise fallback to hide
        if hasattr(styler, 'hide_index'):
            styler = styler.hide_index()
        else:
            styler = styler.hide(axis="index")
        display(styler)
    except ImportError:
        print(tabulate(df, headers='keys', tablefmt='github', showindex=False))
```

______________________________________________________________________

```python
def main():
    """Main entry point: lookup, fetch, build, and display."""
    compound = sys.argv[1] if len(sys.argv) > 1 else 'scopolamine'
    print(f"\n🔍 Looking up '{compound}' in ChEMBL…")
    chembl_id = get_chembl_id(compound)
    print(f"   → Found ChEMBL ID: {chembl_id}\n")

    print("📋 Fetching human bioactivities…")
    acts = fetch_activities(chembl_id)
    df   = build_activity_df(acts)
    print(f"\n🏷  Retrieved {len(df)} records:\n")
    display_df(df)

if __name__ == '__main__':  # pragma: no cover
    main()
```

______________________________________________________________________

# 🧪 Neat ChEMBL Bioactivity Report (Fixed)

**Help**

- Pulls all human bioactivities (IC50, Ki, KA, etc.) for a given compound.
- Looks up each target’s preferred name.
- Builds a DataFrame and prints it as a pretty GitHub-style Markdown table\
  when run in a terminal, or as an HTML table in Jupyter.

**Usage**

```bash
python chembl_bioactivity_fixed.py <compound>
```
