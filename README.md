# Antibody Motif Explorer

Integrated research-use Streamlit application. DRAC is an illustrative amino acid query, not a validated engineering mechanism.

## Local setup (Cursor terminal, Python 3.11 or 3.12)

```sh
git clone https://github.com/Puddin1066/antibody-motif-explorer.git
cd antibody-motif-explorer
python -m venv .venv
```
Activate with `source .venv/bin/activate` (macOS/Linux) or `.venv\\Scripts\\Activate.ps1` (Windows PowerShell). Then:

```sh
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

## Deploy
Use Streamlit Community Cloud with repository main and entrypoint app.py, or build the provided Dockerfile. GitHub Pages only supports the older index.html sequence-only frontend. No API keys needed. Public deployment is not appropriate for confidential data without access controls.

## Workflow
Download a CSV from [Thera-SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/), upload and map columns. Select a motif and a candidate. Enter a known PDB identifier or upload a structure, specify its author chain ID and calculate SASA. Review the highlighted structure and export candidate metadata and manually sourced commercial evidence. RCSB downloads are automated after ID entry; discovery of matching structures and ownership is manual. See STRUCTURE.md for limitations. Full-chain unique exact mapping is required. No CDR proximity scoring or automated commercial rankings are claimed. The integrated app is not runtime validated.
