# Structural screening

Install dependencies: `python -m pip install -r requirements.txt`.
Export hits.json from the browser. Download an experimental PDB or mmCIF from https://www.rcsb.org/. Save the exact imported antibody chain sequence as plain text (no FASTA header) in chain.txt.

```sh
python structure_screen.py --hits hits.json --structure complex.cif --chain H --sequence chain.txt --asset ASSET_NAME --sequence-column 'Heavy chain' --output assessed.json
```
Replace asset, column and chain with actual identifiers. Only a unique exact match of the entire supplied sequence to observed standard residues is accepted; unresolved residues or ambiguous mappings may prevent assessment. SASA is calculated by Biopython Shrake-Rupley with a 1.4 Å water probe and 960 points per atom, using the first supplied model including its other chains. Supply a biologically relevant assembly; crystal contacts, glycans and missing atoms may affect the result. SASA is not proof of reaction feasibility. No binding-site proximity or mechanism-specific pass threshold is asserted. This is a CLI backend, not yet connected to the browser. No runtime validation has been performed.
