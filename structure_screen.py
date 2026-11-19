"""Map sequence motif hits to an experimental PDB/mmCIF chain and estimate SASA.
Requires biopython. Coordinates are sequence indices, not antibody numbering.
Only exact unique sequence mappings are accepted; missing residues and ambiguous
matches are reported rather than silently assigned.
"""
import argparse, json
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser, ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

def assess(hits, structure_path, chain_id, source_sequence, asset, sequence_column):
    parser = MMCIFParser(QUIET=True) if Path(structure_path).suffix.lower() in ('.cif','.mmcif') else PDBParser(QUIET=True)
    structure = parser.get_structure('complex', structure_path)
    model = next(structure.get_models())
    chain = model[chain_id]
    # Calculate accessibility in context of all atoms in the supplied model.
    ShrakeRupley(n_points=960).compute(model, level='R')
    residues = [r for r in chain if is_aa(r, standard=True)]
    observed = ''.join(seq1(r.resname) for r in residues)
    source = ''.join(source_sequence.split()).upper()
    start = observed.find(source)
    unique = start >= 0 and observed.find(source,start+1) < 0
    output=[]
    for h in hits:
        if h.get('asset') != asset or h.get('chain') != sequence_column:
            continue
        record = dict(h)
        record['structure_file'] = str(structure_path)
        record['structure_chain'] = chain_id
        p = int(h['position'])-1
        segment = h['segment']
        if not unique or p < 0 or p+len(segment)>len(source) or source[p:p+len(segment)] != segment:
            record['structure']='Unmapped: exact unique sequence alignment required'
        else:
            rr=residues[start+p:start+p+len(segment)]
            record['structure']='Mapped; accessibility estimated'
            record['residues']=[{'number':r.id[1],'insertion_code':r.id[2].strip(),'aa':seq1(r.resname),'sasa_A2':round(r.sasa,2)} for r in rr]
            record['motif_sasa_A2']=round(sum(r.sasa for r in rr),2)
        output.append(record)
    return output

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('hits','structure','chain','sequence','asset','sequence-column','output'):
        p.add_argument('--'+arg,required=True)
    a=p.parse_args()
    result=assess(json.loads(Path(a.hits).read_text()),a.structure,a.chain,Path(a.sequence).read_text(),a.asset,a.sequence_column)
    Path(a.output).write_text(json.dumps(result,indent=2))
    print(f'{len(result)} selected hits; {sum("residues" in r for r in result)} mapped')
