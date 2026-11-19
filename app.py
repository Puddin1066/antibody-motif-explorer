"""Integrated research-use antibody motif explorer."""
import csv, io, json, re, tempfile, urllib.request
from pathlib import Path
import streamlit as st
import pandas as pd
import py3Dmol
import streamlit.components.v1 as components
from structure_screen import assess

st.set_page_config(page_title="Antibody Motif Explorer",layout="wide")
st.title("Antibody Motif Explorer")
st.caption("Research use only. DRAC is an illustrative motif, not a validated engineering mechanism.")
st.markdown("[Download therapeutic antibody CSV from Thera-SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/)")
upload=st.file_uploader("Therapeutic antibody CSV or TSV",type=["csv","tsv"])
if not upload:
    st.stop()
try:
    df=pd.read_csv(upload,sep="\t" if upload.name.endswith(".tsv") else ",",dtype=str,keep_default_na=False)
except Exception as e:
    st.error(str(e));st.stop()
if df.empty:
    st.warning("No records");st.stop()
cols=list(df.columns)
asset=st.selectbox("Asset name column",cols,index=next((i for i,c in enumerate(cols) if re.search(r"therapeutic|^name$|antibody name",c,re.I)),0))
chains=st.multiselect("Sequence columns",cols,default=[c for c in cols if re.search(r"sequence|^vh$|^vl$|heavy chain|light chain",c,re.I)])
motif=st.text_input("Motif (X = any standard amino acid)","DRAC").strip().upper()
threshold=st.slider("Minimum sequence identity (%)",0,100,100)
aa=set("ACDEFGHIKLMNPQRSTVWY")
if not motif or not set(motif)<=aa|{"X"} or not chains:
    st.info("Select sequence columns and provide a valid motif.");st.stop()
hits=[]
for idx,row in df.iterrows():
    for col in chains:
        seq="".join(row[col].split()).upper()
        if not seq or not set(seq)<=aa|{"X"}:continue
        for p in range(len(seq)-len(motif)+1):
            segment=seq[p:p+len(motif)]
            score=100*sum(b in aa and (a=="X" or a==b) for a,b in zip(motif,segment))/len(motif)
            if score>=threshold:
                hits.append(dict(asset=row[asset],chain=col,position=p+1,segment=segment,identity=score,structure="Not assessed",row_index=int(idx)))
hits.sort(key=lambda h:(-h["identity"],h["asset"],h["position"]))
st.metric("Sequence hits",len(hits))
if not hits:st.stop()
table=pd.DataFrame(hits)
st.dataframe(table,use_container_width=True)
st.download_button("Export sequence hits JSON",json.dumps(hits,indent=2),"hits.json")
def safe_csv(frame):
    return frame.map(lambda x: "'"+x if isinstance(x,str) and x.startswith(("=","+","-","@","\t","\r")) else x).to_csv(index=False)
st.download_button("Export sequence hits CSV",safe_csv(table),"hits.csv")
st.header("Structural assessment")
choice=st.selectbox("Candidate",range(len(hits)),format_func=lambda i:f"{hits[i]['asset']} | {hits[i]['chain']} | {hits[i]['position']} | {hits[i]['segment']}")
hit=hits[choice]
row=df.iloc[hit["row_index"]]
st.write("Source metadata",row.to_dict())
seq="".join(row[hit["chain"]].split()).upper()
pdb=st.text_input("RCSB PDB identifier (from source metadata or external search)")
st.markdown("[Search RCSB PDB](https://www.rcsb.org/)")
local=st.file_uploader("Or upload structure",type=["pdb","cif","mmcif"])
chain=st.text_input("Structure chain identifier","H")
if st.button("Map and calculate accessibility"):
    try:
        if local:
            data=local.getvalue();suffix=Path(local.name).suffix.lower()
        else:
            if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}",pdb.strip()):raise ValueError("Enter a four-character PDB identifier or upload a structure.")
            suffix=".cif"
            with urllib.request.urlopen("https://files.rcsb.org/download/"+pdb.strip().upper()+".cif",timeout=30) as response:
                data=response.read(50_000_001)
        if len(data)>50_000_000:raise ValueError("Structure exceeds 50 MB.")
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/("structure"+suffix);path.write_bytes(data)
            result=assess([hit],path,chain,seq,hit["asset"],hit["chain"])[0]
            result["structure_file"]=local.name if local else pdb.upper()
        st.session_state["assessment"]=(dict(hit),result,data.decode("utf-8"),suffix,chain)
    except Exception as e:st.error(f"Assessment failed: {e}")
saved=st.session_state.get("assessment")
if saved and saved[0]==hit:
    _,result,raw,suffix,selected_chain=saved
    st.write(result)
    st.download_button("Export structural assessment",json.dumps(result,indent=2),"assessed.json")
    if "residues" in result:
        st.metric("Motif SASA (Å²)",result["motif_sasa_A2"])
    view=py3Dmol.view(width=800,height=480)
    view.addModel(raw,"pdb" if suffix==".pdb" else "cif")
    view.setStyle({"cartoon":{"color":"spectrum"}})
    for residue in result.get("residues",[]):
        sel={"chain":selected_chain,"resi":str(residue["number"])+residue["insertion_code"]}
        view.addStyle(sel,{"stick":{"color":"red"}})
    view.zoomTo()
    components.html(view._make_html(),height=500)
st.header("Commercial review")
st.caption("Manual evidence review; sequence identity and SASA do not establish licensing potential.")
company=st.text_input("Verified asset owner")
stage=st.selectbox("Verified development stage",["Unknown","Preclinical","Phase 1","Phase 2","Phase 3","Approved","Discontinued"])
evidence=st.text_area("Source URLs and partnering rationale")
st.download_button("Export candidate review",json.dumps({"candidate":hit,"metadata":row.to_dict(),"owner":company,"stage":stage,"evidence":evidence,"structural_assessment":saved[1] if saved and saved[0]==hit else None},indent=2),"candidate_review.json")
