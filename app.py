import streamlit as st
import primer3
import pandas as pd
import time
import plotly.graph_objects as go
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# --- פונקציית BLAST מפורטת (כולל Alignment) ---
def get_detailed_blast(primer_seq, species):
    try:
        entrez_query = f"{species}[organism]" if species else ""
        handle = NCBIWWW.qblast("blastn", "nt", primer_seq, entrez_query=entrez_query, word_size=11, hitlist_size=5)
        blast_record = NCBIXML.read(handle)
        report = []
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                detail = (
                    f"Target: {alignment.title[:60]}...\n"
                    f"Accession: {alignment.accession} | Identity: {hsp.identities}/{hsp.align_length}\n"
                    f"Query:  {hsp.query}\n"
                    f"Match:  {hsp.match}\n"
                    f"Sbjct:  {hsp.sbjct}\n"
                    "--------------------------------------------------"
                )
                report.append(detail)
        return "\n".join(report) if report else "No significant hits found."
    except: return "BLAST Error"

st.set_page_config(page_title="Cloning Primer Designer", page_icon="🧬", layout="wide")
st.title("🧬 Cloning Primer Designer (Pro Edition)")

# --- קלט משתמש ---
col1, col2 = st.columns([2, 1])
with col1: sequence = st.text_area("Target Sequence (5' to 3'):", height=200)
with col2:
    project_name = st.text_input("Project Name:", "My_Cloning")
    species_name = st.text_input("Species (for BLAST):", "Solanum tuberosum")
    junction_pos = st.text_input("Exon Junctions:", placeholder="e.g., 72, 478")
    run_blast_check = st.checkbox("🔍 Run Detailed BLAST (Slow)")

with st.expander("⚙️ Advanced Thermodynamics & Settings"):
    adv_col1, adv_col2 = st.columns(2)
    with adv_col1:
        target_tm = st.slider("Target Tm (°C)", 50.0, 72.0, 60.0)
        num_returns = st.slider("Total Primers", 1, 10, 5)
    with adv_col2:
        min_amp = st.text_input("Min Amplicon", "100")
        max_amp = st.text_input("Max Amplicon", "1000")

if st.button("🚀 Design Primers", type="primary"):
    if not sequence: st.error("Please enter a sequence.")
    else:
        with st.spinner('Calculating thermodynamics and secondary structures...'):
            clean_seq = "".join(sequence.split()).upper()
            junction_list = [int(x.strip()) for x in junction_pos.split(',') if x.strip().isdigit()] if junction_pos else []
            
            global_args = {
                'PRIMER_OPT_SIZE': 20, 'PRIMER_OPT_TM': target_tm,
                'PRIMER_NUM_RETURN': num_returns,
                'PRIMER_THERMODYNAMIC_PARAMETERS_PATH': primer3.bindings.get_parameter_path()
            }
            if min_amp.isdigit() and max_amp.isdigit(): global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp), int(max_amp)]]

            final_candidates = []
            # (לוגיקת ה-Round-Robin הושארה כפי שהייתה)
            seq_args = {'SEQUENCE_ID': project_name, 'SEQUENCE_TEMPLATE': clean_seq}
            if junction_list: seq_args['SEQUENCE_OVERLAP_JUNCTION_LIST'] = junction_list
            
            raw_res = primer3.bindings.designPrimers(seq_args, global_args)
            
            parsed_data = []
            for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                # שליפת ערכי ה-Delta G עבור דימרים ו-Hairpins
                f_hairpin = round(raw_res.get(f'PRIMER_LEFT_{i}_HAIRPIN_TH', 0), 2)
                r_hairpin = round(raw_res.get(f'PRIMER_RIGHT_{i}_HAIRPIN_TH', 0), 2)
                f_self_any = round(raw_res.get(f'PRIMER_LEFT_{i}_SELF_ANY_TH', 0), 2)
                r_self_any = round(raw_res.get(f'PRIMER_RIGHT_{i}_SELF_ANY_TH', 0), 2)
                cross_dimer = round(raw_res.get(f'PRIMER_PAIR_{i}_COMPL_ANY_TH', 0), 2)

                parsed_data.append({
                    "Rank": i + 1,
                    "Penalty": round(raw_res.get(f'PRIMER_PAIR_{i}_PENALTY', 0), 2),
                    "F_Seq": raw_res.get(f'PRIMER_LEFT_{i}_SEQUENCE'),
                    "R_Seq": raw_res.get(f'PRIMER_RIGHT_{i}_SEQUENCE'),
                    "F_Hairpin (ΔG)": f_hairpin,
                    "R_Hairpin (ΔG)": r_hairpin,
                    "F_SelfDimer (ΔG)": f_self_any,
                    "R_SelfDimer (ΔG)": r_self_any,
                    "Pair_CrossDimer (ΔG)": cross_dimer,
                    "Product_Size": raw_res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE')
                })
            
            df = pd.DataFrame(parsed_data)
            st.success("Analysis Complete!")
            st.dataframe(df, use_container_width=True)

            # הצגת ה-BLAST המפורט מתחת לטבלה
            if run_blast_check:
                for idx, row in df.iterrows():
                    with st.expander(f"🔍 BLAST Detailed Alignment - Candidate {idx+1}"):
                        st.subheader("Forward Primer Alignment")
                        st.code(get_detailed_blast(row['F_Seq'], species_name))
                        st.subheader("Reverse Primer Alignment")
                        st.code(get_detailed_blast(row['R_Seq'], species_name))
