import streamlit as st
import primer3
import pandas as pd
import io
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

# Configuration for NCBI
Entrez.email = "spud_researcher@example.com" 

# --- Helper Functions & Thermodynamics ---
def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def check_and_get_blast(primer_seq, species):
    try:
        entrez_query = f"{species}[organism]" if species else ""
        handle = NCBIWWW.qblast("blastn", "nt", primer_seq, entrez_query=entrez_query, word_size=11, hitlist_size=3)
        blast_record = NCBIXML.read(handle)
        significant_hits = 0
        primer_len = len(primer_seq)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.align_length >= primer_len * 0.9:
                    significant_hits += 1
        return f"{significant_hits} Hits" if significant_hits > 1 else "1 Hit ✅" if significant_hits == 1 else "0 Hits ✅"
    except Exception as e: 
        return "Error"

def parse_target(target_str):
    """Parses the Smart Input into Global, Junction, or ROI"""
    target_str = str(target_str).strip()
    if not target_str or target_str.lower() in ['nan', 'none']:
        return "Global", None
    
    if '-' in target_str:
        try:
            start, end = map(int, target_str.split('-'))
            return "ROI", [start, end]
        except: 
            return "Error", None
            
    if target_str.isdigit():
        return "Junction", int(target_str)
        
    return "Error", None

def validate_7bp_anchor(primer_start, primer_length, junction, is_reverse=False):
    """Validates that a primer spanning a junction has at least 7bp on each side."""
    if is_reverse:
        r_end = primer_start - primer_length + 1
        return (primer_start >= junction + 7) and (r_end <= junction - 6)
    else:
        f_end = primer_start + primer_length - 1
        return (primer_start <= junction - 7) and (f_end >= junction + 6)

# --- Page Setup ---
st.set_page_config(page_title="SPUD - Batch Mode", page_icon="🧬", layout="wide")

# --- Documentation Sidebar ---
with st.sidebar:
    st.title("📖 SPUD Documentation")
    st.markdown("Welcome to **SPUD** (Specific Primer Universal Designer) - **Batch Edition**.")
    
    with st.expander("🎯 1. Smart Target Input"):
        st.write("""
        The **Target** field automatically detects your design strategy:
        * **Single Number (e.g., 452):** Triggers **Junction Mode**. SPUD ensures the primer spans an exon-exon junction with a strict **7-bp anchor** on each side.
        * **Range (e.g., 300-450):** Triggers **ROI Mode**. Forces the entire amplicon to be designed **inside** this specific sequence range.
        * **Empty:** Triggers **Global Mode**. Searches the entire sequence for the best thermodynamic pairs.
        """)

    with st.expander("📦 2. Batch Processing"):
        st.write("""
        1. **Download Template:** Use the provided CSV to format your data.
        2. **Upload:** SPUD will process every row independently without crashing on errors.
        3. **BLAST Strategy:** To save time and avoid NCBI blocking, BLAST is only performed on the **Rank 1** candidate for each gene.
        4. **Quality Flagging:** If a primer's Tm deviates by >1.5°C from the batch average, it will be flagged as an **Outlier**.
        """)

    with st.expander("📊 3. Output Parameters"):
        st.write("""
        * **Rank:** Priority based on penalty score.
        * **Tm_Flag:** Marks primers that may not fit the global plate PCR conditions.
        * **Template_Fold:** (v / !) Assesses if the DNA template forms secondary structures that compete with the primer.
        * **Status:** Shows BLAST results for Rank 1 or error messages for failed designs.
        """)

# --- Main Header ---
st.title("🧬 SPUD: Specific Primer Universal Designer (Batch Edition)")
st.divider()

# --- Top Settings (Global Conditions) ---
with st.expander("⚙️ Lab Conditions & Advanced Thermodynamics (Applied to all)", expanded=True):
    col1, col2, col3, col4 = st.columns(4)
    with col1: 
        species_name = st.text_input("Global Species (for BLAST):", "Solanum tuberosum")
        num_returns = st.slider("Max Candidates per Gene", 1, 10, 3, step=1)
    with col2: 
        primer_conc = st.number_input("Primer Conc. (nM)", value=250.0, step=10.0)
        mg_conc = st.number_input("Mg2+ Conc. (mM)", value=2.5, step=0.1)
    with col3: 
        target_tm = st.slider("Target Tm (°C)", 50.0, 72.0, 60.0, step=0.5)
    with col4: 
        min_amp = st.text_input("Min Amplicon Length", value="80")
        max_amp = st.text_input("Max Amplicon Length", value="150")
        run_blast_check = st.checkbox("🔍 Run BLAST on Rank 1", value=True)

# --- Mode Selection ---
mode = st.radio("Select Processing Mode:", ["Batch Upload (CSV/Excel)", "Single Gene (Quick Test)"], horizontal=True)

genes_data = []

if mode == "Batch Upload (CSV/Excel)":
    st.info("Upload a file with columns: **Gene_ID**, **Sequence**, **Target**.")
    
    # Template Download Generator
    template_df = pd.DataFrame({
        "Gene_ID": ["StPOT1", "StGA2ox", "StActin"], 
        "Sequence": ["ATGC...", "TTGC...", "CCGG..."], 
        "Target": ["452", "300-450", ""]
    })
    csv_buffer = io.BytesIO()
    template_df.to_csv(csv_buffer, index=False)
    st.download_button(label="📥 Download Template", data=csv_buffer.getvalue(), file_name="SPUD_Template.csv", mime="text/csv")
    
    uploaded_file = st.file_uploader("Upload Batch File", type=['csv', 'xlsx'])
    if uploaded_file:
        try:
            if uploaded_file.name.endswith('.csv'): 
                df_input = pd.read_csv(uploaded_file)
            else: 
                df_input = pd.read_excel(uploaded_file)
                
            for index, row in df_input.iterrows():
                genes_data.append({
                    "Gene_ID": str(row.get("Gene_ID", f"Gene_{index}")), 
                    "Sequence": str(row.get("Sequence", "")), 
                    "Target": str(row.get("Target", ""))
                })
        except Exception as e: 
            st.error(f"Error reading file: {e}")

else:
    # Single Gene Input UI
    col_a, col_b = st.columns([2, 1])
    with col_a: 
        s_seq = st.text_area("Target Sequence (5' to 3'):", height=150)
    with col_b:
        s_id = st.text_input("Gene ID:", "MyGene")
        s_target = st.text_input("Target (Smart Input):", placeholder="e.g., 452 OR 300-450")
    if s_seq: 
        genes_data.append({"Gene_ID": s_id, "Sequence": s_seq, "Target": s_target})

# --- Execution Engine ---
if st.button("🚀 Run SPUD Engine", type="primary") and genes_data:
    all_results = []
    
    global_args = {
        'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': target_tm, 'PRIMER_MIN_TM': target_tm - 5.0, 'PRIMER_MAX_TM': target_tm + 5.0,
        'PRIMER_DNA_CONC': primer_conc, 'PRIMER_SALT_DIVALENT': mg_conc, 
        'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_DNTP_CONC': 0.8,
        'PRIMER_TM_FORMULA': 1, 'PRIMER_SALT_CORRECTIONS': 1, 'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
        'PRIMER_NUM_RETURN': max(20, num_returns * 5)
    }
    
    if min_amp.isdigit() and max_amp.isdigit(): 
        global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp), int(max_amp)]]

    progress_bar = st.progress(0)
    status_text = st.empty()

    for idx, gene in enumerate(genes_data):
        status_text.text(f"Processing {gene['Gene_ID']} ({idx+1}/{len(genes_data)})...")
        clean_seq = "".join(gene['Sequence'].split()).upper()
        
        if len(clean_seq) < 50:
            all_results.append({"Gene_ID": gene['Gene_ID'], "Status": "Error: Sequence too short"})
            continue

        target_type, target_val = parse_target(gene['Target'])
        seq_args = {'SEQUENCE_ID': gene['Gene_ID'], 'SEQUENCE_TEMPLATE': clean_seq}
        
        if target_type == "Junction":
            if target_val < 20 or target_val > len(clean_seq) - 20:
                all_results.append({"Gene_ID": gene['Gene_ID'], "Status": "Error: Junction too close to edge"})
                continue
            seq_args['SEQUENCE_OVERLAP_JUNCTION_LIST'] = [target_val]
            
        elif target_type == "ROI":
            roi_start = max(0, target_val[0] - 1)
            roi_len = target_val[1] - target_val[0] + 1
            if roi_len > 0 and roi_start + roi_len <= len(clean_seq):
                seq_args['SEQUENCE_INCLUDED_REGION'] = [roi_start, roi_len]
            else:
                all_results.append({"Gene_ID": gene['Gene_ID'], "Status": "Error: Invalid ROI range"})
                continue

        try:
            raw_res = primer3.bindings.designPrimers(seq_args, global_args)
            valid_candidates = []
            
            for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                f_start = raw_res.get(f'PRIMER_LEFT_{i}')[0]
                f_len = raw_res.get(f'PRIMER_LEFT_{i}')[1]
                r_start = raw_res.get(f'PRIMER_RIGHT_{i}')[0]
                r_len = raw_res.get(f'PRIMER_RIGHT_{i}')[1]
                f_seq = raw_res.get(f'PRIMER_LEFT_{i}_SEQUENCE')
                r_seq = raw_res.get(f'PRIMER_RIGHT_{i}_SEQUENCE')
                f_tm = raw_res.get(f'PRIMER_LEFT_{i}_TM')
                r_tm = raw_res.get(f'PRIMER_RIGHT_{i}_TM')
                penalty = raw_res.get(f'PRIMER_PAIR_{i}_PENALTY')
                
                # Check 7bp rule if in Junction mode
                if target_type == "Junction":
                    if not (validate_7bp_anchor(f_start, f_len, target_val, False) or validate_7bp_anchor(r_start, r_len, target_val, True)):
                        continue
                
                # Secondary Structure Penalty
                f_temp_tm = primer3.calc_hairpin(rev_comp(f_seq)).tm
                r_temp_tm = primer3.calc_hairpin(rev_comp(r_seq)).tm
                sec_struct_flag = "!" if (f_temp_tm >= f_tm - 3.0 or r_temp_tm >= r_tm - 3.0) else "v"
                if sec_struct_flag == "!":
                    penalty += 50.0
                
                valid_candidates.append({
                    "Gene_ID": gene['Gene_ID'], 
                    "Forward_Seq": f_seq, 
                    "Reverse_Seq": r_seq,
                    "F_Tm": round(f_tm, 1), 
                    "R_Tm": round(r_tm, 1),
                    "Amp_Len": raw_res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE'),
                    "Template_Fold": sec_struct_flag, 
                    "Penalty": penalty
                })
            
            valid_candidates.sort(key=lambda x: x['Penalty'])
            
            for rank, cand in enumerate(valid_candidates[:num_returns]):
                cand['Rank'] = rank + 1
                cand['Status'] = "Pending"
                all_results.append(cand)
                
            if not valid_candidates: 
                all_results.append({"Gene_ID": gene['Gene_ID'], "Status": "No candidates found under constraints"})
                
        except Exception as e: 
            all_results.append({"Gene_ID": gene['Gene_ID'], "Status": f"Error: {e}"})
            
        progress_bar.progress((idx + 1) / len(genes_data))

    # --- Final Audit & Formatting ---
    status_text.text("Finalizing Quality Checks & BLAST...")
    
    successful_tms = [c['F_Tm'] for c in all_results if 'F_Tm' in c] + [c['R_Tm'] for c in all_results if 'R_Tm' in c]
    global_mean_tm = sum(successful_tms) / len(successful_tms) if successful_tms else target_tm
    
    for cand in all_results:
        if cand.get('Status') == "Pending":
            mean_pair_tm = (cand['F_Tm'] + cand['R_Tm']) / 2.0
            cand['Tm_Flag'] = "Outlier" if abs(mean_pair_tm - global_mean_tm) > 1.5 else "OK"
            
            if run_blast_check and cand['Rank'] == 1:
                status_text.text(f"BLASTing {cand['Gene_ID']}...")
                cand['Status'] = check_and_get_blast(cand['Forward_Seq'], species_name)
            else: 
                cand['Status'] = "Ready / Skipped BLAST"

    status_text.text("Complete!")
    
    # Organize columns explicitly
    df_results = pd.DataFrame(all_results)
    desired_cols = ['Gene_ID', 'Rank', 'Forward_Seq', 'Reverse_Seq', 'F_Tm', 'R_Tm', 'Tm_Flag', 'Amp_Len', 'Template_Fold', 'Penalty', 'Status']
    final_cols = [c for c in desired_cols if c in df_results.columns]
    
    # Add any leftover columns (like error messages)
    for c in df_results.columns:
        if c not in final_cols:
            final_cols.append(c)
            
    df_results = df_results[final_cols]
    
    st.dataframe(df_results, use_container_width=True)
    st.download_button(
        label="📥 Download Results (CSV)", 
        data=df_results.to_csv(index=False).encode('utf-8'), 
        file_name="SPUD_Results.csv", 
        mime="text/csv",
        type="primary"
    )
