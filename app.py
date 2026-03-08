import streamlit as st
import primer3
import pandas as pd
import time
import re
import plotly.graph_objects as go
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# --- פונקציות עזר ותרמודינמיקה ---
def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def check_and_get_blast(primer_seq, species):
    try:
        entrez_query = f"{species}[organism]" if species else ""
        handle = NCBIWWW.qblast("blastn", "nt", primer_seq, entrez_query=entrez_query, word_size=11, hitlist_size=5)
        blast_record = NCBIXML.read(handle)
        report = []
        significant_hits = 0
        primer_len = len(primer_seq)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.align_length >= primer_len * 0.9:
                    significant_hits += 1
                    tags = []
                    for word in alignment.title.split():
                        clean_word = word.strip("()[],:;")
                        if any(c.isalpha() for c in clean_word) and any(c.isdigit() for c in clean_word) and len(clean_word) > 4:
                            tags.append(clean_word)
                    unique_tags = list(dict.fromkeys(tags))[:3]
                    locus_text = ", ".join(unique_tags) if unique_tags else "N/A"
                    gene_desc = alignment.hit_def[:80] + "..." if len(alignment.hit_def) > 80 else alignment.hit_def
                    detail = (
                        f"Locus / Gene ID: {locus_text}\n"
                        f"Description: {gene_desc}\n"
                        f"NCBI Acc: {alignment.accession} | Identity: {hsp.identities}/{hsp.align_length}\n"
                        f"Query:  {hsp.query}\n"
                        f"Match:  {hsp.match}\n"
                        f"Sbjct:  {hsp.sbjct}\n"
                        "--------------------------------------------------"
                    )
                    report.append(detail)
        detailed_text = "\n".join(report) if report else "No significant hits found (>90% coverage)."
        summary = f"{significant_hits} Hits" if significant_hits > 1 else "1 Hit ✅" if significant_hits == 1 else "0 Hits ✅"
        return summary, detailed_text
    except Exception as e: return "Error", str(e)

# --- פונקציות עיצוב לטבלה ---
def color_negative_red(val, threshold):
    try:
        if float(val) < threshold: return 'background-color: #ffcccc; color: red; font-weight: bold;'
    except: pass
    return ''

def color_blast_red(val):
    try:
        hits = int(str(val).split()[0])
        if hits > 1: return 'background-color: #ffcccc; color: red; font-weight: bold;'
    except: pass
    return ''

def color_sec_struct(val):
    if val == '!': return 'background-color: #ffcccc; color: red; font-weight: bold; text-align: center;'
    elif val == 'v': return 'color: green; font-weight: bold; text-align: center;'
    return ''

# --- הגדרות עמוד ---
st.set_page_config(page_title="SPUD - Primer Designer", page_icon="🧬", layout="wide")

# --- תפריט הסברים בצד (Documentation Sidebar) ---
with st.sidebar:
    st.title("📖 SPUD Documentation")
    st.markdown("Welcome to **SPUD** (Specific Primer Universal Designer).")
    
    with st.expander("🔄 1. The SPUD Workflow"):
        st.write("""
        1. **Design & Distribute:** Generates raw primers, forcibly selecting candidates evenly across all provided Exon-Exon junctions (Round-Robin).
        2. **Thermodynamic Penalty:** Calculates standard Primer3 penalties + strict UNAFold-like penalties for template secondary structures.
        3. **BLAST Gatekeeper (Optional):** Screens the top candidates against NCBI. Any pair with >1 hit is pushed to the bottom of the list.
        4. **Final Sort:** Prioritizes specificity (BLAST) while strictly maintaining the diverse junction distribution.
        """)
        
    with st.expander("🎛️ 2. Input Parameters (What you can tweak)"):
        st.write("""
        * **Target Sequence:** Your full gene/cDNA FASTA sequence.
        * **Species:** Scientific name (e.g., *Solanum tuberosum*). Crucial for accurate BLAST results.
        * **Exon-Exon Junctions:** Base-pair indexes where exons meet. SPUD will span primers across these points.
        * **Primer / Mg2+ Conc:** Calibrated by default to Fast SYBR® Green. Adjust only if using a different Master Mix.
        * **Target Tm:** The ideal melting temperature for your reaction.
        * **Amplicon Length:** Keep between 80-150 for optimal qPCR efficiency.
        """)

    with st.expander("📊 3. Output Parameters (Reading the Table)"):
        st.write("""
        * **Penalty:** Objective score. Lower is better.
        * **ΔG (Gibbs Free Energy):** Measures structure stability in *kcal/mol*. 
          * *Self / Cross Dimers:* Safe if > -5.0. If lower (red), they may form primer-dimers.
          * *Hairpins:* Safe if > -3.0. If lower (red), the primer folds on itself.
        * **BLAST:** Target specificity. "1 Hit ✅" is ideal. >1 Hit (red) means off-target binding risk.
        * **Template_Fold (v / !):** Evaluates the DNA template (Amplicon). If marked with "!", the template forms a stable secondary structure within 3°C of your primer's Tm, risking amplification failure.
        """)

# --- כותרת ראשית ---
st.title("🧬 SPUD: Specific Primer Universal Designer")
st.markdown("""
**SPUD** is a high-precision tool for plant scientists, calibrated for **Fast SYBR® Green** reactions. 
It features smart exon-junction distribution, real-time **$\Delta G$** screening, template UNAFold-like risk assessment, 
and automated **Locus Tag** identification for any plant species via NCBI BLAST.
""")
st.divider()

# --- ממשק ---
col1, col2 = st.columns([2, 1])
with col1: sequence = st.text_area("Target Sequence (5' to 3'):", height=200, placeholder="Paste your FASTA sequence...")
with col2:
    project_name = st.text_input("Project Name:", "My_Gene_Cloning")
    species_name = st.text_input("Species (for BLAST):", "Solanum tuberosum")
    junction_pos = st.text_input("Exon-Exon Junctions:", placeholder="e.g., 72, 478")
    run_blast_check = st.checkbox("🔍 Run NCBI BLAST (Gatekeeper Mode)")

with st.expander("⚙️ Advanced Thermodynamics (Fast SYBR Green)"):
    adv_col1, adv_col2, adv_col3, adv_col4 = st.columns(4)
    with adv_col1: primer_conc = st.number_input("Primer Conc. (nM)", value=250.0, step=10.0)
    with adv_col2: mg_conc = st.number_input("Mg2+ Conc. (mM)", value=2.5, step=0.1)
    with adv_col3: target_tm = st.slider("Target Tm (°C)", 50.0, 72.0, 60.0, step=0.5)
    with adv_col4: num_returns = st.slider("Total Primers to Design", 1, 15, 5, step=1)
    
    st.markdown("<br>", unsafe_allow_html=True)
    amp_col1, amp_col2, _, _ = st.columns(4)
    with amp_col1: min_amp = st.text_input("Min Amplicon Length", value="80")
    with amp_col2: max_amp = st.text_input("Max Amplicon Length", value="150")

# --- לוגיקת הרצה ---
if st.button("🚀 Design Primers", type="primary"):
    if not sequence: st.error("Please enter a DNA sequence to proceed.")
    else:
        with st.spinner('Calculating Folds, Simulating Reactions & Distributing...'):
            try:
                clean_seq = "".join(sequence.split()).upper()
                junction_list = [int(x.strip()) for x in junction_pos.split(',') if x.strip().isdigit()] if junction_pos else []
                global_args = {
                    'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': target_tm, 'PRIMER_MIN_TM': target_tm - 5.0, 'PRIMER_MAX_TM': target_tm + 5.0,
                    'PRIMER_DNA_CONC': primer_conc, 'PRIMER_SALT_DIVALENT': mg_conc, 
                    'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_DNTP_CONC': 0.8,
                    'PRIMER_TM_FORMULA': 1, 'PRIMER_SALT_CORRECTIONS': 1, 'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
                    'PRIMER_NUM_RETURN': max(15, num_returns * 3) # מייצר פול רחב לסינון
                }
                if min_amp.isdigit() and max_amp.isdigit(): global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp), int(max_amp)]]

                pools = [] # רשימה של רשימות (כל פריט הוא רשימת מועמדים מאקסון מסוים)
                
                if junction_list:
                    for j in junction_list:
                        seq_args = {'SEQUENCE_ID': f"J_{j}", 'SEQUENCE_TEMPLATE': clean_seq, 'SEQUENCE_OVERLAP_JUNCTION_LIST': [j]}
                        raw_res = primer3.bindings.designPrimers(seq_args, global_args)
                        j_cands = []
                        for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                            j_cands.append({
                                'penalty': raw_res.get(f'PRIMER_PAIR_{i}_PENALTY', 999.0),
                                'left_start': raw_res.get(f'PRIMER_LEFT_{i}')[0], 'right_start': raw_res.get(f'PRIMER_RIGHT_{i}')[0],
                                'forward_seq': raw_res.get(f'PRIMER_LEFT_{i}_SEQUENCE'), 'reverse_seq': raw_res.get(f'PRIMER_RIGHT_{i}_SEQUENCE'),
                                'forward_tm': raw_res.get(f'PRIMER_LEFT_{i}_TM'), 'reverse_tm': raw_res.get(f'PRIMER_RIGHT_{i}_TM'),
                                'product_size': raw_res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE'), 'junction_used': j
                            })
                        if j_cands: pools.append(j_cands)
                else:
                    seq_args = {'SEQUENCE_ID': project_name, 'SEQUENCE_TEMPLATE': clean_seq}
                    raw_res = primer3.bindings.designPrimers(seq_args, global_args)
                    cands = []
                    for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                        cands.append({
                            'penalty': raw_res.get(f'PRIMER_PAIR_{i}_PENALTY', 999.0),
                            'left_start': raw_res.get(f'PRIMER_LEFT_{i}')[0], 'right_start': raw_res.get(f'PRIMER_RIGHT_{i}')[0],
                            'forward_seq': raw_res.get(f'PRIMER_LEFT_{i}_SEQUENCE'), 'reverse_seq': raw_res.get(f'PRIMER_RIGHT_{i}_SEQUENCE'),
                            'forward_tm': raw_res.get(f'PRIMER_LEFT_{i}_TM'), 'reverse_tm': raw_res.get(f'PRIMER_RIGHT_{i}_TM'),
                            'product_size': raw_res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE'), 'junction_used': "None"
                        })
                    if cands: pools.append(cands)

                if not pools: 
                    st.warning("⚠️ No primers found. Try relaxing the constraints.")
                else:
                    # שלב 1: חישוב Template Penalty לכל המועמדים בכל האקסונים
                    for pool in pools:
                        for cand in pool:
                            f_temp_tm = primer3.calc_hairpin(rev_comp(cand['forward_seq'])).tm
                            r_temp_tm = primer3.calc_hairpin(rev_comp(cand['reverse_seq'])).tm
                            if f_temp_tm >= (cand['forward_tm'] - 3.0) or r_temp_tm >= (cand['reverse_tm'] - 3.0):
                                cand['Template_SecStruct'] = "!"
                                cand['penalty'] += 50.0 
                            else: 
                                cand['Template_SecStruct'] = "v"
                        # מיון פנימי של כל אקסון לפי הציון החדש
                        pool.sort(key=lambda x: x['penalty'])

                    # שלב 2: שליפה בשיטת Round-Robin שומרת על חלוקה שווה בין האקסונים
                    candidates_to_process = []
                    max_cands = max(len(p) for p in pools)
                    rr_counter = 0
                    for i in range(max_cands):
                        for pool in pools:
                            if i < len(pool):
                                cand = pool[i]
                                cand['rr_rank'] = rr_counter # שומר את המיקום היחסי כדי לא לאבד את הפיזור
                                candidates_to_process.append(cand)
                                rr_counter += 1
                    
                    # חותך לכמות סבירה כדי שה-BLAST לא ייקח נצח
                    candidates_to_process = candidates_to_process[:max(10, num_returns * 2)]
                    blast_details = []

                    # שלב 3: BLAST Gatekeeper
                    if run_blast_check:
                        for idx, cand in enumerate(candidates_to_process):
                            st.toast(f"BLASTing Candidate {idx+1}/{len(candidates_to_process)}...")
                            f_sum, f_det = check_and_get_blast(cand['forward_seq'], species_name)
                            r_sum, r_det = check_and_get_blast(cand['reverse_seq'], species_name)
                            cand['f_blast_sum'], cand['f_blast_det'] = f_sum, f_det
                            cand['r_blast_sum'], cand['r_blast_det'] = r_sum, r_det
                            
                            try: f_hits = int(f_sum.split()[0])
                            except: f_hits = 99
                            try: r_hits = int(r_sum.split()[0])
                            except: r_hits = 99
                            
                            cand['blast_bad'] = 1 if (f_hits > 1 or r_hits > 1) else 0
                        
                        # שלב 4: מיון סופי - קודם שומר סף BLAST, ואז שומר על סדר ה-Round-Robin
                        candidates_to_process.sort(key=lambda x: (x.get('blast_bad', 0), x['rr_rank']))
                        final_candidates = candidates_to_process[:num_returns]
                        
                        for cand in final_candidates:
                            blast_details.append({"F": cand['f_blast_det'], "R": cand['r_blast_det']})

                    else:
                        # אם ה-BLAST כבוי, פשוט לוקח את ה-Round Robin הנקי
                        final_candidates = candidates_to_process[:num_returns]
                        for cand in final_candidates:
                            cand['f_blast_sum'] = "Skipped"
                            cand['r_blast_sum'] = "Skipped"

                    # שלב 5: בניית הטבלה הסופית (העברת Template_Fold לסוף)
                    parsed_data = []
                    for idx, cand in enumerate(final_candidates):
                        f_seq = cand['forward_seq']
                        r_seq = cand['reverse_seq']
                        
                        parsed_data.append({
                            "Rank": idx + 1,
                            "Penalty": cand['penalty'],
                            "Junction": cand['junction_used'],
                            "Forward_pos": f"F{cand['left_start']}",
                            "Forward_Seq": f_seq,
                            "F_Self (ΔG)": primer3.calc_homodimer(f_seq).dg / 1000.0,
                            "F_Hairpin (ΔG)": primer3.calc_hairpin(f_seq).dg / 1000.0,
                            "F_BLAST": cand['f_blast_sum'],
                            "F_Tm (°C)": cand['forward_tm'],
                            "Reverse_pos": f"R{cand['right_start']}",
                            "Reverse_Seq": r_seq,
                            "R_Self (ΔG)": primer3.calc_homodimer(r_seq).dg / 1000.0,
                            "R_Hairpin (ΔG)": primer3.calc_hairpin(r_seq).dg / 1000.0,
                            "R_BLAST": cand['r_blast_sum'],
                            "R_Tm (°C)": cand['reverse_tm'],
                            "CrossDimer (ΔG)": primer3.calc_heterodimer(f_seq, r_seq).dg / 1000.0,
                            "Amp_Length": cand['product_size'],
                            "Template_Fold": cand['Template_SecStruct'] # <--- הועבר לסוף!
                        })
                    
                    df = pd.DataFrame(parsed_data)
                    st.success("✅ Analysis Complete!")

                    # --- תצוגת מפה ויזואלית ---
                    st.subheader("🗺️ Amplicon Map")
                    fig = go.Figure()
                    fig.add_shape(type="rect", x0=0, y0=0, x1=len(clean_seq), y1=1, line=dict(color="gray", width=2), fillcolor="lightgray")

                    for idx, cand in enumerate(final_candidates):
                        y_center = 1.5 + (idx * 0.6)
                        fig.add_shape(type="rect", x0=cand['left_start'], y0=y_center-0.2, x1=cand['right_start'], y1=y_center+0.2, 
                                      line=dict(color="DarkGreen", width=1), fillcolor="LimeGreen", opacity=0.8)
                        fig.add_annotation(x=(cand['left_start'] + cand['right_start'])/2, y=y_center, 
                                           text=f"Candidate {idx+1}", showarrow=False, font=dict(color="black", size=11))

                    max_y = 1.5 + (len(final_candidates) * 0.6)
                    if junction_list:
                        for j in junction_list:
                            fig.add_shape(type="line", x0=j, y0=-0.5, x1=j, y1=max_y, line=dict(color="red", width=2, dash="dash"))
                            fig.add_annotation(x=j, y=max_y + 0.2, text=f"Exon {j}", showarrow=False, font=dict(color="red", size=10))

                    fig.update_layout(xaxis=dict(range=[-20, len(clean_seq)+20], title="Position (bp)"), 
                                      yaxis=dict(showticklabels=False, range=[-0.5, max_y + 0.5]), 
                                      height=250 + (len(final_candidates) * 30), margin=dict(l=20, r=20, t=30, b=30), plot_bgcolor="white")
                    st.plotly_chart(fig, use_container_width=True)

                    # --- טבלה סופית מעוצבת ---
                    st.markdown("### Detailed Results Table")
                    try:
                        styled_df = df.style.format(precision=1)
                        if hasattr(styled_df, "map"):
                            styled_df = styled_df.map(lambda x: color_negative_red(x, -5.0), subset=['F_Self (ΔG)', 'R_Self (ΔG)', 'CrossDimer (ΔG)']) \
                                                 .map(lambda x: color_negative_red(x, -3.0), subset=['F_Hairpin (ΔG)', 'R_Hairpin (ΔG)']) \
                                                 .map(color_blast_red, subset=['F_BLAST', 'R_BLAST']) \
                                                 .map(color_sec_struct, subset=['Template_Fold'])
                        else:
                            styled_df = styled_df.applymap(lambda x: color_negative_red(x, -5.0), subset=['F_Self (ΔG)', 'R_Self (ΔG)', 'CrossDimer (ΔG)']) \
                                                 .applymap(lambda x: color_negative_red(x, -3.0), subset=['F_Hairpin (ΔG)', 'R_Hairpin (ΔG)']) \
                                                 .applymap(color_blast_red, subset=['F_BLAST', 'R_BLAST']) \
                                                 .applymap(color_sec_struct, subset=['Template_Fold'])
                    except:
                        styled_df = df

                    st.dataframe(styled_df, use_container_width=True)
                    st.download_button("📥 Download Results (CSV)", df.to_csv(index=False).encode('utf-8'), f"{project_name}.csv", "text/csv")

                    # --- הצגת ה-BLAST המפורט ---
                    if run_blast_check:
                        st.markdown("### 🔍 BLAST Alignments")
                        for idx, detail in enumerate(blast_details):
                            with st.expander(f"View Alignments for Candidate {idx+1}"):
                                st.write("**Forward Primer:**")
                                st.code(detail["F"], language="text")
                                st.write("**Reverse Primer:**")
                                st.code(detail["R"], language="text")

            except Exception as e: st.error(f"Error executing analysis: {e}")
