import streamlit as st
import primer3
import pandas as pd
import time
import re
import plotly.graph_objects as go
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# --- פונקציית Reverse Complement לחישוב Template ---
def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

# --- פונקציית BLAST ---
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
        
        if significant_hits == 0:
            summary = "0 Hits ✅"
        elif significant_hits == 1:
            summary = "1 Hit ✅"
        else:
            summary = f"{significant_hits} Hits"
            
        return summary, detailed_text
    except Exception as e:
        return "Error", f"BLAST Error: {str(e)}"

# --- פונקציות עיצוב לטבלה ---
def color_negative_red(val, threshold):
    try:
        if float(val) < threshold:
            return 'background-color: #ffcccc; color: red; font-weight: bold;'
    except:
        pass
    return ''

def color_blast_red(val):
    try:
        hits = int(str(val).split()[0])
        if hits > 1:
            return 'background-color: #ffcccc; color: red; font-weight: bold;'
    except:
        pass
    return ''

def color_sec_struct(val):
    if val == '!':
        return 'background-color: #ffcccc; color: red; font-weight: bold; text-align: center;'
    elif val == 'v':
        return 'color: green; font-weight: bold; text-align: center;'
    return ''

# --- הגדרות עמוד ---
st.set_page_config(page_title="SPUD - Primer Designer", page_icon="🧬", layout="wide")
st.title("🧬 SPUD: Specific Primer Universal Designer")
st.markdown("""
**SPUD** is a high-precision tool for plant scientists, calibrated for **Fast SYBR® Green** reactions. 
It features smart exon-junction distribution, real-time **$\Delta G$** screening, template UNAFold-like risk assessment, 
and automated **Locus Tag** identification for any plant species via NCBI BLAST.
""")
st.divider()

# --- ממשק ---
col1, col2 = st.columns([2, 1])

with col1:
    sequence = st.text_area("Target Sequence (5' to 3'):", height=200, placeholder="Paste your FASTA sequence...")

with col2:
    project_name = st.text_input("Primer Project Name:", "My_Gene_Cloning")
    species_name = st.text_input("Species (for BLAST):", "Solanum tuberosum")
    junction_pos = st.text_input("Exon-Exon Junctions (comma separated):", placeholder="e.g., 72, 478, 727")
    st.markdown("<br>", unsafe_allow_html=True)
    run_blast_check = st.checkbox("🔍 Run NCBI BLAST Specificity Check (Gatekeeper)")

with st.expander("⚙️ Advanced Thermodynamics (Fast SYBR Green Master Mix)"):
    adv_col1, adv_col2, adv_col3, adv_col4 = st.columns(4)
    with adv_col1:
        primer_conc = st.number_input("Primer Conc. (nM)", value=250.0, step=10.0)
    with adv_col2:
        mg_conc = st.number_input("Mg2+ Conc. (mM)", value=2.5, step=0.1)
    with adv_col3:
        target_tm = st.slider("Target Tm (°C)", min_value=50.0, max_value=72.0, value=60.0, step=0.5)
    with adv_col4:
        num_returns = st.slider("Total Primers to Design", min_value=1, max_value=15, value=5, step=1)
        
    st.markdown("<br>", unsafe_allow_html=True)
    amp_col1, amp_col2, _, _ = st.columns(4)
    with amp_col1:
        min_amp = st.text_input("Min Amplicon Length", value="80")
    with amp_col2:
        max_amp = st.text_input("Max Amplicon Length", value="150")

# --- הרצה ---
if st.button("🚀 Design Primers", type="primary"):
    if not sequence:
        st.error("Please enter a DNA sequence to proceed.")
    else:
        with st.spinner('Simulating Reaction, Calculating Folds & Designing Primers...'):
            try:
                clean_seq = "".join(sequence.split()).upper()
                junction_list = [int(x.strip()) for x in junction_pos.split(',') if x.strip().isdigit()] if junction_pos else []

                global_args = {
                    'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': target_tm, 'PRIMER_MIN_TM': target_tm - 5.0, 'PRIMER_MAX_TM': target_tm + 5.0,
                    'PRIMER_DNA_CONC': primer_conc, 'PRIMER_SALT_DIVALENT': mg_conc, 
                    'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_DNTP_CONC': 0.8,
                    'PRIMER_TM_FORMULA': 1, 'PRIMER_SALT_CORRECTIONS': 1, 'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
                    'PRIMER_NUM_RETURN': max(15, num_returns * 2) # יוצר מאגר גדול לסינון חכם
                }

                if min_amp.strip().isdigit() and max_amp.strip().isdigit(): 
                    global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp.strip()), int(max_amp.strip())]]
                elif min_amp.strip().isdigit(): 
                    global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp.strip()), len(clean_seq)]]
                elif max_amp.strip().isdigit(): 
                    global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[50, int(max_amp.strip())]]

                raw_candidates = []

                if junction_list:
                    all_junction_results = []
                    for j in junction_list:
                        seq_args = {'SEQUENCE_ID': f"J_{j}", 'SEQUENCE_TEMPLATE': clean_seq, 'SEQUENCE_OVERLAP_JUNCTION_LIST': [j]}
                        raw_res = primer3.bindings.designPrimers(seq_args, global_args)
                        j_cands = []
                        for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                            j_cands.append({
                                'penalty': raw_res.get(f'PRIMER_PAIR_{i}_PENALTY', 999.0),
                                'left_start': raw_res.get(f'PRIMER_LEFT_{i}')[0],
                                'right_start': raw_res.get(f'PRIMER_RIGHT_{i}')[0],
                                'forward_seq': raw_res.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
                                'reverse_seq': raw_res.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
                                'forward_tm': raw_res.get(f'PRIMER_LEFT_{i}_TM', 0),
                                'reverse_tm': raw_res.get(f'PRIMER_RIGHT_{i}_TM', 0),
                                'product_size': raw_res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0),
                                'junction_used': j
                            })
                        all_junction_results.append(j_cands)

                    max_cands = max([len(lst) for lst in all_junction_results]) if all_junction_results else 0
                    for i in range(max_cands):
                        for j_cands in all_junction_results:
                            if i < len(j_cands):
                                raw_candidates.append(j_cands[i])

                else:
                    seq_args = {'SEQUENCE_ID': project_name, 'SEQUENCE_TEMPLATE': clean_seq}
                    raw_res = primer3.bindings.designPrimers(seq_args, global_args)
                    for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                        raw_candidates.append({
                            'penalty': raw_res.get(f'PRIMER_PAIR_{i}_PENALTY', 999.0),
                            'left_start': raw_res.get(f'PRIMER_LEFT_{i}')[0],
                            'right_start': raw_res.get(f'PRIMER_RIGHT_{i}')[0],
                            'forward_seq': raw_res.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
                            'reverse_seq': raw_res.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
                            'forward_tm': raw_res.get(f'PRIMER_LEFT_{i}_TM', 0),
                            'reverse_tm': raw_res.get(f'PRIMER_RIGHT_{i}_TM', 0),
                            'product_size': raw_res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0),
                            'junction_used': "None"
                        })

                if not raw_candidates: 
                    st.warning("⚠️ No primers found. Try relaxing the constraints.")
                else:
                    # שלב 1: חישוב Template Secondary Structure Penalty
                    for cand in raw_candidates:
                        f_seq = cand['forward_seq']
                        r_seq = cand['reverse_seq']
                        
                        f_temp_tm = primer3.calc_hairpin(rev_comp(f_seq)).tm
                        r_temp_tm = primer3.calc_hairpin(rev_comp(r_seq)).tm
                        
                        # חוק ה-3 מעלות: האם הקיפול קרוב ל-Tm של הפריימר
                        risk_f = f_temp_tm >= (cand['forward_tm'] - 3.0)
                        risk_r = r_temp_tm >= (cand['reverse_tm'] - 3.0)
                        
                        if risk_f or risk_r:
                            cand['Template_SecondaryStracture'] = "!"
                            cand['penalty'] += 50.0 # עונש תרמודינמי כבד
                        else:
                            cand['Template_SecondaryStracture'] = "v"

                    # ממיין ראשונית לפי הציון התרמודינמי החדש וחותך לכמות סבירה לבדיקה
                    raw_candidates.sort(key=lambda x: x['penalty'])
                    candidates_to_process = raw_candidates[:max(10, num_returns * 2)]
                    
                    blast_details = []

                    # שלב 2: BLAST + מיון דו-שלבי
                    if run_blast_check:
                        for idx, cand in enumerate(candidates_to_process):
                            st.toast(f"Running BLAST for Candidate {idx+1}/{len(candidates_to_process)}...")
                            
                            f_blast_sum, f_blast_det = check_and_get_blast(cand['forward_seq'], species_name)
                            time.sleep(1)
                            r_blast_sum, r_blast_det = check_and_get_blast(cand['reverse_seq'], species_name)
                            time.sleep(1)
                            
                            cand['f_blast_sum'] = f_blast_sum
                            cand['f_blast_det'] = f_blast_det
                            cand['r_blast_sum'] = r_blast_sum
                            cand['r_blast_det'] = r_blast_det
                            
                            try:
                                f_hits = int(f_blast_sum.split()[0])
                            except: f_hits = 99
                            try:
                                r_hits = int(r_blast_sum.split()[0])
                            except: r_hits = 99
                            
                            # סימון מועמדים פסולי BLAST (>1)
                            cand['blast_bad_flag'] = 1 if (f_hits > 1 or r_hits > 1) else 0
                        
                        # מיון סופי: קודם פוסל את אלו עם >1 ב-BLAST, ואז לפי Penalty תרמודינמי
                        candidates_to_process.sort(key=lambda x: (x['blast_bad_flag'], x['penalty']))
                        final_candidates = candidates_to_process[:num_returns]
                        
                        for cand in final_candidates:
                            blast_details.append({"F": cand['f_blast_det'], "R": cand['r_blast_det']})

                    else:
                        # אם דילגנו על BLAST
                        final_candidates = candidates_to_process[:num_returns]
                        for cand in final_candidates:
                            cand['f_blast_sum'] = "Skipped"
                            cand['r_blast_sum'] = "Skipped"

                    # שלב 3: הכנת הטבלה הסופית
                    parsed_data = []
                    for idx, cand in enumerate(final_candidates):
                        f_seq = cand['forward_seq']
                        r_seq = cand['reverse_seq']
                        
                        parsed_data.append({
                            "Rank": idx + 1,
                            "Penalty": cand['penalty'],
                            "Template_SecondaryStracture": cand['Template_SecondaryStracture'],
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
                            "Amp_Length": cand['product_size']
                        })
                    
                    df = pd.DataFrame(parsed_data)
                    st.success("✅ Analysis Complete!")

                    # --- המפה הויזואלית ---
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

                    st.markdown("### Detailed Results Table")
                    
                    try:
                        styled_df = df.style.format(precision=1)
                        if hasattr(styled_df, "map"):
                            styled_df = styled_df.map(lambda x: color_negative_red(x, -5.0), subset=['F_Self (ΔG)', 'R_Self (ΔG)', 'CrossDimer (ΔG)']) \
                                                 .map(lambda x: color_negative_red(x, -3.0), subset=['F_Hairpin (ΔG)', 'R_Hairpin (ΔG)']) \
                                                 .map(lambda x: color_blast_red(x), subset=['F_BLAST', 'R_BLAST']) \
                                                 .map(lambda x: color_sec_struct(x), subset=['Template_SecondaryStracture'])
                        else:
                            styled_df = styled_df.applymap(lambda x: color_negative_red(x, -5.0), subset=['F_Self (ΔG)', 'R_Self (ΔG)', 'CrossDimer (ΔG)']) \
                                                 .applymap(lambda x: color_negative_red(x, -3.0), subset=['F_Hairpin (ΔG)', 'R_Hairpin (ΔG)']) \
                                                 .applymap(lambda x: color_blast_red(x), subset=['F_BLAST', 'R_BLAST']) \
                                                 .applymap(lambda x: color_sec_struct(x), subset=['Template_SecondaryStracture'])
                    except:
                        styled_df = df

                    st.dataframe(styled_df, use_container_width=True)
                    st.download_button("📥 Download Results (CSV)", df.to_csv(index=False).encode('utf-8'), f"{project_name}.csv", "text/csv")

                    if run_blast_check:
                        st.markdown("### 🔍 BLAST Alignments")
                        for idx, detail in enumerate(blast_details):
                            with st.expander(f"View Alignments for Candidate {idx+1}"):
                                st.write("**Forward Primer:**")
                                st.code(detail["F"], language="text")
                                st.write("**Reverse Primer:**")
                                st.code(detail["R"], language="text")

            except Exception as e: st.error(f"Error executing analysis: {e}")
