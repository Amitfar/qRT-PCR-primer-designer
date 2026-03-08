import streamlit as st
import primer3
import pandas as pd
import time
import plotly.graph_objects as go
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# --- פונקציית BLAST מעודכנת (מציגה שם ביולוגי אמיתי) ---
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
                    # שולף את השם הביולוגי האמיתי של הגן במקום רק מספר!
                    gene_desc = alignment.hit_def[:90] + "..." if len(alignment.hit_def) > 90 else alignment.hit_def
                    
                    detail = (
                        f"Gene: {gene_desc}\n"
                        f"Accession: {alignment.accession} | Identity: {hsp.identities}/{hsp.align_length}\n"
                        f"Query:  {hsp.query}\n"
                        f"Match:  {hsp.match}\n"
                        f"Sbjct:  {hsp.sbjct}\n"
                        "--------------------------------------------------"
                    )
                    report.append(detail)
        
        detailed_text = "\n".join(report) if report else "No significant hits found (>90% coverage)."
        summary = f"{significant_hits} Hits" if significant_hits > 0 else "0 Hits ✅"
        return summary, detailed_text
    except Exception as e:
        return "Error", f"BLAST Error: {str(e)}"

# --- הגדרות עיצוב לטבלה ---
def color_negative_red(val, threshold):
    try:
        if float(val) < threshold:
            return 'background-color: #ffcccc; color: red; font-weight: bold;'
    except:
        pass
    return ''

st.set_page_config(page_title="Cloning Primer Designer", page_icon="🧬", layout="wide")
st.title("🧬 Cloning Primer Designer (Volcani Edition)")
st.markdown("תכנון פריימרים להשתלה עם פיזור חכם, סימולציית קיט **Fast SYBR**, ופרמטרים תרמודינמיים מדויקים.")
st.divider()

col1, col2 = st.columns([2, 1])

with col1:
    sequence = st.text_area("Target Sequence (5' to 3'):", height=200, placeholder="Paste your FASTA sequence...")

with col2:
    project_name = st.text_input("Primer Project Name:", "My_Gene_Cloning")
    species_name = st.text_input("Species (for BLAST):", "Solanum tuberosum")
    junction_pos = st.text_input("Exon-Exon Junctions (comma separated):", placeholder="e.g., 72, 478, 727")
    st.markdown("<br>", unsafe_allow_html=True)
    run_blast_check = st.checkbox("🔍 Run NCBI BLAST Specificity Check (Takes 1-3 minutes)")

with st.expander("⚙️ Advanced Thermodynamics (Fast SYBR Green Master Mix)"):
    st.info("💡 **Kit Simulation:** The thermodynamics parameters are explicitly calibrated for Fast SYBR® Green Master Mix AB-4385612.")
    adv_col1, adv_col2, adv_col3, adv_col4 = st.columns(4)
    with adv_col1:
        primer_conc = st.number_input("Primer Conc. (nM)", value=250.0, step=10.0)
    with adv_col2:
        mg_conc = st.number_input("Mg2+ Conc. (mM)", value=3.0, step=0.1)
    with adv_col3:
        target_tm = st.slider("Target Tm (°C)", min_value=50.0, max_value=72.0, value=60.0, step=0.5)
    with adv_col4:
        num_returns = st.slider("Total Primers to Design", min_value=1, max_value=15, value=5, step=1)
        
    st.markdown("<br>", unsafe_allow_html=True)
    amp_col1, amp_col2, _, _ = st.columns(4)
    with amp_col1:
        min_amp = st.text_input("Min Amplicon Length", placeholder="e.g., 100")
    with amp_col2:
        max_amp = st.text_input("Max Amplicon Length", placeholder="e.g., 1500")

if st.button("🚀 Design Primers", type="primary"):
    if not sequence:
        st.error("Please enter a DNA sequence to proceed.")
    else:
        with st.spinner('Simulating Reaction and Designing Primers...'):
            try:
                clean_seq = "".join(sequence.split()).upper()
                junction_list = [int(x.strip()) for x in junction_pos.split(',') if x.strip().isdigit()] if junction_pos else []

                global_args = {
                    'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': target_tm, 'PRIMER_MIN_TM': target_tm - 5.0, 'PRIMER_MAX_TM': target_tm + 5.0,
                    'PRIMER_DNA_CONC': primer_conc, 'PRIMER_SALT_DIVALENT': mg_conc, 
                    'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_DNTP_CONC': 0.8,
                    'PRIMER_TM_FORMULA': 1, 'PRIMER_SALT_CORRECTIONS': 1, 'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
                    'PRIMER_NUM_RETURN': num_returns 
                }

                if min_amp.strip().isdigit() and max_amp.strip().isdigit(): 
                    global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp.strip()), int(max_amp.strip())]]
                elif min_amp.strip().isdigit(): 
                    global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp.strip()), len(clean_seq)]]
                elif max_amp.strip().isdigit(): 
                    global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[50, int(max_amp.strip())]]

                final_candidates = []

                # פונקציית איסוף מועמדים (לוגיקת Round-Robin)
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
                    diverse_pool = []
                    for i in range(max_cands):
                        for j_cands in all_junction_results:
                            if i < len(j_cands):
                                diverse_pool.append(j_cands[i])
                                if len(diverse_pool) == num_returns: break
                        if len(diverse_pool) == num_returns: break
                    final_candidates = sorted(diverse_pool, key=lambda x: x['penalty'])

                else:
                    seq_args = {'SEQUENCE_ID': project_name, 'SEQUENCE_TEMPLATE': clean_seq}
                    raw_res = primer3.bindings.designPrimers(seq_args, global_args)
                    for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                        final_candidates.append({
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

                if not final_candidates: 
                    st.warning("⚠️ No primers found. Try relaxing the constraints.")
                else:
                    parsed_data = []
                    blast_details = [] 
                    
                    for idx, cand in enumerate(final_candidates):
                        f_seq = cand['forward_seq']
                        r_seq = cand['reverse_seq']
                        
                        # חישוב אקטיבי של ה-Delta G ב-kcal/mol כדי לקבל ערכים שליליים אמיתיים
                        f_hairpin_dg = primer3.calc_hairpin(f_seq).dg / 1000.0
                        r_hairpin_dg = primer3.calc_hairpin(r_seq).dg / 1000.0
                        f_self_dg = primer3.calc_homodimer(f_seq).dg / 1000.0
                        r_self_dg = primer3.calc_homodimer(r_seq).dg / 1000.0
                        cross_dg = primer3.calc_heterodimer(f_seq, r_seq).dg / 1000.0

                        f_blast_sum, r_blast_sum = "Skipped", "Skipped"
                        f_blast_det, r_blast_det = "", ""
                        
                        if run_blast_check:
                            st.toast(f"Running Specific BLAST for Pair {idx+1}...")
                            f_blast_sum, f_blast_det = check_and_get_blast(f_seq, species_name)
                            time.sleep(1)
                            r_blast_sum, r_blast_det = check_and_get_blast(r_seq, species_name)
                            time.sleep(1)
                            
                        blast_details.append({"F": f_blast_det, "R": r_blast_det})

                        parsed_data.append({
                            "Rank": idx + 1,
                            "Penalty": cand['penalty'],
                            "Junction": cand['junction_used'],
                            "Forward_pos": f"F{cand['left_start']}",
                            "Forward_Seq": f_seq,
                            "F_Self (ΔG)": f_self_dg,
                            "F_Hairpin (ΔG)": f_hairpin_dg,
                            "F_BLAST": f_blast_sum,
                            "F_Tm (°C)": cand['forward_tm'],
                            "Reverse_pos": f"R{cand['right_start']}",
                            "Reverse_Seq": r_seq,
                            "R_Self (ΔG)": r_self_dg,
                            "R_Hairpin (ΔG)": r_hairpin_dg,
                            "R_BLAST": r_blast_sum,
                            "R_Tm (°C)": cand['reverse_tm'],
                            "CrossDimer (ΔG)": cross_dg,
                            "Amp_Length": cand['product_size']
                        })
                    
                    df = pd.DataFrame(parsed_data)
                    st.success("✅ Analysis Complete!")

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
                        # הפעלת עיצוב דיוק של ספרה אחת (.1f) וצביעה באדום לפי Thresholds אמיתיים
                        styled_df = df.style.format(precision=1)
                        if hasattr(styled_df, "map"):
                            styled_df = styled_df.map(lambda x: color_negative_red(x, -5.0), subset=['F_Self (ΔG)', 'R_Self (ΔG)', 'CrossDimer (ΔG)']) \
                                                 .map(lambda x: color_negative_red(x, -3.0), subset=['F_Hairpin (ΔG)', 'R_Hairpin (ΔG)'])
                        else:
                            styled_df = styled_df.applymap(lambda x: color_negative_red(x, -5.0), subset=['F_Self (ΔG)', 'R_Self (ΔG)', 'CrossDimer (ΔG)']) \
                                                 .applymap(lambda x: color_negative_red(x, -3.0), subset=['F_Hairpin (ΔG)', 'R_Hairpin (ΔG)'])
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
