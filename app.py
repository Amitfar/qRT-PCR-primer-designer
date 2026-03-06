import streamlit as st
import primer3
import pandas as pd
import time
import plotly.graph_objects as go
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def check_primer_blast(primer_seq, species):
    try:
        entrez_query = f"{species}[organism]" if species else ""
        handle = NCBIWWW.qblast("blastn", "nt", primer_seq, entrez_query=entrez_query, word_size=7, hitlist_size=10, expect=100)
        blast_record = NCBIXML.read(handle)
        hits = len(blast_record.alignments)
        if hits == 0: return "0 Hits (Check species)"
        elif hits == 1: return "1 Hit (Specific ✅)"
        else: return f"{hits} Hits ⚠️"
    except Exception as e:
        return "BLAST Error"

st.set_page_config(page_title="Cloning Primer Designer", page_icon="🧬", layout="wide")
st.title("🧬 Cloning Primer Designer (Volcani Edition)")
st.markdown("תכנון פריימרים להשתלה עם פיזור חכם, דירוג אובייקטיבי ומפה ויזואלית.")
st.divider()

col1, col2 = st.columns([2, 1])
with col1: sequence = st.text_area("Target Sequence (5' to 3'):", height=200, placeholder="Paste your FASTA sequence...")
with col2:
    project_name = st.text_input("Primer Project Name:", "My_Gene_Cloning")
    species_name = st.text_input("Species (for BLAST):", "Solanum tuberosum")
    junction_pos = st.text_input("Exon-Exon Junctions (comma separated):", placeholder="e.g., 72, 478, 727")
    st.markdown("<br>", unsafe_allow_html=True)
    run_blast_check = st.checkbox("🔍 Run NCBI BLAST Specificity Check (Takes 1-3 minutes)")

with st.expander("⚙️ Advanced Cloning & Reaction Settings"):
    adv_col1, adv_col2, adv_col3, adv_col4 = st.columns(4)
    with adv_col1: cloning_type = st.selectbox("Cloning Method", ["Standard PCR", "Restriction Digestion", "Gibson Assembly"])
    with adv_col2: mg_conc = st.number_input("Mg2+ Concentration (mM)", value=1.5, step=0.1)
    with adv_col3: target_tm = st.slider("Target Tm (°C)", min_value=50.0, max_value=72.0, value=60.0, step=0.5)
    with adv_col4: num_returns = st.slider("Total Primers to Design", min_value=1, max_value=15, value=5, step=1)
        
    st.markdown("<br>", unsafe_allow_html=True)
    amp_col1, amp_col2, _ = st.columns([1, 1, 2])
    with amp_col1: min_amp = st.text_input("Min Amplicon Length (Optional)", placeholder="e.g., 100")
    with amp_col2: max_amp = st.text_input("Max Amplicon Length (Optional)", placeholder="e.g., 1500")

if st.button("🚀 Design Primers", type="primary"):
    if not sequence: st.error("Please enter a DNA sequence to proceed.")
    else:
        with st.spinner('Designing, distributing, and sorting candidates by objective quality...'):
            try:
                clean_seq = "".join(sequence.split()).upper()
                junction_list = [int(x.strip()) for x in junction_pos.split(',') if x.strip().isdigit()] if junction_pos else []

                global_args = {
                    'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': target_tm, 'PRIMER_MIN_TM': target_tm - 5.0, 'PRIMER_MAX_TM': target_tm + 5.0,
                    'PRIMER_SALT_DIVALENT': mg_conc, 'PRIMER_NUM_RETURN': num_returns 
                }

                if min_amp.strip().isdigit() and max_amp.strip().isdigit(): global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp.strip()), int(max_amp.strip())]]
                elif min_amp.strip().isdigit(): global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp.strip()), len(clean_seq)]]
                elif max_amp.strip().isdigit(): global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[50, int(max_amp.strip())]]

                final_candidates = []

                if junction_list:
                    all_junction_results = []
                    for j in junction_list:
                        seq_args = {'SEQUENCE_ID': f"J_{j}", 'SEQUENCE_TEMPLATE': clean_seq, 'SEQUENCE_OVERLAP_JUNCTION_LIST': [j]}
                        raw_res = primer3.bindings.designPrimers(seq_args, global_args)
                        num_ret = raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)
                        
                        j_candidates = []
                        for i in range(num_ret):
                            cand = {
                                'penalty': raw_res.get(f'PRIMER_PAIR_{i}_PENALTY', 999.0), # שליפת הציון האובייקטיבי!
                                'left_start': raw_res.get(f'PRIMER_LEFT_{i}')[0], 'left_len': raw_res.get(f'PRIMER_LEFT_{i}')[1],
                                'right_start': raw_res.get(f'PRIMER_RIGHT_{i}')[0], 'right_len': raw_res.get(f'PRIMER_RIGHT_{i}')[1],
                                'forward_seq': raw_res.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''), 'reverse_seq': raw_res.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
                                'forward_tm': raw_res.get(f'PRIMER_LEFT_{i}_TM', 0), 'reverse_tm': raw_res.get(f'PRIMER_RIGHT_{i}_TM', 0),
                                'product_size': raw_res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0), 'junction_used': j
                            }
                            j_candidates.append(cand)
                        all_junction_results.append(j_candidates)

                    max_cands = max([len(lst) for lst in all_junction_results]) if all_junction_results else 0
                    diverse_pool = []
                    for i in range(max_cands):
                        for j_cands in all_junction_results:
                            if i < len(j_cands):
                                diverse_pool.append(j_cands[i])
                                if len(diverse_pool) == num_returns: break
                        if len(diverse_pool) == num_returns: break
                    
                    # מיון סופי לפי ציון אובייקטיבי!
                    final_candidates = sorted(diverse_pool, key=lambda x: x['penalty'])

                else:
                    seq_args = {'SEQUENCE_ID': project_name, 'SEQUENCE_TEMPLATE': clean_seq}
                    raw_res = primer3.bindings.designPrimers(seq_args, global_args)
                    for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                        cand = {
                            'penalty': raw_res.get(f'PRIMER_PAIR_{i}_PENALTY', 999.0),
                            'left_start': raw_res.get(f'PRIMER_LEFT_{i}')[0], 'left_len': raw_res.get(f'PRIMER_LEFT_{i}')[1],
                            'right_start': raw_res.get(f'PRIMER_RIGHT_{i}')[0], 'right_len': raw_res.get(f'PRIMER_RIGHT_{i}')[1],
                            'forward_seq': raw_res.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''), 'reverse_seq': raw_res.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
                            'forward_tm': raw_res.get(f'PRIMER_LEFT_{i}_TM', 0), 'reverse_tm': raw_res.get(f'PRIMER_RIGHT_{i}_TM', 0),
                            'product_size': raw_res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0), 'junction_used': "None"
                        }
                        final_candidates.append(cand)

                if not final_candidates: st.warning("⚠️ No primers found. Try relaxing the constraints.")
                else:
                    parsed_data = []
                    overhang_note = "Add Restriction Site + Spacer" if cloning_type == "Restriction Digestion" else "Add 20-30bp Vector Overlap" if cloning_type == "Gibson Assembly" else "None (Standard PCR)"

                    for idx, cand in enumerate(final_candidates):
                        forward_blast_res, reverse_blast_res = "Skipped", "Skipped"
                        if run_blast_check:
                            st.toast(f"Running BLAST for Pair {idx+1}...")
                            forward_blast_res = check_primer_blast(cand['forward_seq'], species_name)
                            time.sleep(2)
                            reverse_blast_res = check_primer_blast(cand['reverse_seq'], species_name)
                            time.sleep(2)

                        parsed_data.append({
                            "Rank": idx + 1, "Penalty Score": round(cand['penalty'], 2), "Junction": cand['junction_used'],
                            "Forward_pos": f"F{cand['left_start']}", "Forward_Seq": cand['forward_seq'], "Forward_BLAST": forward_blast_res, "Forward_T_M (°C)": round(cand['forward_tm'], 1),
                            "Reverse_pos": f"R{cand['right_start']}", "Reverse_Seq": cand['reverse_seq'], "Reverse_BLAST": reverse_blast_res, "Reverse_T_M (°C)": round(cand['reverse_tm'], 1),
                            "Amplicon_Length": cand['product_size'], "5'_Modification": overhang_note
                        })
                    
                    df = pd.DataFrame(parsed_data)
                    st.success("✅ Analysis Complete! Candidates successfully generated and sorted by objective quality.")

                    st.subheader("🗺️ Amplicon Map (Distributed Candidates)")
                    fig = go.Figure()
                    fig.add_shape(type="rect", x0=0, y0=0, x1=len(clean_seq), y1=1, line=dict(color="gray", width=2), fillcolor="lightgray")

                    for idx, cand in enumerate(final_candidates):
                        y_center, y_bottom, y_top = 1.5 + (idx * 0.6), 1.5 + (idx * 0.6) - 0.2, 1.5 + (idx * 0.6) + 0.2
                        fig.add_shape(type="rect", x0=cand['left_start'], y0=y_bottom, x1=cand['right_start'], y1=y_top, line=dict(color="DarkGreen", width=1), fillcolor="LimeGreen", opacity=0.8)
                        fig.add_annotation(x=(cand['left_start'] + cand['right_start'])/2, y=y_center, text=f"Candidate {idx+1}", showarrow=False, font=dict(color="black", size=11))

                    max_y = 1.5 + (len(final_candidates) * 0.6)
                    if junction_list:
                        for j in junction_list:
                            fig.add_shape(type="line", x0=j, y0=-0.5, x1=j, y1=max_y, line=dict(color="red", width=2, dash="dash"))
                            fig.add_annotation(x=j, y=max_y + 0.2, text=f"Exon {j}", showarrow=False, font=dict(color="red", size=10))

                    fig.update_layout(xaxis=dict(range=[-20, len(clean_seq)+20], title="Position (bp)"), yaxis=dict(showticklabels=False, range=[-0.5, max_y + 0.5]), height=250 + (len(final_candidates) * 30), margin=dict(l=20, r=20, t=30, b=30), plot_bgcolor="white")
                    st.plotly_chart(fig, use_container_width=True)

                    st.markdown("### Detailed Results Table")
                    st.dataframe(df, use_container_width=True)
                    st.download_button(label="📥 Download Full Results (CSV)", data=df.to_csv(index=False).encode('utf-8'), file_name=f"{project_name}_full_primers.csv", mime='text/csv')

            except Exception as e: st.error(f"Error executing analysis: {e}")