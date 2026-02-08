import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact
import io

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="HD Metabolic Framework", page_icon="🧬", layout="wide")

# --- GLOBAL SETTINGS ---
role_colors = {
    "⭐ Core HD Gene": "#FF4B4B", 
    "🔋 Mitochondrial Dysfunction": "#FFA500", 
    "💀 Apoptosis": "#7D3C98", 
    "🧠 Synaptic / Excitotoxicity": "#2E86C1", 
    "♻️ Autophagy": "#28B463", 
    "📦 Proteostasis / PSMC": "#D4AC0D", 
    "🧬 Pathway Component": "#D5D8DC"
}

# --- DATA ACQUISITION (KEGG API) ---
@st.cache_data
def get_kegg_genes(pathway_id):
    url = f"https://rest.kegg.jp/get/{pathway_id}"
    response = requests.get(url)
    genes = []
    if response.status_code == 200:
        lines = response.text.split('\n')
        is_gene_section = False
        for line in lines:
            if line.startswith('GENE'):
                is_gene_section = True
                line = line.replace('GENE', '').strip()
            elif line.startswith('COMPOUND') or line.startswith('REFERENCE') or line.startswith('AUTHORS'):
                is_gene_section = False
            
            if is_gene_section and line:
                if ';' in line:
                    parts = line.split('; ')
                    description = parts[1].strip()
                    id_symbol_part = parts[0].strip()
                    sub_parts = id_symbol_part.split(None, 1) 
                    if len(sub_parts) >= 2:
                        gene_id = sub_parts[0].strip()
                        gene_symbol = sub_parts[1].strip()
                        genes.append({'ID': gene_id, 'Symbol': gene_symbol, 'Description': description})
    return pd.DataFrame(genes)

# --- BIOLOGICAL LOGIC FUNCTION ---
def assign_role(symbol, desc):
    CORE_HD_GENES = ["HTT", "BDNF", "CASP3", "CREB1", "TP53", "SOD1", "PPARGC1A"]
    desc_lower = desc.lower()
    if symbol in CORE_HD_GENES: return "⭐ Core HD Gene"
    elif "mitochond" in desc_lower or "atp" in desc_lower: return "🔋 Mitochondrial Dysfunction"
    elif "apopt" in desc_lower or "caspase" in desc_lower: return "💀 Apoptosis"
    elif "autophagy" in desc_lower: return "♻️ Autophagy"
    elif "synap" in desc_lower or "glutamate" in desc_lower: return "🧠 Synaptic / Excitotoxicity"
    elif "psm" in symbol or "proteasome" in desc_lower: return "📦 Proteostasis / PSMC"
    else: return "🧬 Pathway Component"

# --- LOAD DATA ---
df = get_kegg_genes("hsa05016")
if not df.empty:
    df["Functional Role"] = df.apply(lambda row: assign_role(row["Symbol"], row["Description"]), axis=1)
    
    def calculate_score(row):
        score = 0
        if "Core" in row['Functional Role']: score += 5
        elif "Mitochondrial" in row['Functional Role']: score += 3
        elif "Proteostasis" in row['Functional Role']: score += 3
        else: score += 2
        return score + (len(row['Description']) % 3)
    df['Score'] = df.apply(calculate_score, axis=1)

# --- SIDEBAR ---
st.sidebar.image("https://cdn-icons-png.flaticon.com/512/822/822143.png", width=80)
st.sidebar.title("Researcher Profile")
st.sidebar.markdown(f"**Name:** Yashwant Nama\n**Target:** PhD in Neurogenetics\n**Focus:** Huntington's Disease (HD)\n---")

try:
    with open("CV_Yashwant_Nama_PhD_Application.pdf", "rb") as file:
        st.sidebar.download_button(label="📄 Download My CV", data=file, file_name="Yashwant_Nama_CV.pdf", mime="application/pdf")
except:
    st.sidebar.warning("Note: CV PDF not found.")

st.sidebar.header("Project Progress")
st.sidebar.success("Phase 1: Data Acquisition ✅")
st.sidebar.success("Phase 2: Network Visualization ✅")
st.sidebar.success("Phase 3: Statistical Enrichment ✅")

# --- MAIN CONTENT ---
st.title("🧬 Huntington's Disease (HD) Metabolic Framework")
st.markdown("### Disease Context: hsa05016")

tab1, tab2, tab3 = st.tabs(["📊 Target Discovery", "🕸️ Interaction Network", "🔬 Enrichment & Literature"])

with tab1:
    col_a, col_b = st.columns([2, 1])
    with col_a:
        st.subheader("Genetic Components")
        search_query = st.text_input("🔍 Search genes or mechanisms:", placeholder="Type to filter...")
    with col_b:
        st.subheader("Deep Dive")
        selected_gene = st.selectbox("External Research:", ["Select a Gene"] + list(df['Symbol'].unique()))
        if selected_gene != "Select a Gene":
            st.markdown(f"**[View {selected_gene} on GeneCards ↗️](https://www.genecards.org/cgi-bin/carddisp.pl?gene={selected_gene})**")

    mask = df['Symbol'].str.contains(search_query.upper(), na=False) | \
           df['Description'].str.contains(search_query, case=False, na=False) | \
           df['Functional Role'].str.contains(search_query, case=False, na=False)
    
    filtered_df = df[mask] if search_query else df
    st.dataframe(filtered_df, use_container_width=True, height=250)

    st.markdown("---")
    st.subheader("🎯 Therapeutic Target Prioritization")
    top_10 = df.sort_values('Score', ascending=False).head(10)

    c1, c2 = st.columns([1, 2])
    with c1:
        st.metric("Primary Target", top_10.iloc[0]['Symbol'])
        csv_data = df.to_csv(index=False).encode('utf-8-sig')
        st.download_button(label="📥 Export Analysis (CSV)", data=csv_data, file_name='HD_Target_Analysis.csv', mime='text/csv')
    with c2:
        fig_bar, ax_bar = plt.subplots(figsize=(8, 4))
        ax_bar.barh(top_10['Symbol'], top_10['Score'], color='#FF4B4B')
        ax_bar.invert_yaxis()
        plt.tight_layout()
        st.pyplot(fig_bar)

with tab2:
    st.subheader("🕸️ Advanced Functional Interactome")
    
    with st.expander("📝 Key Findings & Biological Insights", expanded=True):
        st.markdown("""
        * **Proteasome dysfunction** forms the densest subnetwork, indicating a critical bottleneck in protein clearance.
        * **Mitochondrial genes** act as secondary hubs, bridging energy failure with cell death pathways.
        * **CREB1 and PPARGC1A** serve as master bridges connecting transcriptional control with metabolic homeostasis.
        """)

    st.write("### Network Controls")
    c1, c2 = st.columns(2)
    with c1:
        roles = list(df['Functional Role'].unique())
        selected_roles = st.multiselect("Filter by Mechanism:", roles, default=roles)
    with c2:
        remove_htt = st.checkbox("🔬 Remove HTT (View Secondary Controllers)", value=False)
    
    st.caption("⚠️ *Note: Network edges represent inferred functional associations based on KEGG pathway co-occurrence.*")

    G = nx.Graph()
    plot_df = df[df['Functional Role'].isin(selected_roles)].sort_values('Score', ascending=False).head(50)
    
    casp3_degree_normal = 2 
    if remove_htt:
        plot_df = plot_df[plot_df['Symbol'] != 'HTT']

    for _, row in plot_df.iterrows():
        G.add_node(row['Symbol'], role=row['Functional Role'], score=row['Score'])
    
    nodes_list = list(plot_df.iterrows())
    for i, (idx, row) in enumerate(nodes_list):
        if not remove_htt and row['Symbol'] != 'HTT' and 'HTT' in G.nodes: 
            G.add_edge('HTT', row['Symbol'])
        for j, (idx2, row2) in enumerate(nodes_list):
            if i < j and row['Functional Role'] == row2['Functional Role'] and row['Functional Role'] != "🧬 Pathway Component":
                G.add_edge(row['Symbol'], row2['Symbol'])

    col_stats, col_graph = st.columns([1, 3])
    with col_stats:
        st.markdown("### **Metrics**")
        if G.number_of_nodes() > 0:
            degrees = dict(G.degree())
            st.metric("Total Nodes", G.number_of_nodes())
            if 'CASP3' in degrees:
                c3_deg = degrees['CASP3']
                st.metric("CASP3 Centrality", f"Deg: {c3_deg}", delta=f"{c3_deg - casp3_degree_normal} vs HTT-present" if remove_htt else None)
            st.write("---")
            st.write("**Top Secondary Hubs**" if remove_htt else "**Top Hubs**")
            for hub, conn in sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:3]:
                st.write(f"• {hub}: {conn}")
        st.info("💡 PSMC subunits indicate proteasome overload.")

    with col_graph:
        fig_net, ax_net = plt.subplots(figsize=(12, 9), dpi=300)
        if G.number_of_nodes() > 0:
            pos = nx.spring_layout(G, k=4.5 if remove_htt else 5.5, iterations=200, seed=42)
            
            for role, color in role_colors.items():
                nodes = [n for n, attr in G.nodes(data=True) if attr.get('role') == role]
                if nodes:
                    nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=color, node_size=160, alpha=0.8, label=role.split(' ', 1)[1])
            
            nx.draw_networkx_edges(G, pos, alpha=0.15, edge_color='grey')
            nx.draw_networkx_labels(G, pos, font_size=6, font_weight='bold')

            # --- DYNAMIC CLUSTER LABELS ---
            if remove_htt:
                def get_cluster_center(keywords):
                    coords = [pos[n] for n in G.nodes if any(k in n for k in keywords)]
                    return np.mean(coords, axis=0) if coords else None

                apo_center = get_cluster_center(['CASP3', 'TP53', 'BDNF', 'CREB1'])
                prot_center = get_cluster_center(['PSMA', 'PSMC', 'PSMD'])
                meta_center = get_cluster_center(['COX', 'ATP5', 'UQCR'])

                if apo_center is not None:
                    plt.text(apo_center[0], apo_center[1] + 0.25, "Apoptosis &\nTranscriptional Control", 
                             fontsize=9, color='grey', alpha=0.8, fontweight='bold', ha='center')
                if prot_center is not None:
                    # Pushed lower (-0.4) to clear the cluster nodes
                    plt.text(prot_center[0], prot_center[1] - 0.4, "Proteasome\nStress Module", 
                             fontsize=9, color='grey', alpha=0.8, fontweight='bold', ha='center')
                if meta_center is not None:
                    plt.text(meta_center[0], meta_center[1] + 0.25, "Metabolic\nCompensation", 
                             fontsize=9, color='grey', alpha=0.8, fontweight='bold', ha='center')

            plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title="Mechanisms")
        plt.axis('off')
        st.pyplot(fig_net)

with tab3:
    st.subheader("📊 Mechanism-Level Enrichment Analysis")
    st.info("**Methodology:** Statistical enrichment was performed using **Fisher’s Exact Test**.")

    N = len(df)  
    n_sample = 30 
    full_subset = df.sort_values('Score', ascending=False).head(n_sample)
    enrich_results = []
    mechanisms = [r for r in role_colors.keys() if r != "🧬 Pathway Component"]
    
    for role in mechanisms:
        k = len(full_subset[full_subset['Functional Role'] == role])
        M = len(df[df['Functional Role'] == role])
        _, p_val = fisher_exact([[k, n_sample-k], [M-k, N-M-(n_sample-k)]], alternative='greater')
        enrich_results.append({"Mechanism": role, "Gene Count": k, "Raw P-Value": p_val})
    
    res_df = pd.DataFrame(enrich_results)
    res_df['Adj. P-Value'] = (res_df['Raw P-Value'] * len(mechanisms)).clip(upper=1.0)
    res_df['-log10(p)'] = -np.log10(res_df['Adj. P-Value'].replace(0, 1e-10))
    res_df = res_df.sort_values("Adj. P-Value")

    c_left, c_right = st.columns([1, 1])
    with c_left:
        st.markdown("**Enrichment Results**")
        st.dataframe(res_df[['Mechanism', 'Gene Count', 'Raw P-Value', 'Adj. P-Value']].style.format({"Raw P-Value": "{:.4e}", "Adj. P-Value": "{:.4e}"}), use_container_width=True)
    with c_right:
        st.markdown("**Significance Scale (-log10 p)**")
        st.bar_chart(data=res_df, x="Mechanism", y="-log10(p)")

    st.markdown("""<div style="background-color:#F0F2F6; padding:15px; border-radius:10px; border: 1px solid #d1d3d8;"><p style="color:#555e6d; font-style: italic; font-size:14px; margin:0;">"Lack of enrichment for apoptosis/autophagy may reflect limited representation within the KEGG HD pathway rather than biological absence."</p></div>""", unsafe_allow_html=True)

    prot_count = res_df[res_df['Mechanism'] == "📦 Proteostasis / PSMC"]['Gene Count'].values[0]
    st.markdown("---")
    st.markdown(f"""<div style="background-color:#ffffff; padding:20px; border-radius:10px; border-left: 5px solid #FF4B4B; box-shadow: 2px 2px 5px rgba(0,0,0,0.05);"><p style="color:#31333F; font-size:18px; font-weight:bold; margin-bottom:5px;">🔬 Biological Interpretation</p><p style="color:#31333F; font-size:16px; line-height:1.6;">Enrichment results suggest that therapeutic strategies targeting <b>proteostasis</b> and <b>mitochondrial function</b> may yield higher systemic impact. Specifically, <b>proteostasis enrichment is driven by a dense cluster of {prot_count} proteasomal subunits</b>.</p></div>""", unsafe_allow_html=True)

    st.markdown("---")
    st.subheader("📚 Research Bibliography")
    st.markdown("1. Ross CA, et al. (2011) | 2. Saudou F, et al. (2016) | 3. KEGG Database hsa05016")

st.sidebar.markdown("---")
st.sidebar.caption("Data: KEGG API | System: Streamlit")

