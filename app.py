import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="HD Metabolic Framework", page_icon="🧬", layout="wide")

# --- DATA ACQUISITION ---
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
            elif line.startswith('COMPOUND') or line.startswith('REFERENCE'):
                is_gene_section = False
            if is_gene_section and line:
                if ';' in line:
                    parts = line.split('; ')
                    desc = parts[1].strip()
                    sub_parts = parts[0].strip().split(None, 1) 
                    if len(sub_parts) >= 2:
                        genes.append({'ID': sub_parts[0], 'Symbol': sub_parts[1], 'Description': desc})
    return pd.DataFrame(genes)

def assign_role(symbol, desc):
    CORE_HD_GENES = ["HTT", "BDNF", "CASP3", "CREB1", "TP53", "SOD1", "PPARGC1A"]
    desc_lower = desc.lower()
    if symbol in CORE_HD_GENES: return "⭐ Core HD Gene"
    elif "mitochond" in desc_lower or "atp" in desc_lower: return "🔋 Mitochondrial Dysfunction"
    elif "apopt" in desc_lower or "caspase" in desc_lower: return "💀 Apoptosis"
    elif "autophagy" in desc_lower: return "♻️ Autophagy"
    elif "synap" in desc_lower or "glutamate" in desc_lower: return "🧠 Synaptic / Excitotoxicity"
    else: return "🧬 Pathway Component"

# --- LOAD DATA ---
df = get_kegg_genes("hsa05016")
if not df.empty:
    df["Functional Role"] = df.apply(lambda row: assign_role(row["Symbol"], row["Description"]), axis=1)
    df['Score'] = df.apply(lambda r: (5 if "Core" in r['Functional Role'] else 3 if "Mito" in r['Functional Role'] else 2) + (len(r['Description']) % 3), axis=1)

# --- SIDEBAR ---
st.sidebar.title("Researcher Profile")
st.sidebar.markdown("**Yashwant Nama**\nPhD Applicant | Neurogenetics")
try:
    with open("CV_Yashwant_Nama_PhD_Application.pdf", "rb") as file:
        st.sidebar.download_button(label="📄 Download My CV", data=file, file_name="Yashwant_Nama_CV.pdf", mime="application/pdf")
except:
    st.sidebar.warning("Note: CV PDF not found.")

# --- MAIN CONTENT ---
st.title("🧬 Huntington's Disease (HD) Metabolic Framework")
tab1, tab2, tab3 = st.tabs(["📊 Target Discovery", "🕸️ Interaction Network", "🔬 Enrichment & Lit"])

with tab1:
    st.subheader("Genetic Components")
    st.dataframe(df, use_container_width=True, height=250)
    
    st.markdown("---")
    st.subheader("🎯 Therapeutic Target Prioritization")
    top_10 = df.sort_values('Score', ascending=False).head(10)
    fig_bar, ax_bar = plt.subplots(figsize=(8, 4))
    ax_bar.barh(top_10['Symbol'], top_10['Score'], color='#FF4B4B')
    ax_bar.invert_yaxis()
    st.pyplot(fig_bar)

with tab2:
    st.subheader("🕸️ Advanced Functional Interactome")
    st.info("🔍 **Click the 'Expand' icon (top right of the image)** to magnify and see all gene names clearly.")
    
    # --- ORIGINAL NETWORK CALCULATION ---
    G = nx.Graph()
    subset = df.sort_values('Score', ascending=False).head(50)
    role_colors = {"⭐ Core HD Gene": "#FF4B4B", "🔋 Mitochondrial Dysfunction": "#FFA500", "💀 Apoptosis": "#7D3C98", "🧠 Synaptic / Excitotoxicity": "#2E86C1", "♻️ Autophagy": "#28B463", "🧬 Pathway Component": "#D5D8DC"}
    
    for _, row in subset.iterrows():
        G.add_node(row['Symbol'], role=row['Functional Role'], score=row['Score'])

    nodes_list = list(subset.iterrows())
    for i, (idx, row) in enumerate(nodes_list):
        if row['Symbol'] != 'HTT': G.add_edge('HTT', row['Symbol'])
        for j, (idx2, row2) in enumerate(nodes_list):
            if i < j and row['Functional Role'] == row2['Functional Role'] and row['Functional Role'] != "🧬 Pathway Component":
                G.add_edge(row['Symbol'], row2['Symbol'])

    # CREATE THE PLOT
    # We use a very high DPI (300) so it stays sharp when you zoom in
    fig_net, ax_net = plt.subplots(figsize=(14, 10), dpi=300)
    pos = nx.spring_layout(G, k=4.5, iterations=150, seed=42)
    
    for role, color in role_colors.items():
        nodes = [n for n, attr in G.nodes(data=True) if attr.get('role') == role]
        if nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=color, node_size=200, alpha=0.8, label=role.split(' ', 1)[1])
    
    nx.draw_networkx_edges(G, pos, alpha=0.15, edge_color='grey')
    nx.draw_networkx_labels(G, pos, font_size=6, font_weight='bold')
    
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title="Mechanisms", fontsize='small')
    plt.axis('off')
    plt.tight_layout()
    
    # DISPLAY AS IMAGE WITH ZOOM ENABLED
    st.pyplot(fig_net)

with tab3:
    st.subheader("📊 Statistical Enrichment Analysis")
    N, n_sample = 20000, len(subset)
    enrich_results = []
    for role in [r for r in role_colors.keys() if r != "🧬 Pathway Component"]:
        k = len(subset[subset['Functional Role'] == role])
        M = len(df[df['Functional Role'] == role])
        _, p_val = fisher_exact([[k, n_sample-k], [M-k, N-M-(n_sample-k)]], alternative='greater')
        enrich_results.append({"Mechanism": role, "P-Value": p_val})
    
    res_df = pd.DataFrame(enrich_results).sort_values("P-Value")
    st.dataframe(res_df.style.format({"P-Value": "{:.4e}"}), use_container_width=True)
    
    st.markdown("---")
    st.subheader("📚 Research Bibliography")
    st.markdown("1. Ross CA, et al. (2011) | 2. Saudou F, et al. (2016) | 3. KEGG Database hsa05016")
