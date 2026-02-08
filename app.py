import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
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

df = get_kegg_genes("hsa05016")
if not df.empty:
    df["Functional Role"] = df.apply(lambda row: assign_role(row["Symbol"], row["Description"]), axis=1)
    df['Score'] = df.apply(lambda r: (5 if "Core" in r['Functional Role'] else 3 if "Mito" in r['Functional Role'] else 2) + (len(r['Description']) % 3), axis=1)

# --- SIDEBAR ---
st.sidebar.title("Researcher Profile")
st.sidebar.markdown("**Yashwant Nama**\nPhD Applicant | Neurogenetics")

# --- MAIN CONTENT ---
st.title("🧬 Huntington's Disease (HD) Metabolic Framework")
tab1, tab2, tab3 = st.tabs(["📊 Target Discovery", "🕸️ Interactive Interactome", "🔬 Enrichment"])

with tab1:
    st.subheader("Genetic Components")
    st.dataframe(df, use_container_width=True, height=300)

with tab2:
    st.subheader("🕸️ Interactive Functional Interactome")
    st.info("💡 **Magnification Active:** Use your mouse wheel to zoom. Click and drag to move. Hover over nodes for details.")
    
    # Network Setup
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

    pos = nx.spring_layout(G, k=0.5, seed=42)

    # Create Edges for Plotly
    edge_x, edge_y = [], []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')

    # Create Nodes for Plotly
    node_traces = []
    for role, color in role_colors.items():
        node_x, node_y, node_text, node_size = [], [], [], []
        for node in G.nodes():
            if G.nodes[node]['role'] == role:
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)
                node_text.append(f"Gene: {node}<br>Role: {role}<br>Score: {G.nodes[node]['score']}")
                node_size.append(G.nodes[node]['score'] * 5)
        
        node_traces.append(go.Scatter(
            x=node_x, y=node_y, mode='markers+text',
            text=[n if G.nodes[n]['role'] == "⭐ Core HD Gene" else "" for n in G.nodes() if G.nodes[n]['role'] == role],
            textposition="top center",
            hoverinfo='text', hovertext=node_text,
            name=role,
            marker=dict(size=node_size, color=color, line_width=2)
        ))

    fig = go.Figure(data=[edge_trace] + node_traces,
                 layout=go.Layout(
                    showlegend=True,
                    hovermode='closest',
                    margin=dict(b=0,l=0,r=0,t=0),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
                ))
    
    st.plotly_chart(fig, use_container_width=True)

with tab3:
    st.subheader("📊 Statistical Enrichment")
    # (Keeping your existing Fisher logic here...)
    N, n_sample = 20000, len(subset)
    enrich_results = []
    for role in [r for r in role_colors.keys() if r != "🧬 Pathway Component"]:
        k = len(subset[subset['Functional Role'] == role])
        M = len(df[df['Functional Role'] == role])
        _, p_val = fisher_exact([[k, n_sample-k], [M-k, N-M-(n_sample-k)]], alternative='greater')
        enrich_results.append({"Mechanism": role, "P-Value": p_val})
    st.table(pd.DataFrame(enrich_results))
