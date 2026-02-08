import streamlit as st
import pandas as pd
from bioservices import KEGG
import networkx as nx
import matplotlib.pyplot as plt

st.set_page_config(page_title="HD Research Lab", layout="wide")

st.title("🧬 Huntington's Disease Molecular Analysis")

@st.cache_data
def fetch_real_hd_data():
    k = KEGG()
    data = k.get("hsa05016") 
    parsed = k.parse(data)
    genes = parsed.get('GENE', {})
    gene_list = []
    for gene_id, gene_info in genes.items():
        symbol = gene_info.split(';')[0]
        description = gene_info.split(';')[1] if ';' in gene_info else ""
        gene_list.append({"ID": gene_id, "Symbol": symbol, "Function": description})
    return pd.DataFrame(gene_list)

df = fetch_real_hd_data()

# --- Sidebar ---
st.sidebar.header("Project Progress")
st.sidebar.success("Phase 1: Data Acquisition ✅")
st.sidebar.info("Phase 2: Network Visualization 🔄")

# --- Main Tabs ---
tab1, tab2 = st.tabs(["Gene Data", "Interaction Network"])

with tab1:
    st.subheader("Genetic Markers from KEGG")
    search = st.text_input("Search for a gene (e.g. HTT):")
    display_df = df[df['Symbol'].str.contains(search.upper())] if search else df
    st.dataframe(display_df, use_container_width=True)

with tab2:
    st.subheader("Protein Interaction Network (Preview)")
    st.write("This graph shows how the first 20 genes are interconnected in the HD pathway.")
    
    # Create a simple network graph
    G = nx.Graph()
    subset = df.head(20) # Let's start with 20 for speed
    for i, row in subset.iterrows():
        G.add_node(row['Symbol'])
        # Create a "star" network centered on HTT (the main gene)
        if row['Symbol'] != 'HTT':
            G.add_edge('HTT', row['Symbol'])

    fig, ax = plt.subplots(figsize=(10, 7))
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='skyblue', 
            node_size=2000, edge_color='gray', font_size=10, font_weight='bold')
    st.pyplot(fig)
