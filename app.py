import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="HD Metabolic Framework", page_icon="🧠", layout="wide")

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
                # KEGG format is usually: ID  SYMBOL; DESCRIPTION
                if ';' in line:
                    parts = line.split('; ')
                    description = parts[1].strip()
                    # Now handle the ID and Symbol part
                    id_symbol_part = parts[0].strip()
                    # Split by multiple spaces to separate ID from Symbol
                    sub_parts = id_symbol_part.split(None, 1) # Split only on the first whitespace
                    if len(sub_parts) >= 2:
                        gene_id = sub_parts[0].strip()
                        gene_symbol = sub_parts[1].strip()
                        genes.append({'ID': gene_id, 'Symbol': gene_symbol, 'Description': description})
    return pd.DataFrame(genes)


# Load the data
df = get_kegg_genes("hsa05016")

# --- NEW LOGIC: ANNOTATE GENE ROLES ---
CORE_HD_GENES = ["HTT", "BDNF", "CASP3", "CREB1", "TP53", "SOD1", "PPARGC1A"]

df["Role"] = df["Symbol"].apply(
    lambda x: "⭐ Core HD Gene" if x in CORE_HD_GENES else "🧬 Pathway Component"
)


# Download Button for CV
try:
    with open("CV_Yashwant_Nama_PhD_Application.pdf", "rb") as file:
        st.sidebar.download_button(
            label="📄 Download My CV",
            data=file,
            file_name="Yashwant_Nama_CV.pdf",
            mime="application/pdf"
        )
except FileNotFoundError:
    st.sidebar.error("Upload CV to GitHub to enable download")

st.sidebar.markdown("---")
st.sidebar.header("Project Progress")
st.sidebar.success("Phase 1: Data Acquisition ✅")
st.sidebar.success("Phase 2: Network Visualization ✅")
st.sidebar.info("Phase 3: Pathway Analysis 🔄")

# --- MAIN CONTENT ---
st.title("🧬 Huntington's Disease (HD) Metabolic Framework")
st.markdown("""
### **Disease Context: hsa05016**
This dashboard provides a computational analysis of the **Huntington's Disease pathway**. 
By mapping gene interactions and metabolic disruptions, we can better understand:
*   **Mitochondrial Dysfunction:** Impaired energy production in neurons.
*   **Proteotoxicity:** Accumulation of mutant Huntingtin (HTT) protein.
*   **Excitotoxicity:** Overstimulation of glutamate receptors leading to cell death.
""")

# --- TABS FOR ANALYSIS ---
tab1, tab2, tab3 = st.tabs(["📊 Gene Data", "🕸️ Interaction Network", "📚 Literature"])

with tab1:
    st.subheader("Genetic Components of hsa05016")
    search_query = st.text_input("🔍 Search for a specific gene (e.g., HTT, BDNF, CASP3):")

    if search_query:
        filtered_df = df[df['Symbol'].str.contains(search_query.upper(), na=False)]
        st.write(f"Showing results for: **{search_query}**")
        st.dataframe(filtered_df, use_container_width=True)
    else:
        st.dataframe(df, use_container_width=True)

with tab2:
    st.subheader("Protein Interaction Network (Preview)")
    st.write("Visualizing how key genes interact with the HTT (Huntingtin) hub.")
    
    # Create the graph
    G = nx.Graph()
    subset = df.head(20) # Focusing on top 20 genes for clarity
    for i, row in subset.iterrows():
        G.add_node(row['Symbol'])
        if row['Symbol'] != 'HTT':
            G.add_edge('HTT', row['Symbol'])

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 7))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, 
            node_color='orange', 
            node_size=2500, 
            edge_color='gray', 
            font_size=10, 
            font_weight='bold')
    st.pyplot(fig)

with tab3:
    st.header("Research Bibliography")
    col1, col2 = st.columns(2)

    with col1:
        with st.expander("📖 Essential HD Reading"):
            st.markdown("""
            *   **Ross CA, Tabrizi SJ (2011)** - *Molecular pathogenesis to clinical treatment.* (Lancet Neurol)
            *   **Saudou F, Humbert S (2016)** - *The Biology of Huntingtin.* (Neuron)
            *   **Bates GP, et al. (2015)** - *Huntington disease.* (Nat Rev Dis Primers)
            """)

    with col2:
        with st.expander("🌐 Databases Used"):
            st.markdown("""
            *   **KEGG:** Kyoto Encyclopedia of Genes and Genomes
            *   **STRING:** Protein-Protein Interaction Networks
            *   **GEO:** Gene Expression Omnibus (NCBI)
            """)

# --- FOOTER ---
st.sidebar.markdown("---")
st.sidebar.caption("Data: KEGG API | Developed for PhD Portfolio")
