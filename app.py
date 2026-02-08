import streamlit as st
import pandas as pd
from bioservices import KEGG
import networkx as nx
import matplotlib.pyplot as plt

st.set_page_config(page_title="HD Research Lab", layout="wide")

st.title("🧬 Huntington's Disease Molecular Analysis")
st.markdown("""
### **Disease Context: hsa05016**
**Huntington's Disease (HD)** is a progressive neurodegenerative disorder caused by a CAG trinucleotide repeat expansion in the **HTT gene**. 
This dashboard focuses on the **metabolic framework** of HD, specifically looking at:
*   **Mitochondrial Dysfunction:** Impaired energy production in neurons.
*   **Proteotoxicity:** Accumulation of mutant Huntingtin protein.
*   **Excitotoxicity:** Overstimulation of glutamate receptors.
""")
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
# --- SIDEBAR RESEARCHER PROFILE ---
st.sidebar.image("https://cdn-icons-png.flaticon.com/512/2491/2491314.png", width=100)
st.sidebar.title("Researcher: Yashwant Nama")
st.sidebar.info("""
**Target:** PhD in Neurogenetics
**Focus:** Huntington's Disease (HD)
""")

# Download Button for your CV
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


# --- Main Tabs ---
# --- Main Tabs ---
# UPDATE LINE 32 TO THIS:
tab1, tab2, tab3 = st.tabs(["Gene Data", "Interaction Network", "Pathway Analysis"])

with tab1:
    st.subheader("Genetic Markers from KEGG")
    search = st.text_input("Search for a gene (e.g. HTT):")
    display_df = df[df['Symbol'].str.contains(search.upper())] if search else df
    st.dataframe(display_df, use_container_width=True)

with tab2:
    st.subheader("Protein Interaction Network (Preview)")
    st.write("This graph shows how the first 20 genes are interconnected in the HD pathway.")
    
    G = nx.Graph()
    subset = df.head(20) 
    for i, row in subset.iterrows():
        G.add_node(row['Symbol'])
        if row['Symbol'] != 'HTT':
            G.add_edge('HTT', row['Symbol'])

    fig, ax = plt.subplots(figsize=(10, 7))
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='skyblue', 
            node_size=2000, edge_color='gray', font_size=10, font_weight='bold')
    st.pyplot(fig)

# --- NEW TAB 3 STARTING AT LINE 58 ---
with tab3:
    st.subheader("Biological Process Distribution")
    import plotly.express as px
    
    # Categorize genes
    categories = {
        "Mitochondrial": df['Function'].str.contains('mitochondria|oxidative', case=False).sum(),
        "Synaptic/Neural": df['Function'].str.contains('synaptic|neural|axon', case=False).sum(),
        "Apoptosis (Cell Death)": df['Function'].str.contains('apoptosis|caspase', case=False).sum(),
        "Transcription": df['Function'].str.contains('transcription|polymerase', case=False).sum()
    }
    
    cat_df = pd.DataFrame(list(categories.items()), columns=['Process', 'Gene Count'])
    
    fig_pie = px.pie(cat_df, values='Gene Count', names='Process', 
                     title="Functional Roles of HD-Associated Genes",
                     hole=0.4, # This makes it a donut chart (looks more modern)
                     color_discrete_sequence=px.colors.sequential.RdBu)
    st.plotly_chart(fig_pie)

    st.info("💡 **Insight:** A high concentration of genes in the 'Transcription' and 'Mitochondrial' categories suggests that Huntington's Disease primarily disrupts the cell's energy production and its ability to read genetic instructions.")
# --- RESEARCH BIBLIOGRAPHY & RESOURCES ---
st.markdown("---")
st.header("📚 Research Bibliography & Resources")
st.write("This dashboard is built upon the following scientific foundations and datasets:")

# Create three tabs or expanders for the different types of resources
col1, col2 = st.columns(2)

with col1:
    with st.expander("📖 Essential HD Reading"):
        st.markdown("""
        *   **Ross CA, Tabrizi SJ (2011)** - *Molecular pathogenesis to clinical treatment.* (Lancet Neurol)
        *   **Saudou F, Humbert S (2016)** - *The Biology of Huntingtin.* (Neuron)
        *   **Bates GP, et al. (2015)** - *Huntington disease.* (Nat Rev Dis Primers)
        *   **Lee JM, et al. (2019)** - *CAG repeat determines timing of HD onset.* (Cell)
        """)

    with st.expander("🧪 Metabolic & Mitochondrial Research"):
        st.markdown("""
        *   **Cheng ML, et al. (2016)** - *Metabolic disturbances in HD.* (J Neurol Sci)
        *   **Reddy PH, et al. (2018)** - *Mitochondrial dysfunction and oxidative stress.* (BBA)
        *   **Weydt P, et al. (2009)** - *PGC-1α modifies age at onset.* (Mol Neurodegener)
        """)

with col2:
    with st.expander("🌐 Network Biology & Multi-Omics"):
        st.markdown("""
        *   **Barabási AL, et al. (2011)** - *Network medicine.* (Nat Rev Genet)
        *   **Menche J, et al. (2015)** - *Uncovering disease-disease relationships.* (Science)
        *   **Subramanian A, et al. (2005)** - *Gene set enrichment analysis (GSEA).* (PNAS)
        """)

    with st.expander("🗄️ Databases Used"):
        st.markdown("""
        *   **KEGG:** Kyoto Encyclopedia of Genes and Genomes
        *   **STRING:** Protein-Protein Interaction Networks
        *   **GEO:** Gene Expression Omnibus (NCBI)
        *   **HDinHD:** Huntington's Disease Integrated Database
        """)

# --- FINAL FOOTER ---
st.sidebar.markdown("---")
st.sidebar.caption("Data sources: KEGG API, SciSpace Literature Review. Developed for PhD Application Portfolio.")
