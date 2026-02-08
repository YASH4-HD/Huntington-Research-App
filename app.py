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
    if symbol in CORE_HD_GENES:
        return "⭐ Core HD Gene"
    elif "mitochond" in desc_lower or "atp" in desc_lower:
        return "🔋 Mitochondrial Dysfunction"
    elif "apopt" in desc_lower or "caspase" in desc_lower:
        return "💀 Apoptosis"
    elif "autophagy" in desc_lower:
        return "♻️ Autophagy"
    elif "synap" in desc_lower or "glutamate" in desc_lower:
        return "🧠 Synaptic / Excitotoxicity"
    else:
        return "🧬 Pathway Component"

# --- LOAD AND PROCESS DATA ---
df = get_kegg_genes("hsa05016")
if not df.empty:
    df["Functional Role"] = df.apply(lambda row: assign_role(row["Symbol"], row["Description"]), axis=1)

# --- SIDEBAR RESEARCHER PROFILE ---
st.sidebar.image("https://cdn-icons-png.flaticon.com/512/822/822143.png", width=100)
st.sidebar.title("Researcher: Yashwant Nama")
st.sidebar.info("""
**Target:** PhD in Neurogenetics
**Focus:** Huntington's Disease (HD)
""")

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
    st.subheader("Genetic Components & Biological Roles")
    
    # Gene Deep Dive Selector
    selected_gene = st.selectbox("Select a gene for deeper analysis:", ["None"] + list(df['Symbol'].unique()))
    if selected_gene != "None":
        st.info(f"🧬 **External Link:** [View {selected_gene} on GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene={selected_gene})")
    
    search_query = st.text_input("🔍 Search for a gene or mechanism (e.g., HTT, Mitochondrial, Apoptosis):")

    if search_query:
        # Search across Symbol, Description, and the new Role column
        mask = (df['Symbol'].str.contains(search_query.upper(), na=False)) | \
               (df['Description'].str.contains(search_query, case=False, na=False)) | \
               (df['Functional Role'].str.contains(search_query, case=False, na=False))
        filtered_df = df[mask]
        st.dataframe(filtered_df, use_container_width=True)
    else:
        st.dataframe(df, use_container_width=True)

with tab2:
    st.subheader("Advanced Functional Interactome")
    st.write("This graph clusters genes by biological mechanism. Nodes are colored by functional pathway.")
    
    # 1. Setup Graph and Colors
    G = nx.Graph()
    color_map = []
    
    # Define color scheme for roles
    role_colors = {
        "⭐ Core HD Gene": "#FF4B4B",           # Red
        "🔋 Mitochondrial Dysfunction": "#FFA500", # Orange
        "💀 Apoptosis": "#7D3C98",              # Purple
        "🧠 Synaptic / Excitotoxicity": "#2E86C1", # Blue
        "♻️ Autophagy": "#28B463",              # Green
        "🧬 Pathway Component": "#D5D8DC"       # Grey
    }

    # 2. Add Nodes and internal connections
    subset = df.head(30)
    for i, row in subset.iterrows():
        G.add_node(row['Symbol'], role=row['Functional Role'])
        
        # Connect everything to HTT (The Hub)
        if row['Symbol'] != 'HTT':
            G.add_edge('HTT', row['Symbol'], weight=1)
        
        # EXTRAORDINARY STEP: Connect genes to EACH OTHER if they share the same role
        # This creates the "Clusters"
        for j, row2 in subset.iterrows():
            if i < j and row['Functional Role'] == row2['Functional Role'] and row['Functional Role'] != "🧬 Pathway Component":
                G.add_edge(row['Symbol'], row2['Symbol'], weight=0.5)

    # 3. Create the Visualization
    fig, ax = plt.subplots(figsize=(12, 9), facecolor='white')
    
    # Use Force-Directed Layout (Spring)
    pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
    
    # Draw nodes by role
    for role, color in role_colors.items():
        node_list = [n for n, attr in G.nodes(data=True) if attr.get('role') == role]
        if node_list:
            size = 3000 if "Core" in role else 1200
            nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_color=color, 
                                   node_size=size, alpha=0.9, label=role)

    # Draw edges and labels
    nx.draw_networkx_edges(G, pos, width=1.5, edge_color='lightgrey', alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', font_family='sans-serif')

    plt.legend(scatterpoints=1, loc='upper left', fontsize='small', title="Functional Mechanisms")
    plt.axis('off')
    st.pyplot(fig)
    
    st.info("💡 **Scientific Insight:** Notice how genes with similar colors cluster together. This represents 'Functional Modules' where specific metabolic failures (like Mitochondrial dysfunction) occur in coordination.")


with tab3:
    st.header("Research Bibliography")
    st.write("Foundational papers used to construct this metabolic framework:")
    
    st.markdown("""
    1. **Ross CA, Tabrizi SJ (2011)** - *Huntington's disease: from molecular pathogenesis to clinical treatment.* (Lancet Neurol)
    2. **Saudou F, Humbert S (2016)** - *The Biology of Huntingtin.* (Neuron)
    3. **Bates GP, et al. (2015)** - *Huntington disease.* (Nat Rev Dis Primers)
    4. **Cheng ML, et al. (2016)** - *Metabolic disturbances in Huntington's disease.* (J Neurol Sci)
    5. **Reddy PH, et al. (2018)** - *Mitochondrial dysfunction and oxidative stress.* (BBA)
    """)
    
    st.info("💡 These resources were used to identify the target genes and metabolic nodes analyzed in this dashboard.")

# --- FOOTER ---
st.sidebar.markdown("---")
st.sidebar.caption("Data: KEGG API | Developed for PhD Portfolio")
