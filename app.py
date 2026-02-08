import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="HD Metabolic Framework", page_icon="🧬", layout="wide")

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
    else: return "🧬 Pathway Component"

# --- LOAD DATA ---
df = get_kegg_genes("hsa05016")
if not df.empty:
    df["Functional Role"] = df.apply(lambda row: assign_role(row["Symbol"], row["Description"]), axis=1)

# --- SIDEBAR ---
st.sidebar.image("https://cdn-icons-png.flaticon.com/512/822/822143.png", width=80)
st.sidebar.title("Researcher Profile")
st.sidebar.markdown(f"""
**Name:** Yashwant Nama  
**Target:** PhD in Neurogenetics  
**Focus:** Huntington's Disease (HD)  
---
""")

try:
    with open("CV_Yashwant_Nama_PhD_Application.pdf", "rb") as file:
        st.sidebar.download_button(label="📄 Download My CV", data=file, file_name="Yashwant_Nama_CV.pdf", mime="application/pdf")
except:
    st.sidebar.warning("Note: CV PDF not found in directory.")

st.sidebar.header("Project Progress")
st.sidebar.success("Phase 1: Data Acquisition ✅")
st.sidebar.success("Phase 2: Network Visualization ✅")
st.sidebar.info("Phase 3: Target Prioritization 🔄")

# --- MAIN CONTENT ---
st.title("🧬 Huntington's Disease (HD) Metabolic Framework")
st.markdown("### Disease Context: hsa05016")

tab1, tab2, tab3 = st.tabs(["📊 Target Discovery", "🕸️ Interaction Network", "📚 Literature"])

with tab1:
    col_a, col_b = st.columns([2, 1])
    with col_a:
        st.subheader("Genetic Components")
        search_query = st.text_input("🔍 Search genes or mechanisms (e.g. HTT, Autophagy):", placeholder="Type to filter...")
    with col_b:
        st.subheader("Deep Dive")
        selected_gene = st.selectbox("External Research:", ["Select a Gene"] + list(df['Symbol'].unique()))
        if selected_gene != "Select a Gene":
            st.markdown(f"**[View {selected_gene} on GeneCards ↗️](https://www.genecards.org/cgi-bin/carddisp.pl?gene={selected_gene})**")

    # Filtering Logic
    mask = df['Symbol'].str.contains(search_query.upper(), na=False) | \
           df['Description'].str.contains(search_query, case=False, na=False) | \
           df['Functional Role'].str.contains(search_query, case=False, na=False)
    
    filtered_df = df[mask] if search_query else df
    st.dataframe(filtered_df, use_container_width=True, height=300)

    # --- GENE PRIORITIZATION SCORING ---
    st.markdown("---")
    st.subheader("🎯 Therapeutic Target Prioritization")
    
    def calculate_score(row):
        score = 0
        if "Core" in row['Functional Role']: score += 5
        elif "Mitochondrial" in row['Functional Role']: score += 3
        elif "Apoptosis" in row['Functional Role']: score += 2
        elif "Synaptic" in row['Functional Role']: score += 2
        elif "Autophagy" in row['Functional Role']: score += 2
        return score + (len(row['Description']) % 3)

    df['Score'] = df.apply(calculate_score, axis=1)
    top_10 = df.sort_values('Score', ascending=False).head(10)

    c1, c2 = st.columns([1, 2])
    with c1:
        st.metric("Primary Target", top_10.iloc[0]['Symbol'])
        st.write("Highest metabolic impact score.")
        
        # EXCEL FIX: Use 'utf-8-sig' to ensure emojis/symbols show correctly in CSV
        csv_data = df.to_csv(index=False).encode('utf-8-sig')
        st.download_button(
            label="📥 Export Analysis (CSV)",
            data=csv_data,
            file_name='HD_Target_Analysis.csv',
            mime='text/csv'
        )
        
    with c2:
        fig_bar, ax_bar = plt.subplots(figsize=(8, 4))
        colors = ['#FF4B4B' if i < 3 else '#ff8a8a' for i in range(len(top_10))]
        ax_bar.barh(top_10['Symbol'], top_10['Score'], color=colors)
        ax_bar.invert_yaxis()
        ax_bar.set_title("Top 10 Ranked Therapeutic Targets")
        ax_bar.set_xlabel("Priority Score")
        plt.tight_layout()
        st.pyplot(fig_bar)

with tab2:
    st.subheader("Advanced Functional Interactome")
    st.write("Clustering genes by metabolic mechanism. Nodes sized by prioritization score.")
    
    G = nx.Graph()
    # Using top 50 genes to ensure a diverse network
    subset = df.sort_values('Score', ascending=False).head(50)
    
    role_colors = {
        "⭐ Core HD Gene": "#FF4B4B", 
        "🔋 Mitochondrial Dysfunction": "#FFA500", 
        "💀 Apoptosis": "#7D3C98", 
        "🧠 Synaptic / Excitotoxicity": "#2E86C1", 
        "♻️ Autophagy": "#28B463", 
        "🧬 Pathway Component": "#D5D8DC"
    }
    
    for _, row in subset.iterrows():
        G.add_node(row['Symbol'], role=row['Functional Role'], score=row['Score'])
        if row['Symbol'] != 'HTT':
            G.add_edge('HTT', row['Symbol'])

    fig_net, ax_net = plt.subplots(figsize=(12, 8))
    pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
    
    for role, color in role_colors.items():
        nodes = [n for n, attr in G.nodes(data=True) if attr.get('role') == role]
        if nodes:
            # Scale node size by score
            node_sizes = [G.nodes[n]['score'] * 200 for n in nodes]
            nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=color, 
                                   node_size=node_sizes, alpha=0.9, label=role.split(' ', 1)[1])

    nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color='grey')
    nx.draw_networkx_labels(G, pos, font_size=7, font_weight='bold')
    
    # Clean Legend
    leg = plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title="Mechanisms", fontsize='small')
    # Use legend_handles or legendHandles based on matplotlib version
    handles = getattr(leg, 'legend_handles', getattr(leg, 'legendHandles', []))
    for h in handles:
        h.set_sizes([100])
    
    plt.axis('off')
    plt.subplots_adjust(right=0.8)
    st.pyplot(fig_net)

with tab3:
    st.header("Research Bibliography")
    st.markdown("""
    1. **Ross CA, et al. (2011)** - *Huntington's disease: molecular pathogenesis to clinical treatment.* (Lancet Neurol)
    2. **Saudou F, et al. (2016)** - *The Biology of Huntingtin.* (Neuron)
    3. **Bates GP, et al. (2015)** - *Huntington disease.* (Nat Rev Dis Primers)
    4. **KEGG Pathway Database** - *hsa05016: Huntington disease.*
    """)
    st.info("💡 Computational analysis developed for PhD Application Portfolio.")

# --- FOOTER ---
st.sidebar.markdown("---")
st.sidebar.caption("Data: KEGG API | System: Streamlit")
