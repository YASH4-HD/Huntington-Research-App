import streamlit as st
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -------------------------------------------------
# PAGE CONFIG
# -------------------------------------------------
st.set_page_config(
    page_title="HD Metabolic Framework",
    page_icon="🧠",
    layout="wide"
)

# -------------------------------------------------
# DATA ACQUISITION (KEGG API)
# -------------------------------------------------
@st.cache_data
def get_kegg_genes(pathway_id):
    url = f"https://rest.kegg.jp/get/{pathway_id}"
    response = requests.get(url)
    genes = []

    if response.status_code == 200:
        lines = response.text.split('\n')
        is_gene_section = False

        for line in lines:
            if line.startswith("GENE"):
                is_gene_section = True
                line = line.replace("GENE", "").strip()
            elif line.startswith(("COMPOUND", "REFERENCE", "AUTHORS")):
                is_gene_section = False

            if is_gene_section and ";" in line:
                left, desc = line.split(";", 1)
                parts = left.strip().split(None, 1)
                if len(parts) == 2:
                    genes.append({
                        "ID": parts[0],
                        "Symbol": parts[1],
                        "Description": desc.strip()
                    })

    return pd.DataFrame(genes)

# -------------------------------------------------
# BIOLOGICAL ROLE ASSIGNMENT
# -------------------------------------------------
CORE_HD_GENES = {"HTT", "BDNF", "CASP3", "CREB1", "TP53", "SOD1", "PPARGC1A"}

def assign_role(symbol, desc):
    d = desc.lower()

    if symbol in CORE_HD_GENES:
        return "Core HD Gene"
    elif any(k in d for k in ["mitochond", "atp", "oxidative", "respiratory"]):
        return "Mitochondrial Dysfunction"
    elif any(k in d for k in ["caspase", "apopt"]):
        return "Apoptosis"
    elif "autophagy" in d:
        return "Autophagy"
    elif any(k in d for k in ["synap", "glutamate", "receptor"]):
        return "Synaptic / Excitotoxicity"
    else:
        return "Pathway Component"

# -------------------------------------------------
# LOAD DATA
# -------------------------------------------------
df = get_kegg_genes("hsa05016")
df["Functional Role"] = df.apply(
    lambda r: assign_role(r["Symbol"], r["Description"]), axis=1
)

# -------------------------------------------------
# SIDEBAR
# -------------------------------------------------
st.sidebar.image(
    "https://cdn-icons-png.flaticon.com/512/822/822143.png", width=90
)
st.sidebar.title("Researcher: Yashwant Nama")
st.sidebar.info(
    "**Target:** PhD in Neurogenetics\n\n"
    "**Focus:** Huntington's Disease (HD)"
)

st.sidebar.markdown("---")
st.sidebar.success("Phase 1: Data Acquisition ✅")
st.sidebar.success("Phase 2: Network Visualization ✅")
st.sidebar.info("Phase 3: Functional Interpretation 🔄")

# -------------------------------------------------
# MAIN HEADER
# -------------------------------------------------
st.title("🧬 Huntington’s Disease (HD) Metabolic Framework")
st.markdown("""
Computational dissection of the **KEGG hsa05016 pathway**, integrating  
gene function, metabolic disruption, and interaction topology.
""")

# -------------------------------------------------
# TABS
# -------------------------------------------------
tab1, tab2, tab3 = st.tabs([
    "📊 Gene Data",
    "🕸️ Functional Interactome",
    "📚 Literature"
])

# -------------------------------------------------
# TAB 1 — DATA TABLE
# -------------------------------------------------
with tab1:
    st.subheader("Genetic Components & Functional Roles")

    gene = st.selectbox(
        "Select a gene for external reference:",
        ["None"] + sorted(df["Symbol"].unique())
    )

    if gene != "None":
        st.info(
            f"🔗 **GeneCards:** "
            f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}"
        )

    query = st.text_input("Search (gene / mechanism / role):")

    if query:
        mask = (
            df["Symbol"].str.contains(query, case=False, na=False)
            | df["Description"].str.contains(query, case=False, na=False)
            | df["Functional Role"].str.contains(query, case=False, na=False)
        )
        st.dataframe(df[mask], use_container_width=True)
    else:
        st.dataframe(df, use_container_width=True)

# -------------------------------------------------
# TAB 2 — ADVANCED INTERACTOME
# -------------------------------------------------
with tab2:
    st.subheader("Advanced Functional Interactome")
    st.caption("Nodes cluster by biological mechanism; edges represent functional proximity.")

    role_colors = {
        "Core HD Gene": "#E74C3C",
        "Mitochondrial Dysfunction": "#F39C12",
        "Apoptosis": "#8E44AD",
        "Synaptic / Excitotoxicity": "#2980B9",
        "Autophagy": "#27AE60",
        "Pathway Component": "#BDC3C7",
    }

    G = nx.Graph()
    subset = df.head(35)

    # Nodes
    for _, r in subset.iterrows():
        G.add_node(r["Symbol"], role=r["Functional Role"])

    # Edges
    for _, r in subset.iterrows():
        if r["Symbol"] != "HTT":
            G.add_edge("HTT", r["Symbol"], weight=1)

    for i, r1 in subset.iterrows():
        for j, r2 in subset.iterrows():
            if i < j and r1["Functional Role"] == r2["Functional Role"]:
                if r1["Functional Role"] != "Pathway Component":
                    G.add_edge(r1["Symbol"], r2["Symbol"], weight=0.4)

    pos = nx.spring_layout(G, seed=42, k=0.5)

    fig, ax = plt.subplots(figsize=(13, 9))

    for role, color in role_colors.items():
        nodes = [n for n, a in G.nodes(data=True) if a["role"] == role]
        if nodes:
            nx.draw_networkx_nodes(
                G, pos,
                nodelist=nodes,
                node_color=color,
                node_size=2600 if role == "Core HD Gene" else 1100,
                alpha=0.9
            )

    nx.draw_networkx_edges(G, pos, alpha=0.4, edge_color="gray")
    nx.draw_networkx_labels(G, pos, font_size=9, font_weight="bold")

    # ---- CLEAN LEGEND (FIXED) ----
    legend_items = [
        Line2D([0], [0], marker='o', color='w',
               label=role,
               markerfacecolor=color,
               markersize=12)
        for role, color in role_colors.items()
    ]

    plt.legend(
        handles=legend_items,
        title="Functional Mechanisms",
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        frameon=True
    )

    plt.axis("off")
    st.pyplot(fig)

    st.success(
        "✔ This visualization now reflects **functional modularity**, "
        "not a simple HTT-centric star."
    )

# -------------------------------------------------
# TAB 3 — LITERATURE
# -------------------------------------------------
with tab3:
    st.subheader("Research Bibliography")

    st.markdown("""
1. **Ross & Tabrizi (2011)** – *Lancet Neurology*  
2. **Saudou & Humbert (2016)** – *Neuron*  
3. **Bates et al. (2015)** – *Nature Reviews Disease Primers*  
4. **Cheng et al. (2016)** – *Journal of Neurological Sciences*  
5. **Reddy et al. (2018)** – *BBA – Molecular Basis of Disease*
""")

    st.info(
        "These studies informed functional annotation and pathway interpretation."
    )

# -------------------------------------------------
# FOOTER
# -------------------------------------------------
st.sidebar.markdown("---")
st.sidebar.caption("Data source: KEGG | Designed for PhD portfolio")
