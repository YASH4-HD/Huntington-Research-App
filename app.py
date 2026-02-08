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

tab1, tab2, tab3 = st.tabs(["📊 Target Discovery", "🕸️ Interaction Network", "🔬 Enrichment & Lit"])

with tab1:
    col_a, col_b = st.columns([2, 1])
    with col_a:
        st.subheader("Genetic Components")
        search_query = st.text_input("🔍 Search genes or mechanisms:", placeholder="Type to filter...")
    with col_b:
        st.subheader("Deep Dive")
        selected_gene = st.selectbox("External Research:", ["Select a Gene"] + list(df['Symbol'].unique()))
        if selected_gene != "Select a Gene":
            st
