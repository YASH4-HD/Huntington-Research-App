import streamlit as st
import pandas as pd
from bioservices import KEGG

st.set_page_config(page_title="HD Research Lab", layout="wide")

st.title("🧬 Huntington's Disease Molecular Analysis")
st.write("Real-time data from KEGG Database")

# Initialize KEGG service
@st.cache_data # This makes the app fast
def fetch_real_hd_data():
    k = KEGG()
    data = k.get("hsa05016") # The ID for Huntington's Disease
    parsed = k.parse(data)
    genes = parsed.get('GENE', {})
    
    gene_list = []
    for gene_id, gene_info in genes.items():
        symbol = gene_info.split(';')[0]
        description = gene_info.split(';')[1] if ';' in gene_info else ""
        gene_list.append({"ID": gene_id, "Symbol": symbol, "Function": description})
    return pd.DataFrame(gene_list)

try:
    with st.spinner('Fetching live genomic data...'):
        df = fetch_real_hd_data()
    
    st.success(f"Successfully retrieved {len(df)} genes associated with HD.")
    
    # Add a search bar
    search = st.text_input("Search for a specific gene (e.g., HTT or CASP3):")
    if search:
        df = df[df['Symbol'].str.contains(search.upper())]

    st.dataframe(df, use_container_width=True)

except Exception as e:
    st.error(f"Database connection error: {e}")

st.sidebar.markdown("### Project Status")
st.sidebar.info("Phase 1: Data Acquisition (Complete)")
st.sidebar.warning("Phase 2: Network Analysis (In Progress)")
