import streamlit as st
import pandas as pd

st.title("Huntington's Disease Research Dashboard")
st.write("Independent Research Project - Yashwant")

st.header("Project Overview")
st.info("This app analyzes the molecular pathways of Huntington's Disease using KEGG and OMIM data.")

# A simple table of key genes involved in HD
data = {
    'Gene Symbol': ['HTT', 'BDNF', 'CASP3', 'CREB1', 'TP53'],
    'Role in Disease': ['Primary Mutation (CAG Repeat)', 'Neurotrophic Factor', 'Apoptosis Executioner', 'Transcription Factor', 'Cell Cycle Regulation'],
    'Significance': ['High', 'High', 'Medium', 'Medium', 'Low']
}
df = pd.DataFrame(data)

st.subheader("Key Genes Identified")
st.table(df)

st.success("Data fetching pipeline: Active")
