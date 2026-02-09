🧬 NeuroMetabolic Framework: Multi-Disease Systems Biology
An advanced computational framework for the systems-level analysis of neurodegenerative and metabolic pathologies. This platform integrates KEGG pathway data with statistical enrichment models and functional interactomes to explore the molecular architecture of 13 target diseases, including Huntington’s, Alzheimer’s, and Type II Diabetes.

🔗 Live Application
👉 Interactive Research Dashboard: https://huntington-research-yash.streamlit.app/

🧠 Project Motivation
Traditional genomic analysis often overlooks the complex functional coupling between metabolic failure and proteotoxicity. This project was developed as an independent systems biology initiative to:

Decipher Pathological Hubs: Transition from single-gene focus (e.g., HTT, APP) to multi-node network analysis.
Quantify Mechanism Enrichment: Use Fisher’s Exact Test to statistically identify whether a pathology is driven by mitochondrial dysfunction, proteostasis failure, or synaptic excitotoxicity.
Bridge Data Gaps: Synthesize KEGG-derived pathway architectures into actionable research hypotheses for experimental prioritization.
🚀 Core Features
Dynamic Multi-Disease Analysis: Real-time data acquisition for 13 distinct pathologies via KEGG REST API.
Functional Role Annotation: Automated classification of genes into biological clusters (e.g., Autophagy, Apoptosis, Mitochondrial Energy Metabolism).
Advanced Functional Interactome: Interactive graph visualizations using NetworkX, mapping inferred functional coupling and hub-node topology.
Statistical Enrichment Engine: Integrated Fisher’s Exact Test with p-value adjustment to identify statistically significant pathological drivers.
Automated Manuscript Generation: One-click generation of scientific summaries for hypothesis documentation.
🛠️ Research Stack
Language: Python 3.9+
Interface: Streamlit (High-performance UI)
Data Engineering: Pandas, NumPy
Network Science: NetworkX
Statistics: SciPy (Fisher’s Exact Test)
Visualization: Matplotlib, Seaborn
📚 Data & Methodology
Primary Source: KEGG Pathway Database (e.g., hsa05016, hsa05010, hsa04930).
Scoring Algorithm: A weighted priority score (60% Functional Role / 40% Literature Prevalence) used to rank candidates for therapeutic targeting.
Network Topology: Nodes represent proteins/genes; edges represent functional co-occurrence within conserved metabolic pathways.
📈 Research Roadmap
Phase 1: Multi-Disease Integration ✅ (Complete)
Phase 2: Statistical Enrichment Modeling ✅ (Complete)
Phase 3: STRING-DB PPI Integration 🔄 (In Progress)
Phase 4: Differential Expression Mapping (GEO Data) ⏳ (Planned)
👤 Author
Yashwant Nama
Prospective PhD Researcher | Neurogenetics & Systems Biology
Research Focus: Deciphering the metabolic-genetic axis in neurodegeneration through computational modeling and wet-lab validation.

💡 How to Cite
If using this framework for hypothesis generation, please credit the NeuroMetabolic Framework (Y. Nama) and the KEGG API.
