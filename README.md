🧬 NeuroMetabolic Framework
A Multi-Disease Systems Biology Platform for Neurodegenerative Disorders

The NeuroMetabolic Framework is an advanced computational platform designed for systems-level analysis of neurodegenerative and metabolic diseases.
It integrates KEGG pathway knowledge, statistical enrichment modeling, and functional interactome analysis to uncover shared and disease-specific molecular mechanisms across 13 neurological and metabolic pathologies, including Huntington’s disease, Alzheimer’s disease, Parkinson’s disease, ALS, and Type II Diabetes.

This framework is intended for research hypothesis generation, mechanism prioritization, and systems-guided experimental planning.

🔗 Live Application

👉 Interactive Research Dashboard
https://huntington-research-yash.streamlit.app/

🧠 Motivation

Traditional single-gene or pathway-isolated analyses often fail to capture the systems-level coupling between metabolic dysfunction, proteostasis failure, and neuronal vulnerability. Neurodegenerative diseases emerge from interacting molecular networks, not isolated genetic defects.

This project was developed as an independent systems biology initiative to:

Move beyond single-gene paradigms (e.g., HTT, APP) toward network-level disease modules

Quantify mechanism enrichment using statistical tests rather than narrative bias

Translate curated pathway knowledge into testable hypotheses for experimental prioritization

🚀 Core Features

Multi-Disease Comparative Analysis
Real-time KEGG pathway integration for 13 neurological and metabolic diseases via the KEGG REST API.

Mechanism-Level Functional Annotation
Automated classification of genes into biological processes such as:

Mitochondrial dysfunction

Proteostasis / proteasome stress

Autophagy

Apoptosis

Synaptic / excitotoxic signaling

Advanced Functional Interactome
Interactive network visualizations built with NetworkX, highlighting:

Functional coupling

Hub genes

Secondary controllers

Disease-specific bottlenecks

Statistical Enrichment Engine
Fisher’s Exact Test with multiple-testing correction to identify statistically overrepresented pathological mechanisms.

Automated Scientific Summaries
One-click generation of manuscript-ready biological interpretations to support hypothesis documentation.

🛠️ Research Stack

Language: Python 3.9+

Interface: Streamlit

Data Handling: Pandas, NumPy

Network Science: NetworkX

Statistics: SciPy (Fisher’s Exact Test)

Visualization: Matplotlib, Seaborn (publication-oriented plots)

📚 Data & Methodology

Primary Data Source:
KEGG Pathway Database (e.g., hsa05016, hsa05010, hsa04930)

Enrichment Strategy:
Over-representation analysis using Fisher’s Exact Test with pathway-level background correction.

Scoring Algorithm:
Weighted priority score combining:

Functional relevance (60%)

Literature prevalence (40%)

Network Construction:
Nodes represent genes/proteins.
Edges represent inferred functional coupling based on KEGG pathway co-occurrence.

⚠️ Note: Network edges do not imply direct physical protein–protein interactions or causal regulation.

📈 Research Roadmap

Phase 1: Multi-Disease Data Integration ✅

Phase 2: Statistical Enrichment Modeling ✅

Phase 3: STRING-DB PPI Integration 🔄 (in progress)

Phase 4: Differential Expression Overlay (GEO datasets) ⏳ (planned)

👤 Author

Yashwant Nama
Prospective PhD Researcher | Neurogenetics & Systems Biology

Research Focus:
Deciphering the metabolic–genetic axis of neurodegeneration through computational modeling and experimental validation.

⚠️ Disclaimer

This tool is intended for research hypothesis generation.
Network edges represent inferred functional associations derived from KEGG pathway co-occurrence and do not imply direct physical, regulatory, or causal protein–protein interactions.

Findings are intended to guide experimental prioritization, not replace wet-lab validation.

📖 Citation

If you use this framework for hypothesis generation or exploratory analysis, please cite:

NeuroMetabolic Framework – Y. Nama
KEGG Pathway Database
Fisher’s Exact Test–based enrichment analysis
