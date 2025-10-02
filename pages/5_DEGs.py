import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

st.title("🧬 Differential Expression (Marker Genes)")

st.markdown("""
In this step, we identify **marker genes** that are differentially expressed between clusters.  

The workflow includes:
1. **Find marker genes** – perform differential expression tests between clusters.  
   - Comparison options:  
     - **One cluster vs all other clusters**  
     - **Two specific clusters**  
   - Test method: **Wilcoxon rank-sum test** (default in Scanpy).  

2. **Inspect marker genes** – view top ranked genes and their statistics (e.g., scores, fold-changes, adjusted p-values).  

👉 Differentially expressed genes help us assign **biological meaning** to clusters.
""")


# --- Check if adata exists from Step 5 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete Step 5 (Clustering & UMAP) first.")
    st.stop()

adata = st.session_state["adata"]

if "leiden" not in adata.obs:
    st.error("No clustering found. Please run Step 5 first.")
    st.stop()

# =========================================================
# --- Step 1: Select comparison mode ---
# =========================================================
st.subheader("📌 Step 1: Choose comparison mode")

mode = st.radio(
    "How would you like to compare clusters?",
    ["Cluster vs all other clusters", "Cluster vs cluster"],
    help="""
    **Cluster vs all other clusters** 🧩  
    - Identifies marker genes that are uniquely upregulated in one cluster compared to **all remaining cells**.  
    - Highlights *cluster-specific* features, making it more suitable for **cell type annotation**.  
    - Example: One cluster shows high expression of *MS4A1* → likely B cells.  
    - Best for answering: *“What cell type does this cluster represent?”*  

    **Cluster vs cluster** ⚖️  
    - Performs differential expression between **two selected clusters**.  
    - Highlights *relative differences*, making it more suitable for **subcluster comparison or functional studies**.  
    - Example: Comparing two T cell subsets to see which genes distinguish them.  
    - Best for answering: *“How are these two clusters different?”*  

    👉 Recommendation:  
    - Use *Cluster vs all others* for **cell type identification**.  
    - Use *Cluster vs cluster* for **subpopulation comparison and functional insights**.
    """
)

clusters = sorted(adata.obs["leiden"].unique())
run_deg = False
cluster = None
cluster1 = None
cluster2 = None

if mode == "Cluster vs all other clusters":
    cluster = st.selectbox("Select cluster:", clusters)
    if st.button("Run DE analysis"):
        run_deg = True

elif mode == "Cluster vs cluster":
    col1, col2 = st.columns(2)
    with col1:
        cluster1 = st.selectbox("Cluster 1:", clusters, index=0)
    with col2:
        cluster2 = st.selectbox("Cluster 2:", clusters, index=1)
    if st.button("Run DE analysis"):
        run_deg = True

elif mode == "All clusters vs rest":
    if st.button("Run DE analysis for all clusters"):
        run_deg = True

# =========================================================
# --- Step 2: Run DE analysis ---
# =========================================================
if run_deg: 
    wait_msg = st.empty()
    wait_msg.info("⏳ Running differential expression analysis...")

    if mode == "Cluster vs all other clusters":
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            groups=[cluster],
            reference="rest",
            method="wilcoxon"
        )
        st.success(f"✅ Marker genes for cluster {cluster} vs all others computed.")

    elif mode == "Cluster vs cluster":
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            groups=[cluster1],
            reference=cluster2,
            method="wilcoxon"
        )
        st.success(f"✅ Marker genes for cluster {cluster1} vs cluster {cluster2} computed.")


    elif mode == "All clusters vs rest":
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            reference="rest",
            method="wilcoxon"
        )
        st.success("✅ Marker genes computed for ALL clusters vs rest.")

    wait_msg.empty()
    st.session_state["adata"] = adata

# =========================================================
# --- Step 2 conti.: Show results ---
# =========================================================
if "rank_genes_groups" in adata.uns:
    st.subheader("📊 Step 2: Inspect marker genes")

    st.markdown("Here are the **top ranked marker genes** per cluster (Scanpy visualization):")
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
    fig = plt.gcf()
    st.pyplot(fig)
    plt.close(fig)

    # Convert results to DataFrame
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    dfs = []
    for g in groups:
        df = pd.DataFrame({
            "names": result["names"][g],
            "scores": result["scores"][g],
            "logfoldchanges": result["logfoldchanges"][g],
            "pvals": result["pvals"][g],
            "pvals_adj": result["pvals_adj"][g],
        })
        df["cluster"] = g
        dfs.append(df)
    df_out = pd.concat(dfs)

    st.download_button(
        label="💾 Download DE results (.csv)",
        data=df_out.to_csv(index=False).encode("utf-8"),
        file_name="DE_results.csv",
        mime="text/csv"
    )

# =========================================================
# --- Step 3: Auto detect marker genes ---
# =========================================================
if "rank_genes_groups" in adata.uns:
    st.subheader("✨ Step 3: Automatic Marker Gene Detection")

    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    # Collect top N marker genes per cluster
    top_n = st.number_input("Number of top genes per cluster", min_value=1, max_value=50, value=5, step=1)

    marker_dict = {}
    for g in groups:
        marker_dict[g] = result["names"][g][:top_n].tolist()

    # Convert to DataFrame for display
    marker_table = []
    for cluster, genes in marker_dict.items():
        marker_table.append({
            "Cluster": cluster,
            "Markers": ", ".join(genes)
        })
    df_markers = pd.DataFrame(marker_table)

    st.dataframe(df_markers)

    # Save marker gene list
    st.download_button(
        label="💾 Download marker genes (.csv)",
        data=df_markers.to_csv(index=False).encode("utf-8"),
        file_name="marker_genes.csv",
        mime="text/csv"
    )

    st.info("💡 These top marker genes are extracted from the DE results automatically. You can use them for cell type annotation in the next step.")


# --- Save top marker per cluster ---
if "rank_genes_groups" in adata.uns:
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    
    top_markers = {}
    for g in groups:
        if len(result["names"][g]) > 0:
            top_markers[g] = result["names"][g][0]  # take top 1 per cluster
    
    st.session_state["top_markers"] = top_markers
    st.info(f"💡 Saved top marker genes per cluster: {list(top_markers.values())}")


#if "rank_genes_groups" in adata.uns:
#    st.subheader("📊 Step 2: Inspect marker genes")
#
#    st.markdown("Here are the **top ranked marker genes** per cluster (Scanpy visualization):")
#   sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
#   fig = plt.gcf()
#    st.pyplot(fig)
#   plt.close(fig)

    # Convert results to DataFrame
#    result = adata.uns["rank_genes_groups"]
#    groups = result["names"].dtype.names
#    dfs = []
#    for g in groups:
#       df = pd.DataFrame({
#           "names": result["names"][g],
#             "scores": result["scores"][g],
#            "logfoldchanges": result["logfoldchanges"][g],
#            "pvals": result["pvals"][g],
#            "pvals_adj": result["pvals_adj"][g],
#       })
#        df["cluster"] = g
#        dfs.append(df)
#    df_out = pd.concat(dfs)

#    st.download_button(
#        label="💾 Download DE results (.csv)",
#        data=df_out.to_csv(index=False).encode("utf-8"),
#        file_name="DE_results.csv",
#        mime="text/csv"
#    )

# # =========================================================
# # --- Step 4: Explore marker genes with violin plots ---
# # =========================================================
# st.subheader("🎻 Step 3: Explore marker genes with violin plots")

# st.markdown("""
# Violin plots show the **distribution of expression levels** of selected genes across clusters:  
# - The **width** of the violin = number of cells at that expression level.  
# - The **height** = range of expression.  
# - This allows you to compare how a gene is expressed in different clusters.  

# 👉 Useful for verifying whether a marker gene truly separates specific cell populations.
# """)

# # 提示框，提供推荐 marker genes
# st.info("""
# 💡 **Tip:** Try common marker genes in PBMC data:
# - **CST3** → dendritic cell / monocyte marker  
# - **NKG7** → NK cell / cytotoxic T cell marker  
# - **PPBP** → Platelet marker 
# - **MS4A1** → B cell marker  
# - **CD3D** → T cell marker  
 
# """)

# # Default marker list
# marker_genes = ["CST3", "NKG7", "PPBP"]

# # Select genes
# violin_genes = st.multiselect(
#     "Select one or more genes for violin plot:",
#     options=adata.var_names.tolist(),
#     default=marker_genes,
#     help="Choose genes to visualize expression distributions across clusters."
# )


# if st.button("Plot violin plots"):
#     for gene in violin_genes:
#         st.subheader(f"Violin plot: {gene}")
#         sc.pl.violin(
#             adata,
#             keys=gene,
#             groupby="leiden",
#             show=False
#         )
#         fig = plt.gcf()
#         st.pyplot(fig)
#         plt.close(fig)


    # # =========================================================
    # # --- Step 4: Visualize marker genes on UMAP ---
    # # =========================================================
    # st.subheader("🗺️ Step 3: Visualize marker genes on UMAP")

    # st.markdown("""
    # You can plot expression of selected genes on the UMAP embedding.  
    # - **Dark cells** → low expression  
    # - **Bright cells** → high expression  

    # 👉 This helps connect clusters to **cell types**.
    # """)

    # # 🔹 提示框
    # st.info("""
    # 💡 **Tip:** Marker genes are genes whose expression highlights specific cell types.  
    # Here are some commonly used marker genes in PBMC data:

    # - **CST3** → dendritic cell / monocyte marker  
    # - **NKG7** → NK cell / cytotoxic T cell marker  
    # - **MS4A1** → B cell marker  
    # - **CD3D** → T cell marker  
    # - **PPBP** → Platelet marker  
    # """)

    # # Simplified default PBMC marker list
    # marker_genes = [
    #     "CST3",   # dendritic cell / monocyte marker
    #     "NKG7",   # NK cell / cytotoxic T cell marker
    #     "MS4A1",  # B cell marker
    #     "CD3D",   # T cell marker
    #     "PPBP"    # Platelet marker
    # ]


    # selected_genes = st.multiselect(
    #     "Select marker genes to visualize on UMAP:",
    #     options=adata.var_names.tolist(),
    #     default=marker_genes,
    #     help="Choose one or more genes to display on UMAP."
    # )

    # if st.button("Plot selected marker genes"):
    #     sc.pl.umap(adata, color=selected_genes, show=False)
    #     fig = plt.gcf()
    #     st.pyplot(fig)
    #     plt.close(fig)

# --- Show "Next" link only after DE is done ---
    st.markdown("""
    <style>
    section[data-testid="stMain"] [data-testid="stPageLink"] a,
    section[data-testid="stMain"] [data-testid="stPageLink"] p {
      font-style: italic !important;
    }
    </style>
    """, unsafe_allow_html=True)

    spacer, right = st.columns([0.45, 0.2], gap="small")
    with right:
        st.page_link("pages/6_Assign_Cell_Type_Identity.py", label="➡️ Next: Assign Cell Identity")
        
