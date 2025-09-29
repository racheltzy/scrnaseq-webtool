import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt

st.title("🔗 Clustering & Non-linear Dimensional Reduction (UMAP)")

st.markdown("""
In this step, we will group cells into **clusters** and visualize them in a **low-dimensional space**.

The workflow includes:
1. **Build the neighborhood graph** – define which cells are similar based on their PCA representation.  
2. **Cluster the cells (Leiden algorithm)** – partition the graph into cell groups.  
3. **Run UMAP** – embed the cells into 2D for visualization.  

*👉 The selected number of PCs from PCA strongly influences clustering and UMAP.  
It is recommended to use the same `n_pcs` as chosen in the previous step.*  
""")

# --- Check if adata exists from Step 4 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete Step 4 (PCA) first.")
    st.stop()

adata = st.session_state["adata"]

# =========================================================
# --- Neighborhood graph ---
# =========================================================
st.subheader("📌 Step 1: Build the neighborhood graph")

if "n_pcs" in st.session_state:
    st.info(f"👉 You selected **{st.session_state['n_pcs']} PCs** in the previous step.")
elif "n_pcs_selected" in adata.uns:
    st.info(f"👉 From saved file: using **{adata.uns['n_pcs_selected']} PCs**.")
else:
    st.warning("⚠️ No PC selection found. Default settings may be used.")

n_pcs = st.number_input(
    "Number of PCs to use for neighbors:",
    min_value=5, max_value=100, value=st.session_state.get("n_pcs", 20), step=5
)

if st.button("Run neighbors"):
    wait_msg = st.empty()
    wait_msg.info("⏳ Computing neighbors...")

    sc.pp.neighbors(adata, n_pcs=int(n_pcs))

    # ✅ Save PC selection so downstream always finds it
    adata.uns["n_pcs_selected"] = int(n_pcs)
    st.session_state["n_pcs"] = int(n_pcs)

    wait_msg.empty()
    st.success(f"✅ Nearest-neighbor graph computed (using {n_pcs} PCs).")
    st.session_state["adata"] = adata


# =========================================================
# --- Clustering + UMAP ---
# =========================================================
st.subheader("📌 Step 2: Leiden clustering and UMAP visualization")

resolution = st.number_input(
    "Leiden resolution (higher = more clusters):",
    min_value=0.1, max_value=5.0, value=1.0, step=0.1, format="%.1f",
    help="""
        The **resolution** parameter controls cluster granularity:

        - Lower values (e.g., 0.4) → fewer, larger clusters.  
        - Higher values (e.g., 1.0–2.0) → more, smaller clusters.  
        - For datasets around ~3k cells, values between 0.4–1.2 are often a good starting point.  

        (Default = 1.0)
        """
)

if st.button("Run clustering and UMAP"):
    try:
        wait_msg = st.empty()
        wait_msg.info("⏳ Running Leiden clustering and UMAP...")

        # Nearest neighbors graph (based on PCA)
        sc.pp.neighbors(adata, n_pcs=int(n_pcs))

        # Leiden clustering
        sc.tl.leiden(
            adata,
            resolution=float(resolution),
            random_state=0,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )

        # UMAP embedding
        sc.tl.umap(adata, random_state=0)

        wait_msg.empty()
        st.success(f"✅ Leiden clustering (resolution={resolution}, PCs={n_pcs}) and UMAP complete.")

        # Show number of clusters
        n_clusters = adata.obs["leiden"].nunique()
        st.write(f"**Number of clusters:** {n_clusters}")
        st.dataframe(adata.obs["leiden"].value_counts().rename("Cell count"))

        # UMAP plot (clusters)
        st.subheader("UMAP with Leiden clusters")
        sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False)
        fig = plt.gcf()
        st.pyplot(fig)
        plt.close(fig)

        # Save back to session
        st.session_state["adata"] = adata

    except ImportError as e:
        wait_msg.empty()
        st.error("Leiden requires the `igraph` and `leidenalg` packages. "
                 "Please install them via conda or pip.")
        st.exception(e)
