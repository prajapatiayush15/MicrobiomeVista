import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------- Streamlit Page Setup -------------------
st.set_page_config(
    page_title="Gut Microbiome Dashboard",
    page_icon="🧬",
    layout="wide"
)

# ----------------- App Title -------------------
st.title("🧬 Gut Microbiome Dashboard")
st.subheader("Explore Taxonomy, Body Sites, and Gene Counts in Human Gut")

# ----------------- Load Dataset -------------------
@st.cache_data
def load_data():
    df = pd.read_csv(r"C:\Users\Ayush\gut_microbiome_project\GUTMICROBIOME\project_catalog.csv")
    df.rename(columns={'HMP Isolation Body Site': 'Body Site'}, inplace=True)  # Standardize column
    return df

df = load_data()

# ----------------- Sidebar Filter -------------------
selected_site = st.sidebar.selectbox("🧭 Choose a Body Site", df['Body Site'].dropna().unique())

# ----------------- Filtered Data -------------------
filtered_df = df[df['Body Site'] == selected_site]

# ----------------- Tabs for Sections -------------------
tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "📊 Top 15 Organisms",
    f"🧬 Organisms in {selected_site}",
    "🧬 Gene Count by Superkingdom",
    "🌡️ Heatmap: Body Site × Superkingdom",
    "🏥 Sequencing Centers",
    "🧬 Gene Count: Top 10 Organisms"
])


# ------- Tab 1: Top 15 Organisms Overall -------
with tab1:
    st.subheader("🌿 Top 15 Most Frequent Organisms (Overall)")

    top_organisms = df['Organism Name'].value_counts().nlargest(15)

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    sns.barplot(x=top_organisms.values, y=top_organisms.index, ax=ax1, palette="viridis")
    ax1.set_title("Top 15 Gut Microbiome Organisms", fontsize=14)
    ax1.set_xlabel("Frequency")
    ax1.set_ylabel("Organism Name")

    st.pyplot(fig1)

# ------- Tab 2: Top 10 Organisms in Selected Body Site -------
with tab2:
    st.subheader(f"🔍 Top 10 Organisms in **{selected_site}**")

    site_organisms = filtered_df['Organism Name'].value_counts().nlargest(10)

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    sns.barplot(x=site_organisms.values, y=site_organisms.index, ax=ax2, palette="crest")
    ax2.set_title(f"Top 10 Organisms in {selected_site}", fontsize=14)
    ax2.set_xlabel("Frequency")
    ax2.set_ylabel("Organism Name")

    st.pyplot(fig2)
    
#"🧬 Gene Count by Superkingdom"
with tab3:
    st.subheader("🧬 Gene Count Distribution Across Superkingdoms")

    # --- Clean and Prepare Data ---
    gene_df = df.copy()

    # Strip spaces, convert to string
    gene_df['NCBI Superkingdom'] = gene_df['NCBI Superkingdom'].astype(str).str.strip()

    # Remove bad entries
    gene_df = gene_df[
        (~gene_df['NCBI Superkingdom'].str.lower().str.contains("error")) &
        (~gene_df['NCBI Superkingdom'].str.lower().isin(['nan'])) &
        (gene_df['NCBI Superkingdom'].notna())
    ]

    # Remove rows with missing gene count
    gene_df = gene_df[gene_df['Gene Count'].notna()]
    gene_df = gene_df[gene_df['Gene Count'] > 0]

    # --- Plot Boxplot ---
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=gene_df,
        x='NCBI Superkingdom',
        y='Gene Count',
        palette="Set2"
    )
    ax3.set_title("Gene Count by NCBI Superkingdom", fontsize=14)
    ax3.set_xlabel("NCBI Superkingdom")
    ax3.set_ylabel("Gene Count")

    st.pyplot(fig3)

    #  Display summary stats
    st.markdown("### 📋 Average Gene Count by Superkingdom")
    avg_gene_counts = gene_df.groupby('NCBI Superkingdom')['Gene Count'].mean().round(2).reset_index()
    st.dataframe(avg_gene_counts)

    
#  Heatmap: Body Site × Superkingdom
with tab4:
    st.subheader("🌡️ Heatmap of Organism Counts by Body Site and Superkingdom")

    heatmap_df = df[['Body Site', 'NCBI Superkingdom']].dropna()

    # Clean the data
    heatmap_df['Body Site'] = heatmap_df['Body Site'].astype(str).str.strip()
    heatmap_df['NCBI Superkingdom'] = heatmap_df['NCBI Superkingdom'].astype(str).str.strip()

    # Filter invalid entries
    heatmap_df = heatmap_df[
        (~heatmap_df['Body Site'].str.lower().isin(['nan', 'error!!!'])) &
        (~heatmap_df['NCBI Superkingdom'].str.lower().isin(['nan', 'error!!!']))
    ]

    # Create pivot table
    heatmap_matrix = heatmap_df.pivot_table(
        index='Body Site',
        columns='NCBI Superkingdom',
        aggfunc='size',
        fill_value=0
    )

    # Plot
    fig4, ax4 = plt.subplots(figsize=(12, 8))
    sns.heatmap(heatmap_matrix, cmap="YlGnBu", annot=True, fmt='d', linewidths=0.5)
    ax4.set_title("Heatmap: Body Site vs NCBI Superkingdom", fontsize=14)
    ax4.set_xlabel("NCBI Superkingdom")
    ax4.set_ylabel("Body Site")

    st.pyplot(fig4)

    # Optional: download matrix
    st.download_button(
        label="📥 Download Heatmap Table (CSV)",
        data=heatmap_matrix.to_csv().encode('utf-8'),
        file_name="body_site_vs_superkingdom.csv",
        mime="text/csv"
    )

#Top Sequencing Centers by Project Count
with tab5:
    st.subheader("🏥 Top Sequencing Centers by Project Count")

    # Clean column
    seq_df = df.copy()
    seq_df['Sequencing Center'] = seq_df['Sequencing Center'].astype(str).str.strip()

    # Filter out invalid entries
    seq_df = seq_df[
        (~seq_df['Sequencing Center'].str.lower().isin(['nan', 'error', 'error!!!'])) &
        (seq_df['Sequencing Center'].notna())
    ]

    top_centers = seq_df['Sequencing Center'].value_counts().nlargest(15)

    fig5, ax5 = plt.subplots(figsize=(10, 6))
    sns.barplot(x=top_centers.values, y=top_centers.index, ax=ax5, palette="mako")
    ax5.set_title("Top 15 Sequencing Centers", fontsize=14)
    ax5.set_xlabel("Number of Projects")
    ax5.set_ylabel("Sequencing Center")

    st.pyplot(fig5)

    # Optional: Show raw table
    with st.expander("📄 Click to view raw center counts"):
        st.dataframe(top_centers.reset_index().rename(columns={
            'index': 'Sequencing Center',
            'Sequencing Center': 'Project Count'
        }))

    # Optional: Download CSV
    st.download_button(
        label="📥 Download Sequencing Center Data",
        data=top_centers.to_csv().encode('utf-8'),
        file_name="top_sequencing_centers.csv",
        mime="text/csv"
    )

# Gene Count comparison of top 10 organisms

with tab6:
    st.subheader("🧬 Gene Count Comparison of Top 10 Organisms")

    gene_df = df[['Organism Name', 'Gene Count']].dropna()

    # Clean names and remove outliers
    gene_df['Organism Name'] = gene_df['Organism Name'].astype(str).str.strip()
    gene_df = gene_df[gene_df['Gene Count'] > 0]

    # Group by Organism and average gene count
    top10_gene = (
        gene_df.groupby('Organism Name')['Gene Count']
        .mean()
        .sort_values(ascending=False)
        .head(10)
    )

    fig6, ax6 = plt.subplots(figsize=(10, 6))
    sns.barplot(x=top10_gene.values, y=top10_gene.index, ax=ax6, palette="coolwarm")
    ax6.set_title("Top 10 Organisms by Average Gene Count", fontsize=14)
    ax6.set_xlabel("Average Gene Count")
    ax6.set_ylabel("Organism")

    st.pyplot(fig6)

    # Optional: Show data table
    with st.expander("📋 View Data Table"):
        st.dataframe(top10_gene.reset_index().rename(columns={"Gene Count": "Avg Gene Count"}))

    # Optional: Download
    st.download_button(
        label="📥 Download Gene Count Data",
        data=top10_gene.to_csv().encode('utf-8'),
        file_name="top10_gene_counts.csv",
        mime="text/csv"
    )

# --------: Download Filtered Data --------
    st.download_button(
        label="📥 Download Filtered Data",
        data=filtered_df.to_csv(index=False).encode('utf-8'),
        file_name=f"gut_microbiome_{selected_site}.csv",
        mime="text/csv"
    )

# --------  Raw Dataset Preview --------
with st.expander("🔍 Click to preview raw dataset"):
    st.dataframe(df.head(20))



