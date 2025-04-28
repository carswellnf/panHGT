import pandas
import matplotlib.pyplot as plot
import seaborn as sns

df = pandas.read_csv("./HGT_evidence/categorized_genes.tsv", sep="\t")

# Draw barplot function
def drawBarplot(df, x, y, title, xlabel, ylabel, filename):
    plot.figure(figsize=(12, 6))
    sns.barplot(x=x, y=y, data=df, palette="viridis", hue=x, legend=False)
    plot.xticks(rotation=45, ha='right')
    plot.title(title)
    plot.xlabel(xlabel)
    plot.ylabel(ylabel)
    plot.tight_layout()
    plot.savefig(filename, dpi=300)
    plot.clf()

# Count categories for all genes
category_counts = df["Functional_Category"].str.split(", ").explode().value_counts().reset_index()
category_counts.columns = ["Functional_Category", "Gene_Count"]

drawBarplot(
    df=category_counts,
    x="Functional_Category",
    y="Gene_Count",
    title="Gene Functional Categories",
    xlabel="Functional Category",
    ylabel="Gene Count",
    filename="./HGT_evidence/charts/bar/all_category_counts_bar.png"
)

# Count categories for strong HGT candidates
strong_df = df[df["Confidence"] == "strong"]
strong_category_counts = strong_df["Functional_Category"].str.split(", ").explode().value_counts().reset_index()
strong_category_counts.columns = ["Functional_Category", "Gene_Count"]

drawBarplot(
    df=strong_category_counts,
    x="Functional_Category",
    y="Gene_Count",
    title="Gene Functional Categories (Strong HGT Candidates)",
    xlabel="Functional Category",
    ylabel="Gene Count",
    filename="./HGT_evidence/charts/bar/strong_category_counts_bar.png"
)

# Count categories for moderate HGT candidates
moderate_df = df[df["Confidence"] == "moderate"]
moderate_category_counts = moderate_df["Functional_Category"].str.split(", ").explode().value_counts().reset_index()
moderate_category_counts.columns = ["Functional_Category", "Gene_Count"]

drawBarplot(
    df=moderate_category_counts,
    x="Functional_Category",
    y="Gene_Count",
    title="Gene Functional Categories (Moderate HGT Candidates)",
    xlabel="Functional Category",
    ylabel="Gene Count",
    filename="./HGT_evidence/charts/bar/moderate_category_counts_bar.png"
)

# Count categories for weak HGT candidates
weak_df = df[df["Confidence"] == "weak"]
weak_category_counts = weak_df["Functional_Category"].str.split(", ").explode().value_counts().reset_index()
weak_category_counts.columns = ["Functional_Category", "Gene_Count"]
drawBarplot(
    df=weak_category_counts,
    x="Functional_Category",
    y="Gene_Count",
    title="Gene Functional Categories (Weak HGT Candidates)",
    xlabel="Functional Category",
    ylabel="Gene Count",
    filename="./HGT_evidence/charts/bar/weak_category_counts_bar.png"
)

def drawPieChartsOnOneFigure(dataframes, titles, filename):
    num_charts = len(dataframes)
    rows = 2
    cols = 2
    fig, axes = plot.subplots(rows, cols, figsize=(cols * 6, rows * 5))
    for i, (df, title) in enumerate(zip(dataframes, titles)):
        ax = axes.flat[i] if num_charts > 1 else axes
        if len(df) > 7:
            top_categories = df.nlargest(7, "Gene_Count")
            other_count = df[~df["Functional_Category"].isin(top_categories["Functional_Category"])]["Gene_Count"].sum()
            top_categories = pandas.concat(
                [top_categories, pandas.DataFrame([{"Functional_Category": "Other", "Gene_Count": other_count}])],
                ignore_index=True
            )
            df = top_categories
        ax.pie(
            df["Gene_Count"],
            labels=df["Functional_Category"],
            autopct='%1.1f%%',
            textprops={'fontsize':12},
            startangle=90,
            colors=sns.color_palette("colorblind", n_colors=len(df["Functional_Category"])),
            wedgeprops=dict(edgecolor='w')
        )
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.axis('equal')  
    plot.tight_layout()
    fig.subplots_adjust(top=0.7)  
    plot.suptitle("Top 7 Gene Functional Categories", fontsize=21, fontweight='bold', y=0.85)
    plot.savefig(filename, dpi=800)
    plot.clf()

# Draw all pie charts on one figure
drawPieChartsOnOneFigure(
    dataframes=[category_counts, strong_category_counts, moderate_category_counts, weak_category_counts],
    titles=[
        "All (n = 11896)",
        "Strong (n = 1691)",
        "Moderate (n = 1137)",
        "Weak (n = 9068)"
    ],
    filename="./HGT_evidence/charts/pie/split_pie_charts.png"
)


# Data from the file
categories = [
    "Core",
    "Soft core",
    "Accessory",
    "Cloud"
]
values = [2039, 571, 5352, 34361] 

# Create a pie chart
plot.figure(figsize=(10, 10))
plot.pie(values, labels=categories, autopct='%1.1f%%', startangle=140, textprops={'fontsize': 18},colors=sns.color_palette("colorblind", n_colors=len(df["Functional_Category"])))
plot.title("Pangenome Gene Distribution", fontsize=24, fontweight='bold')
plot.axis('equal')
plot.savefig("./HGT_evidence/charts/pie/gene_distribution_pie_chart.png", dpi=300)