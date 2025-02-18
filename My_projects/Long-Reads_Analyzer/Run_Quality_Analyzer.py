import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load ONT summary txt
summary_file = ""

def Graphs(summary_file:str):
    df = pd.read_csv(summary_file, sep="\t")

    # Check column name
    assert "read_len" in df.columns, "Comlumn 'read_len' not found"
    assert "mean_qscore" in df.columns, "Column 'mean_qscore' not found"

    # Set Seaborn style to visualize better
    sns.set_style(style="whitegrid")

    # Create a figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Read Length Distribution
    sns.histplot(df["read_len"], bins=100, kde=True, color="blue", ax=axes[0, 0])
    axes[0, 0].set_xlabel("Read Length (bp)")
    axes[0, 0].set_ylabel("Frequency")
    axes[0, 0].set_title("Read Length Distribution")

    # Mean Q-Score Distribution
    sns.histplot(df["mean_qscore"], bins=50, kde=True, color="green", ax=axes[0, 1])
    axes[0, 1].set_xlabel("Mean Q-Score")
    axes[0, 1].set_ylabel("Frequency")
    axes[0, 1].set_title("Mean Q-Score Distribution")

    # Read Length vs Q-Score Scatter Plot
    axes[1, 0].scatter(df["read_len"], df["mean_qscore"], alpha=0.5, s=2, color="red")
    axes[1, 0].set_xlabel("Read Length (bp)")
    axes[1, 0].set_ylabel("Mean Q-Score")
    axes[1, 0].set_title("Read Length vs Mean Q-Score")

    # Cumulative Yield Plot (Total vs Read Count)
    df_sorted = df.sort_values(by="read_len", ascending=False)
    df_sorted["cumulative_bases"] = np.cumsum(df_sorted["read_len"])
    axes[1, 1].plot(range(len(df_sorted)), df_sorted["cumulative_bases"], color="purple")
    axes[1, 1].set_xlabel("Read Count (Ordered by Length)")
    axes[1, 1].set_ylabel("Cumulative Bases")
    axes[1, 1].set_title("Cumulative Yield Plot")

    # Adjust layout
    plt.tight_layout()
    plt.show()
