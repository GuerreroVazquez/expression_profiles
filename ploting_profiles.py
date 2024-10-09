import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def get_gene_and_age(df, gene):
    gene_data = df[['Age', gene]].dropna()
    return gene_data

def plot_expression_over_age(gene_data, gene):
    plt.figure(figsize=(8, 6))
    plt.plot(gene_data['Age'], gene_data[gene], marker='o')
    plt.title(f"Expression of {gene} over Age")
    plt.xlabel("Age")
    plt.ylabel(f"Expression of {gene}")
    plt.grid(True)
    plt.show()

def scatter_plot_expression_over_age(gene_data, gene):
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=gene_data['Age'], y=gene_data[gene], hue=gene_data['Age'], palette='viridis', legend=False)
    plt.title(f"Expression of {gene} across Age")
    plt.xlabel("Age")
    plt.ylabel(f"Expression of {gene}")
    plt.grid(True)
    plt.show()

def violin_plot_grouped_by_age(gene_data, gene):
    # Define age groups
    bins = [18, 35, 65, 100]
    labels = ['Young', 'Middle Age', 'Old']
    gene_data['Age Group'] = pd.cut(gene_data['Age'], bins=bins, labels=labels, right=False)
    
    # Plot violin plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(x='Age Group', y=gene, data=gene_data)
    plt.title(f"Violin plot of {gene} grouped by Age")
    plt.xlabel("Age Group")
    plt.ylabel(f"Expression of {gene}")
    plt.grid(True)
    plt.show()

def gene_profile(df, gene):
    gene_data = get_gene_and_age(df, gene)  # Extract gene data
    scatter_plot_expression_over_age(gene_data, gene)  # Plot expression over age
    violin_plot_grouped_by_age(gene_data, gene)  # Make violin plots


# Function 4: Violin plot grouped by sex
def violin_plot_grouped_by_sex(gene_data, gene):

  
    # Plot violin plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(x='Sex', y=gene, data=gene_data)
    plt.title(f"Violin plot of {gene} grouped by Sex")
    plt.xlabel("Sex")
    plt.ylabel(f"Expression of {gene}")
    plt.grid(True)
    plt.show()

# Function 5: Violin plot grouped by experiment
def violin_plot_grouped_by_experiment(gene_data, gene):
    # Plot violin plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(x='Experiment', y=gene, data=gene_data)
    plt.title(f"Violin plot of {gene} grouped by Experiment")
    plt.xlabel("Experiment")
    plt.ylabel(f"Expression of {gene}")
    plt.grid(True)
    plt.show()

# Function 6: Violin plot grouped by sex and age group
def violin_plot_grouped_by_sex_and_age_group(gene_data, gene):
    
    # Define age groups
    bins = [18, 35, 65, 100]
    labels = ['Young', 'Middle Age', 'Old']
    gene_data['Age Group'] = pd.cut(gene_data['Age'], bins=bins, labels=labels, right=False)
    
    # Plot violin plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(x='Age Group', y=gene, hue='Sex', data=gene_data, split=True)
    plt.title(f"Violin plot of {gene} grouped by Sex and Age Group")
    plt.xlabel("Age Group")
    plt.ylabel(f"Expression of {gene}")
    plt.grid(True)
    plt.show()