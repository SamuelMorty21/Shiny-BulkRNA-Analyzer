# BulkRNA Analysis Toolkit for Shiny

![Project Logo](link_to_logo.png) <!-- 如果有项目标志，请添加 -->

The BulkRNA Analysis Toolkit for Shiny is a comprehensive R-based package designed to simplify and enhance the analysis of bulk RNA sequencing data. It offers a user-friendly Shiny web application with interactive tools and visualizations for researchers and bioinformaticians working with gene expression data.

## Features

- **Count Matrix Generation:** Utilize featureCounts data to create count matrices for downstream analysis.
- **PCA Dimensionality Reduction:** Visualize data and reduce complexity using Principal Component Analysis.
- **Differential Gene Analysis:** Identify statistically significant differentially expressed genes.
- **Heatmap Visualization:** Create expressive heatmaps to visualize gene expression patterns.
- **Customizable Workflow:** Adapt the toolkit to your specific analysis needs with ease.

## Getting Started

To get started with the BulkRNA Analysis Toolkit, follow these steps:

1. **Installation:** Clone this repository and install the required R packages.

   ```bash
   git clone https://github.com/yourusername/your-repo.git
   ```
   
2. **Install the dependencies:** Open the ```BulkRNA.Rmd``` file and load the first section to install the dependencies:

- dplyr
- org.Hs.eg.db
- limma
- ggplot2
- circlize
- tibble
- pheatmap
- gridExtra
- ggrepel
- DESeq2
- ComplexHeatmap
- ggplotify
- clusterProfiler
- cowplot
- shiny

Within the document, you'll find instructions for installing the necessary R packages and setting up any additional configuration required for the analysis. Follow these instructions to ensure that all dependencies are correctly installed.

## Preparing Raw Data and Loading Sample Information

### Organizing Raw Data

1. Place your featureCounts output files in the project's root directory under a folder named `raw_material`. Ensure that the files contain the following standard columns:
   - Geneid
   - Chr
   - Start
   - End
   - Strand
   - Length
   - ./Fine-2-input/Fine-2-inputAligned.sortedByCoord.out.bam

2. Make sure that your gene IDs use `gene_name` as the identifier and do not use Ensembl IDs or any other identifiers.

### Loading Sample Information

1. To fill in sample information, execute the following command in R or RStudio:

   ```R
   # Load the sample information
   load_sample_info()
   ```

## Generating Count Matrix and TPM Matrix

1. Execute the following command to generate the count matrix and TPM matrix:

   ```R
   # Generate the count matrix and TPM matrix
   perform_count_matrix_and_logtpm()
   ```

   After running the command, the system will create a count_matrix folder in the root directory of your project.

   Inside the count_matrix folder, you will find the following files:

   ```count_matrix.csv```: This file contains the count matrix generated from the featureCounts data. It includes information on ```gene IDs```, ```chromosomes```, ```start``` and ```end positions```, ```strand```, ```gene length```, and the ```counts``` for each sample.
   
   ```TPM_matrix.csv```: This file contains the TPM (Transcripts Per Million) matrix generated using the TPM algorithm. It provides normalized expression values for each gene across samples.
   
Now, you have successfully generated the count matrix and TPM matrix necessary for your analysis. These files are located in the count_matrix folder and can be used for downstream analyses.

## Performing Principal Component Analysis (PCA)

To perform Principal Component Analysis (PCA) for dimensionality reduction and visualize sample differences, follow these steps:

1. Open your R or RStudio environment.

2. Execute the following command in the R console:

   ```R
   # Perform Principal Component Analysis (PCA)
   perform_PCA()
   ```
