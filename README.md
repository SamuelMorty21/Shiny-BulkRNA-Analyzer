# Shiny BulkRNA Analysis as a Toolkit

<!-- ![Project Logo](link_to_logo.png)--> 

The BulkRNA Analysis Toolkit is a comprehensive R-based package designed to simplify and enhance the analysis of bulk RNA sequencing data. It offers a user-friendly Shiny web application with interactive tools and visualizations for researchers and bioinformaticians working with gene expression data.

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

![GUI](https://github.com/SamuelMorty22/Shiny-BulkRNA-Analyzer/blob/main/screenshots/GUI/sample_info-2.png)

You can freely choose the root directory, then click 'Load Tag Information' to add the root directory to the environment variables and rename the samples by refilling the sample information.

## Generating Count Matrix and TPM Matrix

### Execute the following command to generate the count matrix and TPM matrix:

   ```R
   # Generate the count matrix and TPM matrix
   perform_count_matrix_and_logtpm()
   ```

   1. After running the command, the system will create a count_matrix folder in the root directory of your project.

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
   
![EXAMPLE](https://github.com/SamuelMorty22/Shiny-BulkRNA-Analyzer/blob/main/screenshots/figure/PCA.png) 


## Performing Differential Expression Analysis (DEGs)

To perform Differential Expression Analysis (DEGs) with different log2fold thresholds, you can use the `perform_DEGs` function. This analysis will identify genes that are differentially expressed compared to your specified log2fold threshold.

### 2-Fold DEGs

To identify genes with at least a 2-fold change in expression, execute the following command:

```R
# Perform 2-Fold DEGs analysis
perform_DEGs(log2fold = 1)
```
### 4-Fold DEGs

To identify genes with at least a 4-fold change in expression, execute the following command:

```R
# Perform 4-Fold DEGs analysis
perform_DEGs(log2fold = 2)
```
### 8-Fold DEGs

To identify genes with at least an 8-fold change in expression, execute the following command:

```R
# Perform 8-Fold DEGs analysis
perform_DEGs(log2fold = 3)
```

After running the analysis, the system will save the generated differential gene expression data in the ./results directory. The output will include the following parameters for each gene:

 - __baseMean__: The average count (average expression level) representing a gene's average expression across all samples.
   
 - __log2FoldChange__: The logarithmic fold change in gene expression, indicating the difference in expression between two sample groups, typically expressed in log2 form.
  
 - __lfcSE__: The standard error of the log fold change, representing the estimate's precision for differential expression.
   
 - __stat__: A statistical measure used to assess the degree of expression difference for a gene, often related to hypothesis testing.
   
 - __pvalue__: The p-value from a hypothesis test, measuring the significance of observed differences. Smaller p-values indicate more significant differences.
   
 - __padj__: The adjusted p-value (usually corrected for multiple testing), accounting for the false positive rate in the context of multiple comparisons. Typically used to determine significance.
  

You can utilize this data for further downstream analysis or visualize the differentially expressed genes.

## Generating a Heatmap

1. ### Input the genes you want to include in the heatmap

    ```R
    # Input genes for the heatmap
    gene_info()
    ```
    If you need to input multiple genes, you can separate them with spaces. For example, if you want to include genes like "DDX6," "DND1," and "NANOG," you can enter them as follows:

    __"DDX6 DND1 NANOG"__

    After executing the gene_info() function with your selected genes, you can proceed to create a heatmap based on the provided gene expression data.
   

3. ### Generate the heatmap object and choose your samples

    ```R
    # Generate a heatmap
    Heatmap_result <- perform_Heatmap()
    ```
    This will create a heatmap object based on the gene expression data you provided.
   
   You can further customize the heatmap by selecting specific columns or features from the heatmap object. To do this, you can use the ht_getcol() function:

   ```R
   # Select specific columns or features from the heatmap object
   ht_getcol()
   ```
   
   By customizing the heatmap, you can focus on the aspects of the data that are most relevant to your analysis or visualization needs.

4. ### Draw a Heatmap
   ```R
   heatmap_fig <- create_custom_heatmap(
     data = Heatmap_result[[1]][, selected_cols],
     col_data = Heatmap_result[[3]],
     genes_info = genes_info,
     show_column_names = FALSE,
     # rect_gp = gpar(col = "white", lwd = 2)
     rect_gp = gpar(col = NA)
   )
   pdf(paste0(file_path, "plots/Heatmap4.pdf"),width = Heatmap_result[[4]], height=Heatmap_result[[5]])
   draw(heatmap_fig)
   dev.off()
   ```
   
   This code will generate a customized heatmap and save it as a PDF file named "Heatmap4.pdf" in the specified directory (file_path/plots/). You can adjust the customization options and file path as needed to meet your visualization requirements.

   ![ht](https://github.com/SamuelMorty22/Shiny-BulkRNA-Analyzer/blob/main/screenshots/figure/heatmap.png)
   
## Generating a Dot Plot with Differential Gene Information

To create a dot plot that visualizes differential gene expression information, you can use the `perform_dotplot()` function. This function allows you to incorporate the data related to differentially expressed genes into your dot plot.

Here's an example of how to use `perform_dotplot()` to visualize differential gene expression:

```R
# Generate a dot plot with differential gene information
perform_dotplot()
```
![dotplot](https://github.com/SamuelMorty22/Shiny-BulkRNA-Analyzer/blob/main/screenshots/figure/dotplot.png)

## Performing Gene Ontology (GO) Term Analysis

To perform Gene Ontology (GO) term analysis with different log2fold thresholds, you can use the `perform_GOs()` function. This analysis helps you identify enriched GO terms associated with differentially expressed genes.

### GO Term Analysis with 2-Fold DEGs

To identify enriched GO terms for genes with at least a 2-fold change in expression, execute the following command:

```R
# Perform GO term analysis for 2-Fold DEGs
perform_GOs(log2fold = 1, mode = 2)
```
After running each analysis, the system will provide information on enriched GO terms associated with differentially expressed genes.

![go](https://github.com/SamuelMorty22/Shiny-BulkRNA-Analyzer/blob/main/screenshots/figure/GO.png)

## Summary
The BulkRNA Analysis Toolkit provides a powerful and user-friendly toolset for handling bulk RNA sequencing data. It offers features such as count matrix generation, PCA dimensionality reduction, differential gene analysis, heatmap visualization, and GO term analysis, enabling researchers to gain deep insights into gene expression data. You can customize workflows as needed and leverage interactive tools and visualizations for valuable insights. If you encounter any issues, refer to the documentation or seek support from the community. 

Wishing you a successful analysis!
