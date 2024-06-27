# DDX3X_SR_2023
Code used for analyses in "Post-transcriptional cross- and auto-regulation buffer expression 
of the human RNA helicases DDX3X and DDX3Y"

kallisto_alignment: kallisto transcriptome alignment used for RNAseq data (shell)
outputTPMS: Generates TPMs for RNAseq batch using tximport (R)
 
Depmap_analysis: Selects XY cell lines from the Cancer Dependency Map expression data and plots DDX3X and DDX3Y TPM (Python).
Fig1Scatterplot: Plots 2D scatterplot for X and Y chromosome dosage-sensitivity metrics (Python).
Fig2AAneuploidyplots: Extracts expression data from sex chromosome aneuploidy expression data and linearly models X and Y chromosome gene expression (ggplot,R). 
Fig3stripplot: Plots DDX3X and DDX3Y TPMs from XY and Azfa expression data (Python).
Fig4DifferentialGeneExpression: DESeq2 approach to normalize data and identify DEGs in DDX3X and DDX3Y knockdowns (R).
Fig6curvefit: Fits expression data from metabolic labeling experiments to obtain RNA half-life (scipy.stats, Python). 

Proteomics_DDX3X_DDX3Y_analysis: MaxFLQ approach to normalize and analyze protein abundance performed by Jason Derks (R). 
proteomics_functions: Functions to analyze proteomic data (Jason Derks)
