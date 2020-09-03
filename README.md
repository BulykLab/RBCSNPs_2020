# RBCSNPs_2020
This GitHub repository contains the code (Python, R) and data files needed to support the PBM-related analyses in "Common Variants in Signaling Transcription Factor Binding Sites Drive the Majority of Red Blood Cell Traits" (Choudhuri*, Trompouki*, Abraham* et al., accepted (*equal contributions)).

Updated 3 September 2020. Kian Hong Kock (kianhongkock@g.harvard.edu), Martha Bulyk Lab (mlbulyk@genetics.med.harvard.edu).

### Code files ###

1) "NatGenet2020_Fig6_FigS6_PythonCode.ipynb" is a Python Jupyter notebook which contains code for analysing RBC SNPs for perturbed TF binding, using PBM data. This uses Python 3.6; you'll need to have Biopython and pandas installed.

2) "NatGenet2020_Fig6_FigS6_RCode.R" is an R script which contains code for generating standardised PBM E-score files from processed 8-mer E-score files (e.g., from UniPROBE and CIS-BP), and for performing bootstrapping analysis (empirical distributions) to evaluate the statistical significance of the frequency of perturbed TF binding events in RBC SNPs, versus a background set of common (allele frequency >10%) SNPs from dbSNP. You'll need to have ggplot2 installed.

3) "NatGenet2020_RBCSNPs_individualPBM_binary_10percentSNPs.py" is a Python script, intended to run on a compute cluster with GRCh38 chromosome FASTA files, PBM files and recoded dbSNP VCF file with allele frequencies, for obtaining the background set of perturbed TF binding events in common (allele frequency >10%) SNPs from dbSNP.

4) "vcftools_commonSNPs.sh" is a shell script calling vcftools for generating a recoded VCF from "00-common_all.vcf", a dbSNP VCF file of common human variants.

### Input and intermediate output files and directories ###

5) "RBC_SNPs_single_30072018.txt", "Signalling_Ctr_SNPs_04092018.txt", "FHS_SNP_SignalingCenters_new.xlsx", "absent_SigCtr_RBC_SNPs_caf.recode.vcf" are input data files for "NatGenet2020_Fig6_FigS6_PythonCode.ipynb". These correspond to the list of RBC SNPs of interest, the list of Signalling Centre SNPs of interest, data from the Framingham Heart Study on cis-eQTLs (genes/exons), and allele frequency information from Signalling Centre SNPs that were not present in the original set of RBC SNPs, respectively.

6) "foreground_3263SNPs.txt" is the foreground output file, obtained using "NatGenet2020_Fig6_FigS6_PythonCode.ipynb" and subsequently processed for the rest of the analysis.

7) The "empirical_background" directory contains the background output files generated using "NatGenet2020_RBCSNPs_individualPBM_binary_10percentSNPs.py".

8) The "PBM_datasets" directory contains the processed PBM data files from UniPROBE and CIS-BP.

9) The "Individual_PBMs_Escores" directory contains the standardised PBM E-score files generated from processed 8-mer E-score files using "NatGenet2020_Fig6_FigS6_RCode.R". Note that the PBM files generated from CIS-BP (HES1, HES5, HES7 and TCF4) contain only one of the 8-mer orientations; the "NatGenet2020_Fig6_FigS6_PythonCode.ipynb" Python Jupyter notebook takes this into consideration and fixes this internally by considering both the forward and reverse 8-mer orientations in generating the dictionary versions of the individual PBM E-score datasets.

10) The "rsID_Escores_sourcedata" directory contains the 8-mers and E-scores underlying the bar charts and boxplots for the individual TF-rsID examples depicted in Figure 6 and Supplementary Figure 6 in the manuscript.

The following files are required for the analysis but NOT included in this directory due to space constraints; these are readily available from journal webpages / public databases and resources:

1) GRCh38/hg38 data files by chromosomes - obtain the "chr<chr number>.fa.gz" files from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

2) dbSNP (Build 151, GRCh38p7) common SNPs - obtain "00-common_all.vcf.gz" from ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/

3) Supplementary Table 1 (41588_2019_362_MOESM3_ESM.xlsx, which contains fine-mapped variants with PP > 0.001 from Ulirsch et al., Nat Genet, 2019) from https://www.nature.com/articles/s41588-019-0362-6#Sec27 (input file for certain code blocks in "NatGenet2020_Fig6_FigS6_PythonCode.ipynb" for examining perturbed TF binding in rsIDs present in both this study and Ulirsch et al., Nat Genet, 2019)

4) Mus musculus PBM E-scores from the CIS-BP database at http://cisbp.ccbr.utoronto.ca/bulk.php

### The workflow of this analysis comprises the following: ###

1) Download all necessary files into the directory containing the code files, and into your compute cluster.

2) Use "NatGenet2020_Fig6_FigS6_RCode.R" to generate standardised PBM E-score files from processed PBM data files from UniPROBE and CIS-BP (present in the "PBM_datasets" directory). The output from this is present in the "Individual_PBMs_Escores" directory.

3) Run code blocks in "NatGenet2020_Fig6_FigS6_PythonCode.ipynb" to obtain individual examples of rsIDs that perturb TF binding (inferred from published PBM datasets), along with general trends observed in the foreground of RBC SNPs with regards to perturbation of binding by TFs of interest. The output from these instructions is present as "foreground_3263SNPs.txt" and individual examples in the "rsID_Escores_sourcedata" directory. There are also code blocks within this notebook for supplementary analyses, of the FHS and fine-mapped SNPs datasets, relevant to perturbed TF binding events considered in this study.

4) Obtain the background distributions for perturbation of binding by TFs of interest in the set of dbSNP common SNPs through your compute cluster. Use "vcftools_commonSNPs.sh" to recode the "00-common_all.vcf" file, and then use "NatGenet2020_RBCSNPs_individualPBM_binary_10percentSNPs.py" to obtain the background statistics. The output from these directions is present as the individual files (one per PBM dataset) in the "empirical_background" directory.

5) Use "NatGenet2020_Fig6_FigS6_RCode.R" to perform the bootstrapping computations and the statistical tests, and to generate the figures.
