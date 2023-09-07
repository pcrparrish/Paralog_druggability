# Paralog_druggability

## Background
Paralog synthetic lethal therapies are an attractive option for cancer treatment. Synthetic lethal therapies have the potential to target cancer cells more specifically with fewer off-target effects than oncogene-targeted therapies, and paralogs (AKA duplicated genes) show relatively high rates of synthetic lethality in human cancer cells ([ref1](https://www.sciencedirect.com/science/article/pii/S2211124721010354), [ref2](https://link.springer.com/article/10.1186/s13059-020-02173-2), [ref3](https://www.nature.com/articles/s41587-020-0437-z), [ref4](https://www.nature.com/articles/s41467-021-21478-9)). For more information on paralog synthetic lethal therapies, see [this recent review](https://www.cell.com/trends/cancer/fulltext/S2405-8033(23)00022-5). 

However, one potential drawback to paralog-targeted therapies is that many small molecule drugs target one or more members of the same gene family. I thus wanted to determine if viability of cancer cells with loss of one paralog gene is reduced when treated with a drug targeting one or more other members of that gene family, as illustrated below. 
![Conceptual schematic of paralog druggability analysis](https://github.com/pcrparrish/Paralog_druggability/blob/main/resources/druggability_analysis_background.png?raw=true)

If loss of only one paralog family-member sensitizes cancer cells to a drug that targets the whole paralog family, that drug could still be used to kill cancer cells while avoiding toxicity to healthy cells. 

I therefore analyzed an [existing dataset](https://www.nature.com/articles/s43018-019-0018-6) generated by treating 578 cancer cell lines with over 4,500 drugs. I focused on 367 gene families that were predicted to be synthetic lethal via [a recent computational approach](https://www.sciencedirect.com/science/article/pii/S240547122100329X). 

## Approach and Results
To identify druggable paralog dependencies, I used a two-tiered approach: first, I applied a discovery-focused linear regression analysis to identify cases where low expression of one paralog was significantly associated with lower viability when cancer cells were treated with an inhibitor targeting that paralog family. This first analysis identified 68 cases where loss of expression of one family-member was significantly associated with increased sensitivity to drugs targeting that gene family (FDR < 0.05 via Benjamini-Hochberg correction). Results for the top six hits from this first analysis are shown below in panels A-F:
![Top hits from linear regression analysis](https://github.com/pcrparrish/Paralog_druggability/blob/main/resources/lm_analysis_top_hits.png?raw=true)
**(A)**	Left panel, scatterplot of Gene A2 expression versus cell viability in drug treatment for *AKT1/AKT3* and MK-2206. Right panel, Gene A1 and Gene A2 expression across PRISM cell lines for *AKT1/AKT3*. Cell viability is measured as the log2(fold-change) in drug versus control (DMSO) treatment. Gene expression is measured as log2(TPM). For each pair, Gene A1 refers to the first gene listed (in this case, *AKT1*) and Gene A2 refers to the second member of the pair (here, *AKT3*). 
**(B)**	As in (A), but for *ABL1/ABL2* and saracatinib. 
**(C)**	As in (A), but for *PTPN11/PTPN6* and BVT-948. 
**(D)**	As in (A), but for *CDK1/CDK2* and NU6027. 
**(E)**	Left panel, scatterplot of Gene A1 expression versus cell viability in drug treatment for FYN/SRC and vandetanib. Right panel, Gene A1 and Gene A2 expression across PRISM cell lines for *FYN/SRC*. Cell viability and gene expression are measured as in (A).
**(F)**	As in (E), but for *ATP1A1/ATP1A3* and k-strophanthidin.


I then used a more stringent outlier-based analysis to confirm my findings for 13 promising paralog pair/drug combinations. This approach validated 5 of the hits from my first analysis. The top hit from our analysis was the use of AKT inhibitors in cancer cell lines showing loss of *AKT3* expression. Results for the *AKT1/AKT2/AKT3* gene family and multiple AKT inhibitors are shown below: 
![Outlier analysis hits](https://github.com/pcrparrish/Paralog_druggability/blob/main/resources/AKT_family_AKTi.png?raw=true)
**(A)**	Grouped scatterplot of *AKT3* (Gene A2) expression versus cell viability in GSK2110183 relative to control (DMSO) treatment for the AKT1/AKT3 paralog pair. Overlaid box-and-whisker plots indicate group median, IQR, and outliers. WRST with Benjamini-Hochberg FDR correction was used to compare the median of low versus normal groups. 
**(B)**	As in (A), but for AZD5363. 
**(C)**	Grouped scatterplot of *AKT1* and *AKT3* expression across PRISM cell lines. Gene expression is measured as log2(TPM). Colors indicate the expression outlier group of each cell line: blue is low, gray is normal, red is high. 
**(D)**	As in (A), but for *AKT2/AKT3*. 
**(E)**	As in (B), but for *AKT2/AKT3*. 
**(F)**	As in (C), but for *AKT2* and *AKT3*.


Importantly, I found that *AKT3* showed a bimodal distribution of expression across multiple cancer cell lineages, as seen below: 
![AKT3 expression across DepMap cell lineages](https://github.com/pcrparrish/Paralog_druggability/blob/main/resources/drug_appendix_fig5.png?raw=true)

Cancer biologists have suggested that bimodal expression may be a biomarker for cancer vulnerability ([ref1](https://aacrjournals.org/cancerres/article/82/13/2378/705034/Bimodal-Gene-Expression-in-Patients-with-Cancer), [ref2](https://www.nature.com/articles/s41573-019-0046-z)).

## Conclusions
These results suggest that drugs targeting multiple paralogs in a synthetic lethal family can be used to selectively harm tumor cells while having limited effects on healthy cell viability. Overall, the results from this computational analysis support the use of paralog synthetic lethal therapies in cancer. 

## To Do:
* Actually remove Finan filtering step from expr and CN DFs
* Build expr and CN DFs in the make_drug_dfs.Rmd
* Add a file.exists() command to make_drug_dfs.Rmd => load DepMap data, etc., if files do not exist; else, read in RDS
* CN as continuous variable => look for outliers (same for expression?)
* look for bimodal expr distribution
* figure out if I am calculating FDR correctly (I am re-testing some things)
* plot all DepMap CN and expression values...or something
  * redo expr stat test for both outliers => 2-sided (right now they are 1-sided which is kind of cheating)
  * plot expression outliers by lineage for top pairs
  * summarize n targets with FDR < 0.1 for each analysis (but make sure I am comparing things correctly)
  * two-sided test => look for resistance too, make volcano plot
  * paralog plots: color dots based on CN/expr balance b/w paralogs
* Add a file.exists() command to make_drug_dfs.Rmd => load DepMap data, etc., if files do not exist; else, read in RDS
* get rid of unused code in make_drug_dfs.Rmd and the CN/expression scripts
* look at data from other drug datasets (Sanger, CellMiner, other Broad sets?)
* look at TCGA paralog mutation data
