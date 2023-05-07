# Paralog_SL_meta_analysis

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
