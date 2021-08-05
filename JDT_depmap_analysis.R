################################################################################
## Paralog analysis in DepMap datasets
################################################################################

## Goals:
## - Use the R functions provided by the bioconductor "depmap" package to test
##     if gene1 (from a paralog pair) is more dependent in cell lines with low
##     expression of gene2. This analysis is actually mostly for the paralog
##     screen paper, but I will also explore the MBNL paralog family here.
##     Therefore, to avoid code redundancy, I'm including these analyses here.

## Conclusions:
## - Many paralogs show evidence of synthetic lethality in this analysis (i.e.,
##     the effect of gene1 knockout is strongest in cells with lowest expression
##     of gene2).
## - Very few paralogs (maybe ~3-10 out of >900 included in the analysis below)
##     show evidence of "reciprocal dependency" -- gene1 depenendence in cells
##     with low gene2 expression *and* gene2 dependence in cells with low gene2
##     expression -- in the CCLE data. Need to think more about what the
##     possible explanations could be. One interesting possibility is that this
##     lack of reciprocity could highlight why doing paralog knockout screens
##     is important -- you might miss some synthetic lethal pairs if you were
##     doing purely computational analysis of published screen data.

## Notes:
## - 

## Load packages required for this analysis
library("depmap")
library("ExperimentHub")

## Parameters
n.bins = 4

## For simplicity/speed, I'm only going to analyze gene (pairs) of interest.
##   First, I'll define all of the paralog pairs included in the Phoebe's paralog
##   screen.
file = file.path(dir.gene_expression, "R", "projects", "2019/paralog_screen",
                 "jdthomas", "data", "PC9_Day6_v_Day21_targetLevel_meanLFC_residual_rank.tsv")
d.paralogs = suppressMessages(read_delim(file, delim = "\t")) %>% select(Gene) %>%
    separate(Gene, c("geneName.1", "geneName.2"))

## Sort lexiographically
d.paralogs = d.paralogs %>%
    mutate("foo" = paste(geneName.1, geneName.2)) %>%
    mutate("foo2" = sapply(strsplit(foo, split = " "), function(i) paste(sort(i), collapse = " "))) %>%
    select(foo2) %>% separate("foo2", c("geneName.1", "geneName.2"))

## Add MBNL1 and MBNL2
d.paralogs = rbind(d.paralogs, tibble("geneName.1" = "MBNL1", "geneName.2" = "MBNL2"))

## Finally, make a vector of paralog geneNames
geneNames = d.paralogs %>% gather() %>% pull(value)

#################################!
## Build data structure: d.depmap
#################################!

## Create ExperimentHub query object
eh = ExperimentHub()
eh.depmap = query(eh, "depmap")

## Look up available data
foo = tibble("eh_id" = names(eh.depmap),
             "title" = eh.depmap$title)
foo %>% print(n = Inf)

## Retrieve data from "2020 quarter 1" and restrict to geneNames of interest
foo %>% filter(grepl("20", title))    
## - CRISPR data
d.crispr = eh.depmap[["EH3290"]] %>% select(-cell_line, -gene, -entrez_id) %>%
    filter(gene_name %in% geneNames)
## - RNA expression data
d.tpm = eh.depmap[["EH3292"]] %>% select(-cell_line, -gene, -ensembl_id) %>%
    filter(gene_name %in% geneNames)
## - Metadata
d.metadata = eh.depmap[["EH3294"]]

## Join d.tpm and d.crispr
d.depmap = d.tpm %>% left_join(d.crispr, by = c("depmap_id", "gene_name"))    

## Dependency data isn't available for every cell line
d.depmap = d.depmap %>% filter(!is.na(dependency))

## Rename columns
d.depmap = d.depmap %>%
    rename("id" = depmap_id,
           "tpm" = expression,
           "geneName" = gene_name,
           "CS" = dependency) %>%
    select(id, geneName, tpm, CS)

## Parse data so that paralogs are "next to each other"
foo = d.paralogs %>% mutate("paralog_id" = paste0(geneName.1, "|", geneName.2))

tmp1 = d.paralogs %>% select(geneName.1) %>% pull()
tmp2 = d.paralogs %>% select(geneName.2) %>% pull()

tmp1 = d.depmap %>% filter(geneName %in% tmp1) %>%
    rename("geneName.1" = geneName, "tpm.1" = tpm, "CS.1" = CS)
tmp2 = d.depmap %>% filter(geneName %in% tmp2) %>%
    rename("geneName.2" = geneName, "tpm.2" = tpm, "CS.2" = CS)

tmp1 = tmp1 %>% left_join(foo %>% select(geneName.1, paralog_id), by = c("geneName.1")) %>%
    select(id, paralog_id, everything()) %>% arrange(paralog_id)                           
tmp2 = tmp2 %>% left_join(foo %>% select(geneName.2, paralog_id), by = c("geneName.2")) %>%
    select(id, paralog_id, everything()) %>% arrange(paralog_id)                           

d.depmap = tmp1 %>% full_join(tmp2, by = c("id", "paralog_id"))

## Only consider paralog pairs that have complete information (i.e., drop
##   rows that contain missing data) for a given cell line.
d.depmap = d.depmap %>% drop_na()

########################################################!
## Box plots: gene2 dependency based on gene1 expression
########################################################!

## Parse d.depmap
d.foo = d.depmap %>% select(id, paralog_id, tpm.1, CS.2)

## For each paralog_pair, "bin" based on gene
d.foo = d.foo %>% group_by(paralog_id) %>%
    mutate("bin.tpm.1" = ntile(tpm.1, n.bins)) %>%
    ungroup()

## Restrict to the "top 100" paralog_ids from the screen
paralog_ids = d.paralogs %>%
    mutate("paralog_id" = paste0(geneName.1, "|", geneName.2)) %>%
    head(n = 100) %>% pull()

## - Add in MBNL1 and MBNL2
paralog_ids = c("MBNL1|MBNL2", paralog_ids)
d.foo = d.foo %>% filter(paralog_id %in% paralog_ids)

## To ensure panels are in the desired order, convert paralog id to a factor
##   and specify the order
d.foo = d.foo %>% mutate("paralog_id" = factor(paralog_id, levels = paralog_ids))

## Convert bin.tpm.1 to factor for plotting (and rearrange order)
d.foo$bin.tpm.1 = factor(d.foo$bin.tpm.1, levels = c("4", "3", "2", "1"))

## Base plot
p = ggplot(d.foo, aes(x = bin.tpm.1, y = CS.2, fill = bin.tpm.1))
## Box plot
p = p + geom_boxplot(outlier.color = "black",
                     outlier.size = 1,
                     outlier.stroke = 0,
                     notch = TRUE)
## Rotate x-axis labels
p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
## Faceting
## - specify number of columns
n.cols = 10
## - specify number of rows
n.rows = round(length(paralog_ids)/n.cols)
p = p + facet_wrap(~ paralog_id,
                   ncol = n.cols,
                   nrow = n.rows,
                   scales = "free")
p = p + theme(strip.background = element_blank(),
              strip.text = element_text(size = 11))
## Remove legend
p = p + theme(legend.position = "none")
## Change colors
p = p + scale_fill_brewer(palette="YlGnBu", direction = -1)
## Define figure dimensions
x.dim = n.cols * 3
y.dim = n.rows * 3
## File name
file = file.path(dir.figz,
                 "box.depmap_gene2_dependency_based_on_gene1_expression.pdf")
## Save figure
save_plot(file,
          plot = p,
          device = "pdf",
          base_width = x.dim,
          base_height = y.dim,
          base_asp = 1,
          limitsize = FALSE,
          useDingbats = FALSE)
message("Saved file: ", file)

##################################!
## Summary scatter plot from above
##################################!

## Notes:
## - x-axis == logFC bin1 -vs- bin4; y-axis == p-value from wilcoxon test

## Parameters:
to_label = c("CDS1|CDS2", "DDX39A|DDX39B",
             "STAG1|STAG2", "IMPDH1|IMPDH2",
             "RBBP4|RBBP7")
max.pval = 0.05

## p-values
tmp1 = d.foo %>% filter(bin.tpm.1 %in% c(1, n.bins)) %>%
    select(id, paralog_id, bin.tpm.1, CS.2) %>%
    group_by(paralog_id) %>%
    summarize("pval" = wilcox.test(CS.2 ~ bin.tpm.1)$p.value) %>%
    ungroup()

## fold-changes
tmp2 = d.foo %>% filter(bin.tpm.1 %in% c(1, n.bins)) %>%
    select(id, paralog_id, bin.tpm.1, CS.2) %>%
    group_by(paralog_id, bin.tpm.1) %>%
    summarize("CS.2.median" = median(CS.2)) %>%
    ungroup() %>%
    spread(bin.tpm.1, CS.2.median) %>%
    mutate("CS.2.diff" = `1` - `4`) %>%
    select(paralog_id, CS.2.diff)

## join
d.foo = left_join(tmp1, tmp2, by = c("paralog_id"))

## convert pval to neg_log10 and fc to logfc
d.foo = d.foo %>% mutate("neglog10pval" = log10(pval) * -1)

## add some annotation columns for labelling points and coloring by significance
d.foo = d.foo %>%
    mutate("to_label" = ifelse(paralog_id %in% to_label, TRUE, FALSE)) %>%
    mutate("significant" = ifelse(pval < max.pval, TRUE, FALSE))

## Base plot
p = ggplot(d.foo, aes(x = CS.2.diff, y = neglog10pval))
## Scatter plot
p = p + geom_point(stroke = 0, aes(color = significant))
## Annotations
p = p + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"))
p = p + geom_text_repel(aes(label = ifelse(to_label == TRUE, as.character(paralog_id),'')),
                        hjust=0, vjust=0, color = "black", size = 3)
## Remove legend
p = p + theme(legend.position = "none")
## Define figure dimensions
x.dim = 4
y.dim = 4
## File name
file = file.path(dir.figz,
                 "scatter.depmap_summary.pdf")
## Save figure
save_plot(file,
          plot = p,
          device = "pdf",
          base_width = x.dim,
          base_height = y.dim,
          base_asp = 1,
          limitsize = FALSE,
          useDingbats = FALSE)
message("Saved file: ", file)

d.foo %>% filter(to_label == TRUE)

###################################!
## Test for "reciprocal dependency"
###################################!

to_label = c(
    ## positive diff.gene2
    "GRHL1|GRHL2", "CREBBP|EP300", "SEPHS1|SEPHS2",
    ## top negative diff.gene2
    "BTRC|FBXW11", "CDS1|CDS2", "ARHGEF6|ARHGEF7",
    ## top negative diff.gene1
    "INTS6|INTS6L", "FAM50A|FAM50B", "SNAP23|SNAP25",
    ## some paralog pairs from previous analysis
    "CDS1|CDS2", "DDX39A|DDX39B", "STAG1|STAG2"
)

max.pval = 0.05

## Parse data: gene1 dependency
## - parsing
d.foo1 = d.depmap %>% select(id, paralog_id, tpm.2, CS.1)
d.foo1 = d.foo1 %>% group_by(paralog_id) %>%
    mutate("bin.tpm.2" = ntile(tpm.2, n.bins)) %>%
    ungroup()
## - compute p-values
tmp1 = d.foo1 %>% filter(bin.tpm.2 %in% c(1, n.bins)) %>%
    select(id, paralog_id, bin.tpm.2, CS.1) %>%
    group_by(paralog_id) %>%
    summarize("pval.gene1" = wilcox.test(CS.1 ~ bin.tpm.2)$p.value) %>%
    ungroup()
## - compute fold-changes
tmp2 = d.foo1 %>% filter(bin.tpm.2 %in% c(1, n.bins)) %>%
    select(id, paralog_id, bin.tpm.2, CS.1) %>%
    group_by(paralog_id, bin.tpm.2) %>%
    summarize("CS.1.median" = median(CS.1)) %>%
    ungroup() %>%
    spread(bin.tpm.2, CS.1.median) %>%
    mutate("diff.gene1" = `1` - `4`) %>%
    select(paralog_id, diff.gene1)
## - join
d.foo1 = left_join(tmp1, tmp2, by = c("paralog_id"))

## Parse data: gene2 dependency
## - parsing
d.foo2 = d.depmap %>% select(id, paralog_id, tpm.1, CS.2)
d.foo2 = d.foo2 %>% group_by(paralog_id) %>%
    mutate("bin.tpm.1" = ntile(tpm.1, n.bins)) %>%
    ungroup()
## - compute p-values
tmp1 = d.foo2 %>% filter(bin.tpm.1 %in% c(1, n.bins)) %>%
    select(id, paralog_id, bin.tpm.1, CS.2) %>%
    group_by(paralog_id) %>%
    summarize("pval.gene2" = wilcox.test(CS.2 ~ bin.tpm.1)$p.value) %>%
    ungroup()
## - compute fold-changes
tmp2 = d.foo2 %>% filter(bin.tpm.1 %in% c(1, n.bins)) %>%
    select(id, paralog_id, bin.tpm.1, CS.2) %>%
    group_by(paralog_id, bin.tpm.1) %>%
    summarize("CS.2.median" = median(CS.2)) %>%
    ungroup() %>%
    spread(bin.tpm.1, CS.2.median) %>%
    mutate("diff.gene2" = `1` - `4`) %>%
    select(paralog_id, diff.gene2)
## - join
d.foo2 = left_join(tmp1, tmp2, by = c("paralog_id"))

## Join
d.foo = d.foo1 %>% left_join(d.foo2, by = c("paralog_id"))

## Add column for labelling
d.foo = d.foo %>% mutate("to_label" = ifelse(paralog_id %in% to_label, TRUE, FALSE))    

## Base plot
p = ggplot(d.foo, aes(x = diff.gene1, y = diff.gene2))
## Scatter plot
p = p + geom_point(stroke = 0)
## Annotations
##p = p + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"))
p = p + geom_text_repel(aes(label = ifelse(to_label == TRUE, as.character(paralog_id),'')),
                        hjust=0, vjust=0, color = "black", size = 3)
## Remove legend
p = p + theme(legend.position = "none")
## Adjust labels
p = p + ylab("difference in gene 2 CRISPR score\n(gene 1 low-vs-high expression groups)")
p = p + xlab("difference in gene 1 CRISPR score\n(gene 2 low-vs-high expression groups)")
## Define figure dimensions
x.dim = 5
y.dim = 5
## File name
file = file.path(dir.figz,
                 "scatter.depmap_reciprocal_dependency.pdf")
## Save figure
save_plot(file,
          plot = p,
          device = "pdf",
          base_width = x.dim,
          base_height = y.dim,
          base_asp = 1,
          limitsize = FALSE,
          useDingbats = FALSE)
message("Saved file: ", file)
