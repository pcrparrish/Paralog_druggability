
## STATS

## adjust p-values
BH_adjust <- function(df){
  
  d.p_vals <- df %>%
    arrange(p_val) %>%
    dplyr::select(target, compound_id, p_val) 
  
  p_vals <- d.p_vals %>%
    arrange(p_val) %>%
    dplyr::select(p_val) %>%
    pull()
  
  fdr_vals <- p.adjust(p_vals, method = "BH")
  
  d.fdr_vals <- tibble(fdr = fdr_vals)
  
  d.stats <- d.p_vals %>%
    add_column(d.fdr_vals)
  d.stats
  
  return(d.stats)
  
}

## summarize stats output
summarize_stats <- function(df){
  df_new <- df %>%
    distinct(target, compound_id, .keep_all = TRUE)
  
  n_p_L0.05 <- df_new %>%
    filter(p_val < 0.05) %>%
    nrow()
  
  n_p_L0.01 <- df_new %>%
    filter(p_val < 0.01) %>%
    nrow()
  
  n_fdr_L0.05 <- df_new %>%
    filter(fdr < 0.05) %>%
    nrow()
  
  n_fdr_L0.1 <- df_new %>%
    filter(fdr < 0.1) %>%
    nrow()
  
  cat(paste("df:", deparse(substitute(df)), "\n", sep = " "))
  cat("# target/compound pairs with\n")
  cat(paste("p < 0.05:", n_p_L0.05, "\n", sep = " "))
  cat(paste("p < 0.01:", n_p_L0.01, "\n", sep = " "))
  cat(paste("FDR < 0.1:", n_fdr_L0.1, "\n", sep = " "))
  cat(paste("FDR < 0.05:", n_fdr_L0.05, "\n\n", sep = " "))
  
}

print_kbl <- function(tbl) {
  kbl(tbl) %>%
  kable_styling(full_width = FALSE, 
                position = "left",
                bootstrap_options = c("striped", "hover", "responsive")) %>%
  scroll_box(width = "100%")
}