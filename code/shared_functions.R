
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


## PRINT NICELY-FORMATTED TABLES
print_kbl <- function(tbl) {
  kbl(tbl) %>%
  kable_styling(full_width = FALSE, 
                position = "left",
                bootstrap_options = c("striped", "hover", "responsive")) %>%
  scroll_box(height = "300px", width = "100%")
}

## MAKE FIGURES

## Linear regression plots
make_lm_plot <- function(df, plot_list, x_val){
  d.plot <- df %>%
    unite("target_compound_label", c("target", "compound_name"), sep = "\n", remove = FALSE) %>%
    filter(target_compound_label %in% plot_list) %>%
    mutate("target_compound_label" = reorder(target_compound_label, fdr))
  
  d.label <- d.plot %>%
    group_by(target, compound_id) %>%
    mutate(group_label_x = min({{x_val}}),
           group_label_y = max(dependency)) %>%
    dplyr::ungroup() %>%
    dplyr::select(target_compound_label, r_squared, p_val, fdr, group_label_x, group_label_y) %>%
    distinct(target_compound_label, .keep_all = TRUE) %>%
    mutate(r_squared = round(r_squared, 3),
           p_val = scientific(p_val, 3),
           fdr = scientific(fdr, 3)) 
  # d.label
  
  plot <- ggplot(d.plot, aes(x = {{x_val}}, y = dependency)) +
    geom_hline(yintercept = 0) +
    rasterize(geom_point(size = 1), dpi = 150) +
    geom_smooth(method = lm, se = FALSE, color = "gray60") +
    geom_richtext(data = d.label,
                  aes(x = group_label_x, y = group_label_y,
                      ## this package can do HTML-formatting for text
                      label = paste0("r<sup>2</sup>=", r_squared, "<br>", 
                                     "FDR =", fdr)),
                  ## adjust label alignment and remove borders
                  hjust = 0, vjust = 1, label.color = NA, size = 3,
                  ## label background = transparent, text is not
                  fill = alpha(c("white"), 0.75),
                  label.margin = unit(c(0, 0, 0, 0), "lines"),
                  label.padding = unit(c(0, 0, 0, 0), "lines")) +
    labs(x = "log2(TPM)", y = "Cell_viability") +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black")) +
    facet_wrap(~target_compound_label, scales = "free", ncol = 3)
  
  return(plot)
}

## Z-scaled outliers
make_outlier_plot <- function(df, plot_list, direction, outlier_var, expr_var){
  
  d.plot <- df %>%
    unite("target_compound_label", c("target", "compound_name"), sep = "\n", remove = FALSE) %>%
    filter(target_compound_label %in% plot_list) %>%
    mutate("target_compound_label" = reorder(target_compound_label, fdr))
  
  d.label <- d.plot %>%
    group_by(target, compound_id) %>%
    mutate(y.position = (max(dependency) + 0.1*(max(dependency) - min(dependency)))) %>%
    dplyr::ungroup() %>%
    dplyr::select(target_compound_label, fdr, y.position) %>%
    distinct(target_compound_label, .keep_all = TRUE) %>%
    mutate(group1 = direction,
           group2 = "neither",
           fdr = scientific(fdr, 3)) 
  
  p <-  ggplot(d.plot, aes(x = {{outlier_var}}, y = dependency)) +
    geom_hline(yintercept = 0) +
    rasterize(geom_jitter(aes(color = {{outlier_var}}), width = 0.3, alpha = 0.6), dpi = 150) +
    geom_boxplot(aes(fill = {{outlier_var}}), outlier.shape = NA, coef = 0) +
    stat_pvalue_manual(d.label, label = paste0("FDR = {fdr}"), tip.length = 0.03, size = 3, vjust = -0.5) +
    ## adds extra space on top & bottom so p-value label will fit
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    scale_color_manual(values = c("low" = "#0080C6", "neither" = "#7F7F7F", "high" = "#B84250")) +
    scale_fill_manual(values = c("low" = "#0080C6", "neither" = "#7F7F7F", "high" = "#B84250")) +
    labs(x = "Cell_viability",
         y = "Expression_group") +
    theme_bw() +
    theme(aspect.ratio = 0.75,
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          legend.position = "none") +
    facet_wrap(~target_compound_label, scales = "free_y",
               ncol = 1)
  
  q <- ggplot(d.plot, aes(x = target, y = {{expr_var}})) +
    rasterize(geom_jitter(aes(color = {{outlier_var}}), width = 0.3), dpi = 150) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, coef = 0) +
    scale_color_manual(values = c("low" = "#0080C6", "neither" = "#7F7F7F", "high" = "#B84250")) +
    labs(x = "Target", y = "log2(TPM)") +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black")) +
    facet_wrap(~target_compound_label, scale = "free", ncol = 1)
  
  ggarrange(p, q, ncol = 2, align = "hv", 
            common.legend = TRUE, legend = "right")
  
}

## Expression across lineages
make_lineage_plot <- function(df, plot_list, y_var, color_flag){
  
  d.plot <- df %>% 
    ungroup() %>%
    unite("target_compound_label", c("target", "compound_name"), sep = "\n", remove = FALSE) %>%
    filter(target_compound_label %in% plot_list) %>%
    mutate("target_compound_label" = reorder(target_compound_label, fdr),
           "target" = reorder(target, fdr)) 
  
  ggplot(d.plot, aes(x = lineage, y = {{y_var}})) +
    rasterize(geom_jitter(aes(color = {{color_flag}}), width = 0.3, size = 1), dpi = 150) +
    geom_boxplot(alpha = 0, outlier.shape = NA, coef = 0) +
    scale_color_manual(values = c("low" = "#0080C6", "neither" = "#7F7F7F", "high" = "#B84250")) +
    labs(x = "Cell_lineage", y = "log2(TPM)") +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 45, hjust=1)) +
    facet_wrap(~target, scales = "free_y",
               ncol = 1)
}

