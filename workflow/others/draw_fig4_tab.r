###
#' @Date: 2022-07-20 13:43:25
#' @LastEditors: Hwrn
#' @LastEditTime: 2022-08-05 14:00:50
#' @FilePath: /2021_09-MT10kSW/workflow/others/draw_fig4_tab.r
#' @Description:
###
source("workflow/utils/RLib.local/R/init.r", chdir = TRUE)


### ######################################################################## ###
#### Preprocessing                                                          ####
### ######################################################################## ###
##### INPUT: file_path, keyword_args, fig_out_path                         #####
tab_out = argv[1]

##### GLOBAL CONST vars                                                    #####
fun.aggregate = "sum"

##### LOAD data AND transform TO basic format                              #####
gene_ko_tpm = load__gene_ko_tpm()
TOTAL_GENE_NUM = gene_ko_tpm$KO %>% unique %>% length %>% {. * 2}

key_genes = load__key_genes()
key_genes_N =
  key_genes %>% .[.$pathway == "N",] %>%
  {.$pathway = .$arrow %>% factor(levels = unique(.)); .}

### ######################################################################## ###
#### Define function AND Calculate data                                     ####
### ######################################################################## ###
##### Define function                                                      #####
report__tpm.key_tab_site <- function(tpm_key, fun.aggregate) {
  tpm_key %>%
    {.$sample = str_glue("{.$group}_{.$Layer}"); .} %>%
    reshape2::acast("sample ~ ko", value.var = "tpm",
                    fun.aggregate = get(fun.aggregate)) %>%
    {
      sorted.ko = key_genes_N$KO %>% levels
      .[, sorted.ko[sorted.ko %in% colnames(.)]]
    } %>%
    {.[is.na(.)] = 0; .}
}

##### Calculate data                                                       #####
tpm.key =
  gene_ko_tpm %>%
  merge(key_genes_N,
        by.x = "ko", by.y = "KO") %>%
  {.$Site = .$Layer %>% gsub("^(.+)\\.\\.(.+)$", "\\1", .); .} %>%
  {.$pathway = factor(.$pathway, levels = levels(key_genes_N$pathway)); .} %>%
  merge(site_group) %>%
  merge(sample_meta[c("location", "group")] %>% unique) %>%
  {.$location = factor(.$location, levels = c("Slope", "Bottom")); .}


### ######################################################################## ###
#### Calulate table and OUTPUT                                              ####
### ######################################################################## ###
d =
  tpm.key %>%
  report__tpm.key_tab_site(fun.aggregate)


write.csv(d, tab_out)
