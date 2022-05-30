library(tidyverse)
library(ggpmisc)
# install.packages("devtools") 
# devtools::install_github("johannesbjork/LaCroixColoR")
library(LaCroixColoR)
library(gridExtra)
library(wesanderson)
# devtools::install_github("ricardo-bion/ggradar")
library(ggradar)
library(knitr)

setwd('/Users/dylansosa/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Yakuba  Long Lab/thesis/3_tRNA_evolution/')


pacbioFiles <- list.files(path = 'pacbio/tRNAscan-SE_output/tRNA_gene_count_information/',pattern = '*.tsv')
pacbioDir <- 'pacbio/tRNAscan-SE_output/tRNA_gene_count_information'
nonpacbioFiles <- list.files(path = 'non-pacbio/tRNAscan-SE_output/tRNA_gene_count_information/',pattern = '*.tsv')
nonpacbioDir <- 'non-pacbio/tRNAscan-SE_output/tRNA_gene_count_information'
# get path and file name

pacbioOutput <- file.path(pacbioDir,pacbioFiles)
nonpacbioOutput <- file.path(nonpacbioDir,nonpacbioFiles)
# combine to make path with filename

for(file in pacbioOutput){
  prefix <- str_remove(str_extract(str_extract(file,'[A-Z]{3}_[0-9].*'),'_[A-Z]?[a-z]+'),'_')
  assign(prefix,read_tsv(file))
}

for(file in nonpacbioOutput){
  prefix <- str_remove(str_extract(str_extract(file,'[A-Z]{3}_[0-9].*'),'_[A-Z]?[a-z]+'),'_')
  assign(prefix,read_tsv(file))
}
# loop through paths and make unique variables based on assembly name


# testisotpyessleb <- Sleb %>%
#   dplyr::select(`Isotypes Used in this Species`,`Number of tRNA Genes per Isotype`) %>%
#   na.omit() %>%
#   pivot_wider(names_from = `Isotypes Used in this Species`, values_from = `Number of tRNA Genes per Isotype`) %>%
#   mutate(species = 'Sleb') %>% dplyr::select(species,everything())
# # process data inta a usable format s
# 
# testisotpyessleb
# eventually I need to combine all species into a dataset

# ggRadar(testisotpyessleb, aes(group=species),
#         rescale = FALSE, legend.position = "none",
#         size = 3, interactive = FALSE, use.label = TRUE) +
#   facet_wrap(~species) +
#   # scale_y_discrete(breaks = T) + # don't show ticks 
#   theme(axis.text.x = element_text(size = 10)) + # larger label sizes
#   # adjust colors of radar charts to uniform colors
#   scale_fill_manual(values = rep(mycolor, nrow(testisotpyessleb))) +
#   scale_color_manual(values = rep(mycolor, nrow(testisotpyessleb))) +
#   ggtitle("Number of tRNA genes per isotype")
# single file case, need to generalize for all
# trying to get all files in correct format 

# speciesList <- ls(pattern = '^[A-Z][a-z]{3}')

widenIsotypes <- function(speciesFile){
  x <- deparse(substitute(speciesFile))
  speciesFile %>%
    dplyr::select(`Isotypes Used in this Species`,`Number of tRNA Genes per Isotype`) %>%
    na.omit() %>%
    pivot_wider(names_from = `Isotypes Used in this Species`, values_from = `Number of tRNA Genes per Isotype`) %>%
    mutate(species = x) %>% dplyr::select(species,everything())
}
# function to do the formatting I want for plotting
 
# for(spp in speciesList){
#   te <-widenIsotypes(get(spp))
# }

Dana.wide <- widenIsotypes(Dana)
Dazt.wide <- widenIsotypes(Dazt)
Dere.wide <- widenIsotypes(Dere)
Dhyd.wide <- widenIsotypes(Dhyd)
Dmel.wide <- widenIsotypes(Dmel)
Dmoj.wide <- widenIsotypes(Dmoj)
Dnov.wide <- widenIsotypes(Dnov)
Dore.wide <- widenIsotypes(Dore)
Dana.wide <- widenIsotypes(Dana)
Dper.wide <- widenIsotypes(Dper)
Dpse.wide <- widenIsotypes(Dpse)
Dsim.wide <- widenIsotypes(Dsim)
Dvir.wide <- widenIsotypes(Dvir)
Dwil.wide <- widenIsotypes(Dwil)
Dyak.wide <- widenIsotypes(Dyak)
Sleb.wide <- widenIsotypes(Sleb)

# cat te | cut -f 1 -d ' ' | sed s'/^/full_join(/' | sed s'/$/) %>%/'
# where te is the rows above
totalWideIstoypeCounts <- 
  Dana.wide %>% 
  full_join(Dazt.wide) %>%
  full_join(Dere.wide) %>%
  full_join(Dhyd.wide) %>%
  full_join(Dmel.wide) %>%
  full_join(Dmoj.wide) %>%
  full_join(Dnov.wide) %>%
  full_join(Dore.wide) %>%
  full_join(Dana.wide) %>%
  full_join(Dper.wide) %>%
  full_join(Dpse.wide) %>%
  full_join(Dsim.wide) %>%
  full_join(Dvir.wide) %>%
  full_join(Dwil.wide) %>%
  full_join(Dyak.wide) %>%
  full_join(Sleb.wide) %>% 
  replace(is.na(.), 0)
# once i have all species data in wide format I can join all to add columns for species that have 0 counts for any isotype
# I hope this would force and NA for that isotype for a species missing it; then I could convert to 0

totalIsotypeCounts <- 
  ggRadar(totalWideIstoypeCounts, aes(group=species),
        rescale = FALSE, legend.position = "none",
        size = 2, interactive = FALSE, use.label = TRUE) +
  theme(axis.text.x = element_text(size = 10)) + 
  # scale_fill_manual(values = rep(mycolor, nrow(totalWideIstoypeCounts))) +
  scale_color_manual(values = lacroix_palette("Apricot", n = 15, type = "continuous")) +
  # LaCroix! :)
  ggtitle("Number of tRNA genes per isotype in 15 Drosophilids")+
  facet_wrap(~species)

ggsave(totalIsotypeCounts, filename = '/Users/dylansosa/Library/Mobile Documents/com~apple~CloudDocs/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/totalIsotypeCounts.png')


save.image(f = 'tRNA_evolution.Rdata')
load('tRNA_evolution.Rdata')
