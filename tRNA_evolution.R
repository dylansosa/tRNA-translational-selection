library(tidyverse)
library(ggpmisc)
# install.packages("devtools")
# devtools::install_github("johannesbjork/LaCroixColoR")
library(LaCroixColoR)
library(gridExtra)
library(ggiraphExtra)
# install.packages('philentropy')
library(philentropy)
# install.packages('ggfortify')
library(ggfortify)
# install.packages("plotly")
library(plotly)

setwd(
  '/Users/dylansosa/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Yakuba  Long Lab/thesis/3_tRNA_evolution/'
)

# setwd('/Users/dylan/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/YakubaLong Lab/thesis/3_tRNA_evolution')

####################
#### load and format count data
####################

pacbioFiles <-
  list.files(path = 'pacbio/tRNAscan-SE_output/tRNA_gene_count_information/', pattern = '*.tsv')
pacbioDir <- 'pacbio/tRNAscan-SE_output/tRNA_gene_count_information'
nonpacbioFiles <-
  list.files(path = 'non-pacbio/tRNAscan-SE_output/tRNA_gene_count_information/', pattern = '*.tsv')
nonpacbioDir <-
  'non-pacbio/tRNAscan-SE_output/tRNA_gene_count_information'
# get path and file name

pacbioOutput <- file.path(pacbioDir, pacbioFiles)
nonpacbioOutput <- file.path(nonpacbioDir, nonpacbioFiles)
# combine to make path to files

for (file in pacbioOutput) {
  prefix <-
    str_remove(str_extract(str_extract(file, '[A-Z]{3}_[0-9].*'), '_[A-Z]?[a-z]+'), '_')
  assign(prefix, read_tsv(file))
}

for (file in nonpacbioOutput) {
  prefix <-
    str_remove(str_extract(str_extract(file, '[A-Z]{3}_[0-9].*'), '_[A-Z]?[a-z]+'), '_')
  assign(prefix, read_tsv(file))
}
# loop through paths and make unique variables based on assembly name

widenIsotypes <- function(speciesFile) {
  x <- deparse(substitute(speciesFile))
  speciesFile %>%
    dplyr::select(`Isotypes Used in this Species`,
                  `Number of tRNA Genes per Isotype`) %>%
    na.omit() %>%
    pivot_wider(names_from = `Isotypes Used in this Species`, values_from = `Number of tRNA Genes per Isotype`) %>%
    mutate(species = x) %>% dplyr::select(species, everything())
}
# function to do the formatting I want for plotting

Dana.wide <- widenIsotypes(Dana)
Dazt.wide <- widenIsotypes(Dazt)
Dere.wide <- widenIsotypes(Dere)
Dhyd.wide <- widenIsotypes(Dhyd)
Dmel.wide <- widenIsotypes(Dmel)
Dmoj.wide <- widenIsotypes(Dmoj)
Dnov.wide <- widenIsotypes(Dnov)
Dore.wide <- widenIsotypes(Dore)
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
  full_join(Dper.wide) %>%
  full_join(Dpse.wide) %>%
  full_join(Dsim.wide) %>%
  full_join(Dvir.wide) %>%
  full_join(Dwil.wide) %>%
  full_join(Dyak.wide) %>%
  full_join(Sleb.wide) %>%
  replace(is.na(.), 0)

####################
#### plot of tRNA gene counts per isotype
####################

totalWideIstoypeCounts$species = factor(
  totalWideIstoypeCounts$species,
  levels = c(
    'Dmel',
    'Dsim',
    'Dore',
    'Dere',
    'Dyak',
    'Dana',
    'Dpse',
    'Dper',
    'Dazt',
    'Dwil',
    'Dvir',
    'Dnov',
    'Dmoj',
    'Dhyd',
    'Sleb'
  )
)
# reorder to match phylogeny

totalIsotypeCounts <-
  ggRadar(
    totalWideIstoypeCounts,
    aes(group = species),
    rescale = FALSE,
    legend.position = "none",
    size = 2,
    interactive = FALSE,
    use.label = TRUE
  ) +
  theme(axis.text.x = element_text(size = 10)) +
  # scale_fill_manual(values = rep(mycolor, nrow(totalWideIstoypeCounts))) +
  scale_color_manual(values = lacroix_palette("Apricot", n = 15, type = "continuous")) +
  # LaCroix! :)
  ggtitle("Number of tRNA genes per isotype in 15 Drosophila species") +
  facet_wrap(~ species, nrow = 3)
totalIsotypeCounts

ggsave(totalIsotypeCounts, filename = '/Users/dylansosa/Library/Mobile Documents/com~apple~CloudDocs/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/totalIsotypeCounts.png')


spp <- totalWideIstoypeCounts %>% 
  select(species)

totalWideIstoypeCounts <- as.data.frame(totalWideIstoypeCounts)
rownames(totalWideIstoypeCounts) <- pull(spp)
    
distanceWideIsotypeCounts <-
  totalWideIstoypeCounts %>%
  select(!(species)) %>% 
  distance(method = 'euclidean',use.row.names = T)
# 15 spp x 15 spp
# I need to make rownames the spp in the original file during the distance computation
# then I will have spp on row

######## 
######## Visualizing the distance matrix of tRNA gene counts per isotype for each spp
######## 
  
pca_res <- prcomp(distanceWideIsotypeCounts, scale. = T)

distanceWideIsotypeCounts2 <- rownames_to_column(as.data.frame(distanceWideIsotypeCounts), var = 'Species')
distanceWideIsotypeCounts2$species <- as.factor(distanceWideIsotypeCounts2$species)

pcaPlotDistanceWideIsotypeCounts <-
  autoplot(
    pca_res,
    data = distanceWideIsotypeCounts2,
    colour = 'Species',
    label = T,
    label.vjust = 1.2,
    frame = T,
    legend = F
  ) +
  theme(legend.position = 'none') +
  scale_color_manual(values = lacroix_palette("Apricot", n = 15, type = "continuous"))+
  ggtitle("PCA plot of Euclidean distance of 15 Drosophila species' tRNA gene counts per isotype") 
pcaPlotDistanceWideIsotypeCounts

ggsave(pcaPlotDistanceWideIsotypeCounts, filename = '/Users/dylansosa/Library/Mobile Documents/com~apple~CloudDocs/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/pcaPlotDistanceWideIsotypeCounts.png')


####################
#### plot of anticodon counts per tRNA isotype
####################

widenAnticodons <- function(speciesFile) {
  x <- deparse(substitute(speciesFile))
  speciesFile %>%
    dplyr::select(`Anticodon`,
                  `Number of tRNA genes per Anticodon`) %>%
    na.omit() %>%
    pivot_wider(names_from = `Anticodon`, values_from = `Number of tRNA genes per Anticodon`,values_fn = sum) %>%
    mutate(species = x) %>% dplyr::select(species, everything())
}

Dana.wideAnticodons <- widenAnticodons(Dana)
Dazt.wideAnticodons <- widenAnticodons(Dazt)
Dere.wideAnticodons <- widenAnticodons(Dere)
Dhyd.wideAnticodons <- widenAnticodons(Dhyd)
Dmel.wideAnticodons <- widenAnticodons(Dmel)
Dmoj.wideAnticodons <- widenAnticodons(Dmoj)
Dnov.wideAnticodons <- widenAnticodons(Dnov)
Dore.wideAnticodons <- widenAnticodons(Dore)
Dper.wideAnticodons <- widenAnticodons(Dper)
Dpse.wideAnticodons <- widenAnticodons(Dpse)
Dsim.wideAnticodons <- widenAnticodons(Dsim)
Dvir.wideAnticodons <- widenAnticodons(Dvir)
Dwil.wideAnticodons <- widenAnticodons(Dwil)
Dyak.wideAnticodons <- widenAnticodons(Dyak)
Sleb.wideAnticodons <- widenAnticodons(Sleb)


totalWideAnticodons <-
  Dana.wideAnticodons %>%
  full_join(Dazt.wideAnticodons) %>%
  full_join(Dere.wideAnticodons) %>%
  full_join(Dhyd.wideAnticodons) %>%
  full_join(Dmel.wideAnticodons) %>%
  full_join(Dmoj.wideAnticodons) %>%
  full_join(Dnov.wideAnticodons) %>%
  full_join(Dore.wideAnticodons) %>%
  full_join(Dper.wideAnticodons) %>%
  full_join(Dpse.wideAnticodons) %>%
  full_join(Dsim.wideAnticodons) %>%
  full_join(Dvir.wideAnticodons) %>%
  full_join(Dwil.wideAnticodons) %>%
  full_join(Dyak.wideAnticodons) %>%
  full_join(Sleb.wideAnticodons) %>%
  replace(is.na(.), 0)

totalAnticodons <-
  ggRadar(
    totalWideAnticodons,
    aes(group = species),
    rescale = FALSE,
    legend.position = "none",
    size = 2,
    interactive = FALSE,
    use.label = TRUE
  ) +
  theme(axis.text.x = element_text(size = 10)) +
  scale_color_manual(values = lacroix_palette("Apricot", n = 15, type = "continuous")) +
  # LaCroix! :)
  ggtitle("Number of anticodons used in tRNA genes per isotype in 15 Drosophila species") +
  facet_wrap( ~ species, nrow=3)
totalAnticodons

ggsave(totalAnticodons, filename = '/Users/dylan/Library/Mobile Documents/com~apple~CloudDocs/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/totalAnticodons.png',width = 30,height = 30)

radarPlotsForAnticodons <- function(speciesFile) {
  x <- deparse(substitute(speciesFile))
  x.anticodons <-
    ggRadar(
      totalWideAnticodons %>% 
        filter(species == x),
      aes(group = species),
      rescale = FALSE,
      legend.position = "none",
      size = 2,
      interactive = FALSE,
      use.label = TRUE
    ) +
    theme(axis.text.x = element_text(size = 10)) +
    scale_color_manual(values = lacroix_palette("Apricot", n = 15, type = "continuous")) +
    ggtitle(paste("Number of anticodons used in",x))
}

Dana.anticodonRadarPlot <- radarPlotsForAnticodons(Dana)
Dazt.anticodonRadarPlot <- radarPlotsForAnticodons(Dazt)
Dere.anticodonRadarPlot <- radarPlotsForAnticodons(Dere)
Dhyd.anticodonRadarPlot <- radarPlotsForAnticodons(Dhyd)
Dmel.anticodonRadarPlot <- radarPlotsForAnticodons(Dmel)
Dmoj.anticodonRadarPlot <- radarPlotsForAnticodons(Dmoj)
Dnov.anticodonRadarPlot <- radarPlotsForAnticodons(Dnov)
Dore.anticodonRadarPlot <- radarPlotsForAnticodons(Dore)
Dper.anticodonRadarPlot <- radarPlotsForAnticodons(Dper)
Dpse.anticodonRadarPlot <- radarPlotsForAnticodons(Dpse)
Dsim.anticodonRadarPlot <- radarPlotsForAnticodons(Dsim)
Dvir.anticodonRadarPlot <- radarPlotsForAnticodons(Dvir)
Dwil.anticodonRadarPlot <- radarPlotsForAnticodons(Dwil)
Dyak.anticodonRadarPlot <- radarPlotsForAnticodons(Dyak)
Sleb.anticodonRadarPlot <- radarPlotsForAnticodons(Sleb)


saveAnticodonRadarPlots <- function(radar) {
  x <- deparse(substitute(radar))
  # print(radar)
  ggsave(radar, filename = paste0('/Users/dylan/Library/Mobile Documents/com~apple~CloudDocs/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/',x,'.png'))
}

# cat te | cut -f 1 -d ' ' | sed s'/^/saveAnticodonRadarPlots(/' | sed s'/$/)/'
saveAnticodonRadarPlots(Dana.anticodonRadarPlot)
saveAnticodonRadarPlots(Dazt.anticodonRadarPlot)
saveAnticodonRadarPlots(Dere.anticodonRadarPlot)
saveAnticodonRadarPlots(Dhyd.anticodonRadarPlot)
saveAnticodonRadarPlots(Dmel.anticodonRadarPlot)
saveAnticodonRadarPlots(Dmoj.anticodonRadarPlot)
saveAnticodonRadarPlots(Dnov.anticodonRadarPlot)
saveAnticodonRadarPlots(Dore.anticodonRadarPlot)
saveAnticodonRadarPlots(Dper.anticodonRadarPlot)
saveAnticodonRadarPlots(Dpse.anticodonRadarPlot)
saveAnticodonRadarPlots(Dsim.anticodonRadarPlot)
saveAnticodonRadarPlots(Dvir.anticodonRadarPlot)
saveAnticodonRadarPlots(Dwil.anticodonRadarPlot)
saveAnticodonRadarPlots(Dyak.anticodonRadarPlot)
saveAnticodonRadarPlots(Sleb.anticodonRadarPlot)

#######
### How many anticodons are used?
#######
anticodonsUsedByDrosophilids <- 
  totalWideAnticodons %>% 
  select(!(species)) %>%
  colnames(.) %>% 
  sort(.) 

  
save.image(f = 'tRNA_evolution.Rdata')
load('tRNA_evolution.Rdata')