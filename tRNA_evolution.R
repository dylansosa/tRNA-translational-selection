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
  # geom_text(label = '10',nudge_y = 5) +
  # geom_text() +
  # geom_label(aes(group = species),label='te')+
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

# library(cluster)

pcaPlotDistanceWideIsotypeCounts <-
  autoplot(
    pca_res,
    data = distanceWideIsotypeCounts2,
    colour = 'Species',
    label = T,
    label.vjust = 1.2,
    # frame = T,
    # frame.type = 'norm',
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

totalWideAnticodons$species = factor(
  totalWideAnticodons$species,
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
  # ylim(0,20) +
  theme(axis.text.x = element_text(size = 10)) +
  scale_color_manual(values = lacroix_palette("Apricot", n = 15, type = "continuous")) +
  # LaCroix! :)
  ggtitle("Number of anticodons used in tRNA genes per isotype in 15 Drosophila species") +
  facet_wrap( ~ species, nrow = 3)
totalAnticodons

ggsave(totalAnticodons, filename = '/Users/dylan/Library/Mobile Documents/com~apple~CloudDocs/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/totalAnticodons.png',width = 30,height = 30)

##########
#testing to improve resolution of figure 

# Dore_changed_AGT_value_wideAnticodo <- totalWideAnticodons
# Dore_changed_AGT_value_wideAnticodo$AGT[8] <- 22
# oops, accidentally overwrote the value 

# Dore_changed_AGT_value_wideAnticodo <-
#   ggRadar(
#     Dore_changed_AGT_value_wideAnticodo,
#     aes(group = species),
#     rescale = FALSE,
#     legend.position = "none",
#     size = 2,
#     interactive = FALSE,
#     use.label = TRUE
#   ) +
#   # ylim(0,20) +
#   theme(axis.text.x = element_text(size = 10)) +
#   scale_color_manual(values = lacroix_palette("Apricot", n = 15, type = "continuous")) +
#   # LaCroix! :)
#   ggtitle("Number of anticodons used in tRNA genes per isotype in 15 Drosophila species") +
#   facet_wrap( ~ species, nrow = 3)
# Dore_changed_AGT_value_wideAnticodo
# 
# # ggsave(Dore_changed_AGT_value_wideAnticodo, filename = '/Users/dylan/Library/Mobile Documents/com~apple~CloudDocs/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/Dore_changed_AGT_value_wideAnticodo.png',width = 30,height = 30)
# ggsave(Dore_changed_AGT_value_wideAnticodo, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/resolution_testing/Dore_changed_AGT_value_wideAnticodo.png',width = 30,height = 30)
# ggsave(totalAnticodons, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/resolution_testing/totalAnticodons.png',width = 30,height = 30)
# # totalAnticodons



#testing to improve resolution of figure 
##########

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
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 16)) +
    theme(title = element_text(size = 16)) +
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
  ggsave(radar, filename = paste0('/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/tRNA figures/individual anticodon plots/',x,'.png'))
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

fourSpeciesAnticodon <- grid.arrange(Dmel.anticodonRadarPlot, Dsim.anticodonRadarPlot,Dore.anticodonRadarPlot,Dere.anticodonRadarPlot,ncol=2)
twoSpeciesAnticodon <- grid.arrange(Dore.anticodonRadarPlot,Dere.anticodonRadarPlot,ncol=2)
# I think this is the best way to display several (not all) anticodon plots
ggsave(twoSpeciesAnticodon, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/tRNA figures/multi_codon_plot/twoSpeciesAnticodon.png',width = 30,height = 30)


# grid.arrange(Dmel.anticodonRadarPlot,
#              Dsim.anticodonRadarPlot,
#              Dore.anticodonRadarPlot,
#              Dere.anticodonRadarPlot,
#              Dyak.anticodonRadarPlot,
#              Dana.anticodonRadarPlot,
#              Dpse.anticodonRadarPlot,
#              Dper.anticodonRadarPlot,
#              Dazt.anticodonRadarPlot,
#              Dwil.anticodonRadarPlot,
#              Dvir.anticodonRadarPlot,
#              Dnov.anticodonRadarPlot,
#              Dmoj.anticodonRadarPlot,
#              Dhyd.anticodonRadarPlot,
#              Sleb.anticodonRadarPlot,
#              ncol=5
#              )

#######
### How many anticodons are used?
#######
anticodonsUsedByDrosophilids <- 
  totalWideAnticodons %>% 
  select(!(species)) %>%
  colnames(.) %>% 
  sort(.) 
anticodonsUsedByDrosophilids
length(anticodonsUsedByDrosophilids)
# 51 anticodons used amongst the 15 drosophila species 


######
### Joining total anticodons with count table to make plot of percentage used 
#####

anticodons_64 <- read.csv('anticodons_64.txt', header = F)
rownames(anticodons_64) <- pull(anticodons_64)
anticodons_64 <- 
  anticodons_64 %>% 
  t() %>% 
  as_tibble()
anticodons_64[,] <- NA
# preparing format for joining 
anticodons_64

wide64Anticodons <- 
  full_join(totalWideAnticodons,anticodons_64) %>% 
  slice(.,1:15) %>% 
  replace(is.na(.), 0) 
  
# now in the same format as the totalWide* variables

wide64Anticodons <-
  wide64Anticodons %>% 
  mutate(totalAnticodons = rowSums(across(where(is.numeric)))) %>% 
  relocate(totalAnticodons,.before = AGC)
# now have column with total anticodons per species
# I can use this during plotting to get percentages per anticodon

# ggplot(wide64Anticodons) +
  # geom_point(aes(x = group(),y = Dmel))

# library(reshape2)
# z <- melt(wide64Anticodons, id.vars="species")%>% filter(species == 'Dmel') 

# totalDmel_anti <- wide64Anticodons[5,"totalAnticodons"]

# zz <- z %>% filter(variable!='totalAnticodons') %>% rowwise() %>% mutate(percentage = ((value/totalDmel_anti)))

# may need a function to do this for each species if I want them to all be plotted together
# currently this is just Dmel
# 
# ggplot(zz,aes(reorder(variable, percentage$totalAnticodons),percentage$totalAnticodons)) +
#   geom_point() +
#   theme(axis.text.x = element_text(angle = 70, vjust = 0.5, hjust=1))

# Dpse.wideAnticodons %>% 
#   select(-species) %>% 
#   full_join(.,anticodons_64) %>% 
#   slice(.,1) %>% 
#   replace(is.na(.),0) %>% 
#   rowwise(.) %>% 
#   sum()
  # mutate(total = rowwise(.) %>% sum())
# this gives me the total anticodons, after transformation the unused ones to 0, make column with total anticodons 'total'

Dpse$`Total tRNA genes`[1]
# I can just use this per species to get the total

# Dpse.wideAnticodons %>%
#   select(-species) %>%
#   full_join(.,anticodons_64) %>%
#   slice(.,1) %>%
#   replace(is.na(.),0) %>%
#   t() %>% 
#   as.data.frame() %>% 
#   mutate(percentage = V1/Dpse$`Total tRNA genes`[1]) %>% 
#   rownames_to_column()
# # change colnames later

# ggplot(Dpse.wideAnticodons %>%
#          select(-species) %>%
#          full_join(.,anticodons_64) %>%
#          slice(.,1) %>%
#          replace(is.na(.),0) %>%
#          t() %>%
#          as.data.frame() %>%
#          mutate(percentage = V1/Dpse$`Total tRNA genes`[1]) %>%
#          rownames_to_column(),
#        aes(x=reorder(rowname,percentage),y=percentage)) +
#   geom_point() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(title = 'Anticodon usage in Dpse') +
#   xlab(label = 'Anticodon') +
#   ylab(label = 'Percentage of total tRNA genes')

scatterPlotsForAnticodonPercentages <- function(speciesFile) {
    x <- deparse(substitute(speciesFile))
      ggplot(eval(parse(text=paste0(x,'.wideAnticodons'))) %>%
               select(-species) %>%
               full_join(.,anticodons_64) %>%
               slice(.,1) %>%
               replace(is.na(.),0) %>%
               t() %>%
               as.data.frame() %>%
               mutate(percentage = V1/Dpse$`Total tRNA genes`[1]) %>%
               rownames_to_column(),
             aes(x=reorder(rowname,percentage),y=percentage)) +
      geom_point() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      # labs(title = 'Anticodon usage') +
      xlab(label = 'Anticodon') +
      ylab(label = 'Percentage of anticodons used in all tRNA genes') +
      # theme(axis.text.x = element_text(size = 12)) +
      # theme(axis.text.y = element_text(size = 12)) +
      # theme(title = element_text(size = 16)) +
      # scale_color_manual(values = lacroix_palette("Apricot", n = 15, type = "continuous")) +
      ggtitle(paste("Percentage of anticodons used in",x, 'tRNA genes'))
}

Dana.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dana)
Dazt.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dazt)
Dere.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dere)
Dhyd.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dhyd)
Dmel.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dmel)
Dmoj.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dmoj)
Dnov.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dnov)
Dore.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dore)
Dper.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dper)
Dpse.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dpse)
Dsim.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dsim)
Dvir.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dvir)
Dwil.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dwil)
Dyak.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Dyak)
Sleb.anticodonScatterPlot <- scatterPlotsForAnticodonPercentages(Sleb)


saveAnticodonScatterPlots <- function(scatter) {
  x <- deparse(substitute(scatter))
  ggsave(scatter, filename = paste0('/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/tRNA figures/individual_percentage_of_anticodons_plots/',x,'.png'))
}

saveAnticodonScatterPlots(Dana.anticodonScatterPlot)
saveAnticodonScatterPlots(Dazt.anticodonScatterPlot)
saveAnticodonScatterPlots(Dere.anticodonScatterPlot)
saveAnticodonScatterPlots(Dhyd.anticodonScatterPlot)
saveAnticodonScatterPlots(Dmel.anticodonScatterPlot)
saveAnticodonScatterPlots(Dmoj.anticodonScatterPlot)
saveAnticodonScatterPlots(Dnov.anticodonScatterPlot)
saveAnticodonScatterPlots(Dore.anticodonScatterPlot)
saveAnticodonScatterPlots(Dper.anticodonScatterPlot)
saveAnticodonScatterPlots(Dpse.anticodonScatterPlot)
saveAnticodonScatterPlots(Dsim.anticodonScatterPlot)
saveAnticodonScatterPlots(Dvir.anticodonScatterPlot)
saveAnticodonScatterPlots(Dwil.anticodonScatterPlot)
saveAnticodonScatterPlots(Dyak.anticodonScatterPlot)
saveAnticodonScatterPlots(Sleb.anticodonScatterPlot)


allSpeciesScatterPlot <- 
  grid.arrange(Dmel.anticodonScatterPlot,
             Dsim.anticodonScatterPlot,
             Dore.anticodonScatterPlot,
             Dere.anticodonScatterPlot,
             Dyak.anticodonScatterPlot,
             Dana.anticodonScatterPlot,
             Dpse.anticodonScatterPlot,
             Dper.anticodonScatterPlot,
             Dazt.anticodonScatterPlot,
             Dwil.anticodonScatterPlot,
             Dvir.anticodonScatterPlot,
             Dnov.anticodonScatterPlot,
             Dmoj.anticodonScatterPlot,
             Dhyd.anticodonScatterPlot,
             Sleb.anticodonScatterPlot,
             ncol=5
             )
ggsave(allSpeciesScatterPlot, filename = paste0('/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA_evolution_pacbio/figures/tRNA figures/total anticodon scatter plot/allSpeciesScatterPlot.png'),width = 40,height = 15)
# maybe I should not order the points based on percentage? It would be easier to compare them if not?

save.image(f = 'tRNA_evolution.Rdata')
load('tRNA_evolution.Rdata')
