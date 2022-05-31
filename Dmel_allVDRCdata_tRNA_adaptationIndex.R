# devtools::install_github("dylan# devtools::install_github("dylansosa/tai")
# devtools::install_github("karthik/wesanderson")
require("tAI")
library('tidyverse')
library(ggpmisc)
library(LaCroixColoR)
library(gridExtra)
library(wesanderson)


setwd('/Users/dylansosa/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Yakuba  Long Lab/thesis/2_tRNA_adaptation')
# setwd('/Users/dylansosa/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/YakubaUChicago-1/thesis/2_tRNA_adaptation')
# setwd('/Users/dylansosa/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/YakubaUChicago/thesis/2_tRNA_adaptation')

#######################################################
#######################################################
# init Dmel data
#######################################################
####################################################### 
dmel.trna <- scan("Dmel.tRNA")
# load tRNA gene copy numbers 
# I computed these data using all CDS of D. mel
# need to truncate it to old/ancient later

dmel.trna.ws <- get.ws(tRNA=dmel.trna, sking = 1)
# compute relative adaptiveness values for each codon in the D. mel genes

#######################################################
#######################################################
# essentiality based tAI
# below is to test on different ages
#######################################################
#######################################################
dmel.youngE.m <- matrix(scan("dmel.vdrc.longestCDS.youngE.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script: 
#   a matrix of codon frequencies per ORF.

# dmel.m <- matrix(scan("Dmel.matrixCodonFreq.txt"), ncol=61, byrow=TRUE)
# this is for all gene data 
# this wasn't working before because of some isoforms having no results...

dmel.youngE.m <- dmel.youngE.m[,-33]
# ignore methionine codons

dmel.youngE.tai <- get.tai(dmel.youngE.m, dmel.trna.ws)
# calculate tRNA adaptation index 

# write.table(dmel.youngE.tai, file = 'tAIoutput/dmel.vdrc.youngE.tai.txt',sep = '\t',quote = F)

hist(dmel.youngE.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.
ggplot()+
  # geom_histogram(aes(dmel.youngE.tai),bins = length(dmel.youngE.tai))
  geom_histogram(aes(dmel.youngE.tai),binwidth = 0.05,)
    # geom_histogram(aes(dmel.youngE.tai))

# df <- read.table("dmel.youngE.cds.w", header=TRUE)
# df <- read.table("dmel.vdrc.longestCDS.youngE.w", header=TRUE)
df.youngE <- read.table("dmel.vdrc.longestCDS.youngE.w", header=TRUE, na.strings = "*****")
# cat dmel.vdrc.longestCDS.youngE.fa | egrep '>' | sed s'/>//g' > dmel.youngE.id.temp
# 
# cat dmel.youngE.tai.txt | sed '1d' | cut -f 2 > dmel.youngE.tai.temp
# 
# paste dmel.youngE.id.temp dmel.youngE.tai.temp > dmel.youngE.tai


plot(dmel.youngE.tai ~ df.youngE$Nc)
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons), 
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
# for these young essential genes I am not seeing that!
# Looks like they have low adaptiveness and have low effective codon 
mean(na.omit(df.youngE$Nc))
# [1] 51.10894


cor(dmel.youngE.tai, df.youngE$Nc, use="p")
# [1] -0.3810451
# [1] -0.3848419
# this is for all vdrc, first is only for sig DE
# [1] -0.3944438 after sking=1

dmel.youngE.s <- get.s(dmel.youngE.tai, df.youngE$Nc, df.youngE$GC3s)
dmel.youngE.s
# [1] 0.4817288
# [1] 0.4767317
# second is all vdrc
# [1] 0.4695291 after sking



# Formally, we can calculate the correlation between tAI, 
# and the corrected Nc, f(GC3s) – Nc, we call this correlation S, 
# because it reflects the intensity of translational selection 
# acting on our sample of genes
#######################################################
#######################################################

dmel.youngNE.m <- matrix(scan("dmel.vdrc.longestCDS.youngNE.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script: 
#   a matrix of codon frequencies per ORF.

dmel.youngNE.m <- dmel.youngNE.m[,-33]
# ignore methionine codons

dmel.youngNE.tai <- get.tai(dmel.youngNE.m, dmel.trna.ws)
# calculate tRNA adaptation index 

# write.table(dmel.youngNE.tai, file = 'tAIoutput/dmel.vdrc.youngNE.tai.txt',sep = '\t',quote = F)

hist(dmel.youngNE.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.
ggplot()+
  geom_histogram(aes(dmel.youngNE.tai))

# df.youngNE <- read.table("dmel.ancientE.cds.w", header=TRUE)
# df.youngNE <- read.table("dmel.cds.young.nonessential.w", header=TRUE)
df.youngNE <- read.table("dmel.vdrc.longestCDS.youngNE.w", header=TRUE, na.strings = "*****")
# this is where errors are,   line 368 did not have 5 elements
# same error when I try to run this script on all Dmel genes
# it just happenes to be that the youngE genes don't have the problem of not computing each col 

plot(dmel.youngNE.tai ~ df.youngNE$Nc)
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons), 
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
mean(na.omit(df.youngNE$Nc))
# [1] 51.34412


cor(dmel.youngNE.tai, df.youngNE$Nc, use="p")
# [1] -0.2091027
# [1] -0.2814465 vdrc
# [1] -0.3084078
dmel.youngNE.s <- get.s(dmel.youngNE.tai, df.youngNE$Nc, df.youngNE$GC3s)
dmel.youngNE.s
# [1] 0.3199346
# [1] 0.350783 vdrc
# 0.3682818

# Formally, we can calculate the correlation between tAI, 
# and the corrected Nc, f(GC3s) – Nc, we call this correlation S, 
# because it reflects the intensity of translational selection 
# acting on our sample of genes:

########################################################
########################################################
dmel.ancientE.m <- matrix(scan("dmel.vdrc.longestCDS.ancientE.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script: 
#   a matrix of codon frequencies per ORF.

dmel.ancientE.m <- dmel.ancientE.m[,-33]
# ignore methionine codons

dmel.ancientE.tai <- get.tai(dmel.ancientE.m, dmel.trna.ws)
# calculate tRNA adaptation index 

# write.table(dmel.ancientE.tai, file = 'tAIoutput/dmel.vdrc.ancientE.tai.txt',sep = '\t',quote = F)

hist(dmel.ancientE.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.

# df.ancientE <- read.table("dmel.ancientE.cds.w", header=TRUE)
# df.ancientE <- read.table("dmel.cds.young.nonessential.w", header=TRUE)
df.ancientE <- read.table("dmel.vdrc.longestCDS.ancientE.w", header=TRUE, na.strings = "*****")
plot(dmel.ancientE.tai ~ df.ancientE$Nc)
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons), 
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
mean(na.omit(df.ancientE$Nc))
# [1] 49.50832


cor(dmel.ancientE.tai, df.ancientE$Nc, use="p")
# [1] -0.240699
# [1] -0.2682363 vdrc
# [1] -0.268946 after realizing i didnt remake codonZ/M files after adding missing data 
# [1] -0.2850244



dmel.ancientE.s <- get.s(dmel.ancientE.tai, df.ancientE$Nc, df.ancientE$GC3s)
dmel.ancientE.s
# [1] 0.3053252
# [1] 0.3556138 vdrc
# [1] 0.3562489 final
# [1] 0.356183 after sking change



# Formally, we can calculate the correlation between tAI, 
# and the corrected Nc, f(GC3s) – Nc, we call this correlation S, 
# because it reflects the intensity of translational selection 
# acting on our sample of genes:

########################################################
########################################################

dmel.ancientNE.m <- matrix(scan("dmel.vdrc.longestCDS.ancientNE.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script: 
#   a matrix of codon frequencies per ORF.

dmel.ancientNE.m <- dmel.ancientNE.m[,-33]
# ignore methionine codon

dmel.ancientNE.tai <- get.tai(dmel.ancientNE.m, dmel.trna.ws)
# calculate tRNA adaptation index 

# write.table(dmel.ancientNE.tai, file = 'tAIoutput/dmel.vdrc.ancientNE.tai.txt',sep = '\t',quote = F)


hist(dmel.ancientNE.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.

# df.ancientNE <- read.table("dmel.ancientE.cds.w", header=TRUE)
# df.ancientNE <- read.table("dmel.cds.young.nonessential.w", header=TRUE)
df.ancientNE <- read.table("dmel.vdrc.longestCDS.ancientNE.w", header=TRUE, na.strings = "*****")
plot(dmel.ancientNE.tai ~ df.ancientNE$Nc)
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons), 
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
mean(na.omit(df.ancientNE$Nc))
# [1] 49.61155



cor(dmel.ancientNE.tai, df.ancientNE$Nc, use="p")
# [1] -0.249064
# [1] -0.2914642 vdrc
# [1] -0.2913454 after fixing the codonZ/M files (remaking them)
# [1] -0.3094459




dmel.ancientNE.s <- get.s(dmel.ancientNE.tai, df.ancientNE$Nc, df.ancientNE$GC3s)
dmel.ancientNE.s
# [1] 0.3483473
# [1] 0.3954093 vdrc
# [1] 0.3946952 final after remaking 
# [1] 0.3959318




# Formally, we can calculate the correlation between tAI, 
# and the corrected Nc, f(GC3s) – Nc, we call this correlation S, 
# because it reflects the intensity of translational selection 
# acting on our sample of genes:


########################################################
#######################################################
#######################################################
# total young or ancient
#######################################################
#######################################################

dmel.young.m <- matrix(scan("dmel.vdrc.longestCDS.young.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script: 
#   a matrix of codon frequencies per ORF.

# dmel.m <- matrix(scan("Dmel.matrixCodonFreq.txt"), ncol=61, byrow=TRUE)
# this is for all gene data 
# this wasn't working before because of some isoforms having no results...

dmel.young.m <- dmel.young.m[,-33]
# ignore methionine codons

dmel.young.tai <- get.tai(dmel.young.m, dmel.trna.ws)
# calculate tRNA adaptation index 

# write.table(dmel.young.tai, file = 'tAIoutput/dmel.vdrc.young.tai.txt',sep = '\t',quote = F)

hist(dmel.young.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.

# df <- read.table("dmel.youngE.cds.w", header=TRUE)
# df <- read.table("dmel.vdrc.longestCDS.youngE.w", header=TRUE)
df.young <- read.table("dmel.vdrc.longestCDS.young.w", header=TRUE, na.strings = "*****")
# cat dmel.vdrc.longestCDS.youngE.fa | egrep '>' | sed s'/>//g' > dmel.youngE.id.temp
# 
# cat dmel.youngE.tai.txt | sed '1d' | cut -f 2 > dmel.youngE.tai.temp
# 
# paste dmel.youngE.id.temp dmel.youngE.tai.temp > dmel.youngE.tai


plot(dmel.young.tai ~ df.young$Nc)
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons), 
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
# for these young essential genes I am not seeing that!
# Looks like they have low adaptiveness and have low effective codon 
# params <- lm(dmel.young.tai ~ df.young$Nc)$coef
# intercept <- params[1]
# slope <- params[2]
young_corr<-ggplot() +
  geom_point(aes(y=dmel.young.tai, x=df.young$Nc))+
  labs(title = 'Correlation between tAI and Nc for 702 young genes')+
  xlab(label = 'Nc (effective number of codons)') +
  ylab(label = 'tAI (tRNA adaptation index)') +
  geom_smooth(aes(y=dmel.young.tai, x=df.young$Nc), method = "lm",se = F)
# +
  # geom_hline(yintercept = .35)

mean(na.omit(df.young$Nc))
# [1] 51.28077

mean(na.omit(dmel.young.tai))
# [1] 0.1941226


cor(dmel.young.tai, df.young$Nc, use="p")
# [1] -0.3305453


dmel.young.s <- get.s(dmel.young.tai, df.young$Nc, df.young$GC3s)
dmel.young.s
# [1] 0.3945509

#######################################################
#######################################################

dmel.ancient.m <- matrix(scan("dmel.vdrc.longestCDS.ancient.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script: 
#   a matrix of codon frequencies per ORF.

# dmel.m <- matrix(scan("Dmel.matrixCodonFreq.txt"), ncol=61, byrow=TRUE)
# this is for all gene data 
# this wasn't working before because of some isoforms having no results...

dmel.ancient.m <- dmel.ancient.m[,-33]
# ignore methionine codons

dmel.ancient.tai <- get.tai(dmel.ancient.m, dmel.trna.ws)
# calculate tRNA adaptation index 

# write.table(dmel.ancient.tai, file = 'tAIoutput/dmel.vdrc.ancient.tai.txt',sep = '\t',quote = F)

hist(dmel.ancient.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.

# df <- read.table("dmel.ancientE.cds.w", header=TRUE)
# df <- read.table("dmel.vdrc.longestCDS.ancientE.w", header=TRUE)
df.ancient <- read.table("dmel.vdrc.longestCDS.ancient.w", header=TRUE, na.strings = "*****")
# cat dmel.vdrc.longestCDS.ancientE.fa | egrep '>' | sed s'/>//g' > dmel.ancientE.id.temp
# 
# cat dmel.ancientE.tai.txt | sed '1d' | cut -f 2 > dmel.ancientE.tai.temp
# 
# paste dmel.ancientE.id.temp dmel.ancientE.tai.temp > dmel.ancientE.tai

plot(dmel.ancient.tai ~ df.ancient$Nc)
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons), 
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
# for these ancient essential genes I am not seeing that!
# Looks like they have low adaptiveness and have low effective codon 
# params <- lm(dmel.ancient.tai ~ df.ancient$Nc)$coef
# intercept <- params[1]
# slope <- params[2]
# ggplot() +
#   geom_point(aes(y=dmel.ancient.tai, x=df.ancient$Nc))+
#   labs(title = 'Correlation between tAI and Nc for 10652 ancient genes')+
#   xlab(label = 'Nc (effective number of codons)') +
#   ylab(label = 'tAI (tRNA adaptation index)') +
#   geom_abline(slope=slope, intercept=intercept, color='red')
# library(ggpmisc)
# formula <- y ~ x
ancient_corr <- ggplot() +
  geom_point(aes(y=dmel.ancient.tai, x=df.ancient$Nc))+
  labs(title = 'Correlation between tAI and Nc for 10652 ancient genes')+
  xlab(label = 'Nc (effective number of codons)') +
  ylab(label = 'tAI (tRNA adaptation index)') +
  geom_smooth(aes(y=dmel.ancient.tai, x=df.ancient$Nc), method = "lm",se = F) 
# +
#   stat_poly_eq(formula = formula, 
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                parse = TRUE) 
  
# +
  # stat_regline_equation(label.y = 0.4, aes(label = ..rr.label..))
  # stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")
  # stat_poly_eq(formula = y ~x, 
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #              parse = TRUE) +  
# geom_smooth(method = "lm")
# stat_smooth(method="lm", se=FALSE)
# geom_line() +
# geom_smooth(method = "lm", se = FALSE)

mean(na.omit(df.ancient$Nc))
# [1] 49.58551

mean(na.omit(dmel.ancient.tai))
# [1] 0.2030969


plot(df.ancient$GC3s,df.ancient$Nc)

cor(dmel.ancient.tai, df.ancient$Nc, use="p")
# [1] -0.3032757

dmel.ancient.s <- get.s(dmel.ancient.tai, df.ancient$Nc, df.ancient$GC3s)
dmel.ancient.s
# [1] 0.3860509

correlationPlots <- grid.arrange(young_corr,ancient_corr,nrow=1)
ggsave(correlationPlots,filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/correlation_ageBased.png')

#######################################################
#######################################################
####### total results 
#######################################################
#######################################################

####### OLD VERSION
#######
################
# dmel.m <- matrix(scan("dmel.vdrc.longestCDS.m"), ncol=61, byrow=TRUE)
# # read output of codonM
# # .m contains the output of the codonM script: 
# #   a matrix of codon frequencies per ORF.
# 
# # dmel.m <- matrix(scan("Dmel.matrixCodonFreq.txt"), ncol=61, byrow=TRUE)
# # this is for all gene data 
# # this wasn't working before because of some isoforms having no results...
# 
# dmel.m <- dmel.m[,-33]
# # ignore methionine codons
# 
# dmel.tai <- get.tai(dmel.m, dmel.trna.ws)
# # calculate tRNA adaptation index 
# write.table(dmel.tai, file = 'tAIoutput/dmel.vdrc.tai.txt',sep = '\t',quote = F)
# 
# hist(dmel.tai)
# # Highly expressed genes present high tAI values (> 0.4), 
# # which means that their codon usage resembles the genomic 
# # structure of tRNA genes.
# 
# # df <- read.table("dmel.ancientE.cds.w", header=TRUE)
# # df <- read.table("dmel.vdrc.longestCDS.ancientE.w", header=TRUE)
# df <- read.table("dmel.vdrc.longestCDS.w", header=TRUE, na.strings = "*****")
# # cat dmel.vdrc.longestCDS.ancientE.fa | egrep '>' | sed s'/>//g' > dmel.ancientE.id.temp
# # 
# # cat dmel.ancientE.tai.txt | sed '1d' | cut -f 2 > dmel.ancientE.tai.temp
# # 
# # paste dmel.ancientE.id.temp dmel.ancientE.tai.temp > dmel.ancientE.tai
# 
# 
# plot(dmel.tai ~ df$Nc)
# # abline()
# # plot of relationship between tAI and effective number of codons
# # genes with very low Nc values (highly biased toward certain codons), 
# # correlate with high tAI values (highly co-adapted to the tRNA gene pool).
# # for these ancient essential genes I am not seeing that!
# # Looks like they have low adaptiveness and have low effective codon 
# 
# 
# 
# # this is the ggplot form of the same thing. Consider 
# # mean(na.omit(df$Nc))
# # mean(dmel.tai)
# # 
# # cor(dmel.tai, df$Nc, use="p")
# # # [1] -0.3067634
# # 
# # dmel.s <- get.s(dmel.tai, df$Nc, df$GC3s)
# # dmel.s
# # # [1] 0.3873286


dmel.m <- matrix(scan("dmel.total.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script:
#   a matrix of codon frequencies per ORF.

dmel.m <- dmel.m[,-33]
# ignore methionine codons

dmel.tai <- get.tai(dmel.m, dmel.trna.ws)
# calculate tRNA adaptation index
write.table(dmel.tai, file = 'tAIoutput/dmel.total.tai.txt',sep = '\t',quote = F)

hist(dmel.tai)
# Highly expressed genes present high tAI values (> 0.4),
# which means that their codon usage resembles the genomic
# structure of tRNA genes.

# df <- read.table("dmel.ancientE.cds.w", header=TRUE)
# df <- read.table("dmel.vdrc.longestCDS.ancientE.w", header=TRUE)
dmel.df <- read.table("dmel.total.w", header=TRUE, na.strings = "*****")
# cat dmel.vdrc.longestCDS.ancientE.fa | egrep '>' | sed s'/>//g' > dmel.ancientE.id.temp
#
# cat dmel.ancientE.tai.txt | sed '1d' | cut -f 2 > dmel.ancientE.tai.temp
#
# paste dmel.ancientE.id.temp dmel.ancientE.tai.temp > dmel.ancientE.tai


plot(dmel.tai ~ dmel.df$Nc)
# abline()
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons),
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
# for these ancient essential genes I am not seeing that!
# Looks like they have low adaptiveness and have low effective codon
dmelTotalCorrelation <- ggplot() +
  geom_point(aes(y=dmel.tai, x=dmel.df$Nc))+
  labs(title = 'Correlation between tAI and Nc for 13,698 D. mel CDS')+
  xlab(label = 'Nc (effective number of codons)') +
  ylab(label = 'tAI (tRNA adaptation index)') +
  geom_smooth(aes(y=dmel.ancient.tai, x=df.ancient$Nc), method = "lm",se = F) 

ggsave(dmelTotalCorrelation,filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/dmelTotalCorrelation.png')


# this is the ggplot form of the same thing. Consider
# mean(na.omit(df$Nc))
# mean(dmel.tai)
#
# cor(dmel.tai, df$Nc, use="p")
# # [1] -0.3067634
#
# dmel.s <- get.s(dmel.tai, df$Nc, df$GC3s)
# dmel.s
# # [1] 0.3873286


#######################################################
#######################################################
# Welch's t-test between groups 
#######################################################
#######################################################
t.test(dmel.tai,dmel.young.tai)
# data:  dmel.tai and dmel.young.tai
# t = 2.6366, df = 768.61, p-value = 0.008542
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.001622126 0.011077118
# sample estimates:
#   mean of x mean of y 
# 0.2004722 0.1941226

t.test(dmel.tai,dmel.youngE.tai)
# data:  dmel.tai and dmel.youngE.tai
# t = 1.4316, df = 195.29, p-value = 0.1539
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.002384383  0.015013506
# sample estimates:
#   mean of x mean of y 
# 0.2004722 0.1941576

t.test(dmel.ancient.tai,dmel.young.tai)
# data:  dmel.ancient.tai and dmel.young.tai
# t = 3.7078, df = 784.19, p-value = 0.0002238
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.004223075 0.013725584
# sample estimates:
#   mean of x mean of y 
# 0.2030969 0.1941226 


t.test(dmel.ancientE.tai,dmel.youngE.tai)
# data:  dmel.ancientE.tai and dmel.youngE.tai
# t = 2.025, df = 216.26, p-value = 0.0441
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.0002444128 0.0180823134
# sample estimates:
#   mean of x mean of y 
# 0.2033210 0.1941576

t.test(dmel.tai,dmel.ancient.tai)
# data:  dmel.tai and dmel.ancient.tai
# t = -3.4193, df = 23234, p-value = 0.000629
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.004129304 -0.001120112
# sample estimates:
#   mean of x mean of y 
# 0.2004722 0.2030969


#######################################################
#######################################################
######## Density plotting of tAI
#######################################################
#######################################################

youngEresults <- read.csv('tAIoutput/dmel.vdrc.youngE.tai',sep = '\t',header = F)
ancientEresults <- read.csv('tAIoutput/dmel.vdrc.ancientE.tai',sep = '\t',header = F)
youngNEresults <-read.csv('tAIoutput/dmel.vdrc.youngNE.tai',sep = '\t',header = F)
ancientNEresults <- read.csv('tAIoutput/dmel.vdrc.ancientNE.tai',sep = '\t',header = F)

x <- ggplot(youngEresults) +
  geom_density(aes(V2)) +
  labs(title = 'young essential gene tAI')+
  xlab(label = 'tAI (tRNA adaptation index) value')

y <- ggplot(ancientEresults) +
  geom_density(aes(V2)) +
  labs(title = 'ancient essential gene tAI')+
  xlab(label = 'tAI (tRNA adaptation index) value')

a <- ggplot(youngNEresults) +
  geom_density(aes(V2)) +
  labs(title = 'young non-essential gene tAI')+
  xlab(label = 'tAI (tRNA adaptation index) value')


b <- ggplot(ancientNEresults) +
  geom_density(aes(V2)) +
  labs(title = 'ancient non-essential gene tAI')+
  xlab(label = 'tAI (tRNA adaptation index) value')



# mean(as.vector(ancientEresults))
# class(ancientEresults)
# mean(ancientNEresults)

groupDensity <- grid.arrange(x,y,a,b,ncol=2)
ggsave(groupDensity, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/vdrc_groupDensity2.png')

totalYoungresults <- read.csv('tAIoutput/dmel.vdrc.young.tai',sep = '\t',header = F)
totalAncientresults <- read.csv('tAIoutput/dmel.vdrc.ancient.tai',sep = '\t',header = F)

n <- ggplot(totalYoungresults) +
  geom_density(aes(V2)) +
  labs(title = 'all young gene tAI',xlab(''))+
  xlab(label = 'tAI (tRNA adaptation index) value')


o <- ggplot(totalAncientresults) +
  geom_density(aes(V2)) +
  labs(title = 'all ancient gene tAI') +
  xlab(label = 'tAI (tRNA adaptation index) value')

  

allDensity <- grid.arrange(x,y,a,b,n,o,nrow=3)
ggsave(allDensity, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/vdrc_all_groupDensity.png')


totalAgeDensity <- grid.arrange(n,o)
ggsave(totalAgeDensity, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/vdrc_total_young_ancientDensity.png')

## de novo in melanogaster (my data and yong zhang's annotations)
## li zhang's oryza data, i want to replciate his findings and also compare plant and insect 
youngDeNovoDensity <- read.csv('tAIoutput/dmel.vdrc.young.deNovo.tai',sep = '\t',header = F)
youngDeNovo <- ggplot(youngDeNovoDensity) +
  geom_density(aes(V2)) +
  labs(title = 'tAI of young genes annotated as de novo by Yong Zhang')+
  xlab(label = 'tAI (tRNA adaptation index) value')

ggsave(youngDeNovo, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/youngDeNovoDensity.png')
# 
oryzaDeNovo <- read.csv('oryza.denovo.tai',sep = '\t',header = F)
oryzaDeNovoDensity <- ggplot(oryzaDeNovo) +
  geom_density(aes(V2)) +
  labs(title = 'tAI of de novo Oryza genes from Li Zhang paper')+
  xlab(label = 'tAI (tRNA adaptation index) value')

ggsave(youngDeNovo, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/oryzaDeNovoDensity.png')

DeNovoDensities <- grid.arrange(youngDeNovo,oryzaDeNovoDensity)
ggsave(DeNovoDensities, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/DeNovoDensities.png')

oryzaOldSingletons <- read.csv('oryza.oldsingletons.tai',sep = '\t',header = F)
oryzaOldSingletonsDensity <- ggplot(oryzaOldSingletons) +
  geom_density(aes(V2)) +
  labs(title = 'tAI of old singletons Oryza genes from Li Zhang paper')+
  xlab(label = 'tAI (tRNA adaptation index) value')

ggsave(oryzaOldSingletonsDensity, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/oryzaOldSingletonsDensity.png')

DeNovoAndAncientSingletonDensities <- grid.arrange(youngDeNovo,oryzaDeNovoDensity,o,oryzaOldSingletonsDensity)
ggsave(DeNovoAndAncientSingletonDensities, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/DeNovoAndAncientSingletonDensities.png')


#######################################################
#######################################################
# more t-testing 
#######################################################
#######################################################
# mean(oryzaDeNovo$V2)
# mean(oryzaOldSingletons$V2)
# t.test(oryzaDeNovo$V2,oryzaOldSingletons$V2)
# # same values as li zhang's paper in which he reported them to be signicantly different  

# mean(youngEresults$V2)
# mean(ancientEresults$V2)
# t.test(youngEresults$V2,ancientEresults$V2)
# data:  youngEresults$V2 and ancientEresults$V2
# t = -2.025, df = 216.26, p-value = 0.0441
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.0180823134 -0.0002444128
# sample estimates:
#   mean of x mean of y 
# 0.1941576 0.2033210 

t.test(youngNEresults$V2,ancientNEresults$V2)
# data:  youngNEresults$V2 and ancientNEresults$V2
# t = -3.1099, df = 567.15, p-value = 0.001966
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.014540400 -0.003283201
# sample estimates:
#   mean of x mean of y 
# 0.1941095 0.2030213 

mean(totalYoungresults$V2)
mean(totalAncientresults$V2)
t.test(totalYoungresults$V2,totalAncientresults$V2)
# data:  totalYoungresults$V2 and totalAncientresults$V2
# t = -3.7078, df = 784.19, p-value = 0.0002238
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.013725584 -0.004223075
# sample estimates:
#   mean of x mean of y 
# 0.1941226 0.2030969 
mean(youngDeNovoDensity$V2)
mean(totalAncientresults$V2)
t.test(youngDeNovoDensity$V2,totalAncientresults$V2)
# data:  youngDeNovoDensity$V2 and totalAncientresults$V2
# t = -0.35679, df = 48.363, p-value = 0.7228
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.02197405  0.01534965
# sample estimates:
#   mean of x mean of y 
# 0.1997847 0.2030969 
mean(youngDeNovoDensity$V2)
mean(ancientEresults$V2)
t.test(youngDeNovoDensity$V2,ancientEresults$V2)
# the ancient total and ancientE are practically the same as young de novo!?
# data:  youngDeNovoDensity$V2 and ancientEresults$V2
# t = -0.37882, df = 49.449, p-value = 0.7064
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.02229151  0.01521892
# sample estimates:
#   mean of x mean of y 
# 0.1997847 0.2033210 

t.test(youngDeNovoDensity$V2,dmel.tai)
# data:  youngDeNovoDensity$V2 and dmel.tai
# t = -0.074082, df = 48.297, p-value = 0.9413
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.01934363  0.01796865
# sample estimates:
#   mean of x mean of y 
# 0.1997847 0.2004722 

##################################################################
##################################################################
### Testing harmonic mean 

# reciprocalSum <- function(x){
#   r <- 1./(x+0.1)
# }
# 
# get.taih <- function(x,w) {
#   # x <- dmel.youngE.m
#   # w <- dmel.trna.ws
#   # x2 <- x
#   # x2[x2==0] <- NA
#   # x2
#   # w = log(w)              #calculate log of w
#   w = 1./w
#   x2 <- x
#   x2[x2==0] <- 0.1
#   x2 = 1./x2
#   # print(x+0.001) # matrix 95 x 60 codons
#   # x = x+0.001
#   # print(w) # vector of 60
#   n = apply(x2,1,'*',w)    #multiply each row of by the weights, 60 x 95
#   # -w?
#   # -w instead
#   # print(n)
#   # reciprocals?
#   # zeros?
#   n = t(n)                #transpose, same up to this point 
#   n
#   # print(n)
#   # n = exp(n)
#   n
#   # n = apply(n,1,function(x) 1/x)      #sum rows
#   n = apply(n,1,sum)      #sum rows
#   n
#   # print(n)
#   # print(n)
#   L = apply(x,1,sum)      #get each ORF length
#   L
#   tAIh = L/n    #get tai harmonic mean
#   # hist(log(tAIh))
#   return(tAIh)
#   # return((n[1,]))
# }
# # ???
# dmel.youngE.taih <- get.taih(dmel.youngE.m, dmel.trna.ws)
# dmel.youngE.taih
# # hist(dmel.youngE.taih)


# test with package data 
# 
# get.tai2 <- function(x,w) {
#   w[w==0] <- NA
#   na.omit(w)
#   # w = log(w)              #calculate log of w
#   n = apply(1/x,1,'*',1/w)    #multiply each row of by the weights
#   n = t(n)                #transpose
#   print(n)
#   n = apply(n,1,sum)      #sum rows
#   L = apply(x,1,sum)      #get each ORF length
#   tAI = L/n          #get tai
#   return(tAI)
# }
# 
# dmel.youngE.tai <- get.tai(dmel.youngE.m, dmel.trna.ws)
# dmel.youngE.taih <- get.tai2(dmel.youngE.m, dmel.trna.ws)
# dmel.youngE.tai
# dmel.youngE.taih

# get.taih <- function(x,w) {
#   x2 <- x
#   x[x==0] <- NA
#   na.omit(x)
#   # w = log(w)              #calculate log of w
#   n = apply(1/x,1,'*',1/w)    #multiply each row of by the weights
#   n = t(n)                #transpose
#   print(n)
#   n = apply(n,1,sum)      #sum rows
#   L = apply(x,1,sum)      #get each ORF length
#   tAI = L/n          #get tai
#   return(tAI)
# }
# 
# get.tai.harmonic <- function(x,w) {
#   # print(x) # 195 genes x 60 sense codon frequencies for 195 open reading frames.
#   # this is where 0 is introduced: some genes do not use all codons
#   # print(w) # a vector of length 60 of relative adaptiveness values for codons.
#   x_init <- x
#   x[x==0] <- 0.1 # assign a small value to replace 0
#   # x <- na.omit(x)
#   # print(min(x))
#   # doing this to avoid division by 0
#   # print(x)
#   # x <- 1/x # get reciprocal
#   # w = log(w)              #calculate log of w
#   # w <- 1/w
#   n = apply(x,1,'*',w)    #multiply each row of by the weights
#   # print(min(n))
#   n = t(n)                #transpose back to x dim; 195 by 60
#   # print(n) # same up to here
#     # I can't do reciprocal! Because of 0? 0 is caused by some codons not appearing in a gene
#   # print(n)
#   n <- 1/n 
#   n = apply(n,1,sum)      #sum rows # sum 1/wi_kg
#   L = apply(x_init,1,sum)      #get each ORF length
#   # print(L)
#   # tAI = exp(L/n)          #get tai
#   tAIh = as.double(L)/n
#   # tAIh.normalized <- tAIh - min(tAIh) / max(tAIh) - min(tAIh)
#   return(tAIh)
# }
# 
# # dmel.youngE.taih <- get.taih(dmel.youngE.m, dmel.trna.ws)
# # dmel.youngE.taih
# 
# dmel.youngE.tai.harmonic <- get.taih(dmel.youngE.m, dmel.trna.ws)
# dmel.youngE.tai.harmonic
# 
# dmel.youngE.tai2 <- get.taih(dmel.youngE.m, dmel.trna.ws)
# 
# #############


# get.tai<- function(x,w) {
#   # x is a matrix of codon counts: 195 genes by 60 sense codons
#   # w is a vector of length 60 of relative 
#     # adaptiveness values for sense codons
#   w = log(w)              #calculate log of w
#   # print(w)
#     # there is a w of 1, so the log is 0
#   # print(x)
#     # there are 0 in these codon counts 
#   n = apply(x,1,'*',w)    #multiply each row of by the weights
#   n = t(n)                #transpose
#   # print(n)
#   # this has row for each gene and cols for sense codons 
#   n = apply(n,1,sum)      #sum rows
#   print(n)
#   L = apply(x,1,sum)      #get each ORF length
#   tAI = exp(n/L)          #get tai
#   return(tAI)
# }

#######################################################
#######################################################
# attempting means
#######################################################
#######################################################


get.taih<- function(x,w) {
  # x is a matrix 195 genes by 60 sense codon frequencies for 195 open reading frames.
    # this is where 0 is introduced: some codon frequencies are 0
  # w is a vector of length 60 of relative adaptiveness values for sense codons
  x[x==0] <- 0.5 # probably doesnt make sesnse here to make GCN mean for those 0 ones lol
  n = apply(x,1,'*',w)    #multiply each row of by the weights
  n = t(n)                #transpose
  # print(n)
  # n[n==0] <- 0.1 # avoid division by 0
  n <- n^-1 # reciprocal of wi_kg
  n = apply(n,1,sum)      #sum reciprocals
  print(n)
  L = apply(x,1,sum)      #get each ORF length
  tAIh = (L/n)          #get tai
  return(tAIh)
}

get.taih2 <- function(x,w) {
  # x is a matrix 195 genes by 60 sense codon frequencies for 195 open reading frames.
  # this is where 0 is introduced: some codon frequencies are 0
  # w is a vector of length 60 of relative adaptiveness values for sense codons
  x[x==0] <- 0.1
  n = apply(x,1,'*',w)    #multiply each row of by the weights
  n = t(n)                #transpose
  # print(n)
  # n[n==0] <- 0.1 # avoid division by 0
  n <- n^-1 # reciprocal of wi_kg
  n = apply(n,1,sum)      #sum reciprocals
  print(n)
  L = apply(x,1,sum)      #get each ORF length
  tAIh = (L/n)          #get tai
  return(tAIh)
}

get.taih_notebook<- function(x,w) {
  # x is a matrix of codon counts: 195 genes by 60 sense codons
  # w is a vector of length 60 of relative adaptiveness values for sense codons
  
  w <- w^-1 # reciprocal of wi_kg
  
  n = apply(x,1,'*',w)    #multiply each row of by the weights
  n = t(n)                #transpose
  n = apply(n,1,sum)      #sum reciprocals
  L = apply(x,1,sum)      #get each ORF length
  tAI = (L/n)             #get tai
  
  return(tAI)
}

get.taia <- function(x,w) {
    w = log(w)
    n = apply(x, 1, "*", w)
    n = t(n)
    n = apply(n, 1, sum)
    L = apply(x, 1, sum)
    tAI = exp(n/L)
    return(tAI)
}

eco.trna <- scan("tai/misc/ecolik12.trna")
eco.ws <- get.ws(tRNA=eco.trna, sking=1)
eco.ws
eco.m <- matrix(scan("tai/misc/ecolik12.m"), ncol=61, byrow=TRUE)
eco.m <- eco.m[,-33]

eco.tai <- get.tai(eco.m, eco.ws)
eco.tai
eco.taih <- get.taih(eco.m, eco.ws)
eco.taih
eco.notebook <- get.taih_notebook(eco.m, eco.ws)
eco.notebook

hist(eco.tai)
ggplot(as.data.frame(eco.tai)) +
  geom_histogram(aes(eco.tai))

hist(eco.taih)
hist(eco.notebook)


youngE.tai <- get.tai(dmel.youngE.m, dmel.trna.ws)
youngE.tai
youngE.taih <- get.taih(dmel.youngE.m, dmel.trna.ws)
youngE.taih
youngE.notebook <- get.taih_notebook(dmel.youngE.m, dmel.trna.ws)
youngE.notebook

hist(youngE.tai)
hist(youngE.taih)
hist(youngE.notebook)


ancientE.tai <- get.tai(dmel.ancientE.m, dmel.trna.ws)
ancientE.tai
ancientE.taih <- get.taih(dmel.ancientE.m, dmel.trna.ws)
ancientE.taih
ancientE.notebook <- get.taih_notebook(dmel.ancientE.m, dmel.trna.ws)
ancientE.notebook

hist(ancientE.tai)
hist(ancientE.taih)
hist(ancientE.notebook)
# 
# 
# dmel.youngE.tai <- get.tai(dmel.youngE.m, dmel.trna.ws)
# dmel.youngE.tai
# dmel.youngE.taih <- get.taih(dmel.youngE.m, dmel.trna.ws)
# dmel.youngE.taih
# 
# eco.m
# eco.ws
# dmel.youngE.m
# dmel.trna.ws

# x <- dmel.youngE.m[1,]
# x
# dmel.youngE.tai.harmonic <- get.taih(dmel.youngE.m, dmel.trna.ws)
# dmel.youngE.tai.harmonic

# eco.m[eco.m==0] <- 0.1
# # sum(eco.m[1,]) # length
# # sum((eco.m[1,])^-1) # sum of reciprocal of wi_kg
# 
# (sum((eco.m[1,])^-1) / sum(eco.m[1,]))^-1


# save.image(f = 'tAI.Rdata')
# load('tAI.Rdata')


################################################
# branch-based analysis
################################################

dmel_branchData <- read.table('branch_based_analysis/dmel_ages',header = T)
names(dmel_branchData)[names(dmel_branchData) == 'g_id'] <- 'FBgn'
dmel_geneAgesIDsEssentiality <- read.table('Dmel_geneAgesIDsEssentiality.tsv', header = T)

dmel_geneAgesIDsEssentialityBranches <- full_join(dmel_geneAgesIDsEssentiality,dmel_branchData)
# combine them

# dmel_branchData_noNA <- dmel_geneAgesIDsEssentialityBranches[!is.na(dmel_geneAgesIDsEssentialityBranches$branch),]
# remove genes without branch data
# 16366 remain
# I should now repeat my correlaiton analysis
# knowing the branches I can take my prior results and add the branch information

dmel.total.final.tai <- read.table('dmel.total.final.tai',sep = '\t')
colnames(dmel.total.final.tai) <- c('FBgn','tAI')

dmel.total.final.w <- read.table("dmel.total.final.w", header=TRUE, na.strings = "*****")
# need to get only FBgn in first column
names(dmel.total.final.w)[names(dmel.total.final.w) == 'title'] <- 'FBgn'
# now i have the information on Nc and others that codonW produces

dmel.tAI.branches.essentiality.age <- full_join(dmel_geneAgesIDsEssentialityBranches,dmel.total.final.tai)
# now I have the data with age, essentiality, branch, tAI, and different IDs
# dmel.tAI.branches.essentiality.age %>% filter(tAI != 'NA') %>% distinct(FBgn,.keep_all = TRUE) %>% nrow()


dmel.BranchData <- full_join(dmel.tAI.branches.essentiality.age,dmel.total.final.w) %>% distinct(FBgn,.keep_all = TRUE)
# now I have all information 
# age, branch, IDs, essentiality, codonW info, tAI, chromosome

# dmel.BranchData %>% filter(tAI != 'NA') %>% select(FBgn,tAI) %>% write_tsv(.,'check_tAIvalues')
# dmel.BranchData%>% distinct(FBgn) %>% nrow()

# write_tsv(dmel.BranchData,'dmel.BranchData.temp')

dmel.BranchData %>% filter(tAI != 'NA') %>% nrow()
# number of genes with tAI values 
# 13968

dmel.BranchData %>% filter(branch != 'NA') %>% nrow()
# number of genes with branch placements
# 16360

dmel.BranchData %>% filter(branch != 'NA') %>% filter(tAI != 'NA') %>% nrow()
# number with branch and tAI
# 13851


formula <- y ~ x

branchCorrelationsR <- ggplot(dmel.BranchData %>% filter(branch != 'NA') %>% filter(tAI != 'NA'),aes(y=tAI, x=Nc)) +
  geom_point() +
  geom_smooth(aes(y=tAI, x=Nc), method = "lm",se = F,formula = formula) +
  stat_poly_eq(formula = formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(title = 'Correlation between tAI and Nc for 13,851 D. mel genes on branches -2 through 6') +
  xlab(label = 'Nc (effective number of codons)') +
  ylab(label = 'tAI (tRNA adaptation index)') +
  facet_wrap(. ~ branch, ncol =3, nrow = 3)
  # facet_grid(. ~ branch) 

branchCorrelations <- ggplot(dmel.BranchData %>% filter(branch != 'NA') %>% filter(tAI != 'NA'),aes(y=tAI, x=Nc)) +
  geom_point() +
  geom_smooth(aes(y=tAI, x=Nc), method = "lm",se = F,formula = formula) +
  # stat_poly_eq(formula = formula, 
  #              aes(label = paste(..rr.label.., sep = "~~~")), 
  #              parse = TRUE) +
  labs(title = 'Correlation between tAI and Nc for 13,851 D. mel genes on branches -2 through 6') +
  xlab(label = 'Nc (effective number of codons)') +
  ylab(label = 'tAI (tRNA adaptation index)') +
  facet_wrap(. ~ branch, ncol =3, nrow = 3)

# ggplot(dmel.BranchData %>% filter(branch != å'NA'),aes(y=tAI, x=Nc)) +
#   geom_point()+
#   labs(title = 'Correlation between tAI and Nc for D. mel genes on branches -2 through 6')+
#   xlab(label = 'Nc (effective number of codons)') +
#   ylab(label = 'tAI (tRNA adaptation index)') +
#   geom_smooth(aes(y=tAI, x=Nc), method = "lm",se = F,formula = my.formula) +
#   stat_poly_eq(formula = my.formula,
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                parse = TRUE) +
#   facet_wrap(. ~ branch, ncol =3, nrow = 3)

ggsave(branchCorrelationsR, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/branchCorrelationsR.png')
ggsave(branchCorrelations, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/branchCorrelations.png')

# # ggplot(dmel.BranchData %>% filter(branch != 'NA'))
# ggplot(dmel.BranchData %>% filter(branch != 'NA'), aes(x=tAI)) + 
#   # geom_vline(aes(xintercept=mean(tAI)),color="red", linetype="dashed", size=1) +
#   # geom_histogram(aes(y=..density..),color="black", fill="forestgreen") +
#   geom_histogram(aes(y=..density..),color="black", fill="forestgreen", binwidth = 0.05) +
#   # geom_histogram(color="black", fill="forestgreen") +
#   geom_density(alpha=.2, fill="red") +
#   facet_wrap(. ~ branch, ncol = 3, nrow = 3)
#   # facet_grid(. ~ branch)


mean_tAI_perBranch <- dmel.BranchData %>% filter(branch != 'NA') %>%
  group_by(branch) %>%
  summarize(mean.tai = mean(na.omit(tAI)),median.tai = median(na.omit(tAI))) 

branchtAIdensity <- ggplot(dmel.BranchData %>% filter(branch != 'NA') %>% filter(tAI != 'NA')) +
  geom_density(aes(tAI)) +
  geom_vline(data = mean_tAI_perBranch, mapping = aes(xintercept = mean.tai),
             color="red", linetype="dashed", size=0.5) +
  labs(title = 'tRNA adaptation index values for 13,851 D. mel genes on branches -2 through 6')+
  xlab(label = 'tAI (tRNA adaptation index) value') +
  ylab(label = 'count') +
  geom_text(data = mean_tAI_perBranch, 
            aes(color = 'red',
                label = paste0('Mean: ',signif(mean.tai,digits = 6)), 
                y = 5.5, x = .4)) +
  # geom_text(data = mean_tAI_perBranch, 
  #           aes(color = 'red',
  #               label = paste0('Median: ',signif(median.tai,digits = 6)), 
  #               y = 4.5, x = .4)) +
  theme(legend.position="none") +
  facet_wrap(. ~ branch, ncol =3, nrow = 3)
ggsave(branchtAIdensity, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/branchTAidensity.png')

###################################
##### origination mechanism-based analyses
###################################

dmel.originationMechnanism <- read.table('origination_mech_yongzhang.txt',sep = '\t')
colnames(dmel.originationMechnanism) <- c('FBgn','parent_id','origMech','movementType','alignmentInfo')
# Add originaiton mechansims based on Yong Zhang's data
# thus I have genes separated based on three mechanisms: orphan (candidate de novo), RNA-dup, DNA-dup
# /Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/data/origniation mechanism data, yong zhang/Dmel-updated original mechnism-yongzhang.xlsx

dmel.BranchData.origMech <- left_join(dmel.BranchData,dmel.originationMechnanism)
dmel.fullData <- dmel.BranchData.origMech %>% 
  mutate(melanogasterSpecific = case_when(branch == 4 | branch == 5 | branch == 6 ~ TRUE, TRUE ~ NA)) %>% 
  mutate(sophophoraSpecific = case_when(branch == 1 | branch == 2 | branch == 3 ~ TRUE, TRUE ~ NA)) %>%
  mutate(drosophilaSpecific = case_when(branch == -1 | branch == 0 ~ TRUE, TRUE ~ NA)) 
# %>% distinct(FBgn,.keep_all = TRUE)
  # %>% 
  # mutate(group = case_when(melanogasterSpecific == TRUE ~ 'm' | sophophoraSpecific == TRUE ~ 's' | drosophilaSpecific == TRUE ~ 'd' , TRUE ~ NA))
# %>% select(branch,melanogasterSpecific)
# conditionally mutate!

dmel.fullData %>% filter(origMech != 'NA') %>% nrow()
# 13766 have origMech info 

dmel.fullData %>% filter(tAI != 'NA') %>% nrow()
# number of genes with tAI values 
# 13968

dmel.fullData %>% filter(tAI != 'NA') %>% filter(origMech != 'NA' ) %>% nrow()
# number with tAI and mech
# 13701
# this makes sense as it is lower than those with just mech or just tAI


dmel.fullData %>% filter(branch != 'NA' ) %>% nrow()
# number of genes with branch placements
# 16360

dmel.fullData %>% filter(branch != 'NA') %>% filter(tAI != 'NA') %>% nrow()
# number with branch and tAI
# 13851

mean_tAI_perMech <- dmel.fullData %>% filter(origMech != 'NA') %>%
  group_by(origMech) %>%
  summarize(mean.tai = mean(na.omit(tAI)),median.tai = median(na.omit(tAI))) 

mech_names <- c(
  `A` = "Lineage-specific (orphan)",
  `D` = "DNA-based duplication",
  `Dl` = "DNA-based duplication-like",
  `R` = "RNA-based duplication",
  `Rl` = "RNA-based duplication-like"
)

correlationOrigMech <- ggplot(dmel.fullData %>% filter(origMech != 'NA') %>% filter(tAI != 'NA'),aes(y=tAI, x=Nc)) +
  geom_point() +
  geom_smooth(aes(y=tAI, x=Nc), method = "lm",se = F,formula = formula) +
  # stat_poly_eq(formula = formula, 
  #              aes(label = paste(..rr.label.., sep = "~~~")), 
  #              parse = TRUE) +
  labs(title = 'Correlation between tAI and Nc for 13,701 D. mel genes of different origination mechanisms') +
  xlab(label = 'Nc (effective number of codons)') +
  ylab(label = 'tAI (tRNA adaptation index)') +
  facet_wrap(. ~ origMech, labeller = as_labeller(mech_names))
ggsave(correlationOrigMech, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/correlationOrigMech.png')


dmel.fullData %>% filter(origMech != 'NA') %>% filter(tAI != 'NA') %>% nrow()
# 13701 have mech and tAIs


densityOrigMech <- ggplot(dmel.fullData %>% filter(origMech != 'NA') %>% filter(tAI != 'NA')) +
  geom_density(aes(tAI)) +
  geom_vline(data = mean_tAI_perMech, mapping = aes(xintercept = mean.tai),
             color="red", linetype="dashed", size=0.5) +
  labs(title = 'tRNA adaptation index values for 13,701 D. mel genes of different origination mechanisms')+
  xlab(label = 'tAI (tRNA adaptation index) value') +
  ylab(label = 'count') +
  geom_text(data = mean_tAI_perMech, 
            aes(color = 'red',
                label = paste0('Mean: ',signif(mean.tai,digits = 6)), 
                y = 5.5, x = .4)) +
  # geom_text(data = mean_tAI_perMech, 
  #           aes(color = 'red',
  #               label = paste0('Median: ',signif(median.tai,digits = 6)), 
  #               y = 4.5, x = .4)) +
  theme(legend.position="none") +
  facet_wrap(. ~ origMech, labeller = as_labeller(mech_names))
ggsave(densityOrigMech, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/densityOrigMech.png')




# start working from here 
# dmel.fullData %>% filter(sophophoraSpecific != 'NA') %>% filter(tAI != 'NA') %>% nrow()

dmel.fullData <- dmel.fullData %>% 
  mutate(group = case_when(sophophoraSpecific == TRUE ~ 'Sophophora subgenus-specific (branches 1,2,3)',
                           melanogasterSpecific == TRUE ~ 'Melanogaster subgroup-specific (branches 4,5,6)',
                           drosophilaSpecific == TRUE ~ 'Drosophila subgenus-specific (branches -1,0)'))

 dmel.fullData %>% filter(group != 'NA') %>% filter(tAI != 'NA') %>% nrow()
# 4373 have a group and tAI

dmel.fullData %>% filter(group != 'NA') %>% nrow()
# 6478 in a group; this does not include branch -2 as it is outgroup

dmel.fullData %>% filter(branch != -2) %>% nrow()
# also 6478 not on branch -2; expected 

dmel.BranchData %>% filter(branch == -2) %>% nrow()
# 9882  are in branch -2 which is not considered for grouping

dmel.BranchData %>% filter(branch == -2) %>% filter(tAI != 'NA') %>% nrow()
# 9478 have tAI and are in branch -2 which is not one of 3 groups of interest

# recall:
dmel.BranchData %>% filter(branch != 'NA') %>% nrow()
# number of genes with branch placements
# my groups should equal this
# 16360

dmel.BranchData %>% filter(branch != 'NA') %>% filter(tAI != 'NA') %>% nrow()
# number with branch and tAI
# 13851

dmel.fullData %>% filter(group != 'NA') %>% filter(tAI != 'NA') %>% nrow() + dmel.BranchData %>% filter(branch == -2) %>% filter(tAI != 'NA') %>% nrow()
# 4373 + 9478 = 13851 
# as expected
# these are those with tAI in a group of interest 
# so the -2 branch is a large portion of those we are not intrested in

# dmel.fullData %>% filter()


mean_tAI_perGroup <- dmel.fullData %>% filter(group != 'NA') %>%
  group_by(group) %>%
  summarize(mean.tai = mean(na.omit(tAI)),median.tai = median(na.omit(tAI))) 

densityGroup <- ggplot(dmel.fullData %>% filter(group != 'NA') %>% filter(tAI != 'NA')) +
  geom_density(aes(tAI)) +
  geom_vline(data = mean_tAI_perGroup, mapping = aes(xintercept = mean.tai),
             color="red", linetype="dashed", size=0.5) +
  labs(title = 'tRNA adaptation index values for 4,373 group-specific genes')+
  xlab(label = 'tAI (tRNA adaptation index) value') +
  ylab(label = 'count') +
  geom_text(data = mean_tAI_perGroup,
            aes(color = 'red',
                label = paste0('Mean: ',signif(mean.tai,digits = 6)),
                y = 5.5, x = .35)) +
  # geom_text(data = mean_tAI_perMech, 
  #           aes(color = 'red',
  #               label = paste0('Median: ',signif(median.tai,digits = 6)), 
  #               y = 4.5, x = .4)) +
  theme(legend.position="none")+
  facet_wrap(. ~ group)
ggsave(densityGroup, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/densityGroup.png')

correlationGroups <- ggplot(dmel.fullData %>% filter(group != 'NA') %>% filter(tAI != 'NA'),aes(y=tAI, x=Nc)) +
  geom_point() +
  geom_smooth(aes(y=tAI, x=Nc), method = "lm",se = F,formula = formula) +
  # stat_poly_eq(formula = formula, 
  #              aes(label = paste(..rr.label.., sep = "~~~")), 
  #              parse = TRUE) +
  labs(title = 'Correlation between tAI and Nc for 4,373 group-specific genes') +
  xlab(label = 'Nc (effective number of codons)') +
  ylab(label = 'tAI (tRNA adaptation index)') +
  facet_wrap(. ~ group)
ggsave(correlationGroups, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/correlationGroups.png')


dmel.fullData %>% filter(group == 'Drosophila subgenus-specific (branches -1,0)') %>% filter(tAI != 'NA') %>% nrow()
dmel.fullData %>% filter(group == 'Sophophora subgenus-specific (branches 1,2,3)') %>% filter(tAI != 'NA') %>% nrow()
dmel.fullData %>% filter(group == 'Melanogaster subgroup-specific (branches 4,5,6)') %>% filter(tAI != 'NA') %>% nrow()

############################################
### T-tests for groups and branch-based analyses 
############################################
# t.test(totalYoungresults$V2,totalAncientresults$V2)
# t.test(dmel.fullData %>% filter(group == 'Melanogaster subgroup-specific (branches 4,5,6)') %>% select(tAI),dmel.fullData %>% filter(group != 'NA') %>% select(tAI))
# t = -2.5081, df = 360.06, p-value = 0.01258
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.017584811 -0.002128193
# sample estimates:
#   mean of x mean of y 
# 0.1828398 0.1926963 
# comparing average tAI for all three group tAIs with just melanogaster-specific
# there is a significant difference between newest branches tAI mean and total tAI mean
# melanogaster group compared to all groups

t.test(dmel.fullData %>% filter(group == 'Melanogaster subgroup-specific (branches 4,5,6)') %>% select(tAI),
       dmel.fullData %>% filter(group != 'Melanogaster subgroup-specific (branches 4,5,6)') %>% select(tAI))
# t = -2.7012, df = 362.69, p-value = 0.007234
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.018377051 -0.002892371
# sample estimates:
#   mean of x mean of y 
# 0.1828398 0.1934745 
# melanogaster group compared to those with groups excluding melanogaster 

# t.test(dmel.fullData %>% filter(group == 'Melanogaster subgroup-specific (branches 4,5,6)') %>% select(tAI),dmel.fullData %>% select(tAI))
# t = -4.5834, df = 330.75, p-value = 6.494e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.02520013 -0.01006472
# sample estimates:
#   mean of x mean of y 
# 0.1828398 0.2004722


# t.test(dmel.fullData %>% filter(group == 'Melanogaster subgroup-specific (branches 4,5,6)') %>% select(tAI),
       # dmel.fullData %>% filter(group != 'Melanogaster subgroup-specific (branches 4,5,6)') %>% select(tAI))
# t = -2.7012, df = 362.69, p-value = 0.007234
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.018377051 -0.002892371
# sample estimates:
#   mean of x mean of y 
# 0.1828398 0.1934745 

t.test(dmel.fullData %>% filter(origMech == 'A') %>% select(tAI),
       dmel.fullData %>% filter(origMech == 'D' | origMech == 'R') %>% select(tAI))
# t = -2.3373, df = 4069.5, p-value = 0.01947
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.0085157226 -0.0007465803
# sample estimates:
#   mean of x mean of y 
# 0.1996253 0.2042564 
# not including DNA-duplicaiton-like or RNA-duplication-like
# this shows a significant difference between orphan gene tAI and RNA or DNA-based duplications



# dmel.fullData %>% filter(group != 'NA') %>% filter(tAI!='NA') %>% nrow()
# dmel.fullData %>% filter(tAI != 'NA') %>% mean()
# dmel.fullData %>% select(tAI) %>% filter(tAI !='NA') 

###########
### modifying scatter plot with color
###########
ggplot(dmel.fullData %>% filter(origMech != 'NA') %>% filter(tAI != 'NA'),aes(y=tAI, x=Nc,color=age)) +
  geom_point(alpha=0.7,position=position_jitter(w = 0.1,h = 0)) +
  geom_smooth(aes(y=tAI, x=Nc), method = "lm",se = F,formula = formula) +
  scale_color_manual(values = lacroix_palette('Lemon', type='discrete')) +
  # stat_poly_eq(formula = formula, 
  #              aes(label = paste(..rr.label.., sep = "~~~")), 
  #              parse = TRUE) +
  labs(title = 'Correlation between tAI and Nc for 13,701 D. mel genes of different origination mechanisms') +
  xlab(label = 'Nc (effective number of codons)') +
  ylab(label = 'tAI (tRNA adaptation index)') +
  facet_wrap(. ~ origMech, labeller = as_labeller(mech_names))

ggplot(dmel.fullData %>% filter(group != 'NA') %>% filter(tAI != 'NA'),aes(y=tAI, x=Nc, color = age)) +
  geom_point(alpha=0.8,position=position_jitter(w = 0.1,h = 0)) +
  geom_smooth(aes(y=tAI, x=Nc), method = "lm",se = F,formula = formula) +
  scale_color_manual(values = wes_palette("Chevalier1",3)) +
  # stat_poly_eq(formula = formula, 
  #              aes(label = paste(..rr.label.., sep = "~~~")), 
  #              parse = TRUE) +
  labs(title = 'Correlation between tAI and Nc for 4,373 group-specific genes') +
  xlab(label = 'Nc (effective number of codons)') +
  ylab(label = 'tAI (tRNA adaptation index)') +
  facet_wrap(. ~ group)

###################################
# Codon bias analysis for longest CDS
###################################
dmel.cds.codonUsage <- read.table('codonUsage/dmel.cds.codonUsage',header = T)
# temp <- left_join(dmel.fullData,dmel.cds.codonUsage)
# not sure if I need this yet 

dmel.tRNA.gcn.0 <- read.table('Dmel.tRNA.gcn.0',header=F)
# tRNA genes with no copy number 

# dmel.cds.codonUsage %>% filter(TTA != 0)
# should do some investigations like this 

colSums(dmel.cds.codonUsage[,-1])
# there are no codons with 0 count??

dmel.denovo.yongzhang.codonUsage <- dmel.fullData %>% filter(age=='young') %>% filter(origMech=='A') %>% select(FBgn) %>% left_join(.,dmel.cds.codonUsage)
# number of Yong Zhang's genes annotates as young de novo

colSums(dmel.denovo.yongzhang.codonUsage[,-1])
# no de novo genes have 0 codons?

dmel.tRNA.gcn.0.wider <- dmel.tRNA.gcn.0 %>% pivot_wider(names_from = V1, values_from = V2)

dmel.tRNA.gcn.0.wider %>% ncol()
#  confirm jsut 20 codons with 0 HCN

dmel.denovo.yongzhang.codonUsage.tRNA.gcn.0 <- dmel.denovo.yongzhang.codonUsage %>% select(FBgn,colnames(dmel.tRNA.gcn.0.wider)) 
# de novo genes and the 20 codons that have 0 corresponding tRNA GCN in melanogaster 

dmel.fullData.tRNA.gcn.0 <- dmel.cds.codonUsage %>% select(FBgn,colnames(dmel.tRNA.gcn.0.wider)) 
# this all protein coding genes and the counts of the 20 codons with 0 GCN

############################################################################################
############################################################################################
####### codon bias based on Marty's chapter 
############################################################################################
############################################################################################
dmel.average.Nc <- dmel.fullData %>% 
  group_by(age) %>% 
  filter(Nc != 'NA') %>% 
  summarise_at(vars(Nc), list(average_Nc = mean))

dmel.fullData %>% filter(Nc != 'NA') %>% select(GENEID,Nc,age) %>% filter(GENEID != 'NA') %>% head(.,200) %>% 
  ggplot(aes(GENEID,Nc,group=age,color=age))+
  geom_line()+
  geom_point()+
  # geom_smooth(se=F)+
  theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
  ggtitle("Comparison of effective num. of codons used, Nc, for genes of D. mel") +
  scale_color_manual(values = lacroix_palette('Pamplemousse', type='discrete')) +
  theme(plot.title = element_text(hjust = 0.5)) 

############################################################################################
############################################################################################
####### looking at Li Zhang's de novo data
############################################################################################
############################################################################################
# 
# liZhangDeNovo <- read.table('liZhangDeNovo')
# colnames(liZhangDeNovo) <- c('FBgn','tAI_lz')
# 
# # dmel.fullData %>% filter(.)
# # 
# # filter(dmel.fullData, FBgn %in% liZhangDeNovo$FBgn) %>% select(FBgn,tAI) 
# # # to look up and extract the 31 de novo genes from my data 
# 
# liZhangDeNovo.both.tAI <- left_join(liZhangDeNovo,dmel.fullData) %>% select(FBgn,tAI_lz,tAI)
# 
# liZhangDeNovo.both.tAI.density <- ggplot(liZhangDeNovo.both.tAI)+
#   geom_density(aes(tAI),color='red')+
#   geom_density(aes(tAI_lz),color='blue')+
#   # labs(title = 'tRNA adaptation index values computed two ways for 31 D. mel de novo genes (annotated by Li Zhang)')+
#   labs(title = expression(paste("tRNA adaptation index values computed two ways for 31 ",italic("D. mel"),' de novo genes (annotated by Li Zhang)')))+
#   xlab(label = 'tAI (tRNA adaptation index) value')
# 
# ggsave(liZhangDeNovo.both.tAI.density, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/liZhangDeNovo.both.tAI.density.png')
# 
# mean(liZhangDeNovo.both.tAI$tAI_lz)
# mean(liZhangDeNovo.both.tAI$tAI)

############################################################################################
############################################################################################
####### Making new age column based on Chuan's branch placement 
############################################################################################
############################################################################################
## adding new column based on Chuan data
## this is the age based on Chuan's branch placement
## I will not use age based on Yong Zhang's data, but I will use the lethality observations.
## Thus, I will keep the essentiality annotation for now.
dmel.fullData_chuan <- dmel.fullData %>% 
  mutate(age_Chuan = case_when(branch >= 3 ~ 'young',
                               branch < 3 ~ 'ancient')
         )
  
dmel.fullData_chuan %>% filter(age_Chuan == 'young') %>% nrow()
# 423
# 504 if including branch 3 as young 
dmel.fullData_chuan %>% filter(age_Chuan == 'ancient') %>% nrow()
# 15937

dmel.fullData_chuan %>% filter(age == 'young') %>% nrow()
# 701
dmel.fullData_chuan %>% filter(age == 'ancient') %>% nrow()
# 10595

# ge <- read.csv(header = F, file = '../3_tRNA_evolution/pacbio/tRNAscan-SE_output/GCA_003285725.2_SlebRS2.trnascan.counts_isotypes_anticodons_numGenes', sep = '\t')

save.image(f = 'tAI.Rdata')
load('tAI.Rdata')
