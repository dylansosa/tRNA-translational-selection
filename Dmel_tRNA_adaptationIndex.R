# devtools::install_github("dylansosa/tai")
require("tAI")
library('tidyverse')

setwd('/Users/dylansosa/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Yakuba  Long Lab/thesis/2_tRNA_adaptation')
# setwd('/Users/dylansosa/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/YakubaUChicago-1/thesis/2_tRNA_adaptation')

dmel.trna <- scan("Dmel.tRNA")
# load tRNA gene copy numbers 
# I computed these data using all CDS of D. mel
# need to truncate it to old/ancient later

dmel.trna.ws <- get.ws(tRNA=dmel.trna, sking = 0)
# compute relative adaptiveness values for each codon in the D. mel genes

dmel.youngE.m <- matrix(scan("dmel.longestCDS.youngE.m"), ncol=61, byrow=TRUE)
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
write.table(dmel.youngE.tai, file = 'dmel.youngE.tai.txt',sep = '\t',quote = F)

hist(dmel.youngE.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.

# df <- read.table("dmel.youngE.cds.w", header=TRUE)
# df <- read.table("dmel.longestCDS.youngE.w", header=TRUE)
df.youngE <- read.table("dmel.longestCDS.youngE.w", header=TRUE, na.strings = "*****")
# cat dmel.longestCDS.youngE.fa | egrep '>' | sed s'/>//g' > dmel.youngE.id.temp
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

cor(dmel.youngE.tai, df.youngE$Nc, use="p")
# [1] -0.3810451
# there is no Nc!
# wow, USING THE na.strings SOLVES THIS
dmel.youngE.s <- get.s(dmel.youngE.tai, df.youngE$Nc, df$GC3s)
dmel.youngE.s
# [1] 0.4817288

# Formally, we can calculate the correlation between tAI, 
# and the corrected Nc, f(GC3s) – Nc, we call this correlation S, 
# because it reflects the intensity of translational selection 
# acting on our sample of genes

#######################################################
#######################################################
# below is to test on different ages
dmel.youngNE.m <- matrix(scan("dmel.longestCDS.youngNE.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script: 
#   a matrix of codon frequencies per ORF.

dmel.youngNE.m <- dmel.youngNE.m[,-33]
# ignore methionine codons

dmel.youngNE.tai <- get.tai(dmel.youngNE.m, dmel.trna.ws)
# calculate tRNA adaptation index 
write.table(dmel.youngNE.tai, file = 'dmel.youngNE.tai.txt',sep = '\t',quote = F)

hist(dmel.youngNE.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.

# df.youngNE <- read.table("dmel.ancientE.cds.w", header=TRUE)
# df.youngNE <- read.table("dmel.cds.young.nonessential.w", header=TRUE)
df.youngNE <- read.table("dmel.longestCDS.youngNE.w", header=TRUE, na.strings = "*****")
# this is where errors are,   line 368 did not have 5 elements
# same error when I try to run this script on all Dmel genes
# it just happenes to be that the youngE genes don't have the problem of not computing each col 

plot(dmel.youngNE.tai ~ df.youngNE$Nc)
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons), 
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
cor(dmel.youngNE.tai, df.youngNE$Nc, use="p")
# [1] -0.2091027
dmel.youngNE.s <- get.s(dmel.youngNE.tai, df.youngNE$Nc, df.youngNE$GC3s)
dmel.youngNE.s
# [1] 0.3199346

# Formally, we can calculate the correlation between tAI, 
# and the corrected Nc, f(GC3s) – Nc, we call this correlation S, 
# because it reflects the intensity of translational selection 
# acting on our sample of genes:

########################################################
########################################################
dmel.ancientE.m <- matrix(scan("dmel.longestCDS.ancientE.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script: 
#   a matrix of codon frequencies per ORF.

dmel.ancientE.m <- dmel.ancientE.m[,-33]
# ignore methionine codons

dmel.ancientE.tai <- get.tai(dmel.ancientE.m, dmel.trna.ws)
# calculate tRNA adaptation index 
write.table(dmel.ancientE.tai, file = 'dmel.ancientE.tai.txt',sep = '\t',quote = F)


hist(dmel.ancientE.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.

# df.ancientE <- read.table("dmel.ancientE.cds.w", header=TRUE)
# df.ancientE <- read.table("dmel.cds.young.nonessential.w", header=TRUE)
df.ancientE <- read.table("dmel.longestCDS.ancientE.w", header=TRUE, na.strings = "*****")
plot(dmel.ancientE.tai ~ df.ancientE$Nc)
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons), 
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
cor(dmel.ancientE.tai, df.ancientE$Nc, use="p")
# [1] -0.240699
dmel.ancientE.s <- get.s(dmel.ancientE.tai, df.ancientE$Nc, df.ancientE$GC3s)
dmel.ancientE.s
# [1] 0.3053252
# Formally, we can calculate the correlation between tAI, 
# and the corrected Nc, f(GC3s) – Nc, we call this correlation S, 
# because it reflects the intensity of translational selection 
# acting on our sample of genes:

########################################################
########################################################

dmel.ancientNE.m <- matrix(scan("dmel.longestCDS.ancientNE.m"), ncol=61, byrow=TRUE)
# read output of codonM
# .m contains the output of the codonM script: 
#   a matrix of codon frequencies per ORF.

dmel.ancientNE.m <- dmel.ancientNE.m[,-33]
# ignore methionine codons

dmel.ancientNE.tai <- get.tai(dmel.ancientNE.m, dmel.trna.ws)
# calculate tRNA adaptation index 
write.table(dmel.ancientNE.tai, file = 'dmel.ancientNE.tai.txt',sep = '\t',quote = F)


hist(dmel.ancientNE.tai)
# Highly expressed genes present high tAI values (> 0.4), 
# which means that their codon usage resembles the genomic 
# structure of tRNA genes.

# df.ancientNE <- read.table("dmel.ancientE.cds.w", header=TRUE)
# df.ancientNE <- read.table("dmel.cds.young.nonessential.w", header=TRUE)
df.ancientNE <- read.table("dmel.longestCDS.ancientNE.w", header=TRUE, na.strings = "*****")
plot(dmel.ancientNE.tai ~ df.ancientNE$Nc)
# plot of relationship between tAI and effective number of codons
# genes with very low Nc values (highly biased toward certain codons), 
# correlate with high tAI values (highly co-adapted to the tRNA gene pool).
cor(dmel.ancientNE.tai, df.ancientNE$Nc, use="p")
# [1] -0.249064

dmel.ancientNE.s <- get.s(dmel.ancientNE.tai, df.ancientNE$Nc, df.ancientNE$GC3s)
dmel.ancientNE.s
# [1] 0.3483473

# Formally, we can calculate the correlation between tAI, 
# and the corrected Nc, f(GC3s) – Nc, we call this correlation S, 
# because it reflects the intensity of translational selection 
# acting on our sample of genes:



library(tidyverse)
library(gridExtra)

youngEresults <- read.csv('dmel.youngE.tai',sep = '\t',header = F)
ancientEresults <- read.csv('dmel.ancientE.tai',sep = '\t',header = F)
youngNEresults <-read.csv('dmel.youngNE.tai',sep = '\t',header = F)
ancientNEresults <- read.csv('dmel.ancientNE.tai',sep = '\t',header = F)

x <- ggplot(youngEresults) +
  geom_density(aes(V2)) +
  labs(title = 'young essential gene tAI')
  
y <- ggplot(ancientEresults) +
  geom_density(aes(V2)) +
  labs(title = 'ancient essential gene tAI')

a <- ggplot(youngNEresults) +
  geom_density(aes(V2)) +
  labs(title = 'young non-essential gene tAI')

b <- ggplot(ancientNEresults) +
  geom_density(aes(V2)) +
  labs(title = 'ancient non-essential gene tAI')


# mean(as.vector(ancientEresults))
# class(ancientEresults)
# mean(ancientNEresults)

groupDensity <- grid.arrange(x,y,a,b,ncol=2)
groupDensity
ggsave(groupDensity, filename = '/Users/dylansosa/Documents/UChicago/Long Lab/Projects/Thesis Projects/prediciton_essential_gene_evolution/code/tRNA/figures/roupDensity.png')

##################################################################
##################################################################
### Testing harmonic mean 

reciprocalSum <- function(x){
  r <- 1./(x+0.1)
}

get.taih <- function(x,w) {
# x <- dmel.youngE.m
# w <- dmel.trna.ws
# x2 <- x
# x2[x2==0] <- NA
# x2
  # w = log(w)              #calculate log of w
  w = 1./w
  x2[x2==0] <- 0.1
  x2 = 1./x2
  # print(x+0.001) # matrix 95 x 60 codons
  # x = x+0.001
  # print(w) # vector of 60
  n = apply(x2,1,'*',w)    #multiply each row of by the weights, 60 x 95
  # -w?
  n = apply(x2,1,'-',10)
  # -w instead
  # print(n)
  
  # reciprocals?
  # zeros?
  n = t(n)                #transpose, same up to this point 
  n
  # print(n)
  # n = exp(n)
  n
  # n = apply(n,1,function(x) 1/x)      #sum rows
  n = apply(n,1,sum)      #sum rows
  n
  # print(n)
  # print(n)
  L = apply(x,1,sum)      #get each ORF length
  L
  tAIh = L/n    #get tai harmonic mean
  # hist(log(tAIh))
  return(tAIh)
  # return((n[1,]))
}
# ???
dmel.youngE.taih <- get.taih(dmel.youngE.m, dmel.trna.ws)
dmel.youngE.taih
# hist(dmel.youngE.taih)


# test with package data 

get.tai2 <- function(x,w) {
  w = log(w)              #calculate log of w
  n = apply(x,1,'*',w)    #multiply each row of by the weights
  n = t(n)                #transpose
  print(n)
  n = apply(n,1,sum)      #sum rows
  L = apply(x,1,sum)      #get each ORF length
  tAI = exp(n/L)          #get tai
  return(tAI)
}

dmel.youngE.tai <- get.tai2(dmel.youngE.m, dmel.trna.ws)
dmel.youngE.tai


#############
eco.trna <- scan("tai/misc/ecolik12.trna")
eco.ws <- get.ws(tRNA=eco.trna, sking=1)
eco.m <- matrix(scan("tai/misc/ecolik12.m"), ncol=61, byrow=TRUE)
eco.m <- eco.m[,-33]
eco.tai <- get.tai(eco.m, eco.ws)
hist(eco.tai)
eco.taih <- get.taih(eco.m, eco.ws)

