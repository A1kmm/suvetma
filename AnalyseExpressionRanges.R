geneProfile <- read.csv('/home/andrew/Documents/hs_gene_level_profiles.csv');
gpcs <- colSums(geneProfile[, 2:19])
barplot(gpcs, names.arg=signif(2^(-2.5:14.5), 3), xlab='Expression', ylab='Frequency')
gpPerUnit <- gpcs[2:17]/2^(-2:13)
barplot(gpPerUnit, names.arg=signif(2^(-1.5:13.5), 3), xlab='Expression', ylab='Width Frequency')
