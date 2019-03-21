
# Import Libraries --------------------------------------------------------

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

library(dplyr) # for %>%
library(reshape2) # for melt
library(ggsci)

# Construct metadata dataframe --------------------------------------------

## get smaple name
samples.out <- rownames(seqtab.nochim)
## get patient id
patid<- substr(samples.out, 1, 5)
## get group information 
group <- substr(patid, 1, 3) 
## get patid
patid <- gsub("[[:digit:]]", "", patid)
##construc samplename and group informations
grouping <- data.frame(sample=samples.out, patid = patid, group=group, stringsAsFactors = FALSE)
## rename grouping into full name 
grouping[which(grouping$group == "HAE"), "group"] <- "HAEC"


## import original metadata
met <- read.csv("Ranalysis.csv", stringsAsFactors = FALSE)
### define mode of gage los
met$gage <- met$gage %>% as.numeric()
met$los <- met$los %>% as.numeric()
### rm 错了的 Sample
met <- met[, -1]

##construct grouping and original metadat 
metadata <- merge(grouping, met, by = "patid", sort = FALSE, all = TRUE)
### remove one "group" column
metadata <- metadata[, -3]
colnames(metadata)[3] <- "group"

## define rownames of the metadata 
rownames(metadata) <- metadata$sample



# construct a phyloseq object ---------------------------------------------

#directl from the dada2 output
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps


##!!!!!!!!!!!!!!!!!!!!
# We are now ready to use phyloseq! ---------------------------------------
##!!!!!!!!!!!!!!!!!!!!

# Visualize alpha-diversity:
plot_richness(ps, x="group", measures=c("Shannon", "Simpson"), color="group")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="group", title="Bray NMDS")

# barplot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="patid", fill="Genus") + facet_wrap(~group, scales="free_x")




# using the .plus output as input for taxonomy barplot input ----------------------

#import data
tax <- read.delim("dada2_counts.taxon.species.txt")
tax <- tax[, -1] #remove seq
tax <- tax[, -65] #remove species
#subset genus "column64"
genus <- cbind(tax[ , 64], tax[ , 1:58])
#define genus colname
colnames(genus)[1] <- "genus"

#calculate proportion ie relative abundance
propgenus <- prop.table(genus[ , 2:59] %>% as.matrix(), 2)
propgenus <- propgenus %>% as.data.frame()
propgenus <- cbind(genus$genus, propgenus)
## convert propgenus into dataframe to summarise
colnames(propgenus)[1] <- "genus"
# sum up same genus relative abundance
propgenus <- propgenus %>% group_by(genus) %>% summarise_all(funs(sum))

# which row contains numbers over 0.05
propgenus0.05 <- as.numeric()
for (i in 1:nrow(propgenus)) {
  propgenus0.05[i] <- which(propgenus[i, 2:59] >= 0.15) %>% length()
}
## cutoff 0.05
propgenus <- propgenus[propgenus0.05 != 0,]

#remove NA genus row
propgenus <- propgenus[-22, ]


# construct genus to merge with metadata
#define genus
genusname <- propgenus$genus %>% as.character()
#transpose
propgenus <- t(propgenus) %>% as.data.frame()
colnames(propgenus) <- genusname
#remove firstrow(历史遗留自由而无用的genusname)
propgenus <- propgenus[-1, ]
#add an sample column to facilitate with merging afterwards
propgenus$sample <- rownames(propgenus)


# merge metadata with propgenus
meta_genus <- merge(metadata, propgenus, by = "sample")



# plot by genus -----------------------------------------------------------

pat_genus <- meta_genus[, c(2, 9:29)]
pat_genus <- sapply(pat_genus, as.character) %>% as.data.frame(., stringsAsFactors = FALSE)
pat_genus[, 2:22] <- sapply(pat_genus[, 2:22], as.numeric)

pat_genus <- pat_genus %>% group_by(patid) %>% summarise_all(funs(mean))
patid <- pat_genus$patid

pat_genus <- t(pat_genus) %>% as.data.frame() #traspose for plot

colnames(pat_genus) <- patid
pat_genus <- pat_genus[-1, ]
pat_genus$genus <- rownames(pat_genus)
rownames(pat_genus) <- NULL
pat_genus <- cbind(pat_genus[,25], pat_genus[,1:24])
colnames(pat_genus)[1] <- "genus"

mmeta_genus <- melt(pat_genus, id.vars = "genus")

ggplot(mmeta_genus[order(mmeta_genus$genus,decreasing=T),], aes(x = variable, y = value, fill = genus), decreasing = TRUE) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_ucscgb("default")(21))+ 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))
