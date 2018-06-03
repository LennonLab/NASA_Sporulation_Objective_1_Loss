library("ggplot2")
library("stringr")
library("reshape2")
library("plyr")

#bring in results from hmm based homology searching
Group1 <- read.table(file="Group1_Spore_Marker_Hits.csv", header = F)
head(Group1[,1:5])
colnames(Group1) <- c("target name", "target_accession",  "query name", "query_accession",    "E-value",  "score",  "bias",   "Domain_E-value",  "Domain_score",  "Domain_bias", "exp", "reg", "clu",  "ov", "env", "dom", "rep", "inc")
#Unique Genome Names
Group1$Genome = str_sub(Group1$`target name`,1,4)
#look at distribution of confidence in gene annotations
ggplot(Group1) + geom_histogram(aes(log10(Group1$score))) + facet_wrap(~Group1$`query name`)
#set confidence value to capture the high confidence mode and distribution 
Group1_conf <- Group1[log10(Group1$score) >= 1.5,]
#Plot distributions of confidence for each gene
ggplot(Group1_conf) + geom_histogram(aes(Group1_conf$score))
ggplot(Group1_conf) + geom_bar(aes(Group1_conf$Genome))
ggplot(Group1_conf) + geom_histogram(aes(log10(Group1_conf$score))) + facet_wrap(Group1_conf$`query name`~Group1_conf$Environment)

#This makes 352*30 something plots
#ggplot(Group1_conf) + geom_histogram(aes(log10(Group1_conf$score))) + facet_wrap(Group1_conf$Genome~Group1_conf$`query name`)

#count up number of hits per genome
total_hits <- melt(table(Group1$Genome))
num_hits <- melt(table(Group1_conf$Genome))
zeros <- as.data.frame(cbind(setdiff(total_hits$Var1, num_hits$Var1), rep(0,76)))

colnames(zeros) <- c("Genome", "value")
colnames(num_hits) <- c("Genome", "value")

num_hits <- rbind2(num_hits,zeros)
num_hits$value <- as.numeric(num_hits$value)

#Rank order number of sporulation specfic genes
num_hits$Genome <- factor(num_hits$Genome, levels=num_hits$Genome[order(as.numeric(num_hits$value))])
ggplot(num_hits) + geom_point(aes(x=num_hits$Genome, (num_hits$value/num_hits$Genome.Size..bp.)*1000000)) + theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1))
ggplot(num_hits_NC) + geom_point(aes(x=num_hits_NC$Genome, num_hits_NC$value)) + theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1))


#Add in detail about where metagenomes were from
genome_attr <- read.table(file="tyson_genome_list_Firmicutes_names_env.txt", sep="\t",header = T)
genome_attr <- read.table(file="tyson_genome_list_Firmicutes_plus.txt", sep="\t",header = T)

genome_attr$Genome = str_sub(genome_attr$DDBJ.ENA.GenBank.Accession,1,4)
#genome_attr <- genome_attr[,2:3]
Group1_conf <- merge(genome_attr, Group1_conf)
Group1 <- merge(genome_attr, Group1)
num_hits <- merge(genome_attr, num_hits)
total_num_hits <- merge(genome_attr, total_hits)
num_hits$NCBI_TAXONOMY <- gsub( "\\s.*", "", num_hits$NCBI.Organism.Name)
total_num_hits$NCBI_TAXONOMY <- gsub( "\\s.*", "", total_num_hits$NCBI.Organism.Name)
head(num_hits)
total_num_hits_NC <- total_num_hits[total_num_hits$Genome.Quality == "Near complete",]
num_hits_NC <- num_hits[num_hits$Genome.Quality == "Near complete",]
ggplot(total_num_hits_NC) + geom_boxplot(aes(x=Environment, y=value)) + theme_classic() + geom_jitter(aes(x=Environment, y=value)) + ylab("Number of Sporulation Specific Genes")
unique(num_hits_NC$Genome)
ggplot(num_hits_NC) + geom_boxplot(aes(x=Environment, y=value)) + theme_classic() + geom_jitter(aes(x=Environment, y=value)) + ylab("Number of Sporulation Specific Genes")
ggplot(num_hits) + geom_boxplot(aes(x=Environment, y=value)) + theme_classic() + geom_jitter(aes(x=Environment, y=value)) + ylab("Number of Sporulation Specific Genes")
ggplot(num_hits) + geom_boxplot(aes(x=NCBI_TAXONOMY, y=value)) + theme_classic() + geom_jitter(aes(color=Environment, x=NCBI_TAXONOMY, y=value)) + ylab("Number of Sporulation Specific Genes") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(num_hits) + geom_boxplot(aes(x=Environment, y=num_hits$value/num_hits$No..Predicted.Genes)) + theme_classic() + geom_jitter(aes(x=Environment, y=num_hits$value/num_hits$No..Predicted.Genes)) + ylab("Number of Sporulation Specific Genes per Gene")
ggplot(num_hits) + geom_boxplot(aes(x=Environment, y=(num_hits$value/num_hits$Genome.Size..bp.)*1000000)) + theme_classic() + geom_jitter(aes(x=Environment, y=(num_hits$value/num_hits$Genome.Size..bp.)*1000000)) + ylab("Number of Sporulation Specific Genes per Mb")
ggplot(num_hits) + geom_boxplot(aes(x=Environment, y=(num_hits$value/num_hits$Genome.Size..bp.)*1000000)) + theme_classic() + geom_jitter(aes(color=num_hits$Genome.Quality, x=Environment, y=(num_hits$value/num_hits$Genome.Size..bp.)*1000000)) + ylab("Number of Sporulation Specific Genes per Mb")
ggplot(num_hits) + geom_boxplot(aes(x=Environment, y=(num_hits$value/num_hits$Genome.Size..bp.)*1000000)) + theme_classic() + geom_jitter(aes(color=num_hits$NCBI.Organism.Name, x=Environment, y=(num_hits$value/num_hits$Genome.Size..bp.)*1000000)) + ylab("Number of Sporulation Specific Genes per Mb")
install.packages("data.tree")

ggplot(num_hits) + geom_boxplot(aes(x=Environment, y=(num_hits$value/num_hits$No..Predicted.Genes))) + theme_classic() + geom_jitter(aes(color=num_hits$NCBI_TAXONOMY, x=Environment, y=num_hits$value/num_hits$No..Predicted.Genes)) + ylab("Number of Sporulation Specific Genes per Mb")
n
p <- ggplot(num_hits[num_hits$NCBI_TAXONOMY == "Firmicutes",]) + geom_boxplot(aes(x=Environment, y=(num_hits[num_hits$NCBI_TAXONOMY == "Firmicutes",]$value/num_hits[num_hits$NCBI_TAXONOMY == "Firmicutes",]$Genome.Size..bp.)*1000000)) + theme_classic() + geom_jitter(aes(color=Genome.Quality,x=Environment, y=(num_hits[num_hits$NCBI_TAXONOMY == "Firmicutes",]$value/num_hits[num_hits$NCBI_TAXONOMY == "Firmicutes",]$Genome.Size..bp.)*1000000)) + ylab("Number of Sporulation Specific Genes per Mb") + facet_wrap(~num_hits[num_hits$NCBI_TAXONOMY == "Firmicutes",]$NCBI_TAXONOMY)

p
p <- ggplot(num_hits[num_hits$Environment == "Gut",]) + geom_boxplot(aes(x=num_hits[num_hits$Environment == "Gut",]$NCBI_TAXONOMY, y=(num_hits[num_hits$Environment == "Gut",]$value/num_hits[num_hits$Environment == "Gut",]$Genome.Size..bp.)*1000000)) + theme_classic() + geom_jitter(aes(x=num_hits[num_hits$Environment == "Gut",]$NCBI_TAXONOMY, y=(num_hits[num_hits$Environment == "Gut",]$value/num_hits[num_hits$Environment == "Gut",]$Genome.Size..bp.)*1000000)) + ylab("Number of Sporulation Specific Genes per Mb") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

p
p <- ggplot(num_hits[num_hits$Environment == "Bioreactor",]) + geom_boxplot(aes(x=num_hits[num_hits$Environment == "Bioreactor",]$NCBI_TAXONOMY, y=(num_hits[num_hits$Environment == "Bioreactor",]$value/num_hits[num_hits$Environment == "Bioreactor",]$Genome.Size..bp.)*1000000)) + theme_classic() + geom_jitter(aes(x=num_hits[num_hits$Environment == "Bioreactor",]$NCBI_TAXONOMY, y=(num_hits[num_hits$Environment == "Bioreactor",]$value/num_hits[num_hits$Environment == "Bioreactor",]$Genome.Size..bp.)*1000000)) + ylab("Number of Sporulation Specific Genes per Mb") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

p
env_aov <- aov(num_hits$value ~ num_hits$Environment)
summary(env_aov)
TukeyHSD(env_aov)
oneway.test(value ~ Environment, num_hits)
ggplot(num_hits) + geom_boxplot(aes(x=Environment, y=value)) + theme_classic() + geom_jitter(aes(x=Environment, y=value)) + ylab("Number of Sporulation Specific Genes")
ggplot(num_hits[num_hits$Environment == "Biofilm" | num_hits$Environment == "Bioreactor" | num_hits$Environment == "Built Environment" | num_hits$Environment == "Gut" | num_hits$Environment == "Sediment" | num_hits$Environment == "Terrestial",]) + geom_violin(aes(x=Environment, y=value)) + theme_classic() + geom_jitter(aes(x=Environment, y=value)) + ylab("Number of Sporulation Specific Genes")

ggplot(num_hits) + geom_jitter(aes(x=Environment, y=value)) + ylab("Number of Sporulation Specific Genes")
ddply(num_hits, .(Environment), summarize, mean_hits=mean(value))
ggplot(num_hits[num_hits$Environment == "Biofilm",]) + geom_jitter(x=num_hits[num_hits$Environment == "Biofilm",]$Environment,y=num_hits[num_hits$Environment == "Biofilm",]value)
num_hits$Pos = as.numeric(num_hits$value > 5)
num_hits$Neg = as.numeric(num_hits$value < 5)
as.numeric(num_hits$value > 5)
ggplot(num_hits[num_hits$Environment == "Biofilm" | num_hits$Environment == "Bioreactor" | num_hits$Environment == "Built Environment" | num_hits$Environment == "Gut" | num_hits$Environment == "Sediment" | num_hits$Environment == "Terrestial",], aes(x=Environment, fill = factor(Pos))) +  geom_bar(position = "fill") + scale_fill_grey() + theme_classic()
head(Group1_conf)
unique(Group1_conf$`query name`)
Group_conf_counts <- table(Group1_conf[,c(1,21)])
Group_conf_counts[Group_conf_counts >= 1] <- 1
Group1_counts <- table(Group1[,c(1,21)])
t <- Group1_counts[setdiff(rownames(Group1_counts), rownames(Group_conf_counts)),]
t[t >= 0] <- 0
t
Group_conf_counts <- rbind(Group_conf_counts,t)
heatmap(Group_conf_counts, scale = "none", Colv = NA, col = c("white", "grey"), main = "Presence/Abscence of Sporulation Specific Genes") 

Group1_counts[Group1_counts >= 1] <- 1
heatmap(Group1_counts, scale = "none", Colv = NA, col = c("white", "grey"), main = "Presence/Abscence of Sporulation Specific Genes") 

head(num_hits)
