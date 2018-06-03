library("ggplot2")
library("stringr")
library("reshape2")
library("plyr")
Group1 <- read.table(file="Group1_Spore_Marker_Hits.csv", header = F)
head(Group1)
head(Group1_conf[,1:5])
colnames(Group1) <- c("target name", "target_accession",  "query name", "query_accession",    "E-value",  "score",  "bias",   "Domain_E-value",  "Domain_score",  "Domain_bias", "exp", "reg", "clu",  "ov", "env", "dom", "rep", "inc")

Group1$Genome = str_sub(Group1$`target name`,1,4)

total_hits <- melt(table(Group1$Genome))
Group1_conf <- Group1[log10(Group1$score) >= 1.5,]
num_hits <- melt(table(Group1_conf$Genome))
ggplot(num_hits) + geom_point(aes(x=num_hits$Var1, num_hits$value))
unique(Group1$`query name`)
ggplot(Group1) + geom_histogram(aes(Group1$score))



head(num_hits)
zeros <- as.data.frame(cbind(setdiff(total_hits$Var1, num_hits$Var1), rep(0,76)))

colnames(zeros) <- c("Genome", "value")
colnames(num_hits) <- c("Genome", "value")
head(num_hits)
head(zeros)
num_hits <- rbind2(num_hits,zeros)
num_hits$Genome <- factor(num_hits$Genome, levels=(num_hits$Genome)[order(num_hits$value)])
ggplot(num_hits) + geom_point(aes(x=num_hits$Genome, num_hits$value)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(Group1_conf) + geom_histogram(aes(Group1_conf$score))
ggplot(Group1_conf) + geom_bar(aes(Group1_conf$Genome))
ggplot(Group1) + geom_histogram(aes(log10(Group1$score))) + facet_wrap(~Group1$`query name`)
ggplot(Group1_conf) + geom_histogram(aes(log10(Group1_conf$score))) + facet_wrap(~Group1_conf$`query name`)
p <- ggplot(Group1_conf) + geom_tile(aes(y=Group1_conf$Genome, x=Group1_conf$`query name`))

x<-as.data.frame(cbind(Group1_conf$Genome,Group1$`query name`))
x$V1 <- factor(x$V1, levels=(x$V1)[order(x$V2)])
cbind(Group1_conf$Genome,Group1_conf$`query name`)

genome_attr <- read.table(file="tyson_genome_list_Firmicutes_names_env.txt", sep="\t",header = T)
genome_attr <- read.table(file="tyson_genome_list_Firmicutes_names_env_plus.txt", sep="\t",header = T)

head(genome_attr)
genome_attr$Genome = str_sub(genome_attr$DDBJ.ENA.GenBank.Accession,1,4)

Group1_conf <- merge(genome_attr, Group1_conf)
num_hits <- merge(genome_attr, num_hits)
head(num_hits)

num_hits$value
ggplot(num_hits) + geom_boxplot(aes(x=Environment, y=value)) + ylab("Number of Sporulation Specific Genes")
ggplot(num_hits) + geom_jitter(aes(x=num_hits$Environment, y=as.numeric(num_hits$value), color=GC)) + ylab("Number of Sporulation Specific Genes")
ggplot(num_hits) + geom_jitter(aes(x=num_hits$Environment, y=(as.numeric(num_hits$value)/num_hits$Genome.Size)*1000000, color=GC)) + ylab("Number of Sporulation Specific Genes per Mb")
ddply(num_hits, .(Environment), summarize, mean_hits=mean(value))
plot(x=num_hits$Genome.Size, y=num_hits$value)
abline(a=cor_geno_size_hits$coefficients[1], b=cor_geno_size_hits$coefficients[2])
plot(x=num_hits$Completeness, y=num_hits$value)
abline(a=cor_geno_com_hits$coefficients[1], b=cor_geno_com_hits$coefficients[2])
cor_geno_size_hits <- lm(num_hits$value ~ num_hits$Genome.Size)
cor_geno_com_hits <- lm(num_hits$value ~ num_hits$Completeness)
cor_geno_size_hits$coefficients

summary(cor_geno_size_hits)
summary(cor_geno_com_hits)
ggplot(num_hits[num_hits$Environment == "Biofilm",]) + geom_jitter(x=num_hits[num_hits$Environment == "Biofilm",]$Environment,y=num_hits[num_hits$Environment == "Biofilm",]value)
num_hits[num_hits$Environment == "Biofilm",]
num_hits[num_hits$value <= 10,]
num_hits[num_hits$value == 0, ]
