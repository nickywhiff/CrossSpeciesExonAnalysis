library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
library(scales)

setwd("~/OneDrive - Nexus365/CrossSpeciesExonAnalysis/gff_files/")

data<-read.table("ExonSizes_perOrganism.txt",sep="\t",header=TRUE)

data<-data %>%
	mutate(TotalExonic=CDS+NonCodingRNA+UTR_5prime+UTR_3prime+Other) %>%
	mutate(propCDS=CDS/TotalExonic) %>%
	mutate(propUTR=(UTR_5prime+UTR_3prime)/TotalExonic) %>%
	mutate(propncRNA=NonCodingRNA/TotalExonic)

extended<-data %>%
	gather(key="ExonType",value="size",CDS,NonCodingRNA,UTR_5prime,UTR_3prime,Other)

extended$Organism <- factor(extended$Organism, levels = c("Homo_sapiens.GRCh38.107","Mus_musculus.GRCm39.107","Danio_rerio.GRCz11.107","Drosophila_melanogaster.BDGP6.32.107","Caenorhabditis_elegans.WBcel235.107","Schizosaccharomyces_pombe.ASM294v2.54"))
extended$ExonType <- factor(extended$ExonType, levels = c("Other","NonCodingRNA","UTR_3prime","UTR_5prime","CDS"))

png("ExonRegions_byOrganism.png",height=500,width=900)
ggplot(extended, aes(x=Organism,y=size,fill=ExonType)) +
	geom_col(colour="black") +
	scale_fill_manual(values=c("#D95F02","#7570B3","#E7298A","#E6AB02","#1B9E77"),labels=c("Other"="Other","NonCodingRNA"="Non-coding RNA","UTR_5prime"="5'UTR","UTR_3prime"="3'UTR","CDS"="CDS")) +
	xlab("") +
	scale_x_discrete(labels=c("Homo_sapiens.GRCh38.107"="Homo sapiens\n(Human)","Mus_musculus.GRCm39.107"="Mus musculus\n(Mouse)","Danio_rerio.GRCz11.107"="Danio rerio\n(Zebrafish)","Drosophila_melanogaster.BDGP6.32.107"="Drosophila melanogaster\n(Fruit fly)","Caenorhabditis_elegans.WBcel235.107"="Caenorhabditis elegans\n(Roundworm)","Schizosaccharomyces_pombe.ASM294v2.54"="Schizosaccharomyces\npombe (Fission yeast)")) +
	ylab("Total size (bps)") +
	theme_classic() +
	coord_flip() +
	scale_y_continuous(labels = comma) +
	theme(legend.title=element_blank(),axis.text=element_text(size=16),axis.title=element_text(size=18),legend.text=element_text(size=18))
dev.off()


wes_data<-read.table("../HumanExomeOverlap/WES_captureOverlap_new.txt",sep="\t",header=TRUE)

wes_data<-wes_data %>%
	mutate(Neither=TotalSize-(InCapture+InBuffer)) %>%
	mutate(Capture=(InCapture/1000000)) %>%
	mutate(Not=(TotalSize-InCapture)/1000000) %>%
	gather(key="Capture",value="size",Capture,Not) %>%
	select(ExonType,Capture,size)

wes_data$ExonType <- factor(wes_data$ExonType, levels = c("CDS","UTR_5prime","UTR_3prime","NonCodingRNA","Other"))
wes_data$Capture <- factor(wes_data$Capture, levels = c("Not","Capture"))

png("HumanRegionsInCapture.png",height=500,width=400)
ggplot(wes_data, aes(x=ExonType,y=size,fill=ExonType,pattern=Capture)) +
	geom_col_pattern(color="black", pattern_angle = 45, pattern_spacing = 0.015, pattern_colour="grey", pattern_fill="grey",pattern_density = 1,) +
	scale_pattern_manual(values = c("stripe","none")) +
	scale_fill_manual(values=c("#1B9E77","#E6AB02","#E7298A","#7570B3","#D95F02")) +
	guides(pattern = guide_legend(override.aes = list(fill = "white")),
		colour = guide_legend(override.aes = list(fill = "#black")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
	scale_x_discrete(labels=c("Pseudogene"="Other","NonCodingRNA"="Non-coding\nRNA","UTR_5prime"="5'UTR","UTR_3prime"="3'UTR","CDS"="CDS")) +
	xlab("") +
	ylab("Total size (Mbps)") +
	theme_classic() +
	theme(axis.text=element_text(size=15),axis.title=element_text(size=17),legend.position="none",axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


