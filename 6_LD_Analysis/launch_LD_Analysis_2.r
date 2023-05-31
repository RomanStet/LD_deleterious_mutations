setwd(".../6_LD_Analysis/")

folder_Sift <- ".../Sift_annotation/"
folder_Neslia <- ".../Neslia_sequence/"

# Writes files SpeciesScaffold_Neslia.txt containing the SNPs for which the
# information about the nucleotide in Neslia is available
# Needed to compute the number of commun SNPs with Neslia

species <- c("Cg","Co")
scaffold <- seq(1,8,1)

for (m in species)
{
	for (l in scaffold)
	{
		SNP <- read.table(paste(m,l,sep=""), header=TRUE, sep="\t") 
				
		Neslia <- read.table(paste(folder_Neslia,"Neslia_",l,sep =""), header=TRUE, sep="\t")
		
		colnames(Neslia) <- c("CHROM", "POS", "Neslia")
		
		merged <- merge(x=SNP , y= Neslia, by="POS", all=FALSE)
		
		merged <- merged[order(merged$POS), ]
		
		write.table(merged, file = paste(m,l,"_Neslia.txt",sep =""), sep = "\t", quote = FALSE, row.names = FALSE)
	}
}

# Writes files Sift_Species_Scaffold.txt containing the Sift 
# annotation of SNPs present in the file SpeciesScaffold and
# writes 

species <- c("Cg","Co")
scaffold <- seq(1,8,1)

for (m in species)
{
	for (l in scaffold)
	{
		SNP <- read.table(paste(m,l,sep =""), header=TRUE, sep="\t") 
		
		v_col <- c("CHROM_POS",	"NUM_CHROM", "POS", "REF","ALT")
		
		for (i in c(1:(ncol(SNP)-5)))
		{
			v_col <- c(v_col, paste("ind_",i,sep=""))
		}
		
		colnames(SNP) <- v_col
		
		write.table(SNP, file = paste(m,l,".txt",sep =""), sep = "\t", quote = FALSE, row.names = FALSE)
		
		sites <- data.frame(POS = SNP$POS)

		Sift <- read.table(paste(folder_Sift,"Sift",l,sep =""), header=TRUE, sep="\t")

		write.table(merge(x= Sift, y=sites, by="POS", all=FALSE), file = paste("Sift_",m,"_",l,".txt",sep =""), sep = "\t", quote = FALSE, row.names = FALSE)
	}
}

####

for (m in species)
{
	for (l in scaffold)
	{
		SNP <- read.table(paste(m,l,".txt",sep=""), header=TRUE, sep="\t") ## on reprend les 1 fichiers
		
		Sift <- read.table(paste("Sift_",m,"_",l,".txt",sep =""), header=TRUE, sep="\t") 

		sites <- data.frame(POS = Sift$POS[!duplicated(Sift$POS)]) ## on extrait la colonne POS du fichier Sift_SNP
		
		write.table(merge(x=SNP , y= sites, by="POS", all=FALSE), file = paste(m,"_",l,"_Sift",".txt",sep =""), sep = "\t", quote = FALSE, row.names = FALSE) ## on rassemble les 2 et on écrit le fichier sortant
	}
}