
#### Computes the overlap rate between blocks and duplications detected by smoove and delly

args <- commandArgs(trailingOnly = TRUE)

species <- as.character(args[1])
scaffold <- as.character(args[2])
start <- as.numeric(args[3])
end <- as.numeric(args[4])
bloc <- as.numeric(args[5])

ind_folder <- ".../Accession_SRA_Number/"
block_folder <- ".../7_SV_Analysis/1_Block_detection/"

delly_calls <- ".../7_SV_Analysis/3_SV_calling_Delly/"
smoove_calls <- ".../7_SV_Analysis/4_SV_calling_Smoove/"

ind_file <- read.table(paste(ind_folder,"ID_ind",sep = ""), header=TRUE, sep="\t")

bloc_ind <- read.table(paste(block_folder,species,scaffold,"_bloc_stats_ind_2.txt",sep = ""), header=TRUE, sep="\t")

v_ind <- bloc_ind[which(bloc_ind$bloc==bloc),]$ind

bloc_ind_write_delly <- data.frame()
bloc_ind_write_smoove <- data.frame()

for (i in v_ind)
{
	ID <- ind_file[which(as.character(ind_file$Species)==species & ind_file$Number == i),]$ID[1]
	line <- bloc_ind[which(bloc_ind$bloc==bloc & bloc_ind$ind == i),]
	line$ID <- ID
	if (file.exists(paste(delly_calls,ID,".vcf",sep = "")))
	{
		dup_delly <- read.table(paste(delly_calls,ID,"_DUP.table",sep = ""), header=FALSE, sep=" ")
		dup_delly <- dup_delly[,seq(1,9,1)]
		colnames(dup_delly) <- c("scaffold","SV","genotype","POS","CIPOS","END","CIEND","QUAL","FILTER")
	
		dup_delly_match <- dup_delly[which(as.character(dup_delly$scaffold)==paste("SCF_",scaffold,sep = "") &
		((((dup_delly$POS-dup_delly$CIPOS)<end) & ((dup_delly$POS-dup_delly$CIPOS)>start)) |
		(((dup_delly$END+dup_delly$CIEND)>start) & ((dup_delly$END+dup_delly$CIEND)<end)) |
		(((dup_delly$POS-dup_delly$CIPOS)<start) & ((dup_delly$END+dup_delly$CIEND)>end))) & ((dup_delly$END+dup_delly$CIEND-(dup_delly$POS-dup_delly$CIPOS))<size_dup_max))
		,]
	
		if(nrow(dup_delly_match)!=0)
		{
			for (j in c(1:nrow(dup_delly_match)))
			{
				line$filter <- "PASS"
			
				if (((dup_delly_match$POS[j]-dup_delly_match$CIPOS[j])<end) & ((dup_delly_match$POS[j]-dup_delly_match$CIPOS[j])>start) & ((dup_delly_match$END[j]+dup_delly_match$CIEND[j])>end))
				{
					line$overlap <- (end-(dup_delly_match$POS[j]-dup_delly_match$CIPOS[j]))/(end-start)
				} else if (((dup_delly_match$END[j]+dup_delly_match$CIEND[j])>start) & ((dup_delly_match$END[j]+dup_delly_match$CIEND[j])<end) & ((dup_delly_match$POS[j]-dup_delly_match$CIPOS[j])<start))
				{
					line$overlap <- (dup_delly_match$END[j]+dup_delly_match$CIEND[j]-start)/(end-start)
				} else {
					line$overlap <- (dup_delly_match$END[j]+dup_delly_match$CIEND[j]-(dup_delly_match$POS[j]-dup_delly_match$CIPOS[j]))/(end-start)
				}
			
				if (as.character(dup_delly_match$genotype[j])=="0/1")
				{
					line$genotype <- "Het"
				} else if (as.character(dup_delly_match$genotype[j])=="1/1") {
					line$genotype <- "Hom"
				} else {
					line$genotype <- "NA"
				}
				
				line$QUAL <- dup_delly_match$QUAL[j]
				
				bloc_ind_write_delly <- rbind(bloc_ind_write_delly,line)
			}
		} else {
			line$filter <- "NA"
			line$overlap <- "NA"
			line$genotype <- "NA"
			line$QUAL <- "NA"
			bloc_ind_write_delly <- rbind(bloc_ind_write_delly,line)
		}
	} else {
		line$filter <- "NA"
		line$overlap <- "NA"
		line$genotype <- "NA"
		line$QUAL <- "NA"
		bloc_ind_write_delly <- rbind(bloc_ind_write_delly,line)
	}
	
	if (file.exists(paste(smoove_calls,ID,"_results/",ID,"-smoove.genotyped.vcf",sep = "")))
	{
		dup_smoove <- read.table(paste(smoove_calls,ID,"_results/",ID,"_DUP.table",sep = ""), header=FALSE, sep=" ")
		dup_smoove <- dup_smoove[,seq(1,9,1)]
		colnames(dup_smoove) <- c("scaffold","SV","genotype","POS","CIPOS","END","CIEND","QUAL","FILTER")
	
		dup_smoove_match <- dup_smoove[which(as.character(dup_smoove$scaffold)==paste("SCF_",scaffold,sep = "") &
		((((dup_smoove$POS-dup_smoove$CIPOS)<end) & ((dup_smoove$POS-dup_smoove$CIPOS)>start)) |
		(((dup_smoove$END+dup_smoove$CIEND)>start) & ((dup_smoove$END+dup_smoove$CIEND)<end)) |
		(((dup_smoove$POS-dup_smoove$CIPOS)<start) & ((dup_smoove$END+dup_smoove$CIEND)>end))) & ((dup_smoove$END+dup_smoove$CIEND-(dup_smoove$POS-dup_smoove$CIPOS))<size_dup_max))
		,]
	
		line <- bloc_ind[which(bloc_ind$bloc==bloc & bloc_ind$ind == i),]
		line$ID <- ID
	
		if(nrow(dup_smoove_match)!=0)
		{
			for (j in c(1:nrow(dup_smoove_match)))
			{
				line$filter <- "PASS"
			
				if (((dup_smoove_match$POS[j]-dup_smoove_match$CIPOS[j])<end) & ((dup_smoove_match$POS[j]-dup_smoove_match$CIPOS[j])>start) & ((dup_smoove_match$END[j]+dup_smoove_match$CIEND[j])>end))
				{
					line$overlap <- (end-(dup_smoove_match$POS[j]-dup_smoove_match$CIPOS[j]))/(end-start)
				} else if (((dup_smoove_match$END[j]+dup_smoove_match$CIEND[j])>start) & ((dup_smoove_match$END[j]+dup_smoove_match$CIEND[j])<end) & ((dup_smoove_match$POS[j]-dup_smoove_match$CIPOS[j])<start))
				{
					line$overlap <- (dup_smoove_match$END[j]+dup_smoove_match$CIEND[j]-start)/(end-start)
				} else {
					line$overlap <- (dup_smoove_match$END[j]+dup_smoove_match$CIEND[j]-(dup_smoove_match$POS[j]-dup_smoove_match$CIPOS[j]))/(end-start)
				}
			
				if (as.character(dup_smoove_match$genotype[j])=="0/1")
				{
					line$genotype <- "Het"
				} else if (as.character(dup_smoove_match$genotype[j])=="1/1") {
					line$genotype <- "Hom"
				} else {
					line$genotype <- "NA"
				}
				
				line$QUAL <- dup_smoove_match$QUAL[j]
				
				bloc_ind_write_smoove <- rbind(bloc_ind_write_smoove,line)
			}
		} else {
			line$filter <- "NA"
			line$overlap <- "NA"
			line$genotype <- "NA"
			line$QUAL <- "NA"
			bloc_ind_write_smoove <- rbind(bloc_ind_write_smoove,line)
		}
	} else {
		line$filter <- "NA"
		line$overlap <- "NA"
		line$genotype <- "NA"
		line$QUAL <- "NA"
		bloc_ind_write_smoove <- rbind(bloc_ind_write_smoove,line)
	}
}

colnames(bloc_ind_write_delly) <- c(colnames(bloc_ind),"ID","filter","overlap","genotype","QUAL") 

write.table(bloc_ind_write_delly, paste(species,scaffold,"_duplic_PASS_bloc_",bloc,"_delly.txt",sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")

colnames(bloc_ind_write_smoove) <- c(colnames(bloc_ind),"ID","filter","overlap","genotype","QUAL") 

write.table(bloc_ind_write_smoove, paste(species,scaffold,"_duplic_PASS_bloc_",bloc,"_smoove.txt",sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")

