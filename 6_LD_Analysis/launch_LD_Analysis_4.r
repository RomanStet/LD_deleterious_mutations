setwd(".../6_LD_Analysis/")

folder_Sift <- ".../Sift_annotation/"
folder_Neslia <- ".../Neslia_sequence/"

species <- c("Cg","Co")
scaffold <- seq(1,8,1)

####

for (m in species)
{
	for (l in scaffold)
	{
		SNP <- read.table(paste(m,l,"Sift",sep="_"), header=TRUE, sep="\t")
		
		v_col <- c("POS", "CHROM_POS", "CHROM", "REF", "ALT")
		
		for (i in c(1:(length(colnames(SNP))-5)))
		{
			if (i %% 2 == 1)
			{
				chr <- 1
			}
			else
			{
				chr <- 2
			}
			
			v_col <- c(v_col, paste("ind_",ceiling(i/2),"_",chr,sep=""))
		}
		
		colnames(SNP) <- v_col
		
		Neslia <- read.table(paste(folder_Neslia,"Neslia_",l,sep =""), header=TRUE, sep="\t")
		
		colnames(Neslia) <- c("CHROM", "POS", "Neslia")
		
		merged <- merge(x=SNP , y= Neslia, by=c("POS","CHROM"), all=FALSE)
		
		merged <- merged[order(merged$POS), ]
		
		write.table(merged, file = paste(m,"_",l,"_Sift_2",".txt",sep =""), sep = "\t", quote = FALSE, row.names = FALSE)
	}
}

####

species <- c("Cg","Co")
scaffold <- seq(1,8,1)

for (m in species)
{
	for (l in scaffold)
	{
		SNP_Sift <- read.table(paste(m,"_",l,"_Sift_2.txt",sep =""), header=TRUE, sep="\t")
		Sift_SNP <- read.table(paste("Sift_",m,"_",l,".txt",sep =""), header=TRUE, sep="\t") 
		
		SNP_Sift <- SNP_Sift[which(SNP_Sift$Neslia != "NA"),]

		nrow_doc <- nrow(SNP_Sift)

		v_POS <- numeric()
		v_CHROM <- character()
		v_ind_chr <- character()
		v_nucl_Neslia <- character()
		v_nucl <- character()
		v_SIFT_SCORE <- numeric()

		headers <- data.frame(POS = v_POS, CHROM = v_CHROM, ind_chr = v_ind_chr, nucl_Neslia = v_nucl_Neslia, nucl = v_nucl, SIFT_SCORE = v_SIFT_SCORE)

		write.table(headers, paste(m,l,"_annotated.txt",sep=""), quote = FALSE, sep = "\t")

		sub_SNP <- SNP_Sift[grepl("ind", colnames(SNP_Sift))]
		ncol <- ncol(sub_SNP)

		for (i in c(1:nrow_doc))
		{	
			POS <- SNP_Sift$POS[i]
			sub_site <- SNP_Sift[which(SNP_Sift$POS == POS),]
			sub_SNP_line <- sub_SNP[i,]
			sub_SNP_line_no_na <- sub_SNP_line[,which(sub_SNP_line[1,]!= ".")]
			ident <- ifelse(sub_site$Neslia==sub_SNP_line_no_na,1,0)
			ident_2 <- ifelse(sub_SNP_line_no_na[1,1]==sub_SNP_line_no_na,1,0)
			
			if (sum(ident_2) == length(sub_SNP_line_no_na) | sum(ident) == 0)
			{
				next
			}
			
			sub_SIFT <- Sift_SNP [which(Sift_SNP $POS == POS),]
			nrow_sub <- nrow(sub_SIFT)

			Nes <- as.character(SNP_Sift$Neslia[i])
			
			for (j in c(1:ncol))
			{
				nucl_ind <- as.character(sub_SNP_line[,j])
				
				for (k in c(1:nrow_sub))
				{	
					if(nucl_ind == ".")
					{				
						line <- data.frame(POS = POS, CHROM = l, ind_chr = colnames(sub_SNP)[j], 
						nucl_Neslia = Nes, nucl = as.character("NA"), SIFT_SCORE = "NA")
										
						write.table(line, paste(m,l,"_annotated.txt",sep=""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
						
						break
					}
					if (nucl_ind == as.character(sub_SIFT$New_allele[k]))
					{			
						line <- data.frame(POS = POS, CHROM = l, ind_chr = colnames(sub_SNP)[j],
						nucl_Neslia = Nes, nucl = nucl_ind,
						SIFT_SCORE = sub_SIFT$SIFT_score[k])
										
						write.table(line, paste(m,l,"_annotated.txt",sep=""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
					
						break
					}
					if (k == nrow_sub)
					{
						line <- data.frame(POS = POS, CHROM = l, ind_chr = colnames(sub_SNP)[j], 
						nucl_Neslia = Nes, nucl = nucl_ind, SIFT_SCORE = "NA")
										
						write.table(line, paste(m,l,"_annotated.txt",sep=""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
					}
				}
			}
		}
	}
}

####

species <- c("Cg","Co")
scaffold <- seq(1,8,1)

for (m in species)
{
	for (l in scaffold)
	{
		data <- read.table(paste(m,l,"_annotated.txt",sep=""), header=TRUE, sep="\t")
				
		write.table(data[0,], paste(m,l,"_neutr.txt",sep = ""), quote = FALSE, sep = "\t", append = FALSE, col.names = TRUE, row.names = FALSE)
		
		write.table(data[0,], paste(m,l,"_m_del.txt",sep = ""), quote = FALSE, sep = "\t", append = FALSE, col.names = TRUE, row.names = FALSE)
		
		write.table(data[0,], paste(m,l,"_del.txt",sep = ""), quote = FALSE, sep = "\t", append = FALSE, col.names = TRUE, row.names = FALSE)
				
		i<-1
		
		while (i <= nrow(data))
		{
			site <- data$POS[i]
			sub_site <- data[which(data$POS == site),]
			sub_site_na_omit <- na.omit(sub_site)
			
			if (!all(is.na(sub_site$nucl)) & !all(is.na(sub_site$SIFT_SCORE)))
			{
				for (j in nrow(sub_site_na_omit))
				{
					ident <- ifelse(sub_site_na_omit$nucl[j]==sub_site_na_omit$nucl,1,0)
				}
			
				if (0 %in% ident)
				{					
					if (1 %in% sub_site_na_omit$SIFT_SCORE & var(sub_site_na_omit$SIFT_SCORE) == 0)
					{
						write.table(sub_site, paste(m,l,"_neutr.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
					}
					
					if (sub_site_na_omit$SIFT_SCORE[which(ident %in% 0)][1]<1 & sub_site_na_omit$SIFT_SCORE[which(ident %in% 1)][1] == 1)	
					{
						if (sub_site_na_omit$SIFT_SCORE[which(ident %in% 0)][1]<=0.05) 
						{
							write.table(sub_site, paste(m,l,"_del.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
						} else {
							write.table(sub_site, paste(m,l,"_m_del.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
						}
					}
					
					if (sub_site_na_omit$SIFT_SCORE[which(ident %in% 1)][1]<1 & sub_site_na_omit$SIFT_SCORE[which(ident %in% 0)][1] == 1)
					{
						if (sub_site_na_omit$SIFT_SCORE[which(ident %in% 1)][1]<=0.05) 
						{
							write.table(sub_site, paste(m,l,"_del.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
						} else {
							write.table(sub_site, paste(m,l,"_m_del.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
						}
					}
				}
			}

			i <- i + nrow(sub_site)
		}
	}
}

####

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr")

exclud_ind_Co <- c(4,9,11,14,23,29,30,31)
exclud_ind_Cg <- c(22,45,46,52,53,64,88,90,109,111,146)

for (s in species)
{
	for (j in scaffold)
	{	
		for (k in mut)
		{
			data_mut <- read.table(paste(s,j,"_",k,".txt",sep = ""), header=TRUE, sep="\t")
			
			data_mut$ind_chr <- as.numeric(gsub(".*?([0-9]+).*", "\\1",data_mut$ind_chr))
			
			if (as.character(s)=="Cg")
			{
				data_mut <- data_mut[which(!(data_mut$ind_chr %in% exclud_ind_Cg)),]
			} else {
				data_mut <- data_mut[which(!(data_mut$ind_chr %in% exclud_ind_Co)),]
			}
			
			v_POS <- numeric()
			v_CHROM <- character()
			v_ind <- character()
			v_genotype <- character()
			v_Freq <- numeric()
			v_n_ind <- numeric()
			
			headers <- data.frame(POS = v_POS, CHROM = v_CHROM, ind = v_ind, genotype = v_genotype)
			write.table(headers, paste(s,j,"_",k,"_genotype.txt",sep = ""), quote = FALSE, sep = "\t")

			headers <- data.frame(POS = v_POS, CHROM = v_CHROM, Freq = v_Freq, n_ind = v_n_ind)
			write.table(headers, paste(s,j,"_",k,"_freq_2.txt",sep = ""), quote = FALSE, sep = "\t")
			
			i <- 1
			
			while (i <= nrow(data_mut))
			{
				site <- data_mut$POS[i]
				sub_site <- data_mut[which(data_mut$POS == site),]
				
				n_chr_tot <- 0
				
				for (l in unique(sub_site$ind_chr))
				{
					if (!any(is.na(sub_site[which(sub_site$ind_chr == l),]$nucl)))
					{
						n_chr_tot <- n_chr_tot + 2
					}
				}
				
				if (as.character(k) == "neutr") 
				{
					sub_site_der <- sub_site[which(sub_site$nucl != sub_site$nucl_Neslia | is.na(sub_site$nucl)),]
				} else {
					sub_site_der <- sub_site[which((sub_site$nucl != sub_site$nucl_Neslia & sub_site$SIFT_SCORE < 1) | is.na(sub_site$nucl)),]
				}
		
				ind <- numeric()
				n_chr <- 0
					
					for (m in unique(sub_site_der$ind_chr))
					{
						if (nrow(sub_site_der[which(sub_site_der$ind_chr == m),]) == 2 & !any(is.na(sub_site_der[which(sub_site_der$ind_chr == m),]$nucl)))
						{
							v_POS <- numeric()
							v_CHROM <- character()
							v_ind <- character()
							v_genotype <- character()
								
							v_POS <- c(v_POS, site)
							v_CHROM <- c(v_CHROM, j)
							v_ind <- c(v_ind, m)
							v_genotype <- c(v_genotype, "Hom")	
								
							line <- data.frame(POS_ = v_POS, CHROM = v_CHROM, ind = v_ind, genotype = v_genotype)
							write.table(line, paste(s,j,"_",k,"_genotype.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
							
							n_chr <- n_chr + 2
						}
						
						if (nrow(sub_site_der[which(sub_site_der$ind_chr == m),]) == 1 & !any(is.na(sub_site_der[which(sub_site_der$ind_chr == m),]$nucl)))
						{
							v_POS <- numeric()
							v_CHROM <- character()
							v_ind <- character()
							v_genotype <- character()
								
							v_POS <- c(v_POS, site)
							v_CHROM <- c(v_CHROM, j)
							v_ind <- c(v_ind, m)
							v_genotype <- c(v_genotype, "Het")	
								
							line <- data.frame(POS_ = v_POS, CHROM = v_CHROM, ind = v_ind, genotype = v_genotype)
							write.table(line, paste(s,j,"_",k,"_genotype.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
							
							n_chr <- n_chr + 1
						}
					}
					
					v_POS <- numeric()
					v_CHROM <- numeric()
					v_Freq <- numeric()
					v_n_ind <- numeric()
					
					v_POS <- c(v_POS, site)
					v_CHROM <- c(v_CHROM, j)
					v_Freq <- c(v_Freq, n_chr/n_chr_tot)
					v_n_ind <- c(v_n_ind, n_chr_tot/2)
					
					line <- data.frame(POS_ = v_POS, CHROM = v_CHROM, Freq = v_Freq, n_ind = v_n_ind)
					write.table(line, paste(s,j,"_",k,"_freq_2.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
					
					i <- i + nrow(sub_site)
			}
		}
	}
}

### Polarized by SIFT score for del and m_del

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del")

exclud_ind_Co <- c(4,9,11,14,23,29,30,31)
exclud_ind_Cg <- c(22,45,46,52,53,64,88,90,109,111,146)

for (s in species)
{
	for (j in scaffold)
	{	
		for (k in mut)
		{
			data_mut <- read.table(paste(s,j,"_",k,".txt",sep = ""), header=TRUE, sep="\t")
			
			data_mut$ind_chr <- as.numeric(gsub(".*?([0-9]+).*", "\\1",data_mut$ind_chr))
			
			if (as.character(s)=="Cg")
			{
				data_mut <- data_mut[which(!(data_mut$ind_chr %in% exclud_ind_Cg)),]
			} else {
				data_mut <- data_mut[which(!(data_mut$ind_chr %in% exclud_ind_Co)),]
			}
			
			v_POS <- numeric()
			v_CHROM <- character()
			v_ind <- character()
			v_genotype <- character()
			v_Freq <- numeric()
			v_n_ind <- numeric()
			
			headers <- data.frame(POS = v_POS, CHROM = v_CHROM, ind = v_ind, genotype = v_genotype)
			write.table(headers, paste(s,j,"_",k,"_genotype_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t")

			headers <- data.frame(POS = v_POS, CHROM = v_CHROM, Freq = v_Freq, n_ind = v_n_ind)
			write.table(headers, paste(s,j,"_",k,"_freq_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t")
			
			i <- 1
			
			while (i <= nrow(data_mut))
			{
				site <- data_mut$POS[i]
				sub_site <- data_mut[which(data_mut$POS == site),]
				
				n_chr_tot <- 0
				
				for (l in unique(sub_site$ind_chr))
				{
					if (!any(is.na(sub_site[which(sub_site$ind_chr == l),]$nucl)))
					{
						n_chr_tot <- n_chr_tot + 2
					}
				}

				sub_site_der <- sub_site[which(sub_site$SIFT_SCORE < 1 | is.na(sub_site$nucl)),]
		
				ind <- numeric()
				n_chr <- 0
					
					for (m in unique(sub_site_der$ind_chr))
					{
						if (nrow(sub_site_der[which(sub_site_der$ind_chr == m),]) == 2 & !any(is.na(sub_site_der[which(sub_site_der$ind_chr == m),]$nucl)))
						{
							v_POS <- numeric()
							v_CHROM <- character()
							v_ind <- character()
							v_genotype <- character()
								
							v_POS <- c(v_POS, site)
							v_CHROM <- c(v_CHROM, j)
							v_ind <- c(v_ind, m)
							v_genotype <- c(v_genotype, "Hom")	
								
							line <- data.frame(POS_ = v_POS, CHROM = v_CHROM, ind = v_ind, genotype = v_genotype)
							write.table(line, paste(s,j,"_",k,"_genotype_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
							
							n_chr <- n_chr + 2
						}
						
						if (nrow(sub_site_der[which(sub_site_der$ind_chr == m),]) == 1 & !any(is.na(sub_site_der[which(sub_site_der$ind_chr == m),]$nucl)))
						{
							v_POS <- numeric()
							v_CHROM <- character()
							v_ind <- character()
							v_genotype <- character()
								
							v_POS <- c(v_POS, site)
							v_CHROM <- c(v_CHROM, j)
							v_ind <- c(v_ind, m)
							v_genotype <- c(v_genotype, "Het")	
								
							line <- data.frame(POS_ = v_POS, CHROM = v_CHROM, ind = v_ind, genotype = v_genotype)
							write.table(line, paste(s,j,"_",k,"_genotype_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
							
							n_chr <- n_chr + 1
						}
					}
					
					v_POS <- numeric()
					v_CHROM <- numeric()
					v_Freq <- numeric()
					v_n_ind <- numeric()
					
					v_POS <- c(v_POS, site)
					v_CHROM <- c(v_CHROM, j)
					v_Freq <- c(v_Freq, n_chr/n_chr_tot)
					v_n_ind <- c(v_n_ind, n_chr_tot/2)
					
					line <- data.frame(POS_ = v_POS, CHROM = v_CHROM, Freq = v_Freq, n_ind = v_n_ind)
					write.table(line, paste(s,j,"_",k,"_freq_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
					
					i <- i + nrow(sub_site)
			}
		}
	}
}

####

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr")

data_tot_freq <- numeric()
data_tot_genotype <- numeric()

for (m in species)
{
	for (j in scaffold)
	{	
		data_scaffold_freq <- numeric()
		data_scaffold_genotype <- numeric()
		
		for (k in mut)
		{
			data_freq <- read.table(paste(m,j,"_",k,"_freq_2.txt",sep = ""), header=TRUE, sep="\t")
			
			data_genotype <- read.table(paste(m,j,"_",k,"_genotype.txt",sep = ""), header=TRUE, sep="\t")
			
			nrow_freq <- nrow(data_freq)
			
			nrow_genotype <- nrow(data_genotype)
		
			data_freq$mut_type <- rep(k,nrow_freq)
			
			data_genotype$mut_type <- rep(k,nrow_genotype)
						
			data_scaffold_freq <- rbind(data_scaffold_freq, data_freq)
			data_scaffold_genotype <- rbind(data_scaffold_genotype, data_genotype)
			
			data_tot_freq <- rbind(data_tot_freq, data_freq)
			data_tot_genotype <- rbind(data_tot_genotype, data_genotype)
		}
		
		data_scaffold_freq <- data_scaffold_freq[which(!(data_scaffold_freq$Freq==0 | data_scaffold_freq$Freq==1)),]
		
		write.table(data_scaffold_freq, paste(m,j,"_freq_2.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
		
		data_scaffold_genotype <- data_scaffold_genotype[which(data_scaffold_genotype$POS %in% data_scaffold_freq$POS),]
		
		write.table(data_scaffold_genotype, paste(m,j,"_genotype.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
	}
	
	data_tot_freq <- data_tot_freq[which(!(data_tot_freq$Freq==0 | data_tot_freq$Freq==1)),]
	
	write.table(data_tot_freq, paste(m,"_freq_2.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
}

###

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del")

data_tot_freq <- numeric()
data_tot_genotype <- numeric()

for (m in species)
{
	for (j in scaffold)
	{	
		data_scaffold_freq <- numeric()
		data_scaffold_genotype <- numeric()
		
		for (k in mut)
		{
			if (as.character(k)=="neutr")
			{
				data_freq <- read.table(paste(m,j,"_",k,"_freq_2.txt",sep = ""), header=TRUE, sep="\t")
			
				data_genotype <- read.table(paste(m,j,"_",k,"_genotype.txt",sep = ""), header=TRUE, sep="\t")
			} else {
				data_freq <- read.table(paste(m,j,"_",k,"_freq_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
			
				data_genotype <- read.table(paste(m,j,"_",k,"_genotype_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
			}
			
			nrow_freq <- nrow(data_freq)
			
			nrow_genotype <- nrow(data_genotype)
		
			data_freq$mut_type <- rep(k,nrow_freq)
			
			data_genotype$mut_type <- rep(k,nrow_genotype)
						
			data_scaffold_freq <- rbind(data_scaffold_freq, data_freq)
			data_scaffold_genotype <- rbind(data_scaffold_genotype, data_genotype)
			
			data_tot_freq <- rbind(data_tot_freq, data_freq)
			data_tot_genotype <- rbind(data_tot_genotype, data_genotype)
		}
		
		data_scaffold_freq <- data_scaffold_freq[which(!(data_scaffold_freq$Freq==0 | data_scaffold_freq$Freq==1)),]
		
		write.table(data_scaffold_freq, paste(m,j,"_freq_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
		
		data_scaffold_genotype <- data_scaffold_genotype[which(data_scaffold_genotype$POS %in% data_scaffold_freq$POS),]
		
		write.table(data_scaffold_genotype, paste(m,j,"_genotype_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
	}
	
	data_tot_freq <- data_tot_freq[which(!(data_tot_freq$Freq==0 | data_tot_freq$Freq==1)),]
	
	write.table(data_tot_freq, paste(m,"_freq_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
}

####

species <- "Co"
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr")

data_tot_freq <- numeric()
data_tot_genotype <- numeric()

for (m in species)
{
	for (j in scaffold)
	{	
		data_scaffold_freq <- numeric()
		data_scaffold_genotype <- numeric()
		
		for (k in mut)
		{
			data_freq <- read.table(paste(m,j,"_",k,"_freq_2.txt",sep = ""), header=TRUE, sep="\t")
			
			data_genotype <- read.table(paste(m,j,"_",k,"_genotype.txt",sep = ""), header=TRUE, sep="\t")
			
			nrow_freq <- nrow(data_freq)
			
			nrow_genotype <- nrow(data_genotype)
		
			data_freq$mut_type <- rep(k,nrow_freq)
			
			data_genotype$mut_type <- rep(k,nrow_genotype)
			
			data_freq <- data_freq[which(!(data_freq$Freq==0 | data_freq$Freq==1)),]
		
			data_genotype <- data_genotype[which(data_genotype$POS %in% data_freq$POS),]
		
			data_sites_het <- unique(data_genotype[which(data_genotype$genotype == "Het"),]$POS)
		
			data_freq <- data_freq[which(!(data_freq$POS %in% data_sites_het)),]
		
			data_genotype <- data_genotype[which(!(data_genotype$POS %in% data_sites_het)),]
						
			data_scaffold_freq <- rbind(data_scaffold_freq, data_freq)
			data_scaffold_genotype <- rbind(data_scaffold_genotype, data_genotype)
			
			data_tot_freq <- rbind(data_tot_freq, data_freq)
			data_tot_genotype <- rbind(data_tot_genotype, data_genotype)
		}
		
		data_scaffold_freq <- data_scaffold_freq[which(!(data_scaffold_freq$Freq==0 | data_scaffold_freq$Freq==1)),]
		
		write.table(data_scaffold_freq, paste(m,j,"_freq_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
		
		data_scaffold_genotype <- data_scaffold_genotype[which(data_scaffold_genotype$POS %in% data_scaffold_freq$POS),]
		
		write.table(data_scaffold_genotype, paste(m,j,"_genotype_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
	}
	
	data_tot_freq <- data_tot_freq[which(!(data_tot_freq$Freq==0 | data_tot_freq$Freq==1)),]
	
	write.table(data_tot_freq, paste(m,"_freq_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
}


####

species <- "Co"
scaffold <- seq(1,8,1)
mut <- c("del","m_del")

data_tot_freq <- numeric()
data_tot_genotype <- numeric()

for (m in species)
{
	for (j in scaffold)
	{	
		data_scaffold_freq <- numeric()
		data_scaffold_genotype <- numeric()
		
		for (k in mut)
		{
			if (as.character(k)=="neutr")
			{
				data_freq <- read.table(paste(m,j,"_",k,"_freq_2.txt",sep = ""), header=TRUE, sep="\t")
			
				data_genotype <- read.table(paste(m,j,"_",k,"_genotype.txt",sep = ""), header=TRUE, sep="\t")
			} else {
				data_freq <- read.table(paste(m,j,"_",k,"_freq_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
			
				data_genotype <- read.table(paste(m,j,"_",k,"_genotype_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
			}
			
			nrow_freq <- nrow(data_freq)
			
			nrow_genotype <- nrow(data_genotype)
		
			data_freq$mut_type <- rep(k,nrow_freq)
			
			data_genotype$mut_type <- rep(k,nrow_genotype)
			
			data_freq <- data_freq[which(!(data_freq$Freq==0 | data_freq$Freq==1)),]
		
			data_genotype <- data_genotype[which(data_genotype$POS %in% data_freq$POS),]
		
			data_sites_het <- unique(data_genotype[which(data_genotype$genotype == "Het"),]$POS)
		
			data_freq <- data_freq[which(!(data_freq$POS %in% data_sites_het)),]
		
			data_genotype <- data_genotype[which(!(data_genotype$POS %in% data_sites_het)),]
						
			data_scaffold_freq <- rbind(data_scaffold_freq, data_freq)
			data_scaffold_genotype <- rbind(data_scaffold_genotype, data_genotype)
			
			data_tot_freq <- rbind(data_tot_freq, data_freq)
			data_tot_genotype <- rbind(data_tot_genotype, data_genotype)
		}
		
		write.table(data_scaffold_freq, paste(m,j,"_freq_SIFT_or_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
		
		write.table(data_scaffold_genotype, paste(m,j,"_genotype_SIFT_or_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
	}
		
	write.table(data_tot_freq, paste(m,"_freq_SIFT_or_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE)
}

####

species <- c("Cg","Co")
scaffold <- seq(1,8,1)

for (i in species)
{
	if (as.character(i) == "Co")
	{
		exclud_ind <- c(4,9,11,14,23,29,30,31)
	} else {
		exclud_ind <- c(22,45,46,52,53,64,88,90,109,111,146)
	}
	
	for (j in scaffold)
	{
		annotated <- read.table(paste(i,j,"_annotated.txt",sep = ""), header=TRUE, sep="\t")
		names <- colnames(annotated)
		sub_annotated <- annotated[which(is.na(annotated$nucl) == TRUE | as.numeric(gsub(".*?([0-9]+).*", "\\1",annotated$ind_chr)) %in% exclud_ind),]
		colnames(sub_annotated) <- names
		write.table(sub_annotated, paste(i,j,"_NA.txt",sep = ""), quote = FALSE, sep = "\t")
	}
}

####

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr")

dist_max <- 1500

window <- 60

size <- dist_max/window

for (i in species)
{
	if (as.character(i)=="Co")
	{
		num_ind <- 33
		ind_excl <- c(4,9,11,14,23,29,30,31)
	} else
	{
		num_ind <- 182
		ind_excl <- c(22,45,46,52,53,64,88,90,109,111,146)
	}
	
	list_ind <- seq(1,num_ind,1)[-ind_excl]
	
	min_ind <- length(list_ind)/2
	
	for (j in scaffold)
	{
		data_NA <- read.table(paste(i,j,"_NA.txt",sep = ""), header=TRUE, sep="\t")
		
		for (k in mut)
		{
			genotype <- read.table(paste(i,j,"_",k,"_genotype.txt",sep = ""), header=TRUE, sep="\t")
			freq <- na.omit(read.table(paste(i,j,"_",k,"_freq_2.txt",sep = ""), header=TRUE, sep="\t"))

			n_sites <- nrow(freq)

			v_Delta_AB <- numeric()
			v_Delta_AB_sqrd <- numeric()
			v_D_AB_A <- numeric()
			v_D_AB_B <- numeric()
			v_Delta_AB_A <- numeric()
			v_D_AB_AB <- numeric()
			v_p_A_p_B_Delta_AB <- numeric()
			v_p_A_Delta_AB <- numeric()
			v_p_B_Delta_AB <- numeric()
			v_p_A_plus_p_B_Delta_AB <- numeric()
			v_pq_AB <- numeric()
			v_p_A <- numeric()
			v_p_B <- numeric()
			v_dist <- numeric()
			v_dist_class <- numeric()
			v_pos_A <- numeric()
			v_pos_B <- numeric()
			v_num_ind <- numeric()
			
			headers <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
						D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
						p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
						Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind)
						
			write.table(headers, paste(i,j,"_",k,"_Delta_AB.txt",sep = ""), quote = FALSE, sep = "\t")
				
				for (n in c(1:n_sites))
				{
					if (freq$Freq[n] == 0 | freq$Freq[n] == 1)
					{
						next
					}
					
					sub_class <- freq[which( (freq$POS>freq$POS[n]) & (freq$POS<= (freq$POS[n] + dist_max)) ), ]
						
					n_sites_sub_class <- nrow(sub_class)
						
					if(n_sites_sub_class == 0)
					{
						next
					}
						
					NA_A <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq$POS[n]),]$ind_chr)))
					
					ind_A <- list_ind[!list_ind %in% NA_A]
					
					sub_A <- genotype[which(genotype$POS == freq$POS[n]), ]
					
					p_A <- freq$Freq[n]
							
					for (l in c((n+1):(n+n_sites_sub_class)))
					{		
						if (freq$Freq[l] == 0 | freq$Freq[l] == 1)
						{
							next
						}
						
						distance <- freq$POS[l] - freq$POS[n]
							
						if (distance == 0 | distance > dist_max)
						{
							next
						}
						
						NA_B <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq$POS[l]),]$ind_chr)))
					
						ind_B <- list_ind[!list_ind %in% NA_B]
						
						comm_ind <- intersect(ind_A, ind_B)
						
						if (length(comm_ind) <= min_ind)
						{
							next
						}
						 
						sub_B <- genotype[which(genotype$POS == freq$POS[l]), ]
						
						p_B <- freq$Freq[l]
							
						v_Delta_AB_ind <- numeric()
						v_D_AB_A_ind <- numeric()
						v_D_AB_B_ind <- numeric()
						v_Delta_AB_A_ind <- numeric()
						v_D_AB_AB_ind <- numeric()	
							
						for (m in comm_ind)
						{
							sum_X_A <- 0
							sum_X_B <- 0
							prod_X_A <- 0
							prod_X_B <- 0
					
							if (!(m %in% sub_A$ind))
							{
								sum_X_A <- 0
								prod_X_A <- 0
							}else
							{
								if (m %in% sub_A[which(as.character(sub_A$genotype) == "Het"),]$ind)
								{
									sum_X_A <- 1
									prod_X_A <- 0
								}
								else
								{
									sum_X_A <- 2
									prod_X_A <- 1
								}
							}
					
							if (!(m %in% sub_B$ind))
							{
								sum_X_B <- 0
								prod_X_B <- 0
							}else 
							{
								if (m %in% sub_B[which(as.character(sub_B$genotype) == "Het"),]$ind)
								{
									sum_X_B <- 1
									prod_X_B <- 0
								}
								else
								{
									sum_X_B <- 2
									prod_X_B <- 1
								}
							}
							
							v_Delta_AB_ind <- c(v_Delta_AB_ind, (sum_X_A - 2*p_A)*(sum_X_B - 2*p_B))
							v_D_AB_A_ind <- c(v_D_AB_A_ind, p_A^2*sum_X_B - 2*p_A^2*p_B + 2*p_A*p_B*sum_X_A - p_A*sum_X_A*sum_X_B + prod_X_A*(-2*p_B+sum_X_B))
							v_D_AB_B_ind <- c(v_D_AB_B_ind, p_B^2*sum_X_A - 2*p_B^2*p_A + 2*p_A*p_B*sum_X_B - p_B*sum_X_A*sum_X_B + prod_X_B*(-2*p_A+sum_X_A))
							v_Delta_AB_A_ind <- c(v_Delta_AB_A_ind, prod_X_A*sum_X_B + prod_X_B*sum_X_A - p_A*(sum_X_A*sum_X_B+2*prod_X_B)
							- p_B*(sum_X_A*sum_X_B+2*prod_X_A) + 2*p_A*p_B*(sum_X_A+sum_X_B-p_A-p_B) + p_A^2*sum_X_B+p_B^2*sum_X_A)
							v_D_AB_AB_ind <- c(v_D_AB_AB_ind, p_A^2*p_B^2 - p_A*p_B^2*sum_X_A - p_A^2*p_B*sum_X_B + p_A*p_B*sum_X_A*sum_X_B
							+ prod_X_A*(p_B^2-p_B*sum_X_B) + prod_X_B*(p_A^2-p_A*sum_X_A) + prod_X_A*prod_X_B)
						}
							
						v_pq_AB <- p_A*(1-p_A)*p_B*(1-p_B)
						v_p_A <- p_A
						v_p_B <- p_B
						v_Delta_AB <- mean(v_Delta_AB_ind)/2
						v_Delta_AB_sqrd <- (mean(v_Delta_AB_ind)/2)^2
						v_D_AB_A <- mean(v_D_AB_A_ind)/2
						v_D_AB_B <- mean(v_D_AB_B_ind)/2
						v_Delta_AB_A <- mean(v_Delta_AB_A_ind)/2
						v_D_AB_AB <- mean(v_D_AB_AB_ind)
						v_p_A_p_B_Delta_AB <- (1-2*p_A)*(1-2*p_B)*(mean(v_Delta_AB_ind)/2)
						v_p_A_Delta_AB <- p_A*(mean(v_Delta_AB_ind)/2)
						v_p_B_Delta_AB <- p_B*(mean(v_Delta_AB_ind)/2)
						v_p_A_plus_p_B_Delta_AB <- (p_A+p_B)*(mean(v_Delta_AB_ind)/2)
						v_dist <- distance
						v_dist_class <- (floor(distance/size)*size+floor(distance/size)*size+size)/2
						v_pos_A <- freq$POS[n]
						v_pos_B <- freq$POS[l]
						
						line <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
						D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
						p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
						Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = length(comm_ind))
						
						write.table(line, paste(i,j,"_",k,"_Delta_AB.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
					}
				}
		}
		
		genotype_del <- read.table(paste(i,j,"_del_genotype.txt",sep = ""), header=TRUE, sep="\t")
		freq_del <- read.table(paste(i,j,"_del_freq_2.txt",sep = ""), header=TRUE, sep="\t")
		
		genotype_neutr <- read.table(paste(i,j,"_neutr_genotype.txt",sep = ""), header=TRUE, sep="\t")
		freq_neutr <- read.table(paste(i,j,"_neutr_freq_2.txt",sep = ""), header=TRUE, sep="\t")

		n_sites_del <- nrow(freq_del)
		n_sites_neutr <- nrow(freq_neutr)

		v_Delta_AB <- numeric()
		v_Delta_AB_sqrd <- numeric()
		v_D_AB_A <- numeric()
		v_D_AB_B <- numeric()
		v_Delta_AB_A <- numeric()
		v_D_AB_AB <- numeric()
		v_p_A_p_B_Delta_AB <- numeric()
		v_p_A_Delta_AB <- numeric()
		v_p_B_Delta_AB <- numeric()
		v_p_A_plus_p_B_Delta_AB <- numeric()
		v_pq_AB <- numeric()
		v_p_A <- numeric()
		v_p_B <- numeric()
		v_dist <- numeric()
		v_dist_class <- numeric()
		v_pos_A <- numeric()
		v_pos_B <- numeric()
		v_nucl_mut <- numeric()
		v_num_ind <- numeric()
			
		headers <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
					D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
					p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind)
						
		write.table(headers, paste(i,j,"_del_neutr_Delta_AB.txt",sep = ""), quote = FALSE, sep = "\t")
				
			for (n in c(1:n_sites_del))
			{
				if (freq_del$Freq[n] == 0 | freq_del$Freq[n] == 1)
				{
					next
				}
				
				sub_class_neutr <- freq_neutr[which( (freq_neutr$POS>=freq_del$POS[n] - dist_max) & (freq_neutr$POS<= (freq_del$POS[n] + dist_max)) ), ]
						
				n_sites_sub_class_neutr <- nrow(sub_class_neutr)
						
				if(n_sites_sub_class_neutr == 0)
				{
					next
				}
				
				NA_A <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq_del$POS[n]),]$ind_chr)))
					
				ind_A <- list_ind[!list_ind %in% NA_A]
				
				sub_A <- genotype_del[which(genotype_del$POS == freq_del$POS[n]), ] 
				
				p_A <- freq_del$Freq[n]
							
				for (l in c(1:n_sites_sub_class_neutr))
				{
					if (sub_class_neutr$Freq[l] == 0 | sub_class_neutr$Freq[l] == 1)
					{
						next
					}
					
					distance <- abs(sub_class_neutr$POS[l] - freq_del$POS[n])
							
					if (distance == 0 | distance > dist_max)
					{
						next
					}
						
					NA_B <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == sub_class_neutr$POS[l]),]$ind_chr)))
					
					ind_B <- list_ind[!list_ind %in% NA_B]
					
					comm_ind <- intersect(ind_A, ind_B)
					
					if (length(comm_ind) <= min_ind)
					{
						next
					}
					
					sub_B <- genotype_neutr[which(genotype_neutr$POS == sub_class_neutr$POS[l]), ]
						
					p_B <- sub_class_neutr$Freq[l]
							
					v_Delta_AB_ind <- numeric()
					v_D_AB_A_ind <- numeric()
					v_D_AB_B_ind <- numeric()
					v_Delta_AB_A_ind <- numeric()
					v_D_AB_AB_ind <- numeric()	
							
					for (m in comm_ind)
					{
						sum_X_A <- 0
						sum_X_B <- 0
						prod_X_A <- 0
						prod_X_B <- 0
					
						if (!(m %in% sub_A$ind))
						{
							sum_X_A <- 0
							prod_X_A <- 0
						}else
						{
							if (m %in% sub_A[which(as.character(sub_A$genotype) == "Het"),]$ind)
							{
								sum_X_A <- 1
								prod_X_A <- 0
							}
							else
							{
								sum_X_A <- 2
								prod_X_A <- 1
							}
						}
					
						if (!(m %in% sub_B$ind))
						{
							sum_X_B <- 0
							prod_X_B <- 0
						}else 
						{
							if (m %in% sub_B[which(as.character(sub_B$genotype) == "Het"),]$ind)
							{
								sum_X_B <- 1
								prod_X_B <- 0
							}
							else
							{
								sum_X_B <- 2
								prod_X_B <- 1
							}
						}
							
						v_Delta_AB_ind <- c(v_Delta_AB_ind, (sum_X_A - 2*p_A)*(sum_X_B - 2*p_B))
						v_D_AB_A_ind <- c(v_D_AB_A_ind, p_A^2*sum_X_B - 2*p_A^2*p_B + 2*p_A*p_B*sum_X_A - p_A*sum_X_A*sum_X_B + prod_X_A*(-2*p_B+sum_X_B))
						v_D_AB_B_ind <- c(v_D_AB_B_ind, p_B^2*sum_X_A - 2*p_B^2*p_A + 2*p_A*p_B*sum_X_B - p_B*sum_X_A*sum_X_B + prod_X_B*(-2*p_A+sum_X_A))
						v_Delta_AB_A_ind <- c(v_Delta_AB_A_ind, prod_X_A*sum_X_B + prod_X_B*sum_X_A - p_A*(sum_X_A*sum_X_B+2*prod_X_B)
						- p_B*(sum_X_A*sum_X_B+2*prod_X_A) + 2*p_A*p_B*(sum_X_A+sum_X_B-p_A-p_B) + p_A^2*sum_X_B+p_B^2*sum_X_A)
						v_D_AB_AB_ind <- c(v_D_AB_AB_ind, p_A^2*p_B^2 - p_A*p_B^2*sum_X_A - p_A^2*p_B*sum_X_B + p_A*p_B*sum_X_A*sum_X_B
						+ prod_X_A*(p_B^2-p_B*sum_X_B) + prod_X_B*(p_A^2-p_A*sum_X_A) + prod_X_A*prod_X_B)
					}
							
					v_pq_AB <- p_A*(1-p_A)*p_B*(1-p_B)
					v_p_A <- p_A
					v_p_B <- p_B
					v_Delta_AB <- mean(v_Delta_AB_ind)/2
					v_Delta_AB_sqrd <- (mean(v_Delta_AB_ind)/2)^2
					v_D_AB_A <- mean(v_D_AB_A_ind)/2
					v_D_AB_B <- mean(v_D_AB_B_ind)/2
					v_Delta_AB_A <- mean(v_Delta_AB_A_ind)/2
					v_D_AB_AB <- mean(v_D_AB_AB_ind)
					v_p_A_p_B_Delta_AB <- (1-2*p_A)*(1-2*p_B)*(mean(v_Delta_AB_ind)/2)
					v_p_A_Delta_AB <- p_A*(mean(v_Delta_AB_ind)/2)
					v_p_B_Delta_AB <- p_B*(mean(v_Delta_AB_ind)/2)
					v_p_A_plus_p_B_Delta_AB <- (p_A+p_B)*(mean(v_Delta_AB_ind)/2)
					v_dist <- distance
					v_dist_class <- (floor(distance/size)*size+floor(distance/size)*size+size)/2
					v_pos_A <- sub_A$POS[1]
					v_pos_B <- sub_B$POS[1]
						
					line <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
					D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
					p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = length(comm_ind))
						
					write.table(line, paste(i,j,"_del_neutr_Delta_AB.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
				}
			}
	}
}


####

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del")

dist_max <- 1500

window <- 60

size <- dist_max/window

for (i in species)
{
	if (as.character(i)=="Co")
	{
		num_ind <- 33
		ind_excl <- c(4,9,11,14,23,29,30,31)
	} else
	{
		num_ind <- 182
		ind_excl <- c(22,45,46,52,53,64,88,90,109,111,146)
	}
	
	list_ind <- seq(1,num_ind,1)[-ind_excl]
	
	min_ind <- length(list_ind)/2
	
	for (j in scaffold)
	{
		data_NA <- read.table(paste(i,j,"_NA.txt",sep = ""), header=TRUE, sep="\t")
		
		for (k in mut)
		{
			genotype <- read.table(paste(i,j,"_",k,"_genotype_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
			freq <- na.omit(read.table(paste(i,j,"_",k,"_freq_SIFT_or.txt",sep = ""), header=TRUE, sep="\t"))

			n_sites <- nrow(freq)

			v_Delta_AB <- numeric()
			v_Delta_AB_sqrd <- numeric()
			v_D_AB_A <- numeric()
			v_D_AB_B <- numeric()
			v_Delta_AB_A <- numeric()
			v_D_AB_AB <- numeric()
			v_p_A_p_B_Delta_AB <- numeric()
			v_p_A_Delta_AB <- numeric()
			v_p_B_Delta_AB <- numeric()
			v_p_A_plus_p_B_Delta_AB <- numeric()
			v_pq_AB <- numeric()
			v_p_A <- numeric()
			v_p_B <- numeric()
			v_dist <- numeric()
			v_dist_class <- numeric()
			v_pos_A <- numeric()
			v_pos_B <- numeric()
			v_num_ind <- numeric()
			
			headers <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
						D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
						p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
						Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind)
						
			write.table(headers, paste(i,j,"_",k,"_Delta_AB_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t")
				
				for (n in c(1:n_sites))
				{
					if (freq$Freq[n] == 0 | freq$Freq[n] == 1)
					{
						next
					}
					
					sub_class <- freq[which( (freq$POS>freq$POS[n]) & (freq$POS<= (freq$POS[n] + dist_max)) ), ]
						
					n_sites_sub_class <- nrow(sub_class)
						
					if(n_sites_sub_class == 0)
					{
						next
					}
						
					NA_A <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq$POS[n]),]$ind_chr)))
					
					ind_A <- list_ind[!list_ind %in% NA_A]
					
					sub_A <- genotype[which(genotype$POS == freq$POS[n]), ]
					
					p_A <- freq$Freq[n]
							
					for (l in c((n+1):(n+n_sites_sub_class)))
					{		
						if (freq$Freq[l] == 0 | freq$Freq[l] == 1)
						{
							next
						}
						
						distance <- freq$POS[l] - freq$POS[n]
							
						if (distance == 0 | distance > dist_max)
						{
							next
						}
						
						NA_B <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq$POS[l]),]$ind_chr)))
					
						ind_B <- list_ind[!list_ind %in% NA_B]
						
						comm_ind <- intersect(ind_A, ind_B)
						
						if (length(comm_ind) <= min_ind)
						{
							next
						}
						 
						sub_B <- genotype[which(genotype$POS == freq$POS[l]), ]
						
						p_B <- freq$Freq[l]
							
						v_Delta_AB_ind <- numeric()
						v_D_AB_A_ind <- numeric()
						v_D_AB_B_ind <- numeric()
						v_Delta_AB_A_ind <- numeric()
						v_D_AB_AB_ind <- numeric()	
							
						for (m in comm_ind)
						{
							sum_X_A <- 0
							sum_X_B <- 0
							prod_X_A <- 0
							prod_X_B <- 0
					
							if (!(m %in% sub_A$ind))
							{
								sum_X_A <- 0
								prod_X_A <- 0
							}else
							{
								if (m %in% sub_A[which(as.character(sub_A$genotype) == "Het"),]$ind)
								{
									sum_X_A <- 1
									prod_X_A <- 0
								}
								else
								{
									sum_X_A <- 2
									prod_X_A <- 1
								}
							}
					
							if (!(m %in% sub_B$ind))
							{
								sum_X_B <- 0
								prod_X_B <- 0
							}else 
							{
								if (m %in% sub_B[which(as.character(sub_B$genotype) == "Het"),]$ind)
								{
									sum_X_B <- 1
									prod_X_B <- 0
								}
								else
								{
									sum_X_B <- 2
									prod_X_B <- 1
								}
							}
							
							v_Delta_AB_ind <- c(v_Delta_AB_ind, (sum_X_A - 2*p_A)*(sum_X_B - 2*p_B))
							v_D_AB_A_ind <- c(v_D_AB_A_ind, p_A^2*sum_X_B - 2*p_A^2*p_B + 2*p_A*p_B*sum_X_A - p_A*sum_X_A*sum_X_B + prod_X_A*(-2*p_B+sum_X_B))
							v_D_AB_B_ind <- c(v_D_AB_B_ind, p_B^2*sum_X_A - 2*p_B^2*p_A + 2*p_A*p_B*sum_X_B - p_B*sum_X_A*sum_X_B + prod_X_B*(-2*p_A+sum_X_A))
							v_Delta_AB_A_ind <- c(v_Delta_AB_A_ind, prod_X_A*sum_X_B + prod_X_B*sum_X_A - p_A*(sum_X_A*sum_X_B+2*prod_X_B)
							- p_B*(sum_X_A*sum_X_B+2*prod_X_A) + 2*p_A*p_B*(sum_X_A+sum_X_B-p_A-p_B) + p_A^2*sum_X_B+p_B^2*sum_X_A)
							v_D_AB_AB_ind <- c(v_D_AB_AB_ind, p_A^2*p_B^2 - p_A*p_B^2*sum_X_A - p_A^2*p_B*sum_X_B + p_A*p_B*sum_X_A*sum_X_B
							+ prod_X_A*(p_B^2-p_B*sum_X_B) + prod_X_B*(p_A^2-p_A*sum_X_A) + prod_X_A*prod_X_B)
						}
							
						v_pq_AB <- p_A*(1-p_A)*p_B*(1-p_B)
						v_p_A <- p_A
						v_p_B <- p_B
						v_Delta_AB <- mean(v_Delta_AB_ind)/2
						v_Delta_AB_sqrd <- (mean(v_Delta_AB_ind)/2)^2
						v_D_AB_A <- mean(v_D_AB_A_ind)/2
						v_D_AB_B <- mean(v_D_AB_B_ind)/2
						v_Delta_AB_A <- mean(v_Delta_AB_A_ind)/2
						v_D_AB_AB <- mean(v_D_AB_AB_ind)
						v_p_A_p_B_Delta_AB <- (1-2*p_A)*(1-2*p_B)*(mean(v_Delta_AB_ind)/2)
						v_p_A_Delta_AB <- p_A*(mean(v_Delta_AB_ind)/2)
						v_p_B_Delta_AB <- p_B*(mean(v_Delta_AB_ind)/2)
						v_p_A_plus_p_B_Delta_AB <- (p_A+p_B)*(mean(v_Delta_AB_ind)/2)
						v_dist <- distance
						v_dist_class <- (floor(distance/size)*size+floor(distance/size)*size+size)/2
						v_pos_A <- freq$POS[n]
						v_pos_B <- freq$POS[l]
						
						line <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
						D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
						p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
						Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = length(comm_ind))
						
						write.table(line, paste(i,j,"_",k,"_Delta_AB_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
					}
				}
		}
		
		genotype_del <- read.table(paste(i,j,"_del_genotype_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
		freq_del <- read.table(paste(i,j,"_del_freq_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
		
		genotype_neutr <- read.table(paste(i,j,"_neutr_genotype.txt",sep = ""), header=TRUE, sep="\t")
		freq_neutr <- read.table(paste(i,j,"_neutr_freq_2.txt",sep = ""), header=TRUE, sep="\t")

		n_sites_del <- nrow(freq_del)
		n_sites_neutr <- nrow(freq_neutr)

		v_Delta_AB <- numeric()
		v_Delta_AB_sqrd <- numeric()
		v_D_AB_A <- numeric()
		v_D_AB_B <- numeric()
		v_Delta_AB_A <- numeric()
		v_D_AB_AB <- numeric()
		v_p_A_p_B_Delta_AB <- numeric()
		v_p_A_Delta_AB <- numeric()
		v_p_B_Delta_AB <- numeric()
		v_p_A_plus_p_B_Delta_AB <- numeric()
		v_pq_AB <- numeric()
		v_p_A <- numeric()
		v_p_B <- numeric()
		v_dist <- numeric()
		v_dist_class <- numeric()
		v_pos_A <- numeric()
		v_pos_B <- numeric()
		v_nucl_mut <- numeric()
		v_num_ind <- numeric()
			
		headers <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
					D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
					p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind)
						
		write.table(headers, paste(i,j,"_del_neutr_Delta_AB_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t")
				
			for (n in c(1:n_sites_del))
			{
				if (freq_del$Freq[n] == 0 | freq_del$Freq[n] == 1)
				{
					next
				}
				
				sub_class_neutr <- freq_neutr[which( (freq_neutr$POS>=freq_del$POS[n] - dist_max) & (freq_neutr$POS<= (freq_del$POS[n] + dist_max)) ), ]
						
				n_sites_sub_class_neutr <- nrow(sub_class_neutr)
						
				if(n_sites_sub_class_neutr == 0)
				{
					next
				}
				
				NA_A <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq_del$POS[n]),]$ind_chr)))
					
				ind_A <- list_ind[!list_ind %in% NA_A]
				
				sub_A <- genotype_del[which(genotype_del$POS == freq_del$POS[n]), ] 
				
				p_A <- freq_del$Freq[n]
							
				for (l in c(1:n_sites_sub_class_neutr))
				{
					if (sub_class_neutr$Freq[l] == 0 | sub_class_neutr$Freq[l] == 1)
					{
						next
					}
					
					distance <- abs(sub_class_neutr$POS[l] - freq_del$POS[n])
							
					if (distance == 0 | distance > dist_max)
					{
						next
					}
						
					NA_B <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == sub_class_neutr$POS[l]),]$ind_chr)))
					
					ind_B <- list_ind[!list_ind %in% NA_B]
					
					comm_ind <- intersect(ind_A, ind_B)
					
					if (length(comm_ind) <= min_ind)
					{
						next
					}
					
					sub_B <- genotype_neutr[which(genotype_neutr$POS == sub_class_neutr$POS[l]), ]
						
					p_B <- sub_class_neutr$Freq[l]
							
					v_Delta_AB_ind <- numeric()
					v_D_AB_A_ind <- numeric()
					v_D_AB_B_ind <- numeric()
					v_Delta_AB_A_ind <- numeric()
					v_D_AB_AB_ind <- numeric()	
							
					for (m in comm_ind)
					{
						sum_X_A <- 0
						sum_X_B <- 0
						prod_X_A <- 0
						prod_X_B <- 0
					
						if (!(m %in% sub_A$ind))
						{
							sum_X_A <- 0
							prod_X_A <- 0
						}else
						{
							if (m %in% sub_A[which(as.character(sub_A$genotype) == "Het"),]$ind)
							{
								sum_X_A <- 1
								prod_X_A <- 0
							}
							else
							{
								sum_X_A <- 2
								prod_X_A <- 1
							}
						}
					
						if (!(m %in% sub_B$ind))
						{
							sum_X_B <- 0
							prod_X_B <- 0
						}else 
						{
							if (m %in% sub_B[which(as.character(sub_B$genotype) == "Het"),]$ind)
							{
								sum_X_B <- 1
								prod_X_B <- 0
							}
							else
							{
								sum_X_B <- 2
								prod_X_B <- 1
							}
						}
							
						v_Delta_AB_ind <- c(v_Delta_AB_ind, (sum_X_A - 2*p_A)*(sum_X_B - 2*p_B))
						v_D_AB_A_ind <- c(v_D_AB_A_ind, p_A^2*sum_X_B - 2*p_A^2*p_B + 2*p_A*p_B*sum_X_A - p_A*sum_X_A*sum_X_B + prod_X_A*(-2*p_B+sum_X_B))
						v_D_AB_B_ind <- c(v_D_AB_B_ind, p_B^2*sum_X_A - 2*p_B^2*p_A + 2*p_A*p_B*sum_X_B - p_B*sum_X_A*sum_X_B + prod_X_B*(-2*p_A+sum_X_A))
						v_Delta_AB_A_ind <- c(v_Delta_AB_A_ind, prod_X_A*sum_X_B + prod_X_B*sum_X_A - p_A*(sum_X_A*sum_X_B+2*prod_X_B)
						- p_B*(sum_X_A*sum_X_B+2*prod_X_A) + 2*p_A*p_B*(sum_X_A+sum_X_B-p_A-p_B) + p_A^2*sum_X_B+p_B^2*sum_X_A)
						v_D_AB_AB_ind <- c(v_D_AB_AB_ind, p_A^2*p_B^2 - p_A*p_B^2*sum_X_A - p_A^2*p_B*sum_X_B + p_A*p_B*sum_X_A*sum_X_B
						+ prod_X_A*(p_B^2-p_B*sum_X_B) + prod_X_B*(p_A^2-p_A*sum_X_A) + prod_X_A*prod_X_B)
					}
							
					v_pq_AB <- p_A*(1-p_A)*p_B*(1-p_B)
					v_p_A <- p_A
					v_p_B <- p_B
					v_Delta_AB <- mean(v_Delta_AB_ind)/2
					v_Delta_AB_sqrd <- (mean(v_Delta_AB_ind)/2)^2
					v_D_AB_A <- mean(v_D_AB_A_ind)/2
					v_D_AB_B <- mean(v_D_AB_B_ind)/2
					v_Delta_AB_A <- mean(v_Delta_AB_A_ind)/2
					v_D_AB_AB <- mean(v_D_AB_AB_ind)
					v_p_A_p_B_Delta_AB <- (1-2*p_A)*(1-2*p_B)*(mean(v_Delta_AB_ind)/2)
					v_p_A_Delta_AB <- p_A*(mean(v_Delta_AB_ind)/2)
					v_p_B_Delta_AB <- p_B*(mean(v_Delta_AB_ind)/2)
					v_p_A_plus_p_B_Delta_AB <- (p_A+p_B)*(mean(v_Delta_AB_ind)/2)
					v_dist <- distance
					v_dist_class <- (floor(distance/size)*size+floor(distance/size)*size+size)/2
					v_pos_A <- sub_A$POS[1]
					v_pos_B <- sub_B$POS[1]
						
					line <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
					D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
					p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = length(comm_ind))
						
					write.table(line, paste(i,j,"_del_neutr_Delta_AB_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
				}
			}
	}
}

####

species <- "Co"
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr")

dist_max <- 3000

window <- 60

size <- dist_max/window

for (i in species)
{
	if (as.character(i)=="Co")
	{
		num_ind <- 33
		ind_excl <- c(4,9,11,14,23,29,30,31)
	} else
	{
		num_ind <- 182
		ind_excl <- c(22,45,46,52,53,64,88,90,109,111,146)
	}
	
	list_ind <- seq(1,num_ind,1)[-ind_excl]
	
	min_ind <- length(list_ind)/2
	
	for (j in scaffold)
	{
		data_NA <- read.table(paste(i,j,"_NA.txt",sep = ""), header=TRUE, sep="\t")
		
		genotype_scaffold <- read.table(paste(i,j,"_genotype_no_het.txt",sep = ""), header=TRUE, sep="\t")
		freq_scaffold <- na.omit(read.table(paste(i,j,"_freq_no_het.txt",sep = ""), header=TRUE, sep="\t"))
		
		for (k in mut)
		{
			genotype <- genotype_scaffold[which(as.character(genotype_scaffold$mut_type)==k),]
			freq <- freq_scaffold[which(as.character(freq_scaffold$mut_type)==k),]

			n_sites <- nrow(freq)

			v_D_AB <- numeric()
			v_D_AB_sqrd <- numeric()
			v_p_A_p_B_D_AB <- numeric()
			v_p_A_D_AB <- numeric()
			v_p_B_D_AB <- numeric()
			v_p_A_plus_p_B_D_AB <- numeric()
			v_pq_AB <- numeric()
			v_p_A <- numeric()
			v_p_B <- numeric()
			v_dist <- numeric()
			v_dist_class <- numeric()
			v_pos_A <- numeric()
			v_pos_B <- numeric()
			v_num_ind <- numeric()
			
			headers <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd, p_A_p_B_D_AB = v_p_A_p_B_D_AB, p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
						p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
						Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind)
						
			write.table(headers, paste(i,j,"_",k,"_D_AB_no_het.txt",sep = ""), quote = FALSE, sep = "\t")
				
				for (n in c(1:n_sites))
				{
					if (freq$Freq[n] == 0 | freq$Freq[n] == 1)
					{
						next
					}
					
					sub_class <- freq[which( (freq$POS>freq$POS[n]) & (freq$POS<= (freq$POS[n] + dist_max)) ), ]
						
					n_sites_sub_class <- nrow(sub_class)
						
					if(n_sites_sub_class == 0)
					{
						next
					}
						
					NA_A <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq$POS[n]),]$ind_chr)))
					
					ind_A <- list_ind[!list_ind %in% NA_A]
					
					sub_A <- genotype[which(genotype$POS == freq$POS[n]), ]
					
					p_A <- freq$Freq[n]
							
					for (l in c((n+1):(n+n_sites_sub_class)))
					{		
						if (freq$Freq[l] == 0 | freq$Freq[l] == 1)
						{
							next
						}
						
						distance <- freq$POS[l] - freq$POS[n]
							
						if (distance == 0 | distance > dist_max)
						{
							next
						}
						
						NA_B <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq$POS[l]),]$ind_chr)))
					
						ind_B <- list_ind[!list_ind %in% NA_B]
						
						comm_ind <- intersect(ind_A, ind_B)
						
						if (length(comm_ind) <= min_ind)
						{
							next
						}
						 
						sub_B <- genotype[which(genotype$POS == freq$POS[l]), ]
						
						p_B <- freq$Freq[l]
							
						v_D_AB_ind <- numeric()
							
						for (m in comm_ind)
						{
							X_A <- 0
							X_B <- 0

							if (m %in% sub_A$ind)
							{
								X_A <- 1
							} else {
								X_A <- 0
							}
						
							if (m %in% sub_B$ind)
							{
								X_B <- 1
							} else {
								X_B <- 0
							}
							
							v_D_AB_ind <- c(v_D_AB_ind, (X_A - p_A)*(X_B - p_B))
						}
							
						v_pq_AB <- p_A*(1-p_A)*p_B*(1-p_B)
						v_p_A <- p_A
						v_p_B <- p_B
						v_D_AB <- mean(v_D_AB_ind)
						v_D_AB_sqrd <- (mean(v_D_AB_ind))^2
						v_p_A_p_B_D_AB <- (1-2*p_A)*(1-2*p_B)*(mean(v_D_AB_ind))
						v_p_A_D_AB <- p_A*(mean(v_D_AB_ind))
						v_p_B_D_AB <- p_B*(mean(v_D_AB_ind))
						v_p_A_plus_p_B_D_AB <- (p_A+p_B)*(mean(v_D_AB_ind))
						v_dist <- distance
						v_dist_class <- (floor(distance/size)*size+floor(distance/size)*size+size)/2
						v_pos_A <- freq$POS[n]
						v_pos_B <- freq$POS[l]
						
						line <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd,
						p_A_p_B_D_AB = v_p_A_p_B_D_AB, p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
						p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
						Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = length(comm_ind))
						
						write.table(line, paste(i,j,"_",k,"_D_AB_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
					}
				}
		}
		
			genotype_scaffold <- read.table(paste(i,j,"_genotype_no_het.txt",sep = ""), header=TRUE, sep="\t")
			freq_scaffold <- na.omit(read.table(paste(i,j,"_freq_no_het.txt",sep = ""), header=TRUE, sep="\t"))
			
			genotype_del <- genotype_scaffold[which(as.character(genotype_scaffold$mut_type)=="del"),]
			freq_del <- freq_scaffold[which(as.character(freq_scaffold$mut_type)=="del"),]
		
			genotype_neutr <- genotype_scaffold[which(as.character(genotype_scaffold$mut_type)=="neutr"),]
			freq_neutr <- freq_scaffold[which(as.character(freq_scaffold$mut_type)=="neutr"),]
			
			list_ind <- sort(unique(genotype_del$ind))

			n_sites_del <- nrow(freq_del)

			v_D_AB <- numeric()
			v_D_AB_sqrd <- numeric()
			v_D_AB_A <- numeric()
			v_D_AB_B <- numeric()
			v_D_AB_A <- numeric()
			v_D_AB_AB <- numeric()
			v_p_A_p_B_D_AB <- numeric()
			v_p_A_D_AB <- numeric()
			v_p_B_D_AB <- numeric()
			v_p_A_plus_p_B_D_AB <- numeric()
			v_pq_AB <- numeric()
			v_p_A <- numeric()
			v_p_B <- numeric()
			v_dist <- numeric()
			v_dist_class <- numeric()
			v_pos_A <- numeric()
			v_pos_B <- numeric()
			v_nucl_mut <- numeric()
			v_num_ind <- numeric()
			
			headers <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd,
					p_A_p_B_D_AB = v_p_A_p_B_D_AB, p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
					p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind)
						
			write.table(headers, paste(i,j,"_del_neutr_D_AB_no_het.txt",sep = ""), quote = FALSE, sep = "\t")
				
			for (n in c(1:n_sites_del))
			{	
				sub_class_neutr <- freq_neutr[which( (freq_neutr$POS>=freq_del$POS[n] - dist_max) & (freq_neutr$POS<= (freq_del$POS[n] + dist_max)) ), ]
						
				n_sites_sub_class_neutr <- nrow(sub_class_neutr)
						
				if(n_sites_sub_class_neutr == 0)
				{
					next
				}
				
				NA_A <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq_del$POS[n]),]$ind_chr)))
					
				ind_A <- list_ind[!list_ind %in% NA_A]
				
				sub_A <- genotype_del[which(genotype_del$POS == freq_del$POS[n]), ] 
				
				p_A <- freq_del$Freq[n]
							
				for (l in c(1:n_sites_sub_class_neutr))
				{
					distance <- abs(sub_class_neutr$POS[l] - freq_del$POS[n])
							
					if (distance == 0 | distance > dist_max)
					{
						next
					}
						
					NA_B <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == sub_class_neutr$POS[l]),]$ind_chr)))
					
					ind_B <- list_ind[!list_ind %in% NA_B]
					
					comm_ind <- intersect(ind_A, ind_B)
					
					if (length(comm_ind) <= min_ind)
					{
						next
					}
					
					sub_B <- genotype_neutr[which(genotype_neutr$POS == sub_class_neutr$POS[l]), ]
						
					p_B <- sub_class_neutr$Freq[l]
					
					v_D_AB_ind <- numeric()
					v_D_AB_A_ind <- numeric()
					v_D_AB_B_ind <- numeric()
					v_D_AB_A_ind <- numeric()
					v_D_AB_AB_ind <- numeric()	
							
					for (m in comm_ind)
					{
						X_A <- 0
						X_B <- 0

						if (m %in% sub_A$ind)
						{
							X_A <- 1
						} else{
							X_A <- 0
						}
						
						if (m %in% sub_B$ind)
						{
							X_B <- 1
						} else{
							X_B <- 0
						}
							
						v_D_AB_ind <- c(v_D_AB_ind, (X_A - p_A)*(X_B - p_B))
					}
							
					v_pq_AB <- p_A*(1-p_A)*p_B*(1-p_B)
					v_p_A <- p_A
					v_p_B <- p_B
					v_D_AB <- mean(v_D_AB_ind)
					v_D_AB_sqrd <- (mean(v_D_AB_ind))^2
					v_p_A_p_B_D_AB <- (1-2*p_A)*(1-2*p_B)*(mean(v_D_AB_ind))
					v_p_A_D_AB <- p_A*(mean(v_D_AB_ind))
					v_p_B_D_AB <- p_B*(mean(v_D_AB_ind))
					v_p_A_plus_p_B_D_AB <- (p_A+p_B)*(mean(v_D_AB_ind))
					v_dist <- distance
					v_dist_class <- (floor(distance/size)*size+floor(distance/size)*size+size)/2
					v_pos_A <- sub_A$POS[1]
					v_pos_B <- sub_B$POS[1]
						
					line <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd, p_A_p_B_D_AB = v_p_A_p_B_D_AB, 
					p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
					p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = length(comm_ind))
						
					write.table(line, paste(i,j,"_del_neutr_D_AB_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
				}
			}
	}
}

####

species <- "Co"
scaffold <- seq(1,8,1)
mut <- c("del","m_del")

dist_max <- 3000

window <- 60

size <- dist_max/window

for (i in species)
{
	if (as.character(i)=="Co")
	{
		num_ind <- 33
		ind_excl <- c(4,9,11,14,23,29,30,31)
	} else
	{
		num_ind <- 182
		ind_excl <- c(22,45,46,52,53,64,88,90,109,111,146)
	}
	
	list_ind <- seq(1,num_ind,1)[-ind_excl]
	
	min_ind <- length(list_ind)/2
	
	for (j in scaffold)
	{
		data_NA <- read.table(paste(i,j,"_NA.txt",sep = ""), header=TRUE, sep="\t")
		
		genotype_scaffold <- read.table(paste(i,j,"_genotype_SIFT_or_no_het.txt",sep = ""), header=TRUE, sep="\t")
		freq_scaffold <- na.omit(read.table(paste(i,j,"_freq_SIFT_or_no_het.txt",sep = ""), header=TRUE, sep="\t"))
		
		for (k in mut)
		{
			genotype <- genotype_scaffold[which(as.character(genotype_scaffold$mut_type)==k),]
			freq <- freq_scaffold[which(as.character(freq_scaffold$mut_type)==k),]

			n_sites <- nrow(freq)

			v_D_AB <- numeric()
			v_D_AB_sqrd <- numeric()
			v_p_A_p_B_D_AB <- numeric()
			v_p_A_D_AB <- numeric()
			v_p_B_D_AB <- numeric()
			v_p_A_plus_p_B_D_AB <- numeric()
			v_pq_AB <- numeric()
			v_p_A <- numeric()
			v_p_B <- numeric()
			v_dist <- numeric()
			v_dist_class <- numeric()
			v_pos_A <- numeric()
			v_pos_B <- numeric()
			v_num_ind <- numeric()
			
			headers <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd, p_A_p_B_D_AB = v_p_A_p_B_D_AB, p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
						p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
						Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind)
						
			write.table(headers, paste(i,j,"_",k,"_D_AB_SIFT_or_no_het.txt",sep = ""), quote = FALSE, sep = "\t")
				
				for (n in c(1:n_sites))
				{
					if (freq$Freq[n] == 0 | freq$Freq[n] == 1)
					{
						next
					}
					
					sub_class <- freq[which( (freq$POS>freq$POS[n]) & (freq$POS<= (freq$POS[n] + dist_max)) ), ]
						
					n_sites_sub_class <- nrow(sub_class)
						
					if(n_sites_sub_class == 0)
					{
						next
					}
						
					NA_A <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq$POS[n]),]$ind_chr)))
					
					ind_A <- list_ind[!list_ind %in% NA_A]
					
					sub_A <- genotype[which(genotype$POS == freq$POS[n]), ]
					
					p_A <- freq$Freq[n]
							
					for (l in c((n+1):(n+n_sites_sub_class)))
					{		
						if (freq$Freq[l] == 0 | freq$Freq[l] == 1)
						{
							next
						}
						
						distance <- freq$POS[l] - freq$POS[n]
							
						if (distance == 0 | distance > dist_max)
						{
							next
						}
						
						NA_B <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq$POS[l]),]$ind_chr)))
					
						ind_B <- list_ind[!list_ind %in% NA_B]
						
						comm_ind <- intersect(ind_A, ind_B)
						
						if (length(comm_ind) <= min_ind)
						{
							next
						}
						 
						sub_B <- genotype[which(genotype$POS == freq$POS[l]), ]
						
						p_B <- freq$Freq[l]
							
						v_D_AB_ind <- numeric()
							
						for (m in comm_ind)
						{
							X_A <- 0
							X_B <- 0

							if (m %in% sub_A$ind)
							{
								X_A <- 1
							} else {
								X_A <- 0
							}
						
							if (m %in% sub_B$ind)
							{
								X_B <- 1
							} else {
								X_B <- 0
							}
							
							v_D_AB_ind <- c(v_D_AB_ind, (X_A - p_A)*(X_B - p_B))
						}
							
						v_pq_AB <- p_A*(1-p_A)*p_B*(1-p_B)
						v_p_A <- p_A
						v_p_B <- p_B
						v_D_AB <- mean(v_D_AB_ind)
						v_D_AB_sqrd <- (mean(v_D_AB_ind))^2
						v_p_A_p_B_D_AB <- (1-2*p_A)*(1-2*p_B)*(mean(v_D_AB_ind))
						v_p_A_D_AB <- p_A*(mean(v_D_AB_ind))
						v_p_B_D_AB <- p_B*(mean(v_D_AB_ind))
						v_p_A_plus_p_B_D_AB <- (p_A+p_B)*(mean(v_D_AB_ind))
						v_dist <- distance
						v_dist_class <- (floor(distance/size)*size+floor(distance/size)*size+size)/2
						v_pos_A <- freq$POS[n]
						v_pos_B <- freq$POS[l]
						
						line <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd,
						p_A_p_B_D_AB = v_p_A_p_B_D_AB, p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
						p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
						Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = length(comm_ind))
						
						write.table(line, paste(i,j,"_",k,"_D_AB_SIFT_or_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
					}
				}
		}
		
			genotype_scaffold <- read.table(paste(i,j,"_genotype_SIFT_or_no_het.txt",sep = ""), header=TRUE, sep="\t")
			freq_scaffold <- na.omit(read.table(paste(i,j,"_freq_SIFT_or_no_het.txt",sep = ""), header=TRUE, sep="\t"))
			
			genotype_del <- genotype_scaffold[which(as.character(genotype_scaffold$mut_type)=="del"),]
			freq_del <- freq_scaffold[which(as.character(freq_scaffold$mut_type)=="del"),]
		
			genotype_neutr <- genotype_scaffold[which(as.character(genotype_scaffold$mut_type)=="neutr"),]
			freq_neutr <- freq_scaffold[which(as.character(freq_scaffold$mut_type)=="neutr"),]
			
			list_ind <- sort(unique(genotype_del$ind))

			n_sites_del <- nrow(freq_del)

			v_D_AB <- numeric()
			v_D_AB_sqrd <- numeric()
			v_D_AB_A <- numeric()
			v_D_AB_B <- numeric()
			v_D_AB_A <- numeric()
			v_D_AB_AB <- numeric()
			v_p_A_p_B_D_AB <- numeric()
			v_p_A_D_AB <- numeric()
			v_p_B_D_AB <- numeric()
			v_p_A_plus_p_B_D_AB <- numeric()
			v_pq_AB <- numeric()
			v_p_A <- numeric()
			v_p_B <- numeric()
			v_dist <- numeric()
			v_dist_class <- numeric()
			v_pos_A <- numeric()
			v_pos_B <- numeric()
			v_nucl_mut <- numeric()
			v_num_ind <- numeric()
			
			headers <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd,
					p_A_p_B_D_AB = v_p_A_p_B_D_AB, p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
					p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind)
						
			write.table(headers, paste(i,j,"_del_neutr_D_AB_SIFT_or_no_het.txt",sep = ""), quote = FALSE, sep = "\t")
				
			for (n in c(1:n_sites_del))
			{	
				sub_class_neutr <- freq_neutr[which( (freq_neutr$POS>=freq_del$POS[n] - dist_max) & (freq_neutr$POS<= (freq_del$POS[n] + dist_max)) ), ]
						
				n_sites_sub_class_neutr <- nrow(sub_class_neutr)
						
				if(n_sites_sub_class_neutr == 0)
				{
					next
				}
				
				NA_A <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == freq_del$POS[n]),]$ind_chr)))
					
				ind_A <- list_ind[!list_ind %in% NA_A]
				
				sub_A <- genotype_del[which(genotype_del$POS == freq_del$POS[n]), ] 
				
				p_A <- freq_del$Freq[n]
							
				for (l in c(1:n_sites_sub_class_neutr))
				{
					distance <- abs(sub_class_neutr$POS[l] - freq_del$POS[n])
							
					if (distance == 0 | distance > dist_max)
					{
						next
					}
						
					NA_B <- unique(as.numeric(gsub(".*?([0-9]+).*", "\\1",data_NA[which(data_NA$POS == sub_class_neutr$POS[l]),]$ind_chr)))
					
					ind_B <- list_ind[!list_ind %in% NA_B]
					
					comm_ind <- intersect(ind_A, ind_B)
					
					if (length(comm_ind) <= min_ind)
					{
						next
					}
					
					sub_B <- genotype_neutr[which(genotype_neutr$POS == sub_class_neutr$POS[l]), ]
						
					p_B <- sub_class_neutr$Freq[l]
					
					v_D_AB_ind <- numeric()
					v_D_AB_A_ind <- numeric()
					v_D_AB_B_ind <- numeric()
					v_D_AB_A_ind <- numeric()
					v_D_AB_AB_ind <- numeric()	
							
					for (m in comm_ind)
					{
						X_A <- 0
						X_B <- 0

						if (m %in% sub_A$ind)
						{
							X_A <- 1
						} else{
							X_A <- 0
						}
						
						if (m %in% sub_B$ind)
						{
							X_B <- 1
						} else{
							X_B <- 0
						}
							
						v_D_AB_ind <- c(v_D_AB_ind, (X_A - p_A)*(X_B - p_B))
					}
							
					v_pq_AB <- p_A*(1-p_A)*p_B*(1-p_B)
					v_p_A <- p_A
					v_p_B <- p_B
					v_D_AB <- mean(v_D_AB_ind)
					v_D_AB_sqrd <- (mean(v_D_AB_ind))^2
					v_p_A_p_B_D_AB <- (1-2*p_A)*(1-2*p_B)*(mean(v_D_AB_ind))
					v_p_A_D_AB <- p_A*(mean(v_D_AB_ind))
					v_p_B_D_AB <- p_B*(mean(v_D_AB_ind))
					v_p_A_plus_p_B_D_AB <- (p_A+p_B)*(mean(v_D_AB_ind))
					v_dist <- distance
					v_dist_class <- (floor(distance/size)*size+floor(distance/size)*size+size)/2
					v_pos_A <- sub_A$POS[1]
					v_pos_B <- sub_B$POS[1]
						
					line <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd, p_A_p_B_D_AB = v_p_A_p_B_D_AB, 
					p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
					p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = length(comm_ind))
						
					write.table(line, paste(i,j,"_del_neutr_D_AB_SIFT_or_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)	
				}
			}
	}
}

###

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr")

for (i in species)
{
		v_Delta_AB <- numeric()
		v_Delta_AB_sqrd <- numeric()
		v_D_AB_A <- numeric()
		v_D_AB_B <- numeric()
		v_Delta_AB_A <- numeric()
		v_D_AB_AB <- numeric()
		v_p_A_p_B_Delta_AB <- numeric()
		v_p_A_Delta_AB <- numeric()
		v_p_B_Delta_AB <- numeric()
		v_p_A_plus_p_B_Delta_AB <- numeric()
		v_pq_AB <- numeric()
		v_p_A <- numeric()
		v_p_B <- numeric()
		v_dist <- numeric()
		v_dist_class <- numeric()
		v_pos_A <- numeric()
		v_pos_B <- numeric()
		v_nucl_mut <- numeric()
		v_mut_type <- numeric()
		v_scaffold <- numeric()
		v_num_ind <- numeric()

	
	data_tot <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
					D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
					p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind, scaffold = v_scaffold, mut_type = v_mut_type)
			
	write.table(data_tot, paste(i,"_Delta_AB.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	write.table(data_tot, paste(i,"_del_neutr_Delta_AB.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	for (k in mut)
	{	
		for (j in scaffold)
		{
			data <- read.table(paste(i,j,"_",k,"_Delta_AB.txt",sep = ""), header=TRUE, sep="\t")
			
			data$scaffold <- rep(j,nrow(data))
			
			data$mut_type <- rep(k,nrow(data))

			write.table(data, paste(i,"_Delta_AB.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
		}
	}
	
	for (j in scaffold)
	{
		data <- read.table(paste(i,j,"_del_neutr_Delta_AB.txt",sep = ""), header=TRUE, sep="\t")
			
		data$scaffold <- rep(j,nrow(data))
			
		data$mut_type <- rep("del_neutr",nrow(data))

		write.table(data, paste(i,"_del_neutr_Delta_AB.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
	}
}

###

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del")

for (i in species)
{
		v_Delta_AB <- numeric()
		v_Delta_AB_sqrd <- numeric()
		v_D_AB_A <- numeric()
		v_D_AB_B <- numeric()
		v_Delta_AB_A <- numeric()
		v_D_AB_AB <- numeric()
		v_p_A_p_B_Delta_AB <- numeric()
		v_p_A_Delta_AB <- numeric()
		v_p_B_Delta_AB <- numeric()
		v_p_A_plus_p_B_Delta_AB <- numeric()
		v_pq_AB <- numeric()
		v_p_A <- numeric()
		v_p_B <- numeric()
		v_dist <- numeric()
		v_dist_class <- numeric()
		v_pos_A <- numeric()
		v_pos_B <- numeric()
		v_nucl_mut <- numeric()
		v_mut_type <- numeric()
		v_scaffold <- numeric()
		v_num_ind <- numeric()

	
	data_tot <- data.frame(Delta_AB = v_Delta_AB, Delta_AB_sqrd = v_Delta_AB_sqrd, D_AB_A = v_D_AB_A, D_AB_B = v_D_AB_B, Delta_AB_A = v_Delta_AB_A,
					D_AB_AB = v_D_AB_AB, p_A_p_B_Delta_AB = v_p_A_p_B_Delta_AB, p_A_Delta_AB = v_p_A_Delta_AB, p_B_Delta_AB = v_p_B_Delta_AB,
					p_A_plus_p_B_Delta_AB = v_p_A_plus_p_B_Delta_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind, scaffold = v_scaffold, mut_type = v_mut_type)
			
	write.table(data_tot, paste(i,"_Delta_AB_SIFT_or.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	write.table(data_tot, paste(i,"_del_neutr_Delta_AB_SIFT_or.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	for (k in mut)
	{	
		for (j in scaffold)
		{
			data <- read.table(paste(i,j,"_",k,"_Delta_AB_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
			
			data$scaffold <- rep(j,nrow(data))
			
			data$mut_type <- rep(k,nrow(data))

			write.table(data, paste(i,"_Delta_AB_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
		}
	}
	
	for (j in scaffold)
	{
		data <- read.table(paste(i,j,"_del_neutr_Delta_AB_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
			
		data$scaffold <- rep(j,nrow(data))
			
		data$mut_type <- rep("del_neutr",nrow(data))

		write.table(data, paste(i,"_del_neutr_Delta_AB_SIFT_or.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
	}
}

####

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr")

for (i in species)
{
		v_D_AB <- numeric()
		v_D_AB_sqrd <- numeric()
		v_p_A_p_B_D_AB <- numeric()
		v_p_A_D_AB <- numeric()
		v_p_B_D_AB <- numeric()
		v_p_A_plus_p_B_D_AB <- numeric()
		v_pq_AB <- numeric()
		v_p_A <- numeric()
		v_p_B <- numeric()
		v_dist <- numeric()
		v_dist_class <- numeric()
		v_pos_A <- numeric()
		v_pos_B <- numeric()
		v_nucl_mut <- numeric()
		v_mut_type <- numeric()
		v_scaffold <- numeric()
		v_num_ind <- numeric()

	
	data_tot <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd, p_A_p_B_D_AB = v_p_A_p_B_D_AB, p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
					p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind, scaffold = v_scaffold, mut_type = v_mut_type)
			
	write.table(data_tot, paste(i,"_D_AB_no_het.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	write.table(data_tot, paste(i,"_del_neutr_D_AB_no_het.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	for (k in mut)
	{	
		for (j in scaffold)
		{
			data <- read.table(paste(i,j,"_",k,"_D_AB_no_het.txt",sep = ""), header=TRUE, sep="\t")
			
			data$scaffold <- rep(j,nrow(data))
			
			data$mut_type <- rep(k,nrow(data))

			write.table(data, paste(i,"_D_AB_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
		}
	}
	
	for (j in scaffold)
	{
		data <- read.table(paste(i,j,"_del_neutr_D_AB_no_het.txt",sep = ""), header=TRUE, sep="\t")
			
		data$scaffold <- rep(j,nrow(data))
			
		data$mut_type <- rep("del_neutr",nrow(data))

		write.table(data, paste(i,"_del_neutr_D_AB_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
	}
}

###

species <- c("Cg","Co")
scaffold <- seq(1,8,1)
mut <- c("del","m_del")

for (i in species)
{
		v_D_AB <- numeric()
		v_D_AB_sqrd <- numeric()
		v_p_A_p_B_D_AB <- numeric()
		v_p_A_D_AB <- numeric()
		v_p_B_D_AB <- numeric()
		v_p_A_plus_p_B_D_AB <- numeric()
		v_pq_AB <- numeric()
		v_p_A <- numeric()
		v_p_B <- numeric()
		v_dist <- numeric()
		v_dist_class <- numeric()
		v_pos_A <- numeric()
		v_pos_B <- numeric()
		v_nucl_mut <- numeric()
		v_mut_type <- numeric()
		v_scaffold <- numeric()
		v_num_ind <- numeric()

	
	data_tot <- data.frame(D_AB = v_D_AB, D_AB_sqrd = v_D_AB_sqrd, p_A_p_B_D_AB = v_p_A_p_B_D_AB, p_A_D_AB = v_p_A_D_AB, p_B_D_AB = v_p_B_D_AB,
					p_A_plus_p_B_D_AB = v_p_A_plus_p_B_D_AB, pq_AB = v_pq_AB, p_A = v_p_A, p_B = v_p_B,					
					Distance = v_dist, Class_dist = v_dist_class, pos_A = v_pos_A, pos_B = v_pos_B, num_ind = v_num_ind, scaffold = v_scaffold, mut_type = v_mut_type)
			
	write.table(data_tot, paste(i,"_D_AB_SIFT_or_no_het.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	write.table(data_tot, paste(i,"_del_neutr_D_AB_SIFT_or_no_het.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	for (k in mut)
	{	
		for (j in scaffold)
		{
			data <- read.table(paste(i,j,"_",k,"_D_AB_SIFT_or_no_het.txt",sep = ""), header=TRUE, sep="\t")
			
			data$scaffold <- rep(j,nrow(data))
			
			data$mut_type <- rep(k,nrow(data))

			write.table(data, paste(i,"_D_AB_SIFT_or_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
		}
	}
	
	for (j in scaffold)
	{
		data <- read.table(paste(i,j,"_del_neutr_D_AB_SIFT_or_no_het.txt",sep = ""), header=TRUE, sep="\t")
			
		data$scaffold <- rep(j,nrow(data))
			
		data$mut_type <- rep("del_neutr",nrow(data))

		write.table(data, paste(i,"_del_neutr_D_AB_SIFT_or_no_het.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
	}
}

####

mut <- c("del","m_del","neutr","del_neutr")

species <- c("Cg","Co")

for (i in species)
{
	data_tot <- data.frame(Delta_AB_scaled_mean = numeric(), Delta_AB_scaled_mean_inf_CI = numeric(), Delta_AB_scaled_mean_sup_CI = numeric(),
					Delta_AB_mean_scaled = numeric(), Delta_AB_mean_scaled_inf_CI = numeric(), Delta_AB_mean_scaled_sup_CI = numeric(),
					Delta_AB_sqrd_scaled_mean = numeric(), Delta_AB_sqrd_scaled_mean_inf_CI = numeric(), Delta_AB_sqrd_scaled_mean_sup_CI = numeric(),
					Delta_AB_sqrd_mean_scaled = numeric(), Delta_AB_sqrd_mean_scaled_inf_CI = numeric(), Delta_AB_sqrd_mean_scaled_sup_CI = numeric(),
					D_AB_A_scaled_mean = numeric(), D_AB_A_scaled_mean_inf_CI = numeric(), D_AB_A_scaled_mean_sup_CI = numeric(),
					D_AB_A_mean_scaled = numeric(), D_AB_A_mean_scaled_inf_CI = numeric(), D_AB_A_mean_scaled_sup_CI = numeric(),
					D_AB_B_scaled_mean = numeric(), D_AB_B_scaled_mean_inf_CI = numeric(), D_AB_B_scaled_mean_sup_CI = numeric(),
					D_AB_B_mean_scaled = numeric(), D_AB_B_mean_scaled_inf_CI = numeric(), D_AB_B_mean_scaled_sup_CI = numeric(),
					Delta_AB_A_scaled_mean = numeric(), Delta_AB_A_scaled_mean_inf_CI = numeric(), Delta_AB_A_scaled_mean_sup_CI = numeric(),
					Delta_AB_A_mean_scaled = numeric(), Delta_AB_A_mean_scaled_inf_CI = numeric(), Delta_AB_A_mean_scaled_sup_CI = numeric(),
					D_AB_AB_scaled_mean = numeric(), D_AB_AB_scaled_mean_inf_CI = numeric(), D_AB_AB_scaled_mean_sup_CI = numeric(),
					D_AB_AB_mean_scaled = numeric(), D_AB_AB_mean_scaled_inf_CI = numeric(), D_AB_AB_mean_scaled_sup_CI = numeric(),
					p_A_p_B_Delta_AB_scaled_mean = numeric(), p_A_p_B_Delta_AB_scaled_mean_inf_CI = numeric(), p_A_p_B_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_A_p_B_Delta_AB_mean_scaled = numeric(), p_A_p_B_Delta_AB_mean_scaled_inf_CI = numeric(), p_A_p_B_Delta_AB_mean_scaled_sup_CI = numeric(),
					p_A_Delta_AB_scaled_mean = numeric(), p_A_Delta_AB_scaled_mean_inf_CI = numeric(), p_A_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_A_Delta_AB_mean_scaled = numeric(), p_A_Delta_AB_mean_scaled_inf_CI = numeric(), p_A_Delta_AB_mean_scaled_sup_CI = numeric(),
					p_B_Delta_AB_scaled_mean = numeric(), p_B_Delta_AB_scaled_mean_inf_CI = numeric(), p_B_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_B_Delta_AB_mean_scaled = numeric(), p_B_Delta_AB_mean_scaled_inf_CI = numeric(), p_B_Delta_AB_mean_scaled_sup_CI = numeric(),
					p_A_plus_p_B_Delta_AB_scaled_mean = numeric(), p_A_plus_p_B_Delta_AB_scaled_mean_inf_CI = numeric(), p_A_plus_p_B_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_A_plus_p_B_Delta_AB_mean_scaled = numeric(), p_A_plus_p_B_Delta_AB_mean_scaled_inf_CI = numeric(), p_A_plus_p_B_Delta_AB_mean_scaled_sup_CI = numeric(),
					Class_dist = numeric(), mut_type = numeric(), n_sample = numeric())
			
	write.table(data_tot, paste(i,"_Delta_AB_mean.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	data <- read.table(paste(i,"_Delta_AB.txt",sep = ""), header=TRUE, sep="\t")
	
	library("dplyr")
	
	library(data.table)
	
	dist_class <- data %>% distinct(sort(data$Class_dist)) 
	
	dist_class <- head(dist_class[,1],-1)
	
	for (k in mut)
	{	
		if (k == "del_neutr")
		{
			data <- read.table(paste(i,"_del_neutr_Delta_AB.txt",sep = ""), header=TRUE, sep="\t")
		
			dist_class <- data %>% distinct(sort(data$Class_dist)) 
		
			dist_class <- head(dist_class[,1],-1)
		}
				for (l in dist_class)
				{
					data_sub <- data[which(data$mut == k & data$Class_dist == l),]
					
					n <- nrow(data_sub)
					
					v_Delta_AB_boot <- numeric()
					v_Delta_AB_2_boot <- numeric()
					v_Delta_AB_sqrd_boot <- numeric()
					v_Delta_AB_sqrd_2_boot <- numeric()
					v_D_AB_A_boot <- numeric()
					v_D_AB_A_2_boot <- numeric()
					v_D_AB_B_boot <- numeric()
					v_D_AB_B_2_boot <- numeric()
					v_Delta_AB_A_boot <- numeric()
					v_Delta_AB_A_2_boot <- numeric()
					v_D_AB_AB_boot <- numeric()
					v_D_AB_AB_2_boot <- numeric()
					v_p_A_p_B_Delta_AB_boot <- numeric()
					v_p_A_p_B_Delta_AB_2_boot <- numeric()
					v_p_A_Delta_AB_boot <- numeric()
					v_p_A_Delta_AB_2_boot <- numeric()
					v_p_B_Delta_AB_boot <- numeric()
					v_p_B_Delta_AB_2_boot <- numeric()
					v_p_A_plus_p_B_Delta_AB_boot <- numeric()
					v_p_A_plus_p_B_Delta_AB_2_boot <- numeric()
					
					for (m in c(1:1000))
					{
						data_sub_samp <- numeric()
						
						data_sub_samp <- data_sub[sample(seq_len(n),n, TRUE),]
						
						
						v_Delta_AB_boot <- c(v_Delta_AB_boot, mean(data_sub_samp$Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_Delta_AB_2_boot <- c(v_Delta_AB_2_boot, mean(data_sub_samp$Delta_AB/data_sub_samp$pq_AB))
						
						
						v_Delta_AB_sqrd_boot <- c(v_Delta_AB_sqrd_boot, mean(data_sub_samp$Delta_AB_sqrd)/mean(data_sub_samp$pq_AB))
						
						v_Delta_AB_sqrd_2_boot <- c(v_Delta_AB_sqrd_2_boot, mean(data_sub_samp$Delta_AB_sqrd/data_sub_samp$pq_AB))
						
						
						v_D_AB_A_boot <- c(v_D_AB_A_boot, mean(data_sub_samp$D_AB_A)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_A_2_boot <- c(v_D_AB_A_2_boot, mean(data_sub_samp$D_AB_A/data_sub_samp$pq_AB))
						
						
						v_D_AB_B_boot <- c(v_D_AB_B_boot, mean(data_sub_samp$D_AB_B)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_B_2_boot <- c(v_D_AB_B_2_boot, mean(data_sub_samp$D_AB_B/data_sub_samp$pq_AB))
						
						
						v_Delta_AB_A_boot <- c(v_Delta_AB_A_boot, mean(data_sub_samp$Delta_AB_A)/mean(data_sub_samp$pq_AB))
						
						v_Delta_AB_A_2_boot <- c(v_Delta_AB_A_2_boot, mean(data_sub_samp$Delta_AB_A/data_sub_samp$pq_AB))
						
						
						v_D_AB_AB_boot <- c(v_D_AB_AB_boot, mean(data_sub_samp$D_AB_AB)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_AB_2_boot <- c(v_D_AB_AB_2_boot, mean(data_sub_samp$D_AB_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_p_B_Delta_AB_boot <- c(v_p_A_p_B_Delta_AB_boot, mean(data_sub_samp$p_A_p_B_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_p_B_Delta_AB_2_boot <- c(v_p_A_p_B_Delta_AB_2_boot, mean(data_sub_samp$p_A_p_B_Delta_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_Delta_AB_boot <- c(v_p_A_Delta_AB_boot, mean(data_sub_samp$p_A_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_Delta_AB_2_boot <- c(v_p_A_Delta_AB_2_boot, mean(data_sub_samp$p_A_Delta_AB/data_sub_samp$pq_AB))
						
						
						v_p_B_Delta_AB_boot <- c(v_p_B_Delta_AB_boot, mean(data_sub_samp$p_B_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_B_Delta_AB_2_boot <- c(v_p_B_Delta_AB_2_boot, mean(data_sub_samp$p_B_Delta_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_plus_p_B_Delta_AB_boot <- c(v_p_A_plus_p_B_Delta_AB_boot, mean(data_sub_samp$p_A_plus_p_B_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_plus_p_B_Delta_AB_2_boot <- c(v_p_A_plus_p_B_Delta_AB_2_boot, mean(data_sub_samp$p_A_plus_p_B_Delta_AB/data_sub_samp$pq_AB))
					
					}
					
					v_Delta_AB <- mean(v_Delta_AB_boot)
					
					v_Delta_AB_inf_CI <- quantile(v_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_sup_CI <- quantile(v_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_2 <- mean(v_Delta_AB_2_boot)
					
					v_Delta_AB_2_inf_CI <- quantile(v_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_2_sup_CI <- quantile(v_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_sqrd <- mean(v_Delta_AB_sqrd_boot)
					
					v_Delta_AB_sqrd_inf_CI <- quantile(v_Delta_AB_sqrd_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_sqrd_sup_CI <- quantile(v_Delta_AB_sqrd_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_sqrd_2 <- mean(v_Delta_AB_sqrd_2_boot)
					
					v_Delta_AB_sqrd_2_inf_CI <- quantile(v_Delta_AB_sqrd_2_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_sqrd_2_sup_CI <- quantile(v_Delta_AB_sqrd_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_A <- mean(v_D_AB_A_boot)
					
					v_D_AB_A_inf_CI <- quantile(v_D_AB_A_boot, probs = 0.025)[[1]]
					
					v_D_AB_A_sup_CI <- quantile(v_D_AB_A_boot, probs = 0.975)[[1]]
					
					v_D_AB_A_2 <- mean(v_D_AB_A_2_boot)
					
					v_D_AB_A_2_inf_CI <- quantile(v_D_AB_A_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_A_2_sup_CI <- quantile(v_D_AB_A_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_B <- mean(v_D_AB_B_boot)
					
					v_D_AB_B_inf_CI <- quantile(v_D_AB_B_boot, probs = 0.025)[[1]]
					
					v_D_AB_B_sup_CI <- quantile(v_D_AB_B_boot, probs = 0.975)[[1]]
					
					v_D_AB_B_2 <- mean(v_D_AB_B_2_boot)
					
					v_D_AB_B_2_inf_CI <- quantile(v_D_AB_B_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_B_2_sup_CI <- quantile(v_D_AB_B_2_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_A <- mean(v_Delta_AB_A_boot)
					
					v_Delta_AB_A_inf_CI <- quantile(v_Delta_AB_A_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_A_sup_CI <- quantile(v_Delta_AB_A_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_A_2 <- mean(v_Delta_AB_A_2_boot)
					
					v_Delta_AB_A_2_inf_CI <- quantile(v_Delta_AB_A_2_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_A_2_sup_CI <- quantile(v_Delta_AB_A_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_AB <- mean(v_D_AB_AB_boot)
					
					v_D_AB_AB_inf_CI <- quantile(v_D_AB_AB_boot, probs = 0.025)[[1]]
					
					v_D_AB_AB_sup_CI <- quantile(v_D_AB_AB_boot, probs = 0.975)[[1]]
					
					v_D_AB_AB_2 <- mean(v_D_AB_AB_2_boot)
					
					v_D_AB_AB_2_inf_CI <- quantile(v_D_AB_AB_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_AB_2_sup_CI <- quantile(v_D_AB_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_p_B_Delta_AB <- mean(v_p_A_p_B_Delta_AB_boot)
					
					v_p_A_p_B_Delta_AB_inf_CI <- quantile(v_p_A_p_B_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_Delta_AB_sup_CI <- quantile(v_p_A_p_B_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_p_B_Delta_AB_2 <- mean(v_p_A_p_B_Delta_AB_2_boot)
					
					v_p_A_p_B_Delta_AB_2_inf_CI <- quantile(v_p_A_p_B_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_Delta_AB_2_sup_CI <- quantile(v_p_A_p_B_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_Delta_AB <- mean(v_p_A_Delta_AB_boot)
					
					v_p_A_Delta_AB_inf_CI <- quantile(v_p_A_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_Delta_AB_sup_CI <- quantile(v_p_A_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_Delta_AB_2 <- mean(v_p_A_Delta_AB_2_boot)
					
					v_p_A_Delta_AB_2_inf_CI <- quantile(v_p_A_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_Delta_AB_2_sup_CI <- quantile(v_p_A_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_B_Delta_AB <- mean(v_p_B_Delta_AB_boot)
					
					v_p_B_Delta_AB_inf_CI <- quantile(v_p_B_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_B_Delta_AB_sup_CI <- quantile(v_p_B_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_B_Delta_AB_2 <- mean(v_p_B_Delta_AB_2_boot)
					
					v_p_B_Delta_AB_2_inf_CI <- quantile(v_p_B_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_B_Delta_AB_2_sup_CI <- quantile(v_p_B_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_Delta_AB <- mean(v_p_A_plus_p_B_Delta_AB_boot)
					
					v_p_A_plus_p_B_Delta_AB_inf_CI <- quantile(v_p_A_plus_p_B_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_Delta_AB_sup_CI <- quantile(v_p_A_plus_p_B_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_Delta_AB_2 <- mean(v_p_A_plus_p_B_Delta_AB_2_boot)
					
					v_p_A_plus_p_B_Delta_AB_2_inf_CI <- quantile(v_p_A_plus_p_B_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_Delta_AB_2_sup_CI <- quantile(v_p_A_plus_p_B_Delta_AB_2_boot, probs = 0.975)[[1]]
					

					write.table(data.frame(Delta_AB_scaled_mean = v_Delta_AB, Delta_AB_scaled_mean_inf_CI = v_Delta_AB_inf_CI, Delta_AB_scaled_mean_sup_CI = v_Delta_AB_sup_CI,
					Delta_AB_mean_scaled = v_Delta_AB_2, Delta_AB_mean_scaled_inf_CI = v_Delta_AB_2_inf_CI, Delta_AB_mean_scaled_sup_CI = v_Delta_AB_2_sup_CI,
					Delta_AB_sqrd_scaled_mean = v_Delta_AB_sqrd, Delta_AB_sqrd_scaled_mean_inf_CI = v_Delta_AB_sqrd_inf_CI, Delta_AB_sqrd_scaled_mean_sup_CI = v_Delta_AB_sqrd_sup_CI,
					Delta_AB_sqrd_mean_scaled = v_Delta_AB_sqrd_2, Delta_AB_sqrd_mean_scaled_inf_CI = v_Delta_AB_sqrd_2_inf_CI, Delta_AB_sqrd_mean_scaled_sup_CI = v_Delta_AB_sqrd_2_sup_CI,
					D_AB_A_scaled_mean = v_D_AB_A, D_AB_A_scaled_mean_inf_CI = v_D_AB_A_inf_CI, D_AB_A_scaled_mean_sup_CI = v_D_AB_A_sup_CI,
					D_AB_A_mean_scaled = v_D_AB_A_2, D_AB_A_mean_scaled_inf_CI = v_D_AB_A_2_inf_CI, D_AB_A_mean_scaled_sup_CI = v_D_AB_A_2_sup_CI,
					D_AB_B_scaled_mean = v_D_AB_B, D_AB_B_scaled_mean_inf_CI = v_D_AB_B_inf_CI, D_AB_B_scaled_mean_sup_CI = v_D_AB_B_sup_CI,
					D_AB_B_mean_scaled = v_D_AB_B_2, D_AB_B_mean_scaled_inf_CI = v_D_AB_B_2_inf_CI, D_AB_B_mean_scaled_sup_CI = v_D_AB_B_2_sup_CI,
					Delta_AB_A_scaled_mean = v_Delta_AB_A, Delta_AB_A_scaled_mean_inf_CI = v_Delta_AB_A_inf_CI, Delta_AB_A_scaled_mean_sup_CI = v_Delta_AB_A_sup_CI,
					Delta_AB_A_mean_scaled = v_Delta_AB_A_2, Delta_AB_A_mean_scaled_inf_CI = v_Delta_AB_A_2_inf_CI, Delta_AB_A_mean_scaled_sup_CI = v_Delta_AB_A_2_sup_CI,
					D_AB_AB_scaled_mean = v_D_AB_AB, D_AB_AB_scaled_mean_inf_CI = v_D_AB_AB_inf_CI, D_AB_AB_scaled_mean_sup_CI = v_D_AB_AB_sup_CI,
					D_AB_AB_mean_scaled = v_D_AB_AB_2, D_AB_AB_mean_scaled_inf_CI = v_D_AB_AB_2_inf_CI, D_AB_AB_mean_scaled_sup_CI = v_D_AB_AB_2_sup_CI,
					p_A_p_B_Delta_AB_scaled_mean = v_p_A_p_B_Delta_AB, p_A_p_B_Delta_AB_scaled_mean_inf_CI = v_p_A_p_B_Delta_AB_inf_CI, p_A_p_B_Delta_AB_scaled_mean_sup_CI = v_p_A_p_B_Delta_AB_sup_CI,
					p_A_p_B_Delta_AB_mean_scaled = v_p_A_p_B_Delta_AB_2, p_A_p_B_Delta_AB_mean_scaled_inf_CI = v_p_A_p_B_Delta_AB_2_inf_CI, p_A_p_B_Delta_AB_mean_scaled_sup_CI = v_p_A_p_B_Delta_AB_2_sup_CI,
					p_A_Delta_AB_scaled_mean = v_p_A_Delta_AB, p_A_Delta_AB_scaled_mean_inf_CI = v_p_A_Delta_AB_inf_CI, p_A_Delta_AB_scaled_mean_sup_CI = v_p_A_Delta_AB_sup_CI,
					p_A_Delta_AB_mean_scaled = v_p_A_Delta_AB_2, p_A_Delta_AB_mean_scaled_inf_CI = v_p_A_Delta_AB_2_inf_CI, p_A_Delta_AB_mean_scaled_sup_CI = v_p_A_Delta_AB_2_sup_CI,
					p_B_Delta_AB_scaled_mean = v_p_B_Delta_AB, p_B_Delta_AB_scaled_mean_inf_CI = v_p_B_Delta_AB_inf_CI, p_B_Delta_AB_scaled_mean_sup_CI = v_p_B_Delta_AB_sup_CI,
					p_B_Delta_AB_mean_scaled = v_p_B_Delta_AB_2, p_B_Delta_AB_mean_scaled_inf_CI = v_p_B_Delta_AB_2_inf_CI, p_B_Delta_AB_mean_scaled_sup_CI = v_p_B_Delta_AB_2_sup_CI,
					p_A_plus_p_B_Delta_AB_scaled_mean = v_p_A_plus_p_B_Delta_AB, p_A_plus_p_B_Delta_AB_scaled_mean_inf_CI = v_p_A_plus_p_B_Delta_AB_inf_CI, p_A_plus_p_B_Delta_AB_scaled_mean_sup_CI = v_p_A_plus_p_B_Delta_AB_sup_CI,
					p_A_plus_p_B_Delta_AB_mean_scaled = v_p_A_plus_p_B_Delta_AB_2, p_A_plus_p_B_Delta_AB_mean_scaled_inf_CI = v_p_A_plus_p_B_Delta_AB_2_inf_CI, p_A_plus_p_B_Delta_AB_mean_scaled_sup_CI = v_p_A_plus_p_B_Delta_AB_2_sup_CI,
					Class_dist = l, mut_type = k, n_sample = n), paste(i,"_Delta_AB_mean.txt",sep = ""), 
					quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
				}
	}
}

####

mut <- c("del","m_del","del_neutr")

species <- c("Cg","Co")

for (i in species)
{
	data_tot <- data.frame(Delta_AB_scaled_mean = numeric(), Delta_AB_scaled_mean_inf_CI = numeric(), Delta_AB_scaled_mean_sup_CI = numeric(),
					Delta_AB_mean_scaled = numeric(), Delta_AB_mean_scaled_inf_CI = numeric(), Delta_AB_mean_scaled_sup_CI = numeric(),
					Delta_AB_sqrd_scaled_mean = numeric(), Delta_AB_sqrd_scaled_mean_inf_CI = numeric(), Delta_AB_sqrd_scaled_mean_sup_CI = numeric(),
					Delta_AB_sqrd_mean_scaled = numeric(), Delta_AB_sqrd_mean_scaled_inf_CI = numeric(), Delta_AB_sqrd_mean_scaled_sup_CI = numeric(),
					D_AB_A_scaled_mean = numeric(), D_AB_A_scaled_mean_inf_CI = numeric(), D_AB_A_scaled_mean_sup_CI = numeric(),
					D_AB_A_mean_scaled = numeric(), D_AB_A_mean_scaled_inf_CI = numeric(), D_AB_A_mean_scaled_sup_CI = numeric(),
					D_AB_B_scaled_mean = numeric(), D_AB_B_scaled_mean_inf_CI = numeric(), D_AB_B_scaled_mean_sup_CI = numeric(),
					D_AB_B_mean_scaled = numeric(), D_AB_B_mean_scaled_inf_CI = numeric(), D_AB_B_mean_scaled_sup_CI = numeric(),
					Delta_AB_A_scaled_mean = numeric(), Delta_AB_A_scaled_mean_inf_CI = numeric(), Delta_AB_A_scaled_mean_sup_CI = numeric(),
					Delta_AB_A_mean_scaled = numeric(), Delta_AB_A_mean_scaled_inf_CI = numeric(), Delta_AB_A_mean_scaled_sup_CI = numeric(),
					D_AB_AB_scaled_mean = numeric(), D_AB_AB_scaled_mean_inf_CI = numeric(), D_AB_AB_scaled_mean_sup_CI = numeric(),
					D_AB_AB_mean_scaled = numeric(), D_AB_AB_mean_scaled_inf_CI = numeric(), D_AB_AB_mean_scaled_sup_CI = numeric(),
					p_A_p_B_Delta_AB_scaled_mean = numeric(), p_A_p_B_Delta_AB_scaled_mean_inf_CI = numeric(), p_A_p_B_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_A_p_B_Delta_AB_mean_scaled = numeric(), p_A_p_B_Delta_AB_mean_scaled_inf_CI = numeric(), p_A_p_B_Delta_AB_mean_scaled_sup_CI = numeric(),
					p_A_Delta_AB_scaled_mean = numeric(), p_A_Delta_AB_scaled_mean_inf_CI = numeric(), p_A_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_A_Delta_AB_mean_scaled = numeric(), p_A_Delta_AB_mean_scaled_inf_CI = numeric(), p_A_Delta_AB_mean_scaled_sup_CI = numeric(),
					p_B_Delta_AB_scaled_mean = numeric(), p_B_Delta_AB_scaled_mean_inf_CI = numeric(), p_B_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_B_Delta_AB_mean_scaled = numeric(), p_B_Delta_AB_mean_scaled_inf_CI = numeric(), p_B_Delta_AB_mean_scaled_sup_CI = numeric(),
					p_A_plus_p_B_Delta_AB_scaled_mean = numeric(), p_A_plus_p_B_Delta_AB_scaled_mean_inf_CI = numeric(), p_A_plus_p_B_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_A_plus_p_B_Delta_AB_mean_scaled = numeric(), p_A_plus_p_B_Delta_AB_mean_scaled_inf_CI = numeric(), p_A_plus_p_B_Delta_AB_mean_scaled_sup_CI = numeric(),
					Class_dist = numeric(), mut_type = numeric(), n_sample = numeric())
			
	write.table(data_tot, paste(i,"_Delta_AB_mean_SIFT_or.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	data <- read.table(paste(i,"_Delta_AB_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
	
	library("dplyr")
	
	library(data.table)
	
	dist_class <- data %>% distinct(sort(data$Class_dist)) 
	
	dist_class <- head(dist_class[,1],-1)
	
	for (k in mut)
	{	
		if (k == "del_neutr")
		{
			data <- read.table(paste(i,"_del_neutr_Delta_AB_SIFT_or.txt",sep = ""), header=TRUE, sep="\t")
		
			dist_class <- data %>% distinct(sort(data$Class_dist)) 
		
			dist_class <- head(dist_class[,1],-1)
		}
				for (l in dist_class)
				{
					data_sub <- data[which(data$mut == k & data$Class_dist == l),]
					
					n <- nrow(data_sub)
					
					v_Delta_AB_boot <- numeric()
					v_Delta_AB_2_boot <- numeric()
					v_Delta_AB_sqrd_boot <- numeric()
					v_Delta_AB_sqrd_2_boot <- numeric()
					v_D_AB_A_boot <- numeric()
					v_D_AB_A_2_boot <- numeric()
					v_D_AB_B_boot <- numeric()
					v_D_AB_B_2_boot <- numeric()
					v_Delta_AB_A_boot <- numeric()
					v_Delta_AB_A_2_boot <- numeric()
					v_D_AB_AB_boot <- numeric()
					v_D_AB_AB_2_boot <- numeric()
					v_p_A_p_B_Delta_AB_boot <- numeric()
					v_p_A_p_B_Delta_AB_2_boot <- numeric()
					v_p_A_Delta_AB_boot <- numeric()
					v_p_A_Delta_AB_2_boot <- numeric()
					v_p_B_Delta_AB_boot <- numeric()
					v_p_B_Delta_AB_2_boot <- numeric()
					v_p_A_plus_p_B_Delta_AB_boot <- numeric()
					v_p_A_plus_p_B_Delta_AB_2_boot <- numeric()
					
					for (m in c(1:1000))
					{
						data_sub_samp <- numeric()
						
						data_sub_samp <- data_sub[sample(seq_len(n),n, TRUE),]
						
						
						v_Delta_AB_boot <- c(v_Delta_AB_boot, mean(data_sub_samp$Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_Delta_AB_2_boot <- c(v_Delta_AB_2_boot, mean(data_sub_samp$Delta_AB/data_sub_samp$pq_AB))
						
						
						v_Delta_AB_sqrd_boot <- c(v_Delta_AB_sqrd_boot, mean(data_sub_samp$Delta_AB_sqrd)/mean(data_sub_samp$pq_AB))
						
						v_Delta_AB_sqrd_2_boot <- c(v_Delta_AB_sqrd_2_boot, mean(data_sub_samp$Delta_AB_sqrd/data_sub_samp$pq_AB))
						
						
						v_D_AB_A_boot <- c(v_D_AB_A_boot, mean(data_sub_samp$D_AB_A)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_A_2_boot <- c(v_D_AB_A_2_boot, mean(data_sub_samp$D_AB_A/data_sub_samp$pq_AB))
						
						
						v_D_AB_B_boot <- c(v_D_AB_B_boot, mean(data_sub_samp$D_AB_B)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_B_2_boot <- c(v_D_AB_B_2_boot, mean(data_sub_samp$D_AB_B/data_sub_samp$pq_AB))
						
						
						v_Delta_AB_A_boot <- c(v_Delta_AB_A_boot, mean(data_sub_samp$Delta_AB_A)/mean(data_sub_samp$pq_AB))
						
						v_Delta_AB_A_2_boot <- c(v_Delta_AB_A_2_boot, mean(data_sub_samp$Delta_AB_A/data_sub_samp$pq_AB))
						
						
						v_D_AB_AB_boot <- c(v_D_AB_AB_boot, mean(data_sub_samp$D_AB_AB)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_AB_2_boot <- c(v_D_AB_AB_2_boot, mean(data_sub_samp$D_AB_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_p_B_Delta_AB_boot <- c(v_p_A_p_B_Delta_AB_boot, mean(data_sub_samp$p_A_p_B_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_p_B_Delta_AB_2_boot <- c(v_p_A_p_B_Delta_AB_2_boot, mean(data_sub_samp$p_A_p_B_Delta_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_Delta_AB_boot <- c(v_p_A_Delta_AB_boot, mean(data_sub_samp$p_A_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_Delta_AB_2_boot <- c(v_p_A_Delta_AB_2_boot, mean(data_sub_samp$p_A_Delta_AB/data_sub_samp$pq_AB))
						
						
						v_p_B_Delta_AB_boot <- c(v_p_B_Delta_AB_boot, mean(data_sub_samp$p_B_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_B_Delta_AB_2_boot <- c(v_p_B_Delta_AB_2_boot, mean(data_sub_samp$p_B_Delta_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_plus_p_B_Delta_AB_boot <- c(v_p_A_plus_p_B_Delta_AB_boot, mean(data_sub_samp$p_A_plus_p_B_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_plus_p_B_Delta_AB_2_boot <- c(v_p_A_plus_p_B_Delta_AB_2_boot, mean(data_sub_samp$p_A_plus_p_B_Delta_AB/data_sub_samp$pq_AB))
					
					}
					
					v_Delta_AB <- mean(v_Delta_AB_boot)
					
					v_Delta_AB_inf_CI <- quantile(v_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_sup_CI <- quantile(v_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_2 <- mean(v_Delta_AB_2_boot)
					
					v_Delta_AB_2_inf_CI <- quantile(v_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_2_sup_CI <- quantile(v_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_sqrd <- mean(v_Delta_AB_sqrd_boot)
					
					v_Delta_AB_sqrd_inf_CI <- quantile(v_Delta_AB_sqrd_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_sqrd_sup_CI <- quantile(v_Delta_AB_sqrd_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_sqrd_2 <- mean(v_Delta_AB_sqrd_2_boot)
					
					v_Delta_AB_sqrd_2_inf_CI <- quantile(v_Delta_AB_sqrd_2_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_sqrd_2_sup_CI <- quantile(v_Delta_AB_sqrd_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_A <- mean(v_D_AB_A_boot)
					
					v_D_AB_A_inf_CI <- quantile(v_D_AB_A_boot, probs = 0.025)[[1]]
					
					v_D_AB_A_sup_CI <- quantile(v_D_AB_A_boot, probs = 0.975)[[1]]
					
					v_D_AB_A_2 <- mean(v_D_AB_A_2_boot)
					
					v_D_AB_A_2_inf_CI <- quantile(v_D_AB_A_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_A_2_sup_CI <- quantile(v_D_AB_A_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_B <- mean(v_D_AB_B_boot)
					
					v_D_AB_B_inf_CI <- quantile(v_D_AB_B_boot, probs = 0.025)[[1]]
					
					v_D_AB_B_sup_CI <- quantile(v_D_AB_B_boot, probs = 0.975)[[1]]
					
					v_D_AB_B_2 <- mean(v_D_AB_B_2_boot)
					
					v_D_AB_B_2_inf_CI <- quantile(v_D_AB_B_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_B_2_sup_CI <- quantile(v_D_AB_B_2_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_A <- mean(v_Delta_AB_A_boot)
					
					v_Delta_AB_A_inf_CI <- quantile(v_Delta_AB_A_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_A_sup_CI <- quantile(v_Delta_AB_A_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_A_2 <- mean(v_Delta_AB_A_2_boot)
					
					v_Delta_AB_A_2_inf_CI <- quantile(v_Delta_AB_A_2_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_A_2_sup_CI <- quantile(v_Delta_AB_A_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_AB <- mean(v_D_AB_AB_boot)
					
					v_D_AB_AB_inf_CI <- quantile(v_D_AB_AB_boot, probs = 0.025)[[1]]
					
					v_D_AB_AB_sup_CI <- quantile(v_D_AB_AB_boot, probs = 0.975)[[1]]
					
					v_D_AB_AB_2 <- mean(v_D_AB_AB_2_boot)
					
					v_D_AB_AB_2_inf_CI <- quantile(v_D_AB_AB_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_AB_2_sup_CI <- quantile(v_D_AB_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_p_B_Delta_AB <- mean(v_p_A_p_B_Delta_AB_boot)
					
					v_p_A_p_B_Delta_AB_inf_CI <- quantile(v_p_A_p_B_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_Delta_AB_sup_CI <- quantile(v_p_A_p_B_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_p_B_Delta_AB_2 <- mean(v_p_A_p_B_Delta_AB_2_boot)
					
					v_p_A_p_B_Delta_AB_2_inf_CI <- quantile(v_p_A_p_B_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_Delta_AB_2_sup_CI <- quantile(v_p_A_p_B_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_Delta_AB <- mean(v_p_A_Delta_AB_boot)
					
					v_p_A_Delta_AB_inf_CI <- quantile(v_p_A_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_Delta_AB_sup_CI <- quantile(v_p_A_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_Delta_AB_2 <- mean(v_p_A_Delta_AB_2_boot)
					
					v_p_A_Delta_AB_2_inf_CI <- quantile(v_p_A_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_Delta_AB_2_sup_CI <- quantile(v_p_A_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_B_Delta_AB <- mean(v_p_B_Delta_AB_boot)
					
					v_p_B_Delta_AB_inf_CI <- quantile(v_p_B_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_B_Delta_AB_sup_CI <- quantile(v_p_B_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_B_Delta_AB_2 <- mean(v_p_B_Delta_AB_2_boot)
					
					v_p_B_Delta_AB_2_inf_CI <- quantile(v_p_B_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_B_Delta_AB_2_sup_CI <- quantile(v_p_B_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_Delta_AB <- mean(v_p_A_plus_p_B_Delta_AB_boot)
					
					v_p_A_plus_p_B_Delta_AB_inf_CI <- quantile(v_p_A_plus_p_B_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_Delta_AB_sup_CI <- quantile(v_p_A_plus_p_B_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_Delta_AB_2 <- mean(v_p_A_plus_p_B_Delta_AB_2_boot)
					
					v_p_A_plus_p_B_Delta_AB_2_inf_CI <- quantile(v_p_A_plus_p_B_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_Delta_AB_2_sup_CI <- quantile(v_p_A_plus_p_B_Delta_AB_2_boot, probs = 0.975)[[1]]
					

					write.table(data.frame(Delta_AB_scaled_mean = v_Delta_AB, Delta_AB_scaled_mean_inf_CI = v_Delta_AB_inf_CI, Delta_AB_scaled_mean_sup_CI = v_Delta_AB_sup_CI,
					Delta_AB_mean_scaled = v_Delta_AB_2, Delta_AB_mean_scaled_inf_CI = v_Delta_AB_2_inf_CI, Delta_AB_mean_scaled_sup_CI = v_Delta_AB_2_sup_CI,
					Delta_AB_sqrd_scaled_mean = v_Delta_AB_sqrd, Delta_AB_sqrd_scaled_mean_inf_CI = v_Delta_AB_sqrd_inf_CI, Delta_AB_sqrd_scaled_mean_sup_CI = v_Delta_AB_sqrd_sup_CI,
					Delta_AB_sqrd_mean_scaled = v_Delta_AB_sqrd_2, Delta_AB_sqrd_mean_scaled_inf_CI = v_Delta_AB_sqrd_2_inf_CI, Delta_AB_sqrd_mean_scaled_sup_CI = v_Delta_AB_sqrd_2_sup_CI,
					D_AB_A_scaled_mean = v_D_AB_A, D_AB_A_scaled_mean_inf_CI = v_D_AB_A_inf_CI, D_AB_A_scaled_mean_sup_CI = v_D_AB_A_sup_CI,
					D_AB_A_mean_scaled = v_D_AB_A_2, D_AB_A_mean_scaled_inf_CI = v_D_AB_A_2_inf_CI, D_AB_A_mean_scaled_sup_CI = v_D_AB_A_2_sup_CI,
					D_AB_B_scaled_mean = v_D_AB_B, D_AB_B_scaled_mean_inf_CI = v_D_AB_B_inf_CI, D_AB_B_scaled_mean_sup_CI = v_D_AB_B_sup_CI,
					D_AB_B_mean_scaled = v_D_AB_B_2, D_AB_B_mean_scaled_inf_CI = v_D_AB_B_2_inf_CI, D_AB_B_mean_scaled_sup_CI = v_D_AB_B_2_sup_CI,
					Delta_AB_A_scaled_mean = v_Delta_AB_A, Delta_AB_A_scaled_mean_inf_CI = v_Delta_AB_A_inf_CI, Delta_AB_A_scaled_mean_sup_CI = v_Delta_AB_A_sup_CI,
					Delta_AB_A_mean_scaled = v_Delta_AB_A_2, Delta_AB_A_mean_scaled_inf_CI = v_Delta_AB_A_2_inf_CI, Delta_AB_A_mean_scaled_sup_CI = v_Delta_AB_A_2_sup_CI,
					D_AB_AB_scaled_mean = v_D_AB_AB, D_AB_AB_scaled_mean_inf_CI = v_D_AB_AB_inf_CI, D_AB_AB_scaled_mean_sup_CI = v_D_AB_AB_sup_CI,
					D_AB_AB_mean_scaled = v_D_AB_AB_2, D_AB_AB_mean_scaled_inf_CI = v_D_AB_AB_2_inf_CI, D_AB_AB_mean_scaled_sup_CI = v_D_AB_AB_2_sup_CI,
					p_A_p_B_Delta_AB_scaled_mean = v_p_A_p_B_Delta_AB, p_A_p_B_Delta_AB_scaled_mean_inf_CI = v_p_A_p_B_Delta_AB_inf_CI, p_A_p_B_Delta_AB_scaled_mean_sup_CI = v_p_A_p_B_Delta_AB_sup_CI,
					p_A_p_B_Delta_AB_mean_scaled = v_p_A_p_B_Delta_AB_2, p_A_p_B_Delta_AB_mean_scaled_inf_CI = v_p_A_p_B_Delta_AB_2_inf_CI, p_A_p_B_Delta_AB_mean_scaled_sup_CI = v_p_A_p_B_Delta_AB_2_sup_CI,
					p_A_Delta_AB_scaled_mean = v_p_A_Delta_AB, p_A_Delta_AB_scaled_mean_inf_CI = v_p_A_Delta_AB_inf_CI, p_A_Delta_AB_scaled_mean_sup_CI = v_p_A_Delta_AB_sup_CI,
					p_A_Delta_AB_mean_scaled = v_p_A_Delta_AB_2, p_A_Delta_AB_mean_scaled_inf_CI = v_p_A_Delta_AB_2_inf_CI, p_A_Delta_AB_mean_scaled_sup_CI = v_p_A_Delta_AB_2_sup_CI,
					p_B_Delta_AB_scaled_mean = v_p_B_Delta_AB, p_B_Delta_AB_scaled_mean_inf_CI = v_p_B_Delta_AB_inf_CI, p_B_Delta_AB_scaled_mean_sup_CI = v_p_B_Delta_AB_sup_CI,
					p_B_Delta_AB_mean_scaled = v_p_B_Delta_AB_2, p_B_Delta_AB_mean_scaled_inf_CI = v_p_B_Delta_AB_2_inf_CI, p_B_Delta_AB_mean_scaled_sup_CI = v_p_B_Delta_AB_2_sup_CI,
					p_A_plus_p_B_Delta_AB_scaled_mean = v_p_A_plus_p_B_Delta_AB, p_A_plus_p_B_Delta_AB_scaled_mean_inf_CI = v_p_A_plus_p_B_Delta_AB_inf_CI, p_A_plus_p_B_Delta_AB_scaled_mean_sup_CI = v_p_A_plus_p_B_Delta_AB_sup_CI,
					p_A_plus_p_B_Delta_AB_mean_scaled = v_p_A_plus_p_B_Delta_AB_2, p_A_plus_p_B_Delta_AB_mean_scaled_inf_CI = v_p_A_plus_p_B_Delta_AB_2_inf_CI, p_A_plus_p_B_Delta_AB_mean_scaled_sup_CI = v_p_A_plus_p_B_Delta_AB_2_sup_CI,
					Class_dist = l, mut_type = k, n_sample = n), paste(i,"_Delta_AB_mean_SIFT_or.txt",sep = ""), 
					quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
				}
	}
}

####

mut <- c("del","m_del","del_neutr")

species <- c("Cg","Co")

for (i in species)
{
	data_tot <- data.frame(Delta_AB_scaled_mean = numeric(), Delta_AB_scaled_mean_inf_CI = numeric(), Delta_AB_scaled_mean_sup_CI = numeric(),
					Delta_AB_mean_scaled = numeric(), Delta_AB_mean_scaled_inf_CI = numeric(), Delta_AB_mean_scaled_sup_CI = numeric(),
					Delta_AB_sqrd_scaled_mean = numeric(), Delta_AB_sqrd_scaled_mean_inf_CI = numeric(), Delta_AB_sqrd_scaled_mean_sup_CI = numeric(),
					Delta_AB_sqrd_mean_scaled = numeric(), Delta_AB_sqrd_mean_scaled_inf_CI = numeric(), Delta_AB_sqrd_mean_scaled_sup_CI = numeric(),
					D_AB_A_scaled_mean = numeric(), D_AB_A_scaled_mean_inf_CI = numeric(), D_AB_A_scaled_mean_sup_CI = numeric(),
					D_AB_A_mean_scaled = numeric(), D_AB_A_mean_scaled_inf_CI = numeric(), D_AB_A_mean_scaled_sup_CI = numeric(),
					D_AB_B_scaled_mean = numeric(), D_AB_B_scaled_mean_inf_CI = numeric(), D_AB_B_scaled_mean_sup_CI = numeric(),
					D_AB_B_mean_scaled = numeric(), D_AB_B_mean_scaled_inf_CI = numeric(), D_AB_B_mean_scaled_sup_CI = numeric(),
					Delta_AB_A_scaled_mean = numeric(), Delta_AB_A_scaled_mean_inf_CI = numeric(), Delta_AB_A_scaled_mean_sup_CI = numeric(),
					Delta_AB_A_mean_scaled = numeric(), Delta_AB_A_mean_scaled_inf_CI = numeric(), Delta_AB_A_mean_scaled_sup_CI = numeric(),
					D_AB_AB_scaled_mean = numeric(), D_AB_AB_scaled_mean_inf_CI = numeric(), D_AB_AB_scaled_mean_sup_CI = numeric(),
					D_AB_AB_mean_scaled = numeric(), D_AB_AB_mean_scaled_inf_CI = numeric(), D_AB_AB_mean_scaled_sup_CI = numeric(),
					p_A_p_B_Delta_AB_scaled_mean = numeric(), p_A_p_B_Delta_AB_scaled_mean_inf_CI = numeric(), p_A_p_B_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_A_p_B_Delta_AB_mean_scaled = numeric(), p_A_p_B_Delta_AB_mean_scaled_inf_CI = numeric(), p_A_p_B_Delta_AB_mean_scaled_sup_CI = numeric(),
					p_A_Delta_AB_scaled_mean = numeric(), p_A_Delta_AB_scaled_mean_inf_CI = numeric(), p_A_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_A_Delta_AB_mean_scaled = numeric(), p_A_Delta_AB_mean_scaled_inf_CI = numeric(), p_A_Delta_AB_mean_scaled_sup_CI = numeric(),
					p_B_Delta_AB_scaled_mean = numeric(), p_B_Delta_AB_scaled_mean_inf_CI = numeric(), p_B_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_B_Delta_AB_mean_scaled = numeric(), p_B_Delta_AB_mean_scaled_inf_CI = numeric(), p_B_Delta_AB_mean_scaled_sup_CI = numeric(),
					p_A_plus_p_B_Delta_AB_scaled_mean = numeric(), p_A_plus_p_B_Delta_AB_scaled_mean_inf_CI = numeric(), p_A_plus_p_B_Delta_AB_scaled_mean_sup_CI = numeric(),
					p_A_plus_p_B_Delta_AB_mean_scaled = numeric(), p_A_plus_p_B_Delta_AB_mean_scaled_inf_CI = numeric(), p_A_plus_p_B_Delta_AB_mean_scaled_sup_CI = numeric(),
					Class_dist = numeric(), mut_type = numeric(), n_sample = numeric())
			
	write.table(data_tot, paste(i,"_Delta_AB_mean_Neslia_or.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	data <- read.table(paste(i,"_Delta_AB_Neslia_or.txt",sep = ""), header=TRUE, sep="\t")
	
	library("dplyr")
	
	library(data.table)
	
	dist_class <- data %>% distinct(sort(data$Class_dist)) 
	
	dist_class <- head(dist_class[,1],-1)
	
	for (k in mut)
	{	
		if (k == "del_neutr")
		{
			data <- read.table(paste(i,"_del_neutr_Delta_AB_Neslia_or.txt",sep = ""), header=TRUE, sep="\t")
		
			dist_class <- data %>% distinct(sort(data$Class_dist)) 
		
			dist_class <- head(dist_class[,1],-1)
		}
				for (l in dist_class)
				{
					data_sub <- data[which(data$mut == k & data$Class_dist == l),]
					
					n <- nrow(data_sub)
					
					v_Delta_AB_boot <- numeric()
					v_Delta_AB_2_boot <- numeric()
					v_Delta_AB_sqrd_boot <- numeric()
					v_Delta_AB_sqrd_2_boot <- numeric()
					v_D_AB_A_boot <- numeric()
					v_D_AB_A_2_boot <- numeric()
					v_D_AB_B_boot <- numeric()
					v_D_AB_B_2_boot <- numeric()
					v_Delta_AB_A_boot <- numeric()
					v_Delta_AB_A_2_boot <- numeric()
					v_D_AB_AB_boot <- numeric()
					v_D_AB_AB_2_boot <- numeric()
					v_p_A_p_B_Delta_AB_boot <- numeric()
					v_p_A_p_B_Delta_AB_2_boot <- numeric()
					v_p_A_Delta_AB_boot <- numeric()
					v_p_A_Delta_AB_2_boot <- numeric()
					v_p_B_Delta_AB_boot <- numeric()
					v_p_B_Delta_AB_2_boot <- numeric()
					v_p_A_plus_p_B_Delta_AB_boot <- numeric()
					v_p_A_plus_p_B_Delta_AB_2_boot <- numeric()
					
					for (m in c(1:1000))
					{
						data_sub_samp <- numeric()
						
						data_sub_samp <- data_sub[sample(seq_len(n),n, TRUE),]
						
						
						v_Delta_AB_boot <- c(v_Delta_AB_boot, mean(data_sub_samp$Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_Delta_AB_2_boot <- c(v_Delta_AB_2_boot, mean(data_sub_samp$Delta_AB/data_sub_samp$pq_AB))
						
						
						v_Delta_AB_sqrd_boot <- c(v_Delta_AB_sqrd_boot, mean(data_sub_samp$Delta_AB_sqrd)/mean(data_sub_samp$pq_AB))
						
						v_Delta_AB_sqrd_2_boot <- c(v_Delta_AB_sqrd_2_boot, mean(data_sub_samp$Delta_AB_sqrd/data_sub_samp$pq_AB))
						
						
						v_D_AB_A_boot <- c(v_D_AB_A_boot, mean(data_sub_samp$D_AB_A)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_A_2_boot <- c(v_D_AB_A_2_boot, mean(data_sub_samp$D_AB_A/data_sub_samp$pq_AB))
						
						
						v_D_AB_B_boot <- c(v_D_AB_B_boot, mean(data_sub_samp$D_AB_B)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_B_2_boot <- c(v_D_AB_B_2_boot, mean(data_sub_samp$D_AB_B/data_sub_samp$pq_AB))
						
						
						v_Delta_AB_A_boot <- c(v_Delta_AB_A_boot, mean(data_sub_samp$Delta_AB_A)/mean(data_sub_samp$pq_AB))
						
						v_Delta_AB_A_2_boot <- c(v_Delta_AB_A_2_boot, mean(data_sub_samp$Delta_AB_A/data_sub_samp$pq_AB))
						
						
						v_D_AB_AB_boot <- c(v_D_AB_AB_boot, mean(data_sub_samp$D_AB_AB)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_AB_2_boot <- c(v_D_AB_AB_2_boot, mean(data_sub_samp$D_AB_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_p_B_Delta_AB_boot <- c(v_p_A_p_B_Delta_AB_boot, mean(data_sub_samp$p_A_p_B_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_p_B_Delta_AB_2_boot <- c(v_p_A_p_B_Delta_AB_2_boot, mean(data_sub_samp$p_A_p_B_Delta_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_Delta_AB_boot <- c(v_p_A_Delta_AB_boot, mean(data_sub_samp$p_A_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_Delta_AB_2_boot <- c(v_p_A_Delta_AB_2_boot, mean(data_sub_samp$p_A_Delta_AB/data_sub_samp$pq_AB))
						
						
						v_p_B_Delta_AB_boot <- c(v_p_B_Delta_AB_boot, mean(data_sub_samp$p_B_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_B_Delta_AB_2_boot <- c(v_p_B_Delta_AB_2_boot, mean(data_sub_samp$p_B_Delta_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_plus_p_B_Delta_AB_boot <- c(v_p_A_plus_p_B_Delta_AB_boot, mean(data_sub_samp$p_A_plus_p_B_Delta_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_plus_p_B_Delta_AB_2_boot <- c(v_p_A_plus_p_B_Delta_AB_2_boot, mean(data_sub_samp$p_A_plus_p_B_Delta_AB/data_sub_samp$pq_AB))
					
					}
					
					v_Delta_AB <- mean(v_Delta_AB_boot)
					
					v_Delta_AB_inf_CI <- quantile(v_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_sup_CI <- quantile(v_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_2 <- mean(v_Delta_AB_2_boot)
					
					v_Delta_AB_2_inf_CI <- quantile(v_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_2_sup_CI <- quantile(v_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_sqrd <- mean(v_Delta_AB_sqrd_boot)
					
					v_Delta_AB_sqrd_inf_CI <- quantile(v_Delta_AB_sqrd_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_sqrd_sup_CI <- quantile(v_Delta_AB_sqrd_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_sqrd_2 <- mean(v_Delta_AB_sqrd_2_boot)
					
					v_Delta_AB_sqrd_2_inf_CI <- quantile(v_Delta_AB_sqrd_2_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_sqrd_2_sup_CI <- quantile(v_Delta_AB_sqrd_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_A <- mean(v_D_AB_A_boot)
					
					v_D_AB_A_inf_CI <- quantile(v_D_AB_A_boot, probs = 0.025)[[1]]
					
					v_D_AB_A_sup_CI <- quantile(v_D_AB_A_boot, probs = 0.975)[[1]]
					
					v_D_AB_A_2 <- mean(v_D_AB_A_2_boot)
					
					v_D_AB_A_2_inf_CI <- quantile(v_D_AB_A_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_A_2_sup_CI <- quantile(v_D_AB_A_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_B <- mean(v_D_AB_B_boot)
					
					v_D_AB_B_inf_CI <- quantile(v_D_AB_B_boot, probs = 0.025)[[1]]
					
					v_D_AB_B_sup_CI <- quantile(v_D_AB_B_boot, probs = 0.975)[[1]]
					
					v_D_AB_B_2 <- mean(v_D_AB_B_2_boot)
					
					v_D_AB_B_2_inf_CI <- quantile(v_D_AB_B_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_B_2_sup_CI <- quantile(v_D_AB_B_2_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_A <- mean(v_Delta_AB_A_boot)
					
					v_Delta_AB_A_inf_CI <- quantile(v_Delta_AB_A_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_A_sup_CI <- quantile(v_Delta_AB_A_boot, probs = 0.975)[[1]]
					
					v_Delta_AB_A_2 <- mean(v_Delta_AB_A_2_boot)
					
					v_Delta_AB_A_2_inf_CI <- quantile(v_Delta_AB_A_2_boot, probs = 0.025)[[1]]
					
					v_Delta_AB_A_2_sup_CI <- quantile(v_Delta_AB_A_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_AB <- mean(v_D_AB_AB_boot)
					
					v_D_AB_AB_inf_CI <- quantile(v_D_AB_AB_boot, probs = 0.025)[[1]]
					
					v_D_AB_AB_sup_CI <- quantile(v_D_AB_AB_boot, probs = 0.975)[[1]]
					
					v_D_AB_AB_2 <- mean(v_D_AB_AB_2_boot)
					
					v_D_AB_AB_2_inf_CI <- quantile(v_D_AB_AB_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_AB_2_sup_CI <- quantile(v_D_AB_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_p_B_Delta_AB <- mean(v_p_A_p_B_Delta_AB_boot)
					
					v_p_A_p_B_Delta_AB_inf_CI <- quantile(v_p_A_p_B_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_Delta_AB_sup_CI <- quantile(v_p_A_p_B_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_p_B_Delta_AB_2 <- mean(v_p_A_p_B_Delta_AB_2_boot)
					
					v_p_A_p_B_Delta_AB_2_inf_CI <- quantile(v_p_A_p_B_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_Delta_AB_2_sup_CI <- quantile(v_p_A_p_B_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_Delta_AB <- mean(v_p_A_Delta_AB_boot)
					
					v_p_A_Delta_AB_inf_CI <- quantile(v_p_A_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_Delta_AB_sup_CI <- quantile(v_p_A_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_Delta_AB_2 <- mean(v_p_A_Delta_AB_2_boot)
					
					v_p_A_Delta_AB_2_inf_CI <- quantile(v_p_A_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_Delta_AB_2_sup_CI <- quantile(v_p_A_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_B_Delta_AB <- mean(v_p_B_Delta_AB_boot)
					
					v_p_B_Delta_AB_inf_CI <- quantile(v_p_B_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_B_Delta_AB_sup_CI <- quantile(v_p_B_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_B_Delta_AB_2 <- mean(v_p_B_Delta_AB_2_boot)
					
					v_p_B_Delta_AB_2_inf_CI <- quantile(v_p_B_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_B_Delta_AB_2_sup_CI <- quantile(v_p_B_Delta_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_Delta_AB <- mean(v_p_A_plus_p_B_Delta_AB_boot)
					
					v_p_A_plus_p_B_Delta_AB_inf_CI <- quantile(v_p_A_plus_p_B_Delta_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_Delta_AB_sup_CI <- quantile(v_p_A_plus_p_B_Delta_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_Delta_AB_2 <- mean(v_p_A_plus_p_B_Delta_AB_2_boot)
					
					v_p_A_plus_p_B_Delta_AB_2_inf_CI <- quantile(v_p_A_plus_p_B_Delta_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_Delta_AB_2_sup_CI <- quantile(v_p_A_plus_p_B_Delta_AB_2_boot, probs = 0.975)[[1]]
					

					write.table(data.frame(Delta_AB_scaled_mean = v_Delta_AB, Delta_AB_scaled_mean_inf_CI = v_Delta_AB_inf_CI, Delta_AB_scaled_mean_sup_CI = v_Delta_AB_sup_CI,
					Delta_AB_mean_scaled = v_Delta_AB_2, Delta_AB_mean_scaled_inf_CI = v_Delta_AB_2_inf_CI, Delta_AB_mean_scaled_sup_CI = v_Delta_AB_2_sup_CI,
					Delta_AB_sqrd_scaled_mean = v_Delta_AB_sqrd, Delta_AB_sqrd_scaled_mean_inf_CI = v_Delta_AB_sqrd_inf_CI, Delta_AB_sqrd_scaled_mean_sup_CI = v_Delta_AB_sqrd_sup_CI,
					Delta_AB_sqrd_mean_scaled = v_Delta_AB_sqrd_2, Delta_AB_sqrd_mean_scaled_inf_CI = v_Delta_AB_sqrd_2_inf_CI, Delta_AB_sqrd_mean_scaled_sup_CI = v_Delta_AB_sqrd_2_sup_CI,
					D_AB_A_scaled_mean = v_D_AB_A, D_AB_A_scaled_mean_inf_CI = v_D_AB_A_inf_CI, D_AB_A_scaled_mean_sup_CI = v_D_AB_A_sup_CI,
					D_AB_A_mean_scaled = v_D_AB_A_2, D_AB_A_mean_scaled_inf_CI = v_D_AB_A_2_inf_CI, D_AB_A_mean_scaled_sup_CI = v_D_AB_A_2_sup_CI,
					D_AB_B_scaled_mean = v_D_AB_B, D_AB_B_scaled_mean_inf_CI = v_D_AB_B_inf_CI, D_AB_B_scaled_mean_sup_CI = v_D_AB_B_sup_CI,
					D_AB_B_mean_scaled = v_D_AB_B_2, D_AB_B_mean_scaled_inf_CI = v_D_AB_B_2_inf_CI, D_AB_B_mean_scaled_sup_CI = v_D_AB_B_2_sup_CI,
					Delta_AB_A_scaled_mean = v_Delta_AB_A, Delta_AB_A_scaled_mean_inf_CI = v_Delta_AB_A_inf_CI, Delta_AB_A_scaled_mean_sup_CI = v_Delta_AB_A_sup_CI,
					Delta_AB_A_mean_scaled = v_Delta_AB_A_2, Delta_AB_A_mean_scaled_inf_CI = v_Delta_AB_A_2_inf_CI, Delta_AB_A_mean_scaled_sup_CI = v_Delta_AB_A_2_sup_CI,
					D_AB_AB_scaled_mean = v_D_AB_AB, D_AB_AB_scaled_mean_inf_CI = v_D_AB_AB_inf_CI, D_AB_AB_scaled_mean_sup_CI = v_D_AB_AB_sup_CI,
					D_AB_AB_mean_scaled = v_D_AB_AB_2, D_AB_AB_mean_scaled_inf_CI = v_D_AB_AB_2_inf_CI, D_AB_AB_mean_scaled_sup_CI = v_D_AB_AB_2_sup_CI,
					p_A_p_B_Delta_AB_scaled_mean = v_p_A_p_B_Delta_AB, p_A_p_B_Delta_AB_scaled_mean_inf_CI = v_p_A_p_B_Delta_AB_inf_CI, p_A_p_B_Delta_AB_scaled_mean_sup_CI = v_p_A_p_B_Delta_AB_sup_CI,
					p_A_p_B_Delta_AB_mean_scaled = v_p_A_p_B_Delta_AB_2, p_A_p_B_Delta_AB_mean_scaled_inf_CI = v_p_A_p_B_Delta_AB_2_inf_CI, p_A_p_B_Delta_AB_mean_scaled_sup_CI = v_p_A_p_B_Delta_AB_2_sup_CI,
					p_A_Delta_AB_scaled_mean = v_p_A_Delta_AB, p_A_Delta_AB_scaled_mean_inf_CI = v_p_A_Delta_AB_inf_CI, p_A_Delta_AB_scaled_mean_sup_CI = v_p_A_Delta_AB_sup_CI,
					p_A_Delta_AB_mean_scaled = v_p_A_Delta_AB_2, p_A_Delta_AB_mean_scaled_inf_CI = v_p_A_Delta_AB_2_inf_CI, p_A_Delta_AB_mean_scaled_sup_CI = v_p_A_Delta_AB_2_sup_CI,
					p_B_Delta_AB_scaled_mean = v_p_B_Delta_AB, p_B_Delta_AB_scaled_mean_inf_CI = v_p_B_Delta_AB_inf_CI, p_B_Delta_AB_scaled_mean_sup_CI = v_p_B_Delta_AB_sup_CI,
					p_B_Delta_AB_mean_scaled = v_p_B_Delta_AB_2, p_B_Delta_AB_mean_scaled_inf_CI = v_p_B_Delta_AB_2_inf_CI, p_B_Delta_AB_mean_scaled_sup_CI = v_p_B_Delta_AB_2_sup_CI,
					p_A_plus_p_B_Delta_AB_scaled_mean = v_p_A_plus_p_B_Delta_AB, p_A_plus_p_B_Delta_AB_scaled_mean_inf_CI = v_p_A_plus_p_B_Delta_AB_inf_CI, p_A_plus_p_B_Delta_AB_scaled_mean_sup_CI = v_p_A_plus_p_B_Delta_AB_sup_CI,
					p_A_plus_p_B_Delta_AB_mean_scaled = v_p_A_plus_p_B_Delta_AB_2, p_A_plus_p_B_Delta_AB_mean_scaled_inf_CI = v_p_A_plus_p_B_Delta_AB_2_inf_CI, p_A_plus_p_B_Delta_AB_mean_scaled_sup_CI = v_p_A_plus_p_B_Delta_AB_2_sup_CI,
					Class_dist = l, mut_type = k, n_sample = n), paste(i,"_Delta_AB_mean_Neslia_or.txt",sep = ""), 
					quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
				}
	}
}

####

mut <- c("del","m_del","neutr","del_neutr")

species <- c("Co")

for (i in species)
{
	data_tot <- data.frame(D_AB_scaled_mean = numeric(), D_AB_scaled_mean_inf_CI = numeric(), D_AB_scaled_mean_sup_CI = numeric(),
					D_AB_mean_scaled = numeric(), D_AB_mean_scaled_inf_CI = numeric(), D_AB_mean_scaled_sup_CI = numeric(),
					D_AB_sqrd_scaled_mean = numeric(), D_AB_sqrd_scaled_mean_inf_CI = numeric(), D_AB_sqrd_scaled_mean_sup_CI = numeric(),
					D_AB_sqrd_mean_scaled = numeric(), D_AB_sqrd_mean_scaled_inf_CI = numeric(), D_AB_sqrd_mean_scaled_sup_CI = numeric(),
					p_A_p_B_D_AB_scaled_mean = numeric(), p_A_p_B_D_AB_scaled_mean_inf_CI = numeric(), p_A_p_B_D_AB_scaled_mean_sup_CI = numeric(),
					p_A_p_B_D_AB_mean_scaled = numeric(), p_A_p_B_D_AB_mean_scaled_inf_CI = numeric(), p_A_p_B_D_AB_mean_scaled_sup_CI = numeric(),
					p_A_D_AB_scaled_mean = numeric(), p_A_D_AB_scaled_mean_inf_CI = numeric(), p_A_D_AB_scaled_mean_sup_CI = numeric(),
					p_A_D_AB_mean_scaled = numeric(), p_A_D_AB_mean_scaled_inf_CI = numeric(), p_A_D_AB_mean_scaled_sup_CI = numeric(),
					p_B_D_AB_scaled_mean = numeric(), p_B_D_AB_scaled_mean_inf_CI = numeric(), p_B_D_AB_scaled_mean_sup_CI = numeric(),
					p_B_D_AB_mean_scaled = numeric(), p_B_D_AB_mean_scaled_inf_CI = numeric(), p_B_D_AB_mean_scaled_sup_CI = numeric(),
					p_A_plus_p_B_D_AB_scaled_mean = numeric(), p_A_plus_p_B_D_AB_scaled_mean_inf_CI = numeric(), p_A_plus_p_B_D_AB_scaled_mean_sup_CI = numeric(),
					p_A_plus_p_B_D_AB_mean_scaled = numeric(), p_A_plus_p_B_D_AB_mean_scaled_inf_CI = numeric(), p_A_plus_p_B_D_AB_mean_scaled_sup_CI = numeric(),
					Class_dist = numeric(), mut_type = numeric(), n_sample = numeric())
			
	write.table(data_tot, paste(i,"_D_AB_mean_no_het.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	data <- read.table(paste(i,"_D_AB_no_het.txt",sep = ""), header=TRUE, sep="\t")
	
	library("dplyr")
	
	library(data.table)
	
	dist_class <- data %>% distinct(sort(data$Class_dist)) 
	
	dist_class <- head(dist_class[,1],-1)
	
	for (k in mut)
	{

		if (k == "del_neutr")
		{
			data <- read.table(paste(i,"_del_neutr_D_AB_no_het.txt",sep = ""), header=TRUE, sep="\t")
		
			dist_class <- data %>% distinct(sort(data$Class_dist)) 
		
			dist_class <- head(dist_class[,1],-1)
		}
				
				for (l in dist_class)
				{
					data_sub <- data[which(data$mut == k & data$Class_dist == l),]
					
					n <- nrow(data_sub)
					
					if (n==0)
					{
						next
					}
					
					v_D_AB_boot <- numeric()
					v_D_AB_2_boot <- numeric()
					v_D_AB_sqrd_boot <- numeric()
					v_D_AB_sqrd_2_boot <- numeric()
					v_p_A_p_B_D_AB_boot <- numeric()
					v_p_A_p_B_D_AB_2_boot <- numeric()
					v_p_A_D_AB_boot <- numeric()
					v_p_A_D_AB_2_boot <- numeric()
					v_p_B_D_AB_boot <- numeric()
					v_p_B_D_AB_2_boot <- numeric()
					v_p_A_plus_p_B_D_AB_boot <- numeric()
					v_p_A_plus_p_B_D_AB_2_boot <- numeric()
					
					for (m in c(1:1000))
					{
						data_sub_samp <- numeric()
						
						data_sub_samp <- data_sub[sample(seq_len(n),n, TRUE),]
						
						
						v_D_AB_boot <- c(v_D_AB_boot, mean(data_sub_samp$D_AB)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_2_boot <- c(v_D_AB_2_boot, mean(data_sub_samp$D_AB/data_sub_samp$pq_AB))
						
						
						v_D_AB_sqrd_boot <- c(v_D_AB_sqrd_boot, mean(data_sub_samp$D_AB_sqrd)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_sqrd_2_boot <- c(v_D_AB_sqrd_2_boot, mean(data_sub_samp$D_AB_sqrd/data_sub_samp$pq_AB))
												
						
						v_p_A_p_B_D_AB_boot <- c(v_p_A_p_B_D_AB_boot, mean(data_sub_samp$p_A_p_B_D_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_p_B_D_AB_2_boot <- c(v_p_A_p_B_D_AB_2_boot, mean(data_sub_samp$p_A_p_B_D_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_D_AB_boot <- c(v_p_A_D_AB_boot, mean(data_sub_samp$p_A_D_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_D_AB_2_boot <- c(v_p_A_D_AB_2_boot, mean(data_sub_samp$p_A_D_AB/data_sub_samp$pq_AB))
						
						
						v_p_B_D_AB_boot <- c(v_p_B_D_AB_boot, mean(data_sub_samp$p_B_D_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_B_D_AB_2_boot <- c(v_p_B_D_AB_2_boot, mean(data_sub_samp$p_B_D_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_plus_p_B_D_AB_boot <- c(v_p_A_plus_p_B_D_AB_boot, mean(data_sub_samp$p_A_plus_p_B_D_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_plus_p_B_D_AB_2_boot <- c(v_p_A_plus_p_B_D_AB_2_boot, mean(data_sub_samp$p_A_plus_p_B_D_AB/data_sub_samp$pq_AB))
					
					}
					
					v_D_AB <- mean(v_D_AB_boot)
					
					v_D_AB_inf_CI <- quantile(v_D_AB_boot, probs = 0.025)[[1]]
					
					v_D_AB_sup_CI <- quantile(v_D_AB_boot, probs = 0.975)[[1]]
					
					v_D_AB_2 <- mean(v_D_AB_2_boot)
					
					v_D_AB_2_inf_CI <- quantile(v_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_2_sup_CI <- quantile(v_D_AB_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_sqrd <- mean(v_D_AB_sqrd_boot)
					
					v_D_AB_sqrd_inf_CI <- quantile(v_D_AB_sqrd_boot, probs = 0.025)[[1]]
					
					v_D_AB_sqrd_sup_CI <- quantile(v_D_AB_sqrd_boot, probs = 0.975)[[1]]
					
					v_D_AB_sqrd_2 <- mean(v_D_AB_sqrd_2_boot)
					
					v_D_AB_sqrd_2_inf_CI <- quantile(v_D_AB_sqrd_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_sqrd_2_sup_CI <- quantile(v_D_AB_sqrd_2_boot, probs = 0.975)[[1]]
										
					v_p_A_p_B_D_AB <- mean(v_p_A_p_B_D_AB_boot)
					
					v_p_A_p_B_D_AB_inf_CI <- quantile(v_p_A_p_B_D_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_D_AB_sup_CI <- quantile(v_p_A_p_B_D_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_p_B_D_AB_2 <- mean(v_p_A_p_B_D_AB_2_boot)
					
					v_p_A_p_B_D_AB_2_inf_CI <- quantile(v_p_A_p_B_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_D_AB_2_sup_CI <- quantile(v_p_A_p_B_D_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_D_AB <- mean(v_p_A_D_AB_boot)
					
					v_p_A_D_AB_inf_CI <- quantile(v_p_A_D_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_D_AB_sup_CI <- quantile(v_p_A_D_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_D_AB_2 <- mean(v_p_A_D_AB_2_boot)
					
					v_p_A_D_AB_2_inf_CI <- quantile(v_p_A_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_D_AB_2_sup_CI <- quantile(v_p_A_D_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_B_D_AB <- mean(v_p_B_D_AB_boot)
					
					v_p_B_D_AB_inf_CI <- quantile(v_p_B_D_AB_boot, probs = 0.025)[[1]]
					
					v_p_B_D_AB_sup_CI <- quantile(v_p_B_D_AB_boot, probs = 0.975)[[1]]
					
					v_p_B_D_AB_2 <- mean(v_p_B_D_AB_2_boot)
					
					v_p_B_D_AB_2_inf_CI <- quantile(v_p_B_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_B_D_AB_2_sup_CI <- quantile(v_p_B_D_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_D_AB <- mean(v_p_A_plus_p_B_D_AB_boot)
					
					v_p_A_plus_p_B_D_AB_inf_CI <- quantile(v_p_A_plus_p_B_D_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_D_AB_sup_CI <- quantile(v_p_A_plus_p_B_D_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_D_AB_2 <- mean(v_p_A_plus_p_B_D_AB_2_boot)
					
					v_p_A_plus_p_B_D_AB_2_inf_CI <- quantile(v_p_A_plus_p_B_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_D_AB_2_sup_CI <- quantile(v_p_A_plus_p_B_D_AB_2_boot, probs = 0.975)[[1]]
					

					write.table(data.frame(D_AB_scaled_mean = v_D_AB, D_AB_scaled_mean_inf_CI = v_D_AB_inf_CI, D_AB_scaled_mean_sup_CI = v_D_AB_sup_CI,
					D_AB_mean_scaled = v_D_AB_2, D_AB_mean_scaled_inf_CI = v_D_AB_2_inf_CI, D_AB_mean_scaled_sup_CI = v_D_AB_2_sup_CI,
					D_AB_sqrd_scaled_mean = v_D_AB_sqrd, D_AB_sqrd_scaled_mean_inf_CI = v_D_AB_sqrd_inf_CI, D_AB_sqrd_scaled_mean_sup_CI = v_D_AB_sqrd_sup_CI,
					D_AB_sqrd_mean_scaled = v_D_AB_sqrd_2, D_AB_sqrd_mean_scaled_inf_CI = v_D_AB_sqrd_2_inf_CI, D_AB_sqrd_mean_scaled_sup_CI = v_D_AB_sqrd_2_sup_CI,
					p_A_p_B_D_AB_scaled_mean = v_p_A_p_B_D_AB, p_A_p_B_D_AB_scaled_mean_inf_CI = v_p_A_p_B_D_AB_inf_CI, p_A_p_B_D_AB_scaled_mean_sup_CI = v_p_A_p_B_D_AB_sup_CI,
					p_A_p_B_D_AB_mean_scaled = v_p_A_p_B_D_AB_2, p_A_p_B_D_AB_mean_scaled_inf_CI = v_p_A_p_B_D_AB_2_inf_CI, p_A_p_B_D_AB_mean_scaled_sup_CI = v_p_A_p_B_D_AB_2_sup_CI,
					p_A_D_AB_scaled_mean = v_p_A_D_AB, p_A_D_AB_scaled_mean_inf_CI = v_p_A_D_AB_inf_CI, p_A_D_AB_scaled_mean_sup_CI = v_p_A_D_AB_sup_CI,
					p_A_D_AB_mean_scaled = v_p_A_D_AB_2, p_A_D_AB_mean_scaled_inf_CI = v_p_A_D_AB_2_inf_CI, p_A_D_AB_mean_scaled_sup_CI = v_p_A_D_AB_2_sup_CI,
					p_B_D_AB_scaled_mean = v_p_B_D_AB, p_B_D_AB_scaled_mean_inf_CI = v_p_B_D_AB_inf_CI, p_B_D_AB_scaled_mean_sup_CI = v_p_B_D_AB_sup_CI,
					p_B_D_AB_mean_scaled = v_p_B_D_AB_2, p_B_D_AB_mean_scaled_inf_CI = v_p_B_D_AB_2_inf_CI, p_B_D_AB_mean_scaled_sup_CI = v_p_B_D_AB_2_sup_CI,
					p_A_plus_p_B_D_AB_scaled_mean = v_p_A_plus_p_B_D_AB, p_A_plus_p_B_D_AB_scaled_mean_inf_CI = v_p_A_plus_p_B_D_AB_inf_CI, p_A_plus_p_B_D_AB_scaled_mean_sup_CI = v_p_A_plus_p_B_D_AB_sup_CI,
					p_A_plus_p_B_D_AB_mean_scaled = v_p_A_plus_p_B_D_AB_2, p_A_plus_p_B_D_AB_mean_scaled_inf_CI = v_p_A_plus_p_B_D_AB_2_inf_CI, p_A_plus_p_B_D_AB_mean_scaled_sup_CI = v_p_A_plus_p_B_D_AB_2_sup_CI,
					Class_dist = l, mut_type = k, n_sample = n), paste(i,"_D_AB_mean_no_het.txt",sep = ""), 
					quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
				}
	}
}

####

mut <- c("del","m_del","del_neutr")

species <- c("Co")

for (i in species)
{
	data_tot <- data.frame(D_AB_scaled_mean = numeric(), D_AB_scaled_mean_inf_CI = numeric(), D_AB_scaled_mean_sup_CI = numeric(),
					D_AB_mean_scaled = numeric(), D_AB_mean_scaled_inf_CI = numeric(), D_AB_mean_scaled_sup_CI = numeric(),
					D_AB_sqrd_scaled_mean = numeric(), D_AB_sqrd_scaled_mean_inf_CI = numeric(), D_AB_sqrd_scaled_mean_sup_CI = numeric(),
					D_AB_sqrd_mean_scaled = numeric(), D_AB_sqrd_mean_scaled_inf_CI = numeric(), D_AB_sqrd_mean_scaled_sup_CI = numeric(),
					p_A_p_B_D_AB_scaled_mean = numeric(), p_A_p_B_D_AB_scaled_mean_inf_CI = numeric(), p_A_p_B_D_AB_scaled_mean_sup_CI = numeric(),
					p_A_p_B_D_AB_mean_scaled = numeric(), p_A_p_B_D_AB_mean_scaled_inf_CI = numeric(), p_A_p_B_D_AB_mean_scaled_sup_CI = numeric(),
					p_A_D_AB_scaled_mean = numeric(), p_A_D_AB_scaled_mean_inf_CI = numeric(), p_A_D_AB_scaled_mean_sup_CI = numeric(),
					p_A_D_AB_mean_scaled = numeric(), p_A_D_AB_mean_scaled_inf_CI = numeric(), p_A_D_AB_mean_scaled_sup_CI = numeric(),
					p_B_D_AB_scaled_mean = numeric(), p_B_D_AB_scaled_mean_inf_CI = numeric(), p_B_D_AB_scaled_mean_sup_CI = numeric(),
					p_B_D_AB_mean_scaled = numeric(), p_B_D_AB_mean_scaled_inf_CI = numeric(), p_B_D_AB_mean_scaled_sup_CI = numeric(),
					p_A_plus_p_B_D_AB_scaled_mean = numeric(), p_A_plus_p_B_D_AB_scaled_mean_inf_CI = numeric(), p_A_plus_p_B_D_AB_scaled_mean_sup_CI = numeric(),
					p_A_plus_p_B_D_AB_mean_scaled = numeric(), p_A_plus_p_B_D_AB_mean_scaled_inf_CI = numeric(), p_A_plus_p_B_D_AB_mean_scaled_sup_CI = numeric(),
					Class_dist = numeric(), mut_type = numeric(), n_sample = numeric())
			
	write.table(data_tot, paste(i,"_D_AB_mean_SIFT_or_no_het.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
	
	data <- read.table(paste(i,"_D_AB_SIFT_or_no_het.txt",sep = ""), header=TRUE, sep="\t")
	
	library("dplyr")
	
	library(data.table)
	
	dist_class <- data %>% distinct(sort(data$Class_dist)) 
	
	dist_class <- head(dist_class[,1],-1)
	
	for (k in mut)
	{

		if (k == "del_neutr")
		{
			data <- read.table(paste(i,"_del_neutr_D_AB_SIFT_or_no_het.txt",sep = ""), header=TRUE, sep="\t")
		
			dist_class <- data %>% distinct(sort(data$Class_dist)) 
		
			dist_class <- head(dist_class[,1],-1)
		}
				
				for (l in dist_class)
				{
					data_sub <- data[which(data$mut == k & data$Class_dist == l),]
					
					n <- nrow(data_sub)
					
					if (n==0)
					{
						next
					}
					
					v_D_AB_boot <- numeric()
					v_D_AB_2_boot <- numeric()
					v_D_AB_sqrd_boot <- numeric()
					v_D_AB_sqrd_2_boot <- numeric()
					v_p_A_p_B_D_AB_boot <- numeric()
					v_p_A_p_B_D_AB_2_boot <- numeric()
					v_p_A_D_AB_boot <- numeric()
					v_p_A_D_AB_2_boot <- numeric()
					v_p_B_D_AB_boot <- numeric()
					v_p_B_D_AB_2_boot <- numeric()
					v_p_A_plus_p_B_D_AB_boot <- numeric()
					v_p_A_plus_p_B_D_AB_2_boot <- numeric()
					
					for (m in c(1:1000))
					{
						data_sub_samp <- numeric()
						
						data_sub_samp <- data_sub[sample(seq_len(n),n, TRUE),]
						
						
						v_D_AB_boot <- c(v_D_AB_boot, mean(data_sub_samp$D_AB)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_2_boot <- c(v_D_AB_2_boot, mean(data_sub_samp$D_AB/data_sub_samp$pq_AB))
						
						
						v_D_AB_sqrd_boot <- c(v_D_AB_sqrd_boot, mean(data_sub_samp$D_AB_sqrd)/mean(data_sub_samp$pq_AB))
						
						v_D_AB_sqrd_2_boot <- c(v_D_AB_sqrd_2_boot, mean(data_sub_samp$D_AB_sqrd/data_sub_samp$pq_AB))
												
						
						v_p_A_p_B_D_AB_boot <- c(v_p_A_p_B_D_AB_boot, mean(data_sub_samp$p_A_p_B_D_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_p_B_D_AB_2_boot <- c(v_p_A_p_B_D_AB_2_boot, mean(data_sub_samp$p_A_p_B_D_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_D_AB_boot <- c(v_p_A_D_AB_boot, mean(data_sub_samp$p_A_D_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_D_AB_2_boot <- c(v_p_A_D_AB_2_boot, mean(data_sub_samp$p_A_D_AB/data_sub_samp$pq_AB))
						
						
						v_p_B_D_AB_boot <- c(v_p_B_D_AB_boot, mean(data_sub_samp$p_B_D_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_B_D_AB_2_boot <- c(v_p_B_D_AB_2_boot, mean(data_sub_samp$p_B_D_AB/data_sub_samp$pq_AB))
						
						
						v_p_A_plus_p_B_D_AB_boot <- c(v_p_A_plus_p_B_D_AB_boot, mean(data_sub_samp$p_A_plus_p_B_D_AB)/mean(data_sub_samp$pq_AB))
						
						v_p_A_plus_p_B_D_AB_2_boot <- c(v_p_A_plus_p_B_D_AB_2_boot, mean(data_sub_samp$p_A_plus_p_B_D_AB/data_sub_samp$pq_AB))
					
					}
					
					v_D_AB <- mean(v_D_AB_boot)
					
					v_D_AB_inf_CI <- quantile(v_D_AB_boot, probs = 0.025)[[1]]
					
					v_D_AB_sup_CI <- quantile(v_D_AB_boot, probs = 0.975)[[1]]
					
					v_D_AB_2 <- mean(v_D_AB_2_boot)
					
					v_D_AB_2_inf_CI <- quantile(v_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_2_sup_CI <- quantile(v_D_AB_2_boot, probs = 0.975)[[1]]
					
					v_D_AB_sqrd <- mean(v_D_AB_sqrd_boot)
					
					v_D_AB_sqrd_inf_CI <- quantile(v_D_AB_sqrd_boot, probs = 0.025)[[1]]
					
					v_D_AB_sqrd_sup_CI <- quantile(v_D_AB_sqrd_boot, probs = 0.975)[[1]]
					
					v_D_AB_sqrd_2 <- mean(v_D_AB_sqrd_2_boot)
					
					v_D_AB_sqrd_2_inf_CI <- quantile(v_D_AB_sqrd_2_boot, probs = 0.025)[[1]]
					
					v_D_AB_sqrd_2_sup_CI <- quantile(v_D_AB_sqrd_2_boot, probs = 0.975)[[1]]
										
					v_p_A_p_B_D_AB <- mean(v_p_A_p_B_D_AB_boot)
					
					v_p_A_p_B_D_AB_inf_CI <- quantile(v_p_A_p_B_D_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_D_AB_sup_CI <- quantile(v_p_A_p_B_D_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_p_B_D_AB_2 <- mean(v_p_A_p_B_D_AB_2_boot)
					
					v_p_A_p_B_D_AB_2_inf_CI <- quantile(v_p_A_p_B_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_p_B_D_AB_2_sup_CI <- quantile(v_p_A_p_B_D_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_D_AB <- mean(v_p_A_D_AB_boot)
					
					v_p_A_D_AB_inf_CI <- quantile(v_p_A_D_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_D_AB_sup_CI <- quantile(v_p_A_D_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_D_AB_2 <- mean(v_p_A_D_AB_2_boot)
					
					v_p_A_D_AB_2_inf_CI <- quantile(v_p_A_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_D_AB_2_sup_CI <- quantile(v_p_A_D_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_B_D_AB <- mean(v_p_B_D_AB_boot)
					
					v_p_B_D_AB_inf_CI <- quantile(v_p_B_D_AB_boot, probs = 0.025)[[1]]
					
					v_p_B_D_AB_sup_CI <- quantile(v_p_B_D_AB_boot, probs = 0.975)[[1]]
					
					v_p_B_D_AB_2 <- mean(v_p_B_D_AB_2_boot)
					
					v_p_B_D_AB_2_inf_CI <- quantile(v_p_B_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_B_D_AB_2_sup_CI <- quantile(v_p_B_D_AB_2_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_D_AB <- mean(v_p_A_plus_p_B_D_AB_boot)
					
					v_p_A_plus_p_B_D_AB_inf_CI <- quantile(v_p_A_plus_p_B_D_AB_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_D_AB_sup_CI <- quantile(v_p_A_plus_p_B_D_AB_boot, probs = 0.975)[[1]]
					
					v_p_A_plus_p_B_D_AB_2 <- mean(v_p_A_plus_p_B_D_AB_2_boot)
					
					v_p_A_plus_p_B_D_AB_2_inf_CI <- quantile(v_p_A_plus_p_B_D_AB_2_boot, probs = 0.025)[[1]]
					
					v_p_A_plus_p_B_D_AB_2_sup_CI <- quantile(v_p_A_plus_p_B_D_AB_2_boot, probs = 0.975)[[1]]
					

					write.table(data.frame(D_AB_scaled_mean = v_D_AB, D_AB_scaled_mean_inf_CI = v_D_AB_inf_CI, D_AB_scaled_mean_sup_CI = v_D_AB_sup_CI,
					D_AB_mean_scaled = v_D_AB_2, D_AB_mean_scaled_inf_CI = v_D_AB_2_inf_CI, D_AB_mean_scaled_sup_CI = v_D_AB_2_sup_CI,
					D_AB_sqrd_scaled_mean = v_D_AB_sqrd, D_AB_sqrd_scaled_mean_inf_CI = v_D_AB_sqrd_inf_CI, D_AB_sqrd_scaled_mean_sup_CI = v_D_AB_sqrd_sup_CI,
					D_AB_sqrd_mean_scaled = v_D_AB_sqrd_2, D_AB_sqrd_mean_scaled_inf_CI = v_D_AB_sqrd_2_inf_CI, D_AB_sqrd_mean_scaled_sup_CI = v_D_AB_sqrd_2_sup_CI,
					p_A_p_B_D_AB_scaled_mean = v_p_A_p_B_D_AB, p_A_p_B_D_AB_scaled_mean_inf_CI = v_p_A_p_B_D_AB_inf_CI, p_A_p_B_D_AB_scaled_mean_sup_CI = v_p_A_p_B_D_AB_sup_CI,
					p_A_p_B_D_AB_mean_scaled = v_p_A_p_B_D_AB_2, p_A_p_B_D_AB_mean_scaled_inf_CI = v_p_A_p_B_D_AB_2_inf_CI, p_A_p_B_D_AB_mean_scaled_sup_CI = v_p_A_p_B_D_AB_2_sup_CI,
					p_A_D_AB_scaled_mean = v_p_A_D_AB, p_A_D_AB_scaled_mean_inf_CI = v_p_A_D_AB_inf_CI, p_A_D_AB_scaled_mean_sup_CI = v_p_A_D_AB_sup_CI,
					p_A_D_AB_mean_scaled = v_p_A_D_AB_2, p_A_D_AB_mean_scaled_inf_CI = v_p_A_D_AB_2_inf_CI, p_A_D_AB_mean_scaled_sup_CI = v_p_A_D_AB_2_sup_CI,
					p_B_D_AB_scaled_mean = v_p_B_D_AB, p_B_D_AB_scaled_mean_inf_CI = v_p_B_D_AB_inf_CI, p_B_D_AB_scaled_mean_sup_CI = v_p_B_D_AB_sup_CI,
					p_B_D_AB_mean_scaled = v_p_B_D_AB_2, p_B_D_AB_mean_scaled_inf_CI = v_p_B_D_AB_2_inf_CI, p_B_D_AB_mean_scaled_sup_CI = v_p_B_D_AB_2_sup_CI,
					p_A_plus_p_B_D_AB_scaled_mean = v_p_A_plus_p_B_D_AB, p_A_plus_p_B_D_AB_scaled_mean_inf_CI = v_p_A_plus_p_B_D_AB_inf_CI, p_A_plus_p_B_D_AB_scaled_mean_sup_CI = v_p_A_plus_p_B_D_AB_sup_CI,
					p_A_plus_p_B_D_AB_mean_scaled = v_p_A_plus_p_B_D_AB_2, p_A_plus_p_B_D_AB_mean_scaled_inf_CI = v_p_A_plus_p_B_D_AB_2_inf_CI, p_A_plus_p_B_D_AB_mean_scaled_sup_CI = v_p_A_plus_p_B_D_AB_2_sup_CI,
					Class_dist = l, mut_type = k, n_sample = n), paste(i,"_D_AB_mean_SIFT_or_no_het.txt",sep = ""), 
					quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
				}
	}
}
