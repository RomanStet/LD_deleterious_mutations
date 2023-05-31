
folder_input_files <- ".../6_LS_Analysis/"

library("dplyr")

species <- "Cg"
scaffold <- seq(1,8,1)

# maximum distance between pairs of sites with identical heterozygosity pattern:
dist_max <- 2000
# tolerated error on total number of individuals:
errind <- 2
# parameter modulating the tolerated difference in heterozygosity pattern between pairs of sites:
fact <- 0.3
# minimal number of sites with similar heterozygosity pattern for the region to be recorded:
nbsites_min <- 2

for (m in species)
{
	nbindtot <- 171
		
	for (j in scaffold)
	{
		listNA <- arrange(read.table(paste(folder_input_files,m,j,"_NA.txt",sep = ""), header=TRUE, sep="\t"), POS)
		genotype <- arrange(read.table(paste(folder_input_files,m,j,"_genotype.txt",sep = ""), header=TRUE, sep="\t"), POS)
		freq <- arrange(read.table(paste(folder_input_files,m,j,"_freq_2.txt",sep = ""), header=TRUE, sep="\t"), POS)
			
		headers <- data.frame(chr = numeric(), bloc = numeric(), nucl = numeric(), 
								  ind = numeric(), genotype = numeric(), mut_type = numeric())
		write.table(headers, paste(m,j,"_blocs.txt",sep = ""), quote = FALSE, sep = "\t")
			
		n_sites <- nrow(freq)
		bloc_id <- 0
		site <- 1
		list_loci <- c()
			
		while (site < n_sites)
		{
		  if (length(list_loci) > 0)
		  {
			if (freq$POS[site] > list_loci[length(list_loci)])
			{
			  list_loci <- c()
			}
		  }
		  if (!(freq$POS[site] %in% list_loci))
		  {
			sub_A <- genotype[which(genotype$POS == freq$POS[site]), ] 
			nbNA <- nrow(listNA[which(listNA$POS == freq$POS[site]),])/2
			sub_het <- sub_A[which(as.character(sub_A$genotype) == "Het"), ]
			num_Het_A <- nrow(sub_het)
			num_Hom_A <- nrow(sub_A[which(as.character(sub_A$genotype) == "Hom"), ])
				
			# if the mutation is only present in the heterozygous state, or present in all individuals and heterozygous in some:
			  
			if ((num_Hom_A == 0) || ((num_Het_A > 0) && ((num_Het_A + num_Hom_A + nbNA) >= (nbindtot - errind))))
			{
			  # list of heterozygous individuals:
			  list_het <- sub_het$ind
			  # number of heterozygous individuals:
			  num_het <- length(list_het)
			  # site number:
			  site_het <- freq$POS[site]
			  # block number:
			  bloc_id <- bloc_id + 1
			  nr <- nrow(sub_A)
				  
			  # bloc data.frame holds scaffold nb, bloc nb, site nb and genotypes of the different individuals:
			  bloc <- data.frame(rep(j, nr), rep(bloc_id, nr), rep(site_het, nr), sub_A$ind, sub_A$genotype, sub_A$mut_type)
			  nbsites <- 1
			  site_bloc <- site + 1
				  
			  # looks for other sites with similar heterozygosity pattern in the next dist_max bp:
				
			  while ((site_bloc <= n_sites) && (freq$POS[site_bloc] - site_het <= dist_max))
			  {
				if (!(freq$POS[site_bloc] %in% list_loci))
				{
				  sub_A <- genotype[which(genotype$POS == freq$POS[site_bloc]), ] 
				  nbNA <- nrow(listNA[which(listNA$POS == freq$POS[site_bloc]),])/2
				  sub_het <- sub_A[which(as.character(sub_A$genotype) == "Het"), ]
				  num_Het_A <- nrow(sub_het)
				  num_Hom_A <- nrow(sub_A[which(as.character(sub_A$genotype) == "Hom"), ])
					  
				  # the mutation must be heterozygous in the same individuals, and present only in the heterozygous state or present in all individuals:
					  
				  if ((length(intersect(sub_het$ind, list_het)) >= ceiling(num_het * (1 - fact))) 
					  && (num_Het_A <= floor(num_het * (1 + fact))) 
					  && ((num_Hom_A <= errind) || ((num_Het_A + num_Hom_A + nbNA) >= (nbindtot - errind))))
				  {
					list_loci <- c(list_loci, freq$POS[site_bloc])
					list_loci <- sort(list_loci, decreasing=FALSE)
						
					# update the list of heterozygous individuals if needed:
						
					if (num_Het_A > num_het)
					{
					  list_het <- sub_het$ind
					  num_het <- length(list_het)
					}
					# site number:
					site_het <- freq$POS[site_bloc]
					nbsites <- nbsites + 1
					nr <- nrow(sub_A)
					# add new site info to bloc data.frame:
					bloc2 <- data.frame(rep(j, nr), rep(bloc_id, nr), rep(site_het, nr), sub_A$ind, sub_A$genotype, sub_A$mut_type)
					bloc <- rbind(bloc, bloc2)
				  }
				}
				site_bloc <- site_bloc + 1
			  }
			  # after exiting from the while loop (no more site matching the pattern), write info in output file:
			  if (nbsites >= nbsites_min)
			  {
				write.table(bloc, paste(m,j,"_blocs.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
			  }
			  # if the number of sites matching the pattern is not high enough, block is discarded:
			  else
			  {
				bloc_id <- bloc_id - 1
			  }
			}
		  }
		  site <- site + 1
		}
	}
}

####

species <- "Cg"
scaffold <- seq(1,8,1)

# minimum number of SNPs in block to record it:
mini <- 2

# minimum number of deleterious mutation in block to record it:
min_del <- seq(0,2,6,12)

for (d in min_del)
{
	for (m in species)
	{
		for (j in scaffold)
		{
		  bloc_stats <- data.frame()
		  
		  blocks <- read.table(paste(m,j,"_blocs.txt",sep = ""), header=TRUE, sep="\t")
		  
		  nbloc <- blocks[nrow(blocks),]$bloc
		  for (b in c(1:nbloc))
		  {
			sub_bloc <- blocks[which(blocks$bloc == b), ] 
			site_begin <- sub_bloc[1,]$nucl
			site_end <- sub_bloc[nrow(sub_bloc),]$nucl
			sites <- unique(sub_bloc$nucl)
			nbhet <- numeric()
			mtype <- numeric()
			list_ind <- numeric()
			for (site in sites)
			{
			  sub_site <- sub_bloc[which((sub_bloc$nucl == site) & (sub_bloc$genotype == "Het")), ]
			  nbhet <- c(nbhet, nrow(sub_site))
			  mtype <- c(mtype, sub_bloc[match(site, sub_bloc[,"nucl"]),]$mut_type)
			  list_ind <- c(list_ind, sub_site$ind)
			}
			nbhet_max <- max(nbhet)
			ndel <- length(mtype[mtype == "del"])
			list_ind <- sort(unique(list_ind), decreasing=FALSE)
			sz <- length(list_ind)
			if ((length(sites) >= mini) & (ndel >= d))
			{
			  # duplic_stats holds scaffold nb, block nb, start and en sites, nb of SNPs matching the pattern, max nb of heterozygous individuals and list of those individuals:
			  bloc_b <- data.frame(rep(j, sz), rep(b, sz), rep(site_begin, sz), rep(site_end, sz), rep(length(sites), sz), rep(nbhet_max, sz), list_ind)
			  bloc_stats <- rbind(bloc_stats, bloc_b)
			}
		  }
		  colnames(bloc_stats) <- c("chr", "bloc", "start", "end", "nb_sites", "nb_het", "ind")
		  write.table(bloc_stats, paste(m,j,"_bloc_stats_ind_",d,".txt",sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
		}
	}
}

###

species <- "Cg"
scaffold <- seq(1,8,1)

# minimum number of SNPs in block to record it:
mini <- 2

# minimum number deleterious mutations in block to record it:
min_del <- seq(0,2,6,12)

for (d in min_del)
{
	for (m in species)
	{
		for (j in scaffold)
		{
		  bloc_stats <- data.frame()
		  
		  blocks <- read.table(paste(m,j,"_blocs.txt",sep = ""), header=TRUE, sep="\t")

		  nbloc <- blocks[nrow(blocks),]$bloc
		  for (b in c(1:nbloc))
		  {
			sub_bloc <- blocks[which(blocks$bloc == b), ] 
			site_begin <- sub_bloc[1,]$nucl
			site_end <- sub_bloc[nrow(sub_bloc),]$nucl
			sites <- unique(sub_bloc$nucl)
			nbhet <- numeric()
			mtype <- numeric()
			ninv <- 0
			for (site in sites)
			{
			  nbhet <- c(nbhet, nrow(sub_bloc[which((sub_bloc$nucl == site) & (sub_bloc$genotype == "Het")), ]))
			  mtype <- c(mtype, sub_bloc[match(site, sub_bloc[,"nucl"]),]$mut_type)
			  if (nrow(sub_bloc[which((sub_bloc$nucl == site) & (sub_bloc$genotype == "Hom")), ]) > 2)
			  {
				ninv <- ninv + 1
			  }
			}
			nbhet_max <- max(nbhet)
			ndel <- length(mtype[mtype == "del"])
			if ((length(sites) >= mini) & (ndel >= d))
			{
			  # bloc_stats holds scaffold nb, block nb, start and end sites, nb of mutations (total and del) matching the pattern and max nb of heterozygous individuals:
			  bloc_b <- data.frame(j, b, site_begin, site_end, length(sites), ninv, ndel, nbhet_max)
			  bloc_stats <- rbind(bloc_stats, bloc_b)
			}
		  }
		  colnames(bloc_stats) <- c("chr", "bloc", "start", "end", "nb_sites", "nb_inv", "nb_del", "nb_het")
		  write.table(bloc_stats, paste(m,j,"_bloc_stats_",d,".txt",sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
		}
	}
}
