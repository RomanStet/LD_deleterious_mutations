

#### Different overlap rates between blocks and duplication calls ccording to quality

library("dplyr")

species <- "Cg"

scaffold <- seq(1,8,1)

min_del <- c(2,6,12)

folder_overlap <- paste(".../7_SV_Analysis/5_SV_Analysis/match_duplic_",species,sep="")

quant_delly <- c(1,260,449,1011)

quant_smoove <- c(1,158.35,269.84,544.8)

quant_delly_smoove <- data.frame(threshold=c(100,50,30,10),quant_smoove,quant_delly)

for (m in species)
{
	for (d in min_del)
	{
		data_blocks <- numeric()
		
		for (j in scaffold)
		{
			blocks <- read.table(paste(species,j,"_bloc_stats_",d,".txt",sep = ""), header=TRUE, sep="\t")
			
			data_blocks <- rbind(data_blocks, blocks)
		}
		
		data_blocks <- data_blocks[which(data_blocks$nb_het != 0),]
			
		v_PASS_smoove <- rep(0,nrow(data_blocks))
		v_PASS_delly <- rep(0,nrow(data_blocks))
		v_chr <- numeric()
		v_block <- numeric()
		v_nb_del <- numeric()
		
			for (i in c(1:nrow(data_blocks)))
			{
				overlap_smoove <- read.table(paste(folder_overlap,data_blocks$chr[i],"/",m,data_blocks$chr[i],"_duplic_PASS_bloc_",data_blocks$bloc[i],"_smoove.txt",sep = ""), header=TRUE, sep="\t")
				overlap_delly <- read.table(paste(folder_overlap,data_blocks$chr[i],"/",m,data_blocks$chr[i],"_duplic_PASS_bloc_",data_blocks$bloc[i],"_delly.txt",sep = ""), header=TRUE, sep="\t")
				
				ind_het_smoove <- unique(overlap_smoove$ind)
				ind_het_delly <- unique(overlap_delly$ind)
				
				PASS_smoove <- rep(0,length(ind_het_smoove))
				PASS_delly <- rep(0,length(ind_het_delly))
				
				for (l in c(1:nrow(quant_delly_smoove)))
				{
					for (k in c(1:length(ind_het_smoove)))
					{
						if (nrow(na.omit(overlap_smoove)[which((na.omit(overlap_smoove)$ind==ind_het_smoove[k]) & (na.omit(overlap_smoove)$QUAL>quant_delly_smoove$quant_smoove[l])),])!=0)
						{
							PASS_smoove[k] <- quant_delly_smoove$threshold[l]
						}
					}
				
					for (k in c(1:length(ind_het_delly)))
					{
						if (nrow(na.omit(overlap_delly)[which((na.omit(overlap_delly)$ind==ind_het_delly[k]) & (na.omit(overlap_delly)$QUAL>quant_delly_smoove$quant_delly[l])),])!=0)
						{
							PASS_delly[k] <- quant_delly_smoove$threshold[l]
						}
					}
				}
				
				for (l in c(1:nrow(quant_delly_smoove)))
				{
					if (length(PASS_smoove[which((PASS_smoove<=quant_delly_smoove$threshold[l]) & (PASS_smoove!= 0))])/length(ind_het_smoove) > 0.5)
					{
						v_PASS_smoove[i] <- quant_delly_smoove$threshold[l]
					} 
					
					if (length(PASS_delly[which((PASS_delly<=quant_delly_smoove$threshold[l]) & (PASS_delly!= 0))])/length(ind_het_delly) > 0.5)
					{
						v_PASS_delly[i] <- quant_delly_smoove$threshold[l]
					} 
				}
				
				v_chr <- c(v_chr,data_blocks$chr[i])
				v_block <- c(v_block,data_blocks$bloc[i])
				v_nb_del <- c(v_nb_del,data_blocks$nb_del[i])
			}
		
		
		write.table(data.frame(scaffold=v_chr, block=v_block, nb_del=v_nb_del, filter_delly_PASS=v_PASS_delly,filter_smoove_PASS=v_PASS_smoove), paste(m,"_duplic_overlap_block.txt",sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
	}
}	


#### Density of blocks

annotated_folder <- ".../6_LD_Analysis/"

species <- "Cg"

scaffold <- seq(1,8,1)

for (m in species)
{
	for (j in scaffold)
	{
		data_POS <- read.table(paste(annotated_folder,m,j,"_annotated.txt",sep = ""), header=TRUE, sep="\t")
		
		data_POS <- unique(data_POS$POS)
		
		for (d in min_del)
		{
			data_blocks <- read.table(paste(m,j,"_bloc_stats_",d,".txt",sep = ""), header=TRUE, sep="\t")
			
			data_density <- data.frame(POS=numeric(), density=numeric(), num_del=numeric())
			
			v_density <- rep(0,length(data_POS))
			
			v_num_del <- rep(0,length(data_POS))
			
			for (i in c(1:nrow(data_blocks)))
			{
				v_density[match(data_POS[which((data_POS>=data_blocks$start[i]) & (data_POS<=data_blocks$end[i]))],data_POS)]<-1
				v_num_del[match(data_POS[which((data_POS>=data_blocks$start[i]) & (data_POS<=data_blocks$end[i]))],data_POS)]<-data_blocks$nb_del[i]
			}
			
			write.table(data.frame(POS=data_POS,density=v_density, num_del=v_num_del), paste(m,j,"_block_density_",d,"_del.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
		}
	}
}

#### Proportion of deleterious mutations inside blocks

annotated_folder <- ".../6_LD_Analysis/"

species <- "Cg"

scaffold <- seq(1,8,1)

	for (m in species)
	{
		v_scaffold <- numeric()
		v_ratio_del <- numeric()
		v_nb_del <- numeric()
		v_Delta_AB_del <- numeric()
		v_size <- numeric()
		
		write.table(data.frame(scaffold = v_scaffold, nb_del = v_nb_del, ratio_del_annotated = v_ratio_del, Delta_AB_del=v_Delta_AB_del, size= v_size), 
		paste(m,"_corr_del_blocks.txt",sep = ""), quote = FALSE, append = FALSE, row.names = FALSE, sep = "\t")

		
		for (j in scaffold)
		{
			data_POS_m_del <- read.table(paste(annotated_folder,m,j,"_m_del_freq_2.txt",sep = ""), header=TRUE, sep="\t")
			data_POS_m_del <- unique(data_POS_m_del$POS)
			
			data_POS_neutr <- read.table(paste(annotated_folder,m,j,"_neutr_freq_2.txt",sep = ""), header=TRUE, sep="\t")
			data_POS_neutr <- unique(data_POS_neutr$POS)
				
			data_blocks <- read.table(paste(m,j,"_bloc_stats_2.txt",sep = ""), header=TRUE, sep="\t")
			
			data_Delta_AB_del <- read.table(paste(annotated_folder,m,j,"_del_Delta_AB.txt",sep = ""), header=TRUE, sep="\t")
										
			v_ratio_del <- numeric()
			v_nb_del <- numeric()
			v_Delta_AB_del <- numeric()
			v_size <- numeric()
						
			for (i in c(1:nrow(data_blocks)))
			{
				nb_m_del <- length(data_POS_m_del[which((data_POS_m_del>=data_blocks$start[i]) & (data_POS_m_del<=data_blocks$end[i]))])
				
				nb_neutr <- length(data_POS_neutr[which((data_POS_neutr>=data_blocks$start[i]) & (data_POS_neutr<=data_blocks$end[i]))])
				
				sub_data_Delta_AB_del <- data_Delta_AB_del[which(((data_Delta_AB_del$pos_A>=data_blocks$start[i]) & (data_Delta_AB_del$pos_A<=data_blocks$end[i])) & 
				((data_Delta_AB_del$pos_B>=data_blocks$start[i]) & (data_Delta_AB_del$pos_B<=data_blocks$end[i]))),]
				
				v_Delta_AB_del <- c(v_Delta_AB_del, mean(sub_data_Delta_AB_del$Delta_AB)/mean(sub_data_Delta_AB_del$pq_AB))
				v_ratio_del <- c(v_ratio_del, data_blocks$nb_del[i]/(data_blocks$nb_del[i]+nb_m_del+nb_neutr))
				v_nb_del <- c(v_nb_del, data_blocks$nb_del[i])
				v_size <- c(v_size, data_blocks$end[i]-data_blocks$start[i])
			}
					
			write.table(data.frame(scaffold = rep(j, length(v_ratio_del)), nb_del=v_nb_del, ratio_del_annotated=v_ratio_del, Delta_AB_del=v_Delta_AB_del, size= v_size), paste(m,"_corr_del_blocks.txt",sep = "") 
						, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
		}
	}
	
#### Proportion of deleterious mutations inside blocks compared to outside blocks

annotated_folder <- ".../6_LD_Analysis/"

species <- "Cg"

scaffold <- seq(1,8,1)

	for (m in species)
	{
		v_scaffold <- numeric()
		v_ratio_del <- numeric()
		v_block <- numeric()
		v_nb_del <- numeric()
		v_nb_m_del <- numeric()
		v_nb_neutr <- numeric()
		
		for (j in scaffold)
		{
			data_POS_del <- read.table(paste(annotated_folder,m,j,"_del_freq_2.txt",sep = ""), header=TRUE, sep="\t")
			data_POS_del <- unique(data_POS_del$POS)
		
			data_POS_m_del <- read.table(paste(annotated_folder,m,j,"_m_del_freq_2.txt",sep = ""), header=TRUE, sep="\t")
			data_POS_m_del <- unique(data_POS_m_del$POS)
			
			data_POS_neutr <- read.table(paste(annotated_folder,m,j,"_neutr_freq_2.txt",sep = ""), header=TRUE, sep="\t")
			data_POS_neutr <- unique(data_POS_neutr$POS)
				
			data_blocks <- read.table(paste(m,j,"_bloc_stats_0.txt",sep = ""), header=TRUE, sep="\t")
			
			list_POS_del_block <- numeric()
			list_POS_m_del_block <- numeric()
			list_POS_neutr_block <- numeric()
			list_POS_del_no_block <- numeric()
			list_POS_m_del_no_block <- numeric()
			list_POS_neutr_no_block <- numeric()

			for (i in c(1:nrow(data_blocks)))
			{				
				list_POS_del_block <- c(list_POS_del_block, data_POS_del[which((data_POS_del>=data_blocks$start[i]) & (data_POS_del<=data_blocks$end[i]))])
				
				list_POS_m_del_block <- c(list_POS_m_del_block, data_POS_m_del[which((data_POS_m_del>=data_blocks$start[i]) & (data_POS_m_del<=data_blocks$end[i]))])
				
				list_POS_neutr_block <- c(list_POS_neutr_block, data_POS_neutr[which((data_POS_neutr>=data_blocks$start[i]) & (data_POS_neutr<=data_blocks$end[i]))])
			}

			list_POS_del_block <- unique(list_POS_del_block)
			list_POS_m_del_block <- unique(list_POS_m_del_block)
			list_POS_neutr_block <- unique(list_POS_neutr_block)
			
			v_ratio_del <- c(v_ratio_del, length(list_POS_del_block)/(length(list_POS_del_block)+length(list_POS_m_del_block)+length(list_POS_neutr_block)))
			v_nb_del <- c(v_nb_del, length(list_POS_del_block))
			v_nb_m_del <- c(v_nb_m_del, length(list_POS_m_del_block))
			v_nb_neutr <- c(v_nb_neutr, length(list_POS_neutr_block))
			v_scaffold <- c(v_scaffold,j)
			v_block <- c(v_block,1)
			
			v_ratio_del <- c(v_ratio_del, (length(data_POS_del)-length(list_POS_del_block))/((length(data_POS_del)-length(list_POS_del_block))+
			(length(data_POS_m_del)-length(list_POS_m_del_block))+(length(data_POS_neutr)-length(list_POS_neutr_block))))
			v_nb_del <- c(v_nb_del, length(data_POS_del)-length(list_POS_del_block))
			v_nb_m_del <- c(v_nb_m_del, length(data_POS_m_del)-length(list_POS_m_del_block))
			v_nb_neutr <- c(v_nb_neutr, length(data_POS_neutr)-length(list_POS_neutr_block))
			v_scaffold <- c(v_scaffold,j)
			v_block <- c(v_block,0)
			
		}			
			write.table(data.frame(scaffold = v_scaffold, ratio_del_annotated=v_ratio_del, nb_del=v_nb_del, nb_m_del=v_nb_m_del, 
			nb_neutr=v_nb_neutr, block=v_block), paste(m,"_ratio_del_block_no_block.txt",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
	}
	