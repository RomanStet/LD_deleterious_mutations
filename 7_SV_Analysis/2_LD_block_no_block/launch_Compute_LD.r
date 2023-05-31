
setwd(".../7_SV_Analysis/2_LD_block_no_block/")

### Keeps only sites outside of blocks

folder_block <- ".../7_SV_Analysis/1_Block_detection/"

min_del <- c(2,6,12)
species <- "Cg"
scaffold <- seq(1,8,1)
mut <- ("del","m_del","neutr","del_neutr")

for (d in min_del)
{
	for (m in species)
	{
		for (j in scaffold)
		{
			blocks <- read.table(paste(folder_block,m,j,"_bloc_stats_",d,".txt",sep = ""), header=TRUE, sep="\t")
			
			for (k in mut)
			{
				data_Delta <- read.table(paste(m,j,"_",k,"_Delta_AB.txt",sep = ""), header=TRUE, sep="\t")
				
				for (i in c(1:nrow(blocks)))
				{
					data_Delta <- data_Delta[which(!((data_Delta$pos_A>=blocks$start[i] & data_Delta$pos_A<=blocks$end[i]) | (data_Delta$pos_B>=blocks$start[i] & data_Delta$pos_B<=blocks$end[i]))),]
				}
				
				write.table(data_Delta, paste(m,j,"_",k,"_Delta_AB_no_bloc_",d,"_del.txt",sep = ""), quote = FALSE, row.names=FALSE, sep = "\t")
			}
		}
	}
}

### Keeps only sites inside blocks

library("dplyr")

folder_block <- ".../7_SV_Analysis/1_Block_detection/"

min_del <- c(2,6,12)
species <- "Cg"
scaffold <- seq(1,8,1)
mut <- ("del","m_del","neutr","del_neutr")

for (d in min_del)
{
	for (m in species)
	{
		for (j in scaffold)
		{
			blocks <- read.table(paste(folder_block,m,j,"_bloc_stats_",d,".txt",sep = ""), header=TRUE, sep="\t")
			
			for (k in mut)
			{
				data_Delta_tot <- read.table(paste(m,j,"_",k,"_Delta_AB.txt",sep = ""), header=TRUE, sep="\t")
				
				data_Delta_tot$pair <- c(1:nrow(data_Delta_tot))
				
				data_Delta_no_bloc <- data_Delta_tot
				
				for (i in c(1:nrow(blocks)))
				{
					data_Delta_no_bloc <- data_Delta_no_bloc[which(!((data_Delta_no_bloc$pos_A>=blocks$start[i] & data_Delta_no_bloc$pos_A<=blocks$end[i]) | (data_Delta_no_bloc$pos_B>=blocks$start[i] & data_Delta_no_bloc$pos_B<=blocks$end[i]))),]
				}
				
				data_Delta_bloc <- data_Delta_tot[which(!(data_Delta_tot$pair %in% data_Delta_no_bloc$pair)),]
				
				write.table(data_Delta_bloc[,-ncol(data_Delta_bloc)], paste(m,j,"_",k,"_Delta_AB_bloc_",d,"_del.txt",sep = ""), quote = FALSE, row.names=FALSE, sep = "\t")
			}
		}
	}
}

####

species <- "Cg"
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr")
min_del <- c(0,2,6,12)

for (d in min_del)
{
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
				
		write.table(data_tot, paste(i,"_Delta_AB_no_bloc_",d,"_del.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
		
		write.table(data_tot, paste(i,"_del_neutr_Delta_AB_no_bloc_",d,"_del.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
		
		for (k in mut)
		{	
			for (j in scaffold)
			{
				data <- read.table(paste(i,j,"_",k,"_Delta_AB_no_bloc_",d,"_del.txt",sep = ""), header=TRUE, sep="\t")
				
				data$scaffold <- rep(j,nrow(data))
				
				data$mut_type <- rep(k,nrow(data))

				write.table(data, paste(i,"_Delta_AB_no_bloc_",d,"_del.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
			}
		}
		
		for (j in scaffold)
		{
			data <- read.table(paste(i,j,"_del_neutr_Delta_AB_no_bloc_",d,"_del.txt",sep = ""), header=TRUE, sep="\t")
				
			data$scaffold <- rep(j,nrow(data))
				
			data$mut_type <- rep("del_neutr",nrow(data))

			write.table(data, paste(i,"_del_neutr_Delta_AB_no_bloc_",d,"_del.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
		}
	}
}

####

species <- "Cg"
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr")
min_del <- c(0,2,6,12)

for (d in min_del)
{
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
				
		write.table(data_tot, paste(i,"_Delta_AB_bloc_",d,"_del.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
		
		write.table(data_tot, paste(i,"_del_neutr_Delta_AB_bloc_",d,"_del.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
		
		for (k in mut)
		{	
			for (j in scaffold)
			{
				data <- read.table(paste(i,j,"_",k,"_Delta_AB_bloc_",d,"_del.txt",sep = ""), header=TRUE, sep="\t")
				
				data$scaffold <- rep(j,nrow(data))
				
				data$mut_type <- rep(k,nrow(data))

				write.table(data, paste(i,"_Delta_AB_bloc_",d,"_del.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
			}
		}
		
		for (j in scaffold)
		{
			data <- read.table(paste(i,j,"_del_neutr_Delta_AB_bloc_",d,"_del.txt",sep = ""), header=TRUE, sep="\t")
				
			data$scaffold <- rep(j,nrow(data))
				
			data$mut_type <- rep("del_neutr",nrow(data))

			write.table(data, paste(i,"_del_neutr_Delta_AB_bloc_",d,"_del.txt",sep = ""), quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
		}
	}
}

####

species <- "Cg"
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr","del_neutr")
min_del <- c(0,2,6,12)

for (i in species)
{
	for (d in min_del)
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
				
		write.table(data_tot, paste(i,"_Delta_AB_mean_no_bloc_",d,"_del.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
		
		data <- read.table(paste(i,"_Delta_AB_no_bloc_",d,"_del.txt",sep = ""), header=TRUE, sep="\t")
		
		library("dplyr")
		
		library(data.table)
		
		dist_class <- data %>% distinct(sort(data$Class_dist)) 
		
		dist_class <- head(dist_class[,1],-1)
		
		for (k in mut)
		{	
			if (k == "del_neutr")
			{
				data <- read.table(paste(i,"_del_neutr_Delta_AB_no_bloc_",d,"_del.txt",sep = ""), header=TRUE, sep="\t")
			
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
						Class_dist = l, mut_type = k, n_sample = n), paste(i,"_Delta_AB_mean_no_bloc_",d,"_del.txt",sep = ""), 
						quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
					}
		}
	}
}

####

species <- "Cg"
scaffold <- seq(1,8,1)
mut <- c("del","m_del","neutr","del_neutr")
min_del <- c(0,2,6,12)

for (i in species)
{
	for (d in min_del)
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
				
		write.table(data_tot, paste(i,"_Delta_AB_mean_bloc_",d,"_del.txt",sep = ""), quote = FALSE, append = FALSE, sep = "\t")
		
		data <- read.table(paste(i,"_Delta_AB_bloc_",d,"_del.txt",sep = ""), header=TRUE, sep="\t")
		
		library("dplyr")
		
		library(data.table)
		
		dist_class <- data %>% distinct(sort(data$Class_dist)) 
		
		dist_class <- head(dist_class[,1],-1)
		
		for (k in mut)
		{	
			if (k == "del_neutr")
			{
				data <- read.table(paste(i,"_del_neutr_Delta_AB_bloc_",d,"_del.txt",sep = ""), header=TRUE, sep="\t")
			
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
						Class_dist = l, mut_type = k, n_sample = n), paste(i,"_Delta_AB_mean_bloc_",d,"_del.txt",sep = ""), 
						quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
					}
		}
	}
}