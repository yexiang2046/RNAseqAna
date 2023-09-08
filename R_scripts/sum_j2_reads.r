library(tidyverse)

files <- list.files(path = "/home/xiang/New_disk2/vacv_j2_ip_yz/Read1/", pattern = "tx.type", include.dirs = TRUE)

files <- paste0("/home/xiang/New_disk2/vacv_j2_ip_yz/Read1/", files)


sample_count <- lapply(files, function(x){read.table(file = x)})

for (i in 2:length(sample_count)){
    if (i == 2){
        y <- sample_count[[1]]
    }
    y <- full_join(y, sample_count[[i]], by = "V2")
}

y <- y[, c(2, 1, 3:21)]
colnames(y) <- c("type", paste0("9221", "-S", 1:12), paste0("9814", "-S", 1:8))
colnames(y) <- c("Cat", "WT_uninf_inputA", "WT_uninf_inputB", "WT_inf24h_inputA", "WT_inf24h_inputB", "KO_inf24h_inputA", "KO_inf24h_inputB",
                 "WT_uninf_ipA", "WT_uninf_ipB", "WT_inf24h_ipA", "WT_inf24_ipB", "KO_inf24h_ipA", "KO_inf24h_ipB",
                 "WT_inf6h_inputA", "WT_inf6h_inputB", "KO_inf6h_inputA", "KO_inf6h_inputB", 
                 "WT_inf6h_ipA", "WT_inf6h_ipB", "KO_inf6h_ipA", "KO_inf6h_ipB")
# add mtRNA reads
mtRNA <- read.table(file = "VACV_bin500_counts")
colnames(mtRNA) <- c("Bins", "chr", "start", "end", "strand", "length", "WT_uninf_inputA", "WT_uninf_inputB", "WT_inf24h_inputA", "WT_inf24h_inputB", "KO_inf24h_inputA", "KO_inf24h_inputB",
                 "WT_uninf_ipA", "WT_uninf_ipB", "WT_inf24h_ipA", "WT_inf24_ipB", "KO_inf24h_ipA", "KO_inf24h_ipB",
                 "WT_inf6h_inputA", "WT_inf6h_inputB", "KO_inf6h_inputA", "KO_inf6h_inputB", 
                 "WT_inf6h_ipA", "WT_inf6h_ipB", "KO_inf6h_ipA", "KO_inf6h_ipB")

mtRNA <- mtRNA[mtRNA$chr == "chrM",]
mtRNA <- mtRNA[, 7:26]
mtRNA_df <- as.vector(sapply(mtRNA, function(x){sum(as.numeric(x))}))
mtRNA_df <- c("mtRNA", mtRNA_df)

# add VACV reads
VACV_read <- read.table(file = "VACV_bin500_counts")
colnames(VACV_read) <- c("Bins", "chr", "start", "end", "strand", "length", "WT_uninf_inputA", "WT_uninf_inputB", "WT_inf24h_inputA", "WT_inf24h_inputB", "KO_inf24h_inputA", "KO_inf24h_inputB",
                 "WT_uninf_ipA", "WT_uninf_ipB", "WT_inf24h_ipA", "WT_inf24_ipB", "KO_inf24h_ipA", "KO_inf24h_ipB",
                 "WT_inf6h_inputA", "WT_inf6h_inputB", "KO_inf6h_inputA", "KO_inf6h_inputB", 
                 "WT_inf6h_ipA", "WT_inf6h_ipB", "KO_inf6h_ipA", "KO_inf6h_ipB")
VACV_read <- VACV_read[VACV_read$chr == "NC_006998.1",]
VACV_read <- VACV_read[, 7:26]
VACV_df <- as.vector(sapply(VACV_read, function(x){sum(as.numeric(x))}))
VACV_df <- c("VACV_reads", VACV_df)

y <- rbind(y, mtRNA_df)
y[is.na(y)] <- 0

y <- rbind(y, VACV_df)


# write.csv(y, file = "J2_IP_catergory_reads.csv")
y[, 2:21] <- lapply(y[, 2:21], function(x){ z <- as.numeric(x)
                                            100*z/sum(z)})






write.csv(y, file = "J2_IP_category_percentage.csv")


# y[, 2:21] <- round(y[, 2:21])


y %>% pivot_longer(names_to = "Samples", values_to = "Percentage", cols = 2:21) %>%
    ggplot(aes(x = Samples, y = Percentage, fill = Cat)) +
    geom_bar(stat = "identity") +
    coord_flip()



