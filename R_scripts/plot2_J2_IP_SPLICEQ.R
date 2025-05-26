library(tidyverse)
library(tidyplots)
library(ggridges)

data_clean <- read.csv(file = "VACV_j2_ip_spliceQ_clean.csv")


data_clean_ip <- data_clean[data_clean$treatment %in% 
                                            c("mean_WT_Noninf_ip",
                                            "mean_WT_VACV24h_ip",
                                            "mean_WT_deltaE3L24h_ip"),]

data_clean_ip$treatment <- factor(data_clean_ip$treatment, 
                                levels = c("mean_WT_Noninf_ip",
                                            "mean_WT_VACV24h_ip",
                                            "mean_WT_deltaE3L24h_ip"))

data_clean_ip %>%
    tidyplot(x = score, y = treatment, color = treatment) +
    geom_density_ridges()
ggsave("j2_ip_spliceQ_ridges.pdf")

data_clean_ip %>% 
    tidyplot(x = treatment, y = score, color = treatment) %>%
    add_boxplot(show_outliers = FALSE) %>%
    #add_test_pvalue(ref.group = c("mean_WT_Noninf_ip")) %>%
    #add_data_points_beeswarm() %>%
    adjust_x_axis(labels = c("Mock", "VACV24h", "dE3L24h"))
ggsave("j2_spliceQ_plots.pdf")

data_clean_ip %>%
    #filter(is.na(score) == FALSE) %>%
    #filter(score != Inf) %>%
    tidyplot(x = treatment, y = score, color = treatment) %>%
    add_violin() %>%
    #add_test_pvalue(ref.group = c("mean_WT_Noninf_ip")) %>%
    #add_data_points_beeswarm() %>%
    adjust_x_axis(labels = c("Mock", "VACV24h", "dE3L24h")) 
    #geom_hline(yintercept = 1.0, linetype = "dashed", color = "red")
ggsave("j2_spliceQ_vlnplots.pdf")




data_clean_input <- data_clean[data_clean$treatment %in% 
                                            c("mean_WT_Noninf_input",
                                            "mean_WT_VACV24h_input",
                                            "mean_WT_deltaE3L24h_input"),]

data_clean_input$treatment <- factor(data_clean_input$treatment, 
                                levels = c("mean_WT_Noninf_input",
                                            "mean_WT_VACV24h_input",
                                            "mean_WT_deltaE3L24h_input"))

data_clean_input%>% 
    tidyplot(x = score, y = treatment, color = treatment) +
    geom_density_ridges()
ggsave("j2_input_spliceQ_ridges.pdf")

data_clean_input %>% 
    tidyplot(x = treatment, y = score, color = treatment) %>%
    add_boxplot(show_outliers = FALSE) %>%
    #add_test_pvalue(ref.group = c("mean_WT_Noninf_ip")) %>%
    #add_data_points_beeswarm() %>%
    adjust_x_axis(labels = c("Mock", "VACV24h", "dE3L24h"))
ggsave("j2_input_spliceQ_plots.pdf")

data_clean_input %>% 
    tidyplot(x = treatment, y = score, color = treatment) %>%
    add_violin() %>%
    #add_test_pvalue(ref.group = c("mean_WT_Noninf_ip")) %>%
    #add_data_points_beeswarm() %>%
    adjust_x_axis(labels = c("Mock", "VACV24h", "dE3L24h"))
    
ggsave("j2_input_spliceQ_vlnplots.pdf")

