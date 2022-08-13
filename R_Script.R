setwd("~/path to working directory/")
library(tidyverse)
library(ggplot2)

# 1. Generate genetic architectures
# here, we generate 100 simulations for each of the 192 genetic architectures (19200 rows)
N <- c(100000, 200000)
NGE <- 500
MODEL <- c(1, "1pct", "10pct")
H2G <- 0.1
H2GE <- c(0.00000, 0.00005, 0.00010, 0.00050, 0.00100, 0.00250, 0.00500, 0.01000)
LINEAR_MODEL <- c("lasso", "enet", "ridge", "trueqtl")

df_1 <- crossing(N, NGE, MODEL, H2G, H2GE, LINEAR_MODEL)
df_2 <- df_1 %>%
  mutate(sim = row_number())
df_3 <- crossing(sim = df_2$sim, n_sim = 1:100) %>%
  left_join(df_2, by = "sim") %>%
  mutate(id = row_number()) 
df_4 <- subset(df_3, select = -n_sim) %>%
  select(sim, id, everything())
write_tsv(df_4,"param_twas_sim.tsv")

# 2. Merge output summaries into one summary file
# here, we merge a total of 38400 (19200 summaries from fast GWAS mode + 19200 summaries from standard GWAS mode) 
# output summaries into one large file.
df_list_fast <- list.files(path = "/output directory/",
                           recursive = TRUE,
                           pattern = "*.fast.summary.tsv",
                           full.names = TRUE)
df_list_std <- list.files(path = "/output directory/",
                          recursive = TRUE,
                          pattern = "*.std.summary.tsv",
                          full.names = TRUE)
df_fast = df_list_fast %>% map_df(read_tsv, col_types = cols())
df_std = df_list_std %>% map_df(read_tsv, col_types = cols())
df_all <- rbind(df_fast, df_std)
write_tsv(df_all, "summary_raw.tsv")

# 3. Data cleansing
# Format output summary
df_summary_raw <- read_tsv("summary_raw.tsv") # Import file
df_summary <- df_summary_raw %>% mutate(twas.chi2 = twas.z^2) # Calculate chi2
df_summary$snp_model <- gsub('NumCausals := 1 SNPs', '1snp', df_summary$snp_model) # Rename parameters in certain columns
df_summary$snp_model <- gsub('NumCausals := 1.0% of observed SNPs', '1pct', df_summary$snp_model)
df_summary$snp_model <- gsub('NumCausals := 10.0% of observed SNPs', '10pct', df_summary$snp_model)
df_summary$linear_model <- gsub('ridge', 'GBLUP', df_summary$linear_model)
df_summary$linear_model <- gsub('lasso', 'LASSO', df_summary$linear_model)
df_summary$linear_model <- gsub('enet', 'Elastic Net', df_summary$linear_model)
df_summary$linear_model <- gsub('trueqtl', 'True eQTL', df_summary$linear_model)
df_summary <- df_summary %>% arrange(df_summary$sim, df_summary$id, df_summary$gwas) # Sort by sim and id
write_tsv(df_summary,"summary.tsv")

# 4. Line plot for power
df_summary <- read_tsv("summary.tsv")
# Bootstrap allows us to estimate distribution of population from a single sample
# Create a data frame to summarize confidence limits for twas.p < 0.05 / 22000 from nonparametric bootstrap
df_lineplot <- df_summary %>%
  group_by(gwas, h2ge, linear_model, snp_model) %>% 
  filter(ngwas == 1e5) %>%
  summarise(ci = list(mean_cl_boot(twas.p < 0.05 / 22000) %>% 
                        rename(meanval=y, lwr=ymin, upr=ymax))) %>%
  unnest (cols = c(ci))

# Plot labeller 
gwas_names <- as_labeller(c("fast" = "Fast GWAS", "std" = "Standard GWAS"))
snp_model_names <- as_labeller(c("1snp" = "1 SNP", "1pct" = "1% SNP Model", "10pct" = "10% SNP Model"))

# Plot TWAS Power vs. h2ge
plt_lineplot <- ggplot(df_lineplot, aes(x= h2ge, y=meanval, ymin=lwr, ymax=upr, color=linear_model)) +
  geom_point() + 
  geom_errorbar() + 
  geom_line() +
  scale_x_continuous(labels = function(x) format(x, scientific=TRUE)) +
  labs(title = bquote("TWAS Power"~~(alpha==2.27e-6)~~""),
       subtitle = bquote("Increasing"~~ h[GE]^2 ~~"from 0.00 to 0.01, GWAS sample size = 100,000"),
       x = bquote(h[GE]^2),
       y = "TWAS Power") +
  theme(axis.text = element_text(size = 10), 
        axis.text.x = element_text(angle =45, vjust=0.75, hjust=0.5)) +
  guides(color = guide_legend(title="Linear Model")) +
  facet_grid( gwas ~ snp_model, 
              scales = "free", 
              labeller = labeller(gwas = as_labeller(gwas_names), 
                                  snp_model = as_labeller(snp_model_names))) +
  theme(axis.text = element_text(size = 8)) +
  guides(fill = guide_legend(title="Linear Model"),
         color = guide_legend(title="Linear Model"))
ggsave(file = "TWAS_Power_h2ge.pdf", plot = plt_lineplot, width = 6.85, height=4.0, units = "in")