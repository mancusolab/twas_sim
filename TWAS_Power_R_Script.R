library(tidyverse)
library(ggplot2)

df_summary <- read_tsv("summary.tsv")

df <- df_summary %>%
  group_by(gwas, h2ge, linear_model, snp_model) %>% filter(ngwas == 1e5) %>%
  summarise(ci = list(mean_cl_boot(twas.p < 0.05 / 22000) %>% rename(meanval=y, lwr=ymin, upr=ymax))) %>%
  unnest (cols = c(ci))

# Facet labels
gwas_names <- as_labeller(c("fast" = "Fast GWAS", "std" = "Standard GWAS"))
snp_model_names <- as_labeller(c("1snp" = "1 SNP", "1pct" = "1% SNP Model", "10pct" = "10% SNP Model"))

# Line Plot: TWAS Power vs. h2ge
plt <- ggplot(df, aes(x= h2ge, y=meanval, ymin=lwr, ymax=upr, color=linear_model)) +
  geom_point() + geom_errorbar() + geom_line() +
  scale_x_continuous(labels = function(x) format(x, scientific=TRUE)) +
  labs(title = bquote("TWAS Power"~~(alpha==0.05 / 22000)~~""),
       subtitle = bquote("Increasing"~~ h[GE]^2 ~~"from 0.00 to 0.01, GWAS sample size = 100,000"),
       x = bquote(h[GE]^2), y = "TWAS Power") +
  theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle =45, vjust=0.75, hjust=0.5)) +
  guides(color=guide_legend(title="Linear Model")) +
  facet_grid( gwas ~ snp_model, scales = "free", 
              labeller = labeller(gwas = as_labeller(gwas_names), snp_model = as_labeller(snp_model_names))) +
  theme(axis.text = element_text(size = 8)) + 
  guides(fill=guide_legend(title="Linear Model"), color=guide_legend(title="Linear Model"))
ggsave(file = "TWAS_Power_h2ge.pdf", plot = plt, width = 6.85, height=6.0, units = "in")
