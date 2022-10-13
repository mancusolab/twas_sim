library(tidyverse)
library(ggplot2)
library(rsample)
library(cowplot)
library(gridExtra)
library(boot)
library(broom)
library(dplyr)
library(rlang)

# Canonical Parameters 
  # GWAS simulation mode: fast, 
  # linear model: enet,
  # SNP model: 1pct, 
  # eQTL sample size: 250,
  # GWAS sample size: 200000, 
  # h2ge: 0.0001, 
  # h2g: 0.1. 

### Main Paper
# Import data frame
setwd("~/Local Documents/Software_Output/R_Studio_Output/TWAS/Results0914")
df_mem <- read_tsv("mem.tsv")
df_barchart_mem <- df_mem %>%
  group_by(ngwas, gwas) %>%
  summarise(ngwas = mean(ngwas), 
            n=n(), 
            mean.MaxRSS = mean(MaxRSS),
            sd = sd(MaxRSS),
            se = sd / sqrt(n))
df_barchart_mem$ngwas <- factor(df_barchart_mem$ngwas, levels = c(50000, 100000, 500000), labels = c("FiftyK", "OneHundK", "FiveHundK"))
df_barchart_mem$gwas <- gsub('Fast mode', 'fast', df_barchart_mem$gwas)
df_barchart_mem$gwas <- gsub('Standard mode', 'std', df_barchart_mem$gwas)
df_barchart_cputime$ngwas <- factor(df_barchart_cputime$ngwas, levels = c(2500, 5000, 7000),  labels = c("2500", "5000", "7000"))
df_summary <- read_tsv("all.summary.tsv")
df_qq <- df_summary %>% filter(gwas.sim == "fast", ngwas == 200000 , nqtl == 250, snp_model == "1pct", h2ge == 0) 
df_h2ge_power <- df_summary %>% filter(gwas.sim == "fast", ngwas == 200000, nqtl == 250, snp_model == "1pct")

# Plot preset
# Color Preset
lm.group.colors <- c(enet="#F8766D", lasso="#7CAE00", ridge="#00BFC4", trueqtl="#C77CFF")
gwas.mode.group.colors <- c(fast="darkblue", std="steelblue3")
# Param Preset
thold <- 0.05/22000
# Plot Preset
font.size <- 10
text.angel <- 45
text.vjust <- 0.75
text.hjust <- 0.50
cowplot.margin <- margin(20, 0, 40, 0)
# Legend Names
ngwas_names1 <- c('FiftyK', 'OneHundK', 'FiveHundK')
ngwas_names2 <- c("50K", "100K", "500K")
lm_names1 <- c('enet', 'lasso', 'ridge', 'trueqtl')
lm_names2 <- c("Elastic Net", "LASSO", "GBLUP", "True eQTL")
gwas_mode_names1 <- c('fast', 'std')
gwas_mode_names2 <- c("Fast", "Standard")

# Figure 1
plot_1 <- ggplot(df_qq, aes(sample = twas.z^2, color= linear_model, fill= linear_model, group= linear_model)) +
  geom_qq(distribution = stats::qchisq, dparams=list(df=1)) +
  geom_abline(slope=1, intercept=0, color="black") +
  labs(x = "Expected",
       y = "Observed") +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="Linear Model"),
         color=guide_legend(title="Linear Model", nrow = 1, byrow = TRUE)) + 
  scale_color_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) + 
  scale_fill_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) 
plot_2 <- ggplot(df_h2ge_power %>% mutate(power = as.numeric(twas.p < thold)), aes(x=as.factor(h2ge))) + 
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", geom="point") + 
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", geom="errorbar", fun.args = list(conf.int = 0.95), width = 0.2) + 
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", geom="line") +
  ylim(0,1) +
  labs(x = bquote(""~~ h[GE]^2 ~~""), y = "Power") +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="Linear Model"),
         color=guide_legend(title="Linear Model", nrow = 1, byrow = TRUE)) + 
  scale_color_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) + 
  scale_fill_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) 
plot_3 <- ggplot(df_barchart_mem) +
  geom_bar(aes(x=ngwas, y=mean.MaxRSS/1000000, color = gwas, fill = gwas, group = gwas), 
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(x=ngwas, ymin=mean.MaxRSS-sd, ymax=mean.MaxRSS+sd, group = gwas), 
                colour="orange", position = position_dodge(width = 0.9), width=0.2, size=0.5) +
  labs(x = bquote("GWAS Sample Size"),
       y = bquote("Memory (MB)")) +
  ylim(0,9000) +
  scale_x_discrete(breaks=ngwas_names1, labels=ngwas_names2) +
  theme(axis.text = element_text(size = font.size), 
        axis.text.x = element_text(angle = 0, vjust=text.vjust, hjust=text.hjust), 
        legend.direction = "horizontal", 
        legend.title=element_text(size=font.size), 
        legend.text=element_text(size=font.size), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="GWAS Mode"),
         color=guide_legend(title="GWAS Mode", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=gwas.mode.group.colors, breaks=gwas_mode_names1, labels=gwas_mode_names2) + 
  scale_fill_manual(values=gwas.mode.group.colors, breaks=gwas_mode_names1, labels=gwas_mode_names2)
prow <- plot_grid(
  plot_1 + theme(plot.margin = margin(6, 0, 40, 0), legend.position="none"),
  plot_2 + theme(plot.margin = margin(6, 0, 40, 0), legend.position="none"),
  plot_3 + theme(plot.margin = margin(6, 0, 40, 0), legend.position="none"),#legend.position=c(0.4, 0.8)),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)
prow
legend1 <- get_legend(plot_1 + theme(legend.position = c(0.72, 2.4)))
legend3 <- get_legend(plot_3 +theme(legend.position = c(-0.6, 0.9)))
combineLegend <- plot_grid(legend1, legend3)
cow_plot <- plot_grid(prow, combineLegend, ncol = 1, rel_heights = c(1, .1))
save_plot(file = "mainpaper_plot.pdf", plot = cow_plot, base_width = 17.4, base_height=6.5, units = "cm")

### Supplements 
# Canonical Parameters 
canonical.gwas.sim <- "fast"
canonical.linear.model <- "enet"
canonical.snp.model <- "1pct"
canonical.nqtl <- 250
canonical.ngwas <- 200000
canonical.h2ge <- 0.0005

# FWER
# Plot 2: FWER vs Linear Model
df_lm_fwer <- df_summary %>% filter(gwas.sim == canonical.gwas.sim, snp_model == canonical.snp.model, nqtl == canonical.nqtl, ngwas == canonical.ngwas, h2ge == 0)
# Plot 3: FWER vs SNP Model
df_snp_fwer <- df_summary %>% filter(gwas.sim == canonical.gwas.sim, linear_model == canonical.linear.model, nqtl == canonical.nqtl, ngwas == canonical.ngwas, h2ge == 0)
df_snp_fwer$snp_model <- factor(df_snp_fwer$snp_model, levels = c("1snp", "1pct", "10pct"), labels = c("OneSNP", "OnePCT", "TenPCT"))
# Plot 4: FWER vs NGWAS
df_ngwas_fwer <- df_summary %>% filter(gwas.sim == canonical.gwas.sim, linear_model == canonical.linear.model, snp_model == canonical.snp.model, nqtl == canonical.nqtl, h2ge == 0)
df_ngwas_fwer$ngwas <- factor(df_ngwas_fwer$ngwas, levels = c(50000, 100000, 200000, 500000), labels = c("FiftyK", "OneHundK", "TwoHundK", "FiveHundK"))
# Plot 5: FWER vs NQTL
df_nqtl_fwer <- df_summary %>% filter(gwas.sim == canonical.gwas.sim, linear_model == canonical.linear.model, snp_model == canonical.snp.model, ngwas == canonical.ngwas, h2ge == 0)
df_nqtl_fwer$nqtl <- factor(df_nqtl_fwer$nqtl, levels = c(100, 250, 500), labels = c("OneHund", "TwoHundFifty", "FiveHund"))
# Plot 6: FWER vs GWAS Mode
df_gwas_sim_fwer <- df_summary %>% filter(linear_model == canonical.linear.model, snp_model == canonical.snp.model, ngwas == canonical.ngwas, nqtl == canonical.nqtl, h2ge == 0)

# Inflation
get_lgc <- function(data, B=1000) {
  tibble(est = map_dbl(bootstraps(data, B)$splits, 
                       function(x) median(as.data.frame(x)$twas.z^2)/.455)) %>%
    summarize(Inflation = mean(est), L95 = quantile(est, .025), U95 = quantile(est, 0.975))
}
inf_lm <- df_summary %>% 
  filter(gwas.sim == canonical.gwas.sim, snp_model == canonical.snp.model, nqtl == canonical.nqtl, ngwas == canonical.ngwas, h2ge == 0) %>%
  group_by(linear_model) %>%
  do(est = get_lgc(.)) %>%
  unnest(est) %>% 
  rename(parameter = linear_model)
inf_snp <- df_summary %>% 
  filter(gwas.sim == canonical.gwas.sim, linear_model == canonical.linear.model, nqtl == canonical.nqtl, ngwas == canonical.ngwas, h2ge == 0) %>%
  group_by(snp_model) %>%
  do(est = get_lgc(.)) %>%
  unnest(est) %>% 
  rename(parameter = snp_model)
inf_snp$parameter <- factor(inf_snp$parameter, levels = c("1snp", "1pct", "10pct"), labels = c("OneSNP", "OnePCT", "TenPCT"))
inf_gwas_sim <- df_summary %>% 
  filter(snp_model == canonical.snp.model, linear_model == canonical.linear.model, nqtl == canonical.nqtl, ngwas == canonical.ngwas, h2ge == 0) %>%
  group_by(gwas.sim) %>%
  do(est = get_lgc(.)) %>%
  unnest(est) %>% 
  rename(parameter = gwas.sim)
inf_nqtl <- df_summary %>% 
  filter(gwas.sim == canonical.gwas.sim, snp_model == canonical.snp.model, linear_model == canonical.linear.model, ngwas == canonical.ngwas, h2ge == 0) %>%
  group_by(nqtl) %>%
  do(est = get_lgc(.)) %>%
  unnest(est) %>% 
  rename(parameter = nqtl)
inf_nqtl$parameter <- as.character(inf_nqtl$parameter)
inf_nqtl$parameter <- factor(inf_nqtl$parameter, levels = c(100, 250, 500), labels = c("OneHund", "TwoHundFifty", "FiveHund"))
inf_ngwas <- df_summary %>% 
  filter(gwas.sim == canonical.gwas.sim, snp_model == canonical.snp.model, linear_model == canonical.linear.model, nqtl == canonical.nqtl, h2ge == 0) %>%
  group_by(ngwas) %>%
  do(est = get_lgc(.)) %>%
  unnest(est) %>% 
  rename(parameter = ngwas)
inf_ngwas$parameter <- factor(inf_ngwas$parameter, levels = c(50000, 100000, 200000, 500000), labels = c("FiftyK", "OneHundK", "TwoHundK", "FiveHundK"))
inf_ngwas$parameter <- as.character(inf_ngwas$parameter)
inf <- rbind(inf_lm, inf_snp, inf_gwas_sim, inf_nqtl, inf_ngwas)
inf$parameter <- gsub('enet', 'Elastic Net', inf$parameter)
inf$parameter <- gsub('lasso', 'LASSO', inf$parameter)
inf$parameter <- gsub('ridge', 'GBLUP', inf$parameter)
inf$parameter <- gsub('trueqtl', 'True eQTL', inf$parameter)
inf$parameter <- gsub('fast', 'Fast', inf$parameter)
inf$parameter <- gsub('std', 'Standard', inf$parameter)
inf$parameter <- gsub('1snp', '1 SNP', inf$parameter)
inf$parameter <- gsub('1pct', '1% SNPs', inf$parameter)
inf$parameter <- gsub('10pct', '10% SNPs', inf$parameter)

# Power
# Plot 1: Power vs H2GE
df_h2ge_power <- df_summary %>% filter(gwas.sim == canonical.gwas.sim, linear_model == canonical.linear.model, snp_model == canonical.snp.model, nqtl == canonical.nqtl, ngwas == canonical.ngwas)
# Plot 2: Power vs Linear Model
df_lm_power <- df_summary %>% filter(gwas.sim == canonical.gwas.sim, snp_model == canonical.snp.model, nqtl == canonical.nqtl, ngwas == canonical.ngwas, h2ge == canonical.h2ge)
# Plot 3: Power vs SNP Model
df_snp_power <- df_summary %>% filter(gwas.sim == canonical.gwas.sim, linear_model == canonical.linear.model, nqtl == canonical.nqtl, ngwas == canonical.ngwas, h2ge == canonical.h2ge)
df_snp_power$snp_model <- factor(df_snp_power$snp_model, levels = c("1snp", "1pct", "10pct"), labels = c("OneSNP", "OnePCT", "TenPCT"))
# Plot 4: Power vs NGWAS
df_ngwas_power <- df_summary %>% filter(gwas.sim == canonical.gwas.sim, linear_model == canonical.linear.model, snp_model == canonical.snp.model, nqtl == canonical.nqtl, h2ge == canonical.h2ge)
df_ngwas_power$ngwas <- factor(df_ngwas_power$ngwas, levels = c(50000, 100000, 200000, 500000), labels = c("FiftyK", "OneHundK", "TwoHundK", "FiveHundK"))
# Plot 5: Power vs NQTL
df_nqtl_power <- df_summary %>% filter(gwas.sim == canonical.gwas.sim, linear_model == canonical.linear.model, snp_model == canonical.snp.model, ngwas == canonical.ngwas, h2ge == canonical.h2ge)
df_nqtl_power$nqtl <- factor(df_nqtl_power$nqtl, levels = c(100, 250, 500), labels = c("OneHund", "TwoHundFifty", "FiveHund"))
# Plot 6: Power vs GWAS Mode
df_gwas_sim_power <- df_summary %>% filter(linear_model == canonical.linear.model, snp_model == canonical.snp.model, ngwas == canonical.ngwas, nqtl == canonical.nqtl, h2ge == canonical.h2ge)

# CPU Time
df_cpu <- df_summary %>% filter(linear_model == canonical.linear.model, snp_model == canonical.snp.model, nqtl == canonical.nqtl, h2ge == canonical.h2ge)
df_barchart_cputime <- df_cpu %>% # Create a data frame to summarize confidence limits for twas.p < 0.05 / 22000 obtained from nonparametric bootstrap
  group_by(ngwas, gwas.sim) %>%
  summarise(sim = mean(id), ngwas = mean(ngwas), mean.cpu.time = list(mean(cpu.time))) %>% 
  unnest (cols = c(mean.cpu.time))
df_barchart_cputime$mean.cpu.time <- round(df_barchart_cputime$mean.cpu.time, 2)
df_barchart_cputime$ngwas <- factor(df_barchart_cputime$ngwas, levels = c(50000, 100000, 200000, 500000), labels = c("FiftyK", "OneHundK", "TwoHundK", "FiveHundK"))

# Plot preset
# Color Preset
lm.group.colors <- c(enet="#F8766D", lasso="#7CAE00", ridge="#00BFC4", trueqtl="#C77CFF")
snp.group.colors <- c(OneSNP="#AED0EE", OnePCT="#BBA1CB", TenPCT="#FFFBC7")
ngwas.group.colors <- c(FiftyK="#faee91", OneHundK="#edd832", TwoHundK="#bdb362", FiveHundK="#827722")
nqtl.group.colors <- c(OneHund="#5DA39D", TwoHundFifty="#32788A", FiveHund="#535164")
gwas.mode.group.colors <- c(fast="darkblue", std="steelblue3")
# Param Preset
thold <- 0.05/22000
# Plot Preset
font.size <- 10
text.angel <- 30
text.vjust <- 0.75
text.hjust <- 0.50
cowplot.margin <- margin(10, 10, 40, 10)
plot.margin <- margin(10, 10, 10, 10)
ylim.power <- c(0, 1)
ylim.fwer <- c(0, 0.15)
ylim.inf.ab <- c(0, 2)
ylim.inf.c <- c(0, 3)
ylim.inf.de <- c(0, 2.5)
# Plot Label
gwas_names <- as_labeller(c("fast" = "Fast GWAS", "std" = "Standard GWAS"))
# Legend Names
lm_names1 <- c('enet', 'lasso', 'ridge', 'trueqtl')
lm_names2 <- c("Elastic Net", "LASSO", "GBLUP", "True eQTL")
snp_names1 <- c("OneSNP", "OnePCT", "TenPCT")
snp_names2 <- c("1 SNP", "1% SNPs", "10% SNPs")
ngwas_names1 <- c('FiftyK', 'OneHundK', 'TwoHundK', 'FiveHundK')
ngwas_names2 <- c("50K", "100K", "200K", "500K")
nqtl_names1 <- c("OneHund", "TwoHundFifty", "FiveHund")
nqtl_names2 <- c("100", "250", "500")
gwas_mode_names1 <- c("fast", "std")
gwas_mode_names2 <- c("Fast", "Standard")

# Plot
# Supplementary figure 2: FWER
plot_2 <- ggplot(df_lm_fwer %>% mutate(fwer = as.numeric(twas.p < thold)), aes(x=linear_model)) +
  stat_summary(aes(y=fwer, group=linear_model, color=linear_model, fill=linear_model), fun.data = "mean_cl_boot", geom ="bar", fun.args = list(B=1000)) +
  stat_summary(aes(y=fwer, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.fwer) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.2) +
  labs(x = "Linear Model", 
       y = "FWER") +
  scale_x_discrete(breaks=lm_names1, labels=lm_names2) +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="Linear Model"),
         color=guide_legend(title="Linear Model", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) + 
  scale_fill_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2)
plot_3 <- ggplot(df_snp_fwer %>% mutate(fwer = as.numeric(twas.p < thold)), aes(x=snp_model)) +
  stat_summary(aes(y=fwer, group=snp_model, color=snp_model, fill=snp_model), fun.data = "mean_cl_boot", geom ="bar", fun.args = list(B=1000)) +
  stat_summary(aes(y=fwer, group=snp_model, color=snp_model), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.fwer) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.2) +
  labs(x = "SNP Model", 
       y = "FWER") +
  scale_x_discrete(breaks=snp_names1, labels=snp_names2) +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="SNP Model"),
         color=guide_legend(title="SNP Model", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=snp.group.colors, breaks = snp_names1, labels = snp_names2) + 
  scale_fill_manual(values=snp.group.colors, breaks = snp_names1, labels = snp_names2)
plot_4 <- ggplot(df_ngwas_fwer %>% mutate(fwer = as.numeric(twas.p < thold)), aes(x=ngwas)) +
  stat_summary(aes(y=fwer, group=ngwas, color=ngwas, fill=ngwas), fun.data = "mean_cl_boot", geom ="bar", fun.args = list(B=1000)) +
  stat_summary(aes(y=fwer, group=ngwas, color=ngwas), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.fwer) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.2) +
  labs(x = "GWAS Sample Size", y = "FWER") +
  scale_x_discrete(breaks=ngwas_names1, labels=ngwas_names2) +
  theme(axis.text = element_text(size = font.size), 
        axis.text.x = element_text(angle = 0, vjust=text.vjust, hjust=text.hjust), 
        legend.direction = "horizontal", 
        legend.title=element_text(size=font.size), 
        legend.text=element_text(size=font.size), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="GWAS Sample Size"),
         color=guide_legend(title="GWAS Sample Size", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=ngwas.group.colors, breaks=ngwas_names1, labels=ngwas_names2) + 
  scale_fill_manual(values=ngwas.group.colors, breaks=ngwas_names1, labels=ngwas_names2)
plot_5 <- ggplot(df_nqtl_fwer %>% mutate(fwer = as.numeric(twas.p < thold)), aes(x=nqtl)) +
  stat_summary(aes(y=fwer, group=nqtl, color=nqtl, fill=nqtl), fun.data = "mean_cl_boot", geom ="bar", fun.args = list(B=1000)) +
  stat_summary(aes(y=fwer, group=nqtl, color=nqtl), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.fwer) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.2) +
  labs(x = "eQTL Sample Size", y = "FWER") +
  scale_x_discrete(breaks=nqtl_names1, labels=nqtl_names2) +
  theme(axis.text = element_text(size = font.size), 
        axis.text.x = element_text(angle = 0, vjust=text.vjust, hjust=text.hjust), 
        legend.direction = "horizontal", 
        legend.title=element_text(size=font.size), 
        legend.text=element_text(size=font.size), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="eQTL Sample Size"),
         color=guide_legend(title="eQTL Sample Size", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=nqtl.group.colors, breaks=nqtl_names1, labels=nqtl_names2) + 
  scale_fill_manual(values=nqtl.group.colors, breaks=nqtl_names1, labels=nqtl_names2) 
plot_6 <- ggplot(df_gwas_sim_fwer %>% mutate(fwer = as.numeric(twas.p < thold)), aes(x=gwas.sim)) +
  stat_summary(aes(y=fwer, group=gwas.sim, color=gwas.sim, fill=gwas.sim), fun.data = "mean_cl_boot", geom ="bar", fun.args = list(B=1000)) +
  stat_summary(aes(y=fwer, group=gwas.sim, color=gwas.sim), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.fwer) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.2) +
  labs(x = "GWAS Mode", 
       y = "FWER") +
  scale_x_discrete(breaks=gwas_mode_names1, labels=gwas_mode_names2) +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =0, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="GWAS Mode"),
         color=guide_legend(title="GWAS Mode", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=gwas.mode.group.colors, breaks = gwas_mode_names1, labels = gwas_mode_names2) + 
  scale_fill_manual(values=gwas.mode.group.colors, breaks = gwas_mode_names1, labels = gwas_mode_names2)
prow <- plot_grid(
  plot_2 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_3 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_4 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_5 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_6 + theme(plot.margin = cowplot.margin, legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C", "D", "E"),
  hjust = -1,
  nrow = 2
)
prow
legend2 <- get_legend(plot_2 + theme(legend.position = c(1.685, 2), legend.justification = c("right", "top"), legend.box.just = "right",))
legend3 <- get_legend(plot_3 + theme(legend.position = c(1.95, 2), legend.justification = c("right", "top"), legend.box.just = "right",))
legend4 <- get_legend(plot_4 + theme(legend.position = c(-0.53, 1), legend.justification = c("right", "top"), legend.box.just = "right",))
legend5 <- get_legend(plot_5 + theme(legend.position = c(2.79, 2), legend.justification = c("right", "top"), legend.box.just = "right",))
legend6 <- get_legend(plot_6 + theme(legend.position = c(-.0875, 1), legend.justification = c("right", "top"), legend.box.just = "right",))
combineLegend <- plot_grid(legend2, legend3, legend4, legend5, legend6)
cow_plot <- plot_grid(prow, combineLegend, ncol = 1, rel_heights = c(1, .1))
cow_plot
save_plot(file = "fwer_strict_cowplot.pdf", plot = cow_plot, base_width = 8.5, base_height=8, units = "in")
# extract data
plot2_fwer <- ggplot_build(plot_2)$data[[1]]
plot3_fwer <- ggplot_build(plot_3)$data[[1]]
plot4_fwer <- ggplot_build(plot_4)$data[[1]]
plot5_fwer <- ggplot_build(plot_5)$data[[1]]
plot6_fwer <- ggplot_build(plot_6)$data[[1]]

# Supplementary figure 3: Inflation
plot_2 <- ggplot(inf_lm, aes(x=parameter, y=Inflation, ymin=L95, ymax=U95, group=parameter, color=parameter, fill=parameter)) +
  geom_bar(stat="identity", width=0.5) +
  geom_errorbar(width = 0.1, color="grey37") +
  geom_point(color="grey37", size=3) +
  coord_cartesian(ylim = ylim.inf.ab) +
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.2) +
  labs(x = "Linear Model",
       y = "Inflation") +
  scale_x_discrete(breaks=lm_names1, labels=lm_names2)+
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="Linear Model"),
         color=guide_legend(title="Linear Model", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) + 
  scale_fill_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2)
plot_3 <- ggplot(inf_snp, aes(x=parameter, y=Inflation, ymin=L95, ymax=U95, group=parameter, color=parameter, fill=parameter)) +
  geom_bar(stat="identity", width=0.5) +
  geom_errorbar(width = 0.1, color="grey37") +
  geom_point(color="grey37", size=3) +
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.2) +
  coord_cartesian(ylim = ylim.inf.ab) +
  labs(x = "SNP Model",
       y = "Inflation") +
  scale_x_discrete(breaks=snp_names1, labels=snp_names2)+
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="SNP Model"),
         color=guide_legend(title="SNP Model", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=snp.group.colors, breaks = snp_names1, labels = snp_names2) + 
  scale_fill_manual(values=snp.group.colors, breaks = snp_names1, labels = snp_names2)
plot_4 <- ggplot(inf_ngwas, aes(x=parameter, y=Inflation, ymin=L95, ymax=U95, group=parameter, color=parameter, fill=parameter)) +
  geom_bar(stat="identity", width=0.5) +
  geom_errorbar(width = 0.1, color="grey37") +
  geom_point(color="grey37", size=3) +
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.2) +
  coord_cartesian(ylim = ylim.inf.c) +
  labs(x = "GWAS Sample Size",
       y = "Inflation") +
  scale_x_discrete(breaks=ngwas_names1, labels=ngwas_names2)+
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =0, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="GWAS Sample Size"),
         color=guide_legend(title="GWAS Sample Size", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=ngwas.group.colors, breaks = ngwas_names1, labels = ngwas_names2) + 
  scale_fill_manual(values=ngwas.group.colors, breaks = ngwas_names1, labels = ngwas_names2)
plot_5 <- ggplot(inf_nqtl, aes(x=parameter, y=Inflation, ymin=L95, ymax=U95, group=parameter, color=parameter, fill=parameter)) +
  geom_bar(stat="identity", width=0.5) +
  geom_errorbar(width = 0.1, color="grey37") +
  geom_point(color="grey37", size=3) +  
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.2) +
  coord_cartesian(ylim = ylim.inf.de) +
  labs(x = "eQTL Sample Size",
       y = "Inflation") +
  scale_x_discrete(breaks=nqtl_names1, labels=nqtl_names2)+
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =0, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="eQTL Sample Size"),
         color=guide_legend(title="eQTL Sample Size", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=nqtl.group.colors, breaks = nqtl_names1, labels = nqtl_names2) + 
  scale_fill_manual(values=nqtl.group.colors, breaks = nqtl_names1, labels = nqtl_names2)
plot_6 <- ggplot(inf_gwas_sim, aes(x=parameter, y=Inflation, ymin=L95, ymax=U95, group=parameter, color=parameter, fill=parameter)) +
  geom_bar(stat="identity", width=0.5) +
  geom_errorbar(width = 0.1, color="grey37") +
  geom_point(color="grey37", size=3) +
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.2) +
  coord_cartesian(ylim = ylim.inf.de) +
  labs(x = "GWAS Mode",
       y = "Inflation") +
  scale_x_discrete(breaks=gwas_mode_names1, labels=gwas_mode_names2)+
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =0, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="GWAS Mode"),
         color=guide_legend(title="GWAS Mode", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=gwas.mode.group.colors, breaks = gwas_mode_names1, labels = gwas_mode_names2) + 
  scale_fill_manual(values=gwas.mode.group.colors, breaks = gwas_mode_names1, labels = gwas_mode_names2)
prow <- plot_grid(
  plot_2 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_3 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_4 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_5 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_6 + theme(plot.margin = cowplot.margin, legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C", "D", "E", "F"),
  hjust = -1,
  nrow = 2
)
prow
legend2 <- get_legend(plot_2 + theme(legend.position = c(1.685, 2), legend.justification = c("right", "top"), legend.box.just = "right",))
legend3 <- get_legend(plot_3 + theme(legend.position = c(1.95, 2),legend.justification = c("right", "top"), legend.box.just = "right",))
legend4 <- get_legend(plot_4 + theme(legend.position = c(-0.53, 1), legend.justification = c("right", "top"), legend.box.just = "right",))
legend5 <- get_legend(plot_5 + theme(legend.position = c(2.79, 2), legend.justification = c("right", "top"), legend.box.just = "right",))
legend6 <- get_legend(plot_6 + theme(legend.position = c(-.0875, 1), legend.justification = c("right", "top"), legend.box.just = "right",))
combineLegend <- plot_grid(legend2, legend3, legend4, legend5, legend6)
cow_plot <- plot_grid(prow, combineLegend, ncol = 1, rel_heights = c(1, .1))
cow_plot
save_plot(file = "inflation.pdf", plot = cow_plot, base_width = 8.5, base_height=8, units = "in")

# Supplementary figure 4: Power
plot_1 <- ggplot(df_h2ge_power %>% mutate(power = as.numeric(twas.p < thold)), aes(x=as.factor(h2ge))) + 
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", geom="point", fun.args = list(B=1000)) + 
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) + 
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", geom="line", fun.args = list(B=1000)) +
  coord_cartesian(ylim = ylim.power) +
  labs(x = bquote(""~~ h[GE]^2 ~~""), y = "Power") +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="Linear Model"),
         color=guide_legend(title="Linear Model", nrow = 1, byrow = TRUE)) + 
  scale_color_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) + 
  scale_fill_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) 
plot_2 <- ggplot(df_lm_power %>% mutate(power = as.numeric(twas.p < thold)), aes(x=linear_model)) +
  stat_summary(aes(y=power, group=linear_model, color=linear_model, fill=linear_model), fun.data = "mean_cl_boot", geom ="bar",fun.args = list(B=1000)) +
  stat_summary(aes(y=power, group=linear_model, color=linear_model), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.power) +
  labs(x = "Linear Model", 
       y = "Power") +
  scale_x_discrete(breaks=lm_names1, labels=lm_names2) +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="Linear Model"),
         color=guide_legend(title="Linear Model", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2) + 
  scale_fill_manual(values=lm.group.colors, breaks = lm_names1, labels = lm_names2)
plot_3 <- ggplot(df_snp_power %>% mutate(power = as.numeric(twas.p < thold)), aes(x=snp_model)) +
  stat_summary(aes(y=power, group=snp_model, color=snp_model, fill=snp_model), fun.data = "mean_cl_boot", geom ="bar", fun.args = list(B=1000)) +
  stat_summary(aes(y=power, group=snp_model, color=snp_model), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.power) +
  labs(x = "SNP Model", 
       y = "Power") +
  scale_x_discrete(breaks=snp_names1, labels=snp_names2) +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =text.angel, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="SNP Model"),
         color=guide_legend(title="SNP Model", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=snp.group.colors, breaks = snp_names1, labels = snp_names2) + 
  scale_fill_manual(values=snp.group.colors, breaks = snp_names1, labels = snp_names2)
plot_4 <- ggplot(df_ngwas_power %>% mutate(power = as.numeric(twas.p < thold)), aes(x=ngwas)) +
  stat_summary(aes(y=power, group=ngwas, color=ngwas, fill=ngwas), fun.data = "mean_cl_boot", geom ="bar", fun.args = list(B=1000)) +
  stat_summary(aes(y=power, group=ngwas, color=ngwas), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.power) +
  labs(x = "GWAS Sample Size", y = "Power") +
  scale_x_discrete(breaks=ngwas_names1, labels=ngwas_names2) +
  theme(axis.text = element_text(size = font.size), 
        axis.text.x = element_text(angle = 0, vjust=text.vjust, hjust=text.hjust), 
        legend.direction = "horizontal", 
        legend.title=element_text(size=font.size), 
        legend.text=element_text(size=font.size), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="GWAS Sample Size"),
         color=guide_legend(title="GWAS Sample Size", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=ngwas.group.colors, breaks=ngwas_names1, labels=ngwas_names2) + 
  scale_fill_manual(values=ngwas.group.colors, breaks=ngwas_names1, labels=ngwas_names2)
plot_5 <- ggplot(df_nqtl_power %>% mutate(power = as.numeric(twas.p < thold)), aes(x=nqtl)) +
  stat_summary(aes(y=power, group=nqtl, color=nqtl, fill=nqtl), fun.data = "mean_cl_boot", geom ="bar", fun.args = list(B=1000)) +
  stat_summary(aes(y=power, group=nqtl, color=nqtl), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.power) +
  labs(x = "eQTL Sample Size", y = "Power") +
  scale_x_discrete(breaks=nqtl_names1, labels=nqtl_names2) +
  theme(axis.text = element_text(size = font.size), 
        axis.text.x = element_text(angle = 0, vjust=text.vjust, hjust=text.hjust), 
        legend.direction = "horizontal", 
        legend.title=element_text(size=font.size), 
        legend.text=element_text(size=font.size), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="eQTL Sample Size"),
         color=guide_legend(title="eQTL Sample Size", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=nqtl.group.colors, breaks=nqtl_names1, labels=nqtl_names2) + 
  scale_fill_manual(values=nqtl.group.colors, breaks=nqtl_names1, labels=nqtl_names2) 
plot_6 <- ggplot(df_gwas_sim_power %>% mutate(power = as.numeric(twas.p < thold)), aes(x=gwas.sim)) +
  stat_summary(aes(y=power, group=gwas.sim, color=gwas.sim, fill=gwas.sim), fun.data = "mean_cl_boot", geom ="bar", fun.args = list(B=1000)) +
  stat_summary(aes(y=power, group=gwas.sim, color=gwas.sim), fun.data = "mean_cl_boot", color="gray37", geom="errorbar", fun.args = list(conf.int = 0.95, B=1000), width = 0.2) +
  coord_cartesian(ylim = ylim.power) +
  labs(x = "GWAS Mode", 
       y = "Power") +
  scale_x_discrete(breaks=gwas_mode_names1, labels=gwas_mode_names2) +
  theme(axis.text = element_text(size = font.size),
        legend.direction = "horizontal",
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        axis.text.x = element_text(angle =0, vjust=text.vjust, hjust=text.hjust), 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="GWAS Mode"),
         color=guide_legend(title="GWAS Mode", nrow = 1, byrow = TRUE)) +
  scale_color_manual(values=gwas.mode.group.colors, breaks = gwas_mode_names1, labels = gwas_mode_names2) + 
  scale_fill_manual(values=gwas.mode.group.colors, breaks = gwas_mode_names1, labels = gwas_mode_names2)
prow <- plot_grid(
  plot_1 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_2 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_3 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_4 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_5 + theme(plot.margin = cowplot.margin, legend.position="none"),
  plot_6 + theme(plot.margin = cowplot.margin, legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C", "D", "E", "F"),
  hjust = -1,
  nrow = 2
)
prow
legend2 <- get_legend(plot_2 + theme(legend.position = c(1.685, 2), legend.justification = c("right", "top"), legend.box.just = "right",))
legend3 <- get_legend(plot_3 + theme(legend.position = c(1.95, 2),legend.justification = c("right", "top"), legend.box.just = "right",))
legend4 <- get_legend(plot_4 + theme(legend.position = c(-0.53, 1), legend.justification = c("right", "top"), legend.box.just = "right",))
legend5 <- get_legend(plot_5 + theme(legend.position = c(2.79, 2), legend.justification = c("right", "top"), legend.box.just = "right",))
legend6 <- get_legend(plot_6 + theme(legend.position = c(-.0875, 1), legend.justification = c("right", "top"), legend.box.just = "right",))
combineLegend <- plot_grid(legend2, legend3, legend4, legend5, legend6)
cow_plot <- plot_grid(prow, combineLegend, ncol = 1, rel_heights = c(1, .1))
cow_plot
save_plot(file = "power_strict_cowplot.pdf", plot = cow_plot, base_width = 8.5, base_height=8, units = "in")
# extract power
plot1_power <- ggplot_build(plot_1)$data[[1]]
plot2_power <- ggplot_build(plot_2)$data[[1]]
plot3_power <- ggplot_build(plot_3)$data[[1]]
plot4_power <- ggplot_build(plot_4)$data[[1]]
plot5_power <- ggplot_build(plot_5)$data[[1]]
plot6_power <- ggplot_build(plot_6)$data[[1]]

# Supplementary figure 5: CPU Time
plot_time <- ggplot(df_barchart_cputime, aes(x=ngwas, y=mean.cpu.time, color=ngwas, fill=ngwas)) +
  geom_bar(stat = "identity", position=position_dodge())  + 
  facet_wrap(~gwas.sim, labeller = labeller(gwas.sim = as_labeller(gwas_names))) +
  labs(x = "GWAS Sample Size", 
       y = "Mean CPU Time (in seconds)") +
  scale_x_discrete(breaks=ngwas_names1, labels=ngwas_names2) +
  theme(axis.text = element_text(size = font.size),
        legend.title=element_text(size=font.size),
        legend.text=element_text(size=font.size), 
        plot.margin = plot.margin, 
        panel.background  = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  guides(fill=guide_legend(title="GWAS Sample Size"),
         color=guide_legend(title="GWAS Sample Size")) +
  scale_color_manual(values=ngwas.group.colors, breaks=ngwas_names1, labels=ngwas_names2) + 
  scale_fill_manual(values=ngwas.group.colors, breaks=ngwas_names1, labels=ngwas_names2)
ggsave(file = "CPU_time.pdf", plot = plot_time, width = 8.5, height=5.0, units = "in")