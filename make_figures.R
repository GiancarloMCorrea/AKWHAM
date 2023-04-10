source('aux_fun.R')
theme_set(theme_bw())
require(ggplot2)

# Make diagram ------------------------------------------------------------

require(DiagrammeR)
require(htmltools)

diag1 = DiagrammeR::grViz("digraph {

graph [layout = dot, rankdir = LB]

# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = Linen, fontsize=25]

LAA1 [label = 'LAA: \n Parametric', fillcolor = Pink]
LAA2 [label = 'LAA: \n Non-parametric', fillcolor = Pink]
WAA1 [label =  'WAA: \n Empirical WAA', fillcolor = lightskyblue2]
WAA2 [label =  'WAA: \n Through L-W relationship', fillcolor = lightskyblue2]
WAA3 [label =  'WAA: \n Non-parametric', fillcolor = lightskyblue2]
SSB [label= 'SSB and reference points \n calculation', fillcolor = Beige]
DAT1 [label = 'Marginal length comps \n Conditional age-at-length', fillcolor = White, style = rounded]
DAT2 [label = 'Observed mean \n weight-at-age', fillcolor = White, style = rounded]

{rank = min; DAT1 DAT2}
{rank = same; LAA1 LAA2}
{rank = same; WAA1 WAA2 WAA3}
{rank = max; SSB}

# edge definitions with the node IDs
DAT1 -> {LAA1 LAA2}
DAT2 -> {WAA2 WAA3}
edge [arrowhead = none]
{WAA1 WAA2 WAA3}  -> SSB
{LAA1 LAA2}  -> WAA2

}")


# 2. Convert to SVG, then save as png
tmp = DiagrammeRsvg::export_svg(diag1)
html_print(HTML(tmp))
tmp = charToRaw(tmp) # flatten
rsvg::rsvg_svg(tmp, "diagram2.svg", width = 1700) # saved graph as png in current working directory

bitmap <- rsvg::rsvg("diagram2.svg")
dim(bitmap)

jpeg::writeJPEG(bitmap, "diagram2.jpg", quality = 1)


# Make data plots ---------------------------------------------------------

require(r4ss)
# GOA pollock:
mydat = readRDS('aux_data/datfile.RDS')

ebscod = r4ss::SS_output('SS_models/EBS_pcod')
goacod = r4ss::SS_output('SS_models/GOA_pcod')

png(filename = 'Data_cases.png', width = 190, height = 160, units = 'mm', res = 500)
par(mfrow = c(2,2))
plot_data_overview(datlist = mydat, sectionCex = 0.9)
r4ss::SSplotData(replist = goacod, subplots = 2, margins = c(1.7,1,2,4.5))
title(main = 'GOA Pacific cod')
r4ss::SSplotData(replist = ebscod, subplots = 2, margins = c(1.7,1,2,4.5))
title(main = 'EBS Pacific cod')
dev.off()


# GOA pollock plots -------------------------------------------------------
model_names = c('ADMB', 'wham_ewaa', 'wham_iid', 'wham_2dar1')

# ADMB model:
asdrep = readRDS('aux_data/adsdrep_selex_fixed.RDS') # ADMB pollock model output
SSBdata0 = cbind(model = model_names[1],
                 filter(asdrep, name=='log_ssb') %>% select(year, est, sd=se))
SSBdata0$est = log(exp(SSBdata0$est)*1e+06) # in tons
# Model a:
load('GOA_pollock/fit_a.RData')
this_model = fit_a
model_name = model_names[2]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata1 = cbind(model = model_name, year=1970:2021,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
# Model b:
load('GOA_pollock/fit_b.RData')
this_model = fit_b
model_name = model_names[3]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata2 = cbind(model = model_name, year=1970:2021,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
# Model c:
load('GOA_pollock/fit_c.RData')
this_model = fit_c
model_name = model_names[4]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata3 = cbind(model = model_name, year=1970:2021,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))

# Compare SSB among models:
plot_data = rbind(SSBdata0, SSBdata1, SSBdata2, SSBdata3)
plot_data$model = factor(plot_data$model, levels = model_names)
plot_data$ssb_min = exp(plot_data$est-1.96*plot_data$sd)
plot_data$ssb_max = exp(plot_data$est+1.96*plot_data$sd)
plot_data$est = exp(plot_data$est)

p1 = ggplot(plot_data, aes(year, est*1e-06, ymin=ssb_min*1e-06, ymax=ssb_max*1e-06,
                           fill=model, color=model)) +
  labs(y='SSB (million tons)', x = NULL) +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  coord_cartesian(xlim = c(1970,2021), ylim = c(0,1.2)) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  annotate("text", label = 'A', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 2), 
         color=guide_legend(title=NULL,nrow = 2),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.9), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())

# Plot projections:
# Model a:
load('GOA_pollock/proj_a.RData')
this_model = proj_a
model_name = model_names[2]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBproj1 = cbind(model = model_name, year=1970:2024,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
# Model b:
load('GOA_pollock/proj_b.RData')
this_model = proj_b
model_name = model_names[3]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBproj2 = cbind(model = model_name, year=1970:2024,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
# Model c:
load('GOA_pollock/proj_c.RData')
this_model = proj_c
model_name = model_names[4]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBproj3 = cbind(model = model_name, year=1970:2024,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))

# Merge data:
plot_data = rbind(SSBproj1, SSBproj2, SSBproj3)
plot_data$model = factor(plot_data$model, levels = model_names)
plot_data$ssb_min = exp(plot_data$est-1.96*plot_data$sd)
plot_data$ssb_max = exp(plot_data$est+1.96*plot_data$sd)
plot_data$est = exp(plot_data$est)

thisCols = RColorBrewer::brewer.pal(n = 9, name = 'Set1')[2:4]
p2 = ggplot(plot_data, aes(year, est*1e-06, ymin=ssb_min*1e-06, ymax=ssb_max*1e-06,
                           fill=model, color=model)) +
  labs(y='SSB (million tons)', x = NULL) +
  coord_cartesian(xlim = c(2020,2024), ylim = c(0,0.7)) +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  geom_vline(xintercept = 2021, linetype = 'dashed') +
  scale_fill_manual(values = thisCols) +
  scale_color_manual(values = thisCols) +
  theme_bw() +
  annotate("text", label = 'B', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 2), 
         color=guide_legend(title=NULL,nrow = 2),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.9), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())

# Plot random effects iid and 2DAR1:
# iid model:
WAA_data = fit_b$rep$WAA_re
colnames(WAA_data) = fit_b$input$ages.lab
rownames(WAA_data) = fit_b$input$years
WAA_plot = reshape2::melt(data = WAA_data, varnames = c('year', 'age'))
p3 = ggplot(WAA_plot, aes(x = year, y = value, color=age)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  labs(y='WAA random effects (wham_iid)', x= NULL) +
  labs( color=NULL) +
  scale_color_brewer(palette = 'Spectral') +
  ylim(-0.75, 1) +
  theme_bw() +
  annotate("text", label = 'C', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 2), 
         color=guide_legend(title=NULL,nrow = 2),
         shape = guide_legend(override.aes = list(size = 0.1))) +
  theme(legend.position=c(0.5,0.9), 
        legend.text = element_text(size = 7.5),
        legend.background = element_blank()) 
# 2dAR1 model:
WAA_data = fit_c$rep$WAA_re
colnames(WAA_data) = fit_c$input$ages.lab
rownames(WAA_data) = fit_c$input$years
WAA_plot = reshape2::melt(data = WAA_data, varnames = c('year', 'age'))
p4 = ggplot(WAA_plot, aes(x = year, y = value, color=age)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  labs(y='WAA random effects (wham_2dar1)', x= NULL) +
  labs( color=NULL) +
  scale_color_brewer(palette = 'Spectral') +
  ylim(-0.75, 1) +
  theme_bw() +
  annotate("text", label = 'D', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 2), 
         color=guide_legend(title=NULL,nrow = 2),
         shape = guide_legend(override.aes = list(size = 0.1))) +
  theme(legend.position=c(0.5,0.9), 
        legend.text = element_text(size = 7.5),
        legend.background = element_blank())

# Merge plots and save:
png(filename = 'GOA_pollock/main_GOApollock.png', width = 190, height = 160, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

# Plot WAA fits:
plot_waa_fit(fit = fit_b, minyr=1990, maxyr=2009)
ggsave(filename = 'GOA_pollock/summary_WAA_fit_b.png', width = 190, height = 140, units = 'mm', dpi = 500)
plot_waa_fit(fit = fit_c, minyr=1990, maxyr=2009)
ggsave(filename = 'GOA_pollock/summary_WAA_fit_c.png', width = 190, height = 140, units = 'mm', dpi = 500)
plot_waa_fit(fit = fit_b, minyr=1977, maxyr=2021, by.cohort = FALSE)
ggsave(filename = 'GOA_pollock/summary_WAA_year_fit_b.png', width = 190, height = 140, units = 'mm', dpi = 500)
plot_waa_fit(fit = fit_c, minyr=1977, maxyr=2021, by.cohort = FALSE)
ggsave(filename = 'GOA_pollock/summary_WAA_year_fit_c.png', width = 190, height = 140, units = 'mm', dpi = 500)

# Plot WAA proj:
plot_waa_proj(mods = list(proj_b, proj_c), minyr=2020, maxyr=2024, myCols = thisCols[2:3], 
              modNames = model_names[3:4], projYear = 2021)
ggsave(filename = 'GOA_pollock/summary_WAA_proj.png', width = 190, height = 140, units = 'mm', dpi = 500)


# GOA pcod plots -------------------------------------------------------

n_ages = 10
# Model names in plot:
model_names = c('SS', 'wham')

SS_report = r4ss::SS_output(dir = 'SS_models/GOA_pcod') # from OM
# WHAM model:
load('GOA_pcod/fit_a.RData')

# Get SSB from SS:
SS_SSB = SS_report$derived_quants[grep(pattern = 'SSB_', x = SS_report$derived_quants$Label),]
data1 = cbind(data.frame(model = model_names[1], year = fit_a$years), 
              SS_SSB[3:48, c('Value', 'StdDev')])
colnames(data1)[3:4] = c('est', 'sd')
data1$ssb_min = data1$est - 1.96*data1$sd
data1$ssb_max = data1$est + 1.96*data1$sd

this_model = fit_a
model_name = model_names[2]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
data2 = cbind(model = model_name, year=1977:2022,
              filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
data2$ssb_min = exp(data2$est-1.96*data2$sd)
data2$ssb_max = exp(data2$est+1.96*data2$sd)
data2$est = exp(data2$est)

# Make plot:
plot_data =rbind(data1, data2)
plot_data$model = factor(plot_data$model, levels = model_names)

p1 = ggplot(plot_data, aes(year, est*1e-06, ymin=ssb_min*1e-06, ymax=ssb_max*1e-06,
                           fill=model, color=model)) +
  ylim(0,NA) + labs(y='SSB (million tons)', x = NULL) +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  annotate("text", label = 'A', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 1), 
         color=guide_legend(title=NULL,nrow = 1),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.95), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())

# Get LAA from SS:
data1 = data.frame(ages = SS_report$endgrowth$int_Age[2:11], 
                   LAA = SS_report$endgrowth$Len_Beg[2:11], 
                   SD = SS_report$endgrowth$SD_Beg[2:11], 
                   model = model_names[1])
data1$len_min = data1$LAA - 1.96*data1$SD
data1$len_max = data1$LAA + 1.96*data1$SD

# Get LAA from WHAM:
data2 = data.frame(ages = 1:n_ages, LAA = fit_a$rep$LAA[1,], 
                   SD = fit_a$rep$SDAA[1,], model = model_names[2])
data2$len_min = data2$LAA - 1.96*data2$SD
data2$len_max = data2$LAA + 1.96*data2$SD

# Make plot:
plot_data =rbind(data1, data2)
plot_data$model = factor(plot_data$model, levels = model_names)

p2 = ggplot(plot_data, aes(ages, LAA, ymin=len_min, ymax=len_max,
                           fill=model, color=model)) +
  ylim(0,NA) + labs(y='Mean length (cm)', x = NULL) +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  annotate("text", label = 'B', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  guides(fill=guide_legend(title=NULL,nrow = 1), 
         color=guide_legend(title=NULL,nrow = 1),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.95), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())

# Merge plots:
png(filename = 'GOA_pcod/main_GOApcod.png', width = 190, height = 70, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()


# EBS pcod plots ----------------------------------------------------------

n_ages = 20
years = 1977:2022
SS_report = r4ss::SS_output(dir = 'SS_models/EBS_pcod') # from OM

# Model names in plot:
model_names = c('SS', 'wham_ar1', 'wham_ecov')

# Get SSB estimates:
load('EBS_pcod/fit_a.RData')
this_model = fit_a
model_name = model_names[2]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata1 = cbind(model = model_name, year=years,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
SSBdata1$ssb_min = exp(SSBdata1$est-1.96*SSBdata1$sd)
SSBdata1$ssb_max = exp(SSBdata1$est+1.96*SSBdata1$sd)
SSBdata1$est = exp(SSBdata1$est)

# Get SSB estimates:
load('EBS_pcod/fit_b.RData')
this_model = fit_b
model_name = model_names[3]
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
SSBdata2 = cbind(model = model_name, year=years,
                 filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
SSBdata2$ssb_min = exp(SSBdata2$est-1.96*SSBdata2$sd)
SSBdata2$ssb_max = exp(SSBdata2$est+1.96*SSBdata2$sd)
SSBdata2$est = exp(SSBdata2$est)

# Get SSB from SS:
SS_SSB = SS_report$derived_quants[grep(pattern = 'SSB_', x = SS_report$derived_quants$Label),]
SSBdata0 = cbind(data.frame(model = model_names[1], year = years), 
                 SS_SSB[3:48, c('Value', 'StdDev')])
colnames(SSBdata0)[3:4] = c('est', 'sd')
SSBdata0$ssb_min = SSBdata0$est - 1.96*SSBdata0$sd
SSBdata0$ssb_max = SSBdata0$est + 1.96*SSBdata0$sd

# Merge data:
plot_data = rbind(SSBdata0, SSBdata1, SSBdata2)
plot_data$est = plot_data$est*1E-06
plot_data$ssb_min = plot_data$ssb_min*1E-06
plot_data$ssb_max = plot_data$ssb_max*1E-06

plot_data$model = factor(plot_data$model, levels = model_names)
# Make plot mean SSB:
p1 = ggplot(plot_data, aes(year, est, ymin=ssb_min, ymax=ssb_max,
                           fill=model, color=model)) +
  ylim(0,NA) + labs(y='SSB (million tons)', x = NULL) +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  labs( color=NULL, fill=NULL) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  annotate("text", label = 'A', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 1), 
         color=guide_legend(title=NULL,nrow = 1),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.075), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())

# Plot L1 temporal variability:
LAA_a = data.frame(year = years, L1 = fit_a$rep$LAA[,1], model = model_names[2])
LAA_b = data.frame(year = years, L1 = fit_b$rep$LAA[,1], model = model_names[3])
LAA_SS = data.frame(year = years[-length(years)], L1 = SS_report$growthseries[,6], model = model_names[1])

# Merge datasets:
plot_data = rbind(LAA_a, LAA_b, LAA_SS)
plot_data$model = factor(plot_data$model, levels = model_names)

# Make plot:
p2 = ggplot(plot_data, aes(year, L1, color=model)) +
  geom_line(lwd=1) +
  labs(y='Mean length-at-age 1 (cm)', x = NULL) +
  labs( color=NULL, fill=NULL) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  coord_cartesian(ylim=c(6, 18)) +
  annotate("text", label = 'B', x = -Inf, y = Inf, hjust = -1, vjust = 1.5) +
  guides(fill=guide_legend(title=NULL,nrow = 1), 
         color=guide_legend(title=NULL,nrow = 1),
         shape = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme(legend.position=c(0.5,0.075), 
        legend.text = element_text(size = 10),
        legend.background = element_blank())
# ggsave(filename = 'EBS_pcod/compare_L1.png', width = 190, height = 120, units = 'mm', dpi = 500)

# Plot Ecov:
p3 = plot_ecov_fit(fit_a, label = 'C', myCol = "#377EB8")
p4 = plot_ecov_fit(fit_b, label = 'D', myCol = "#4DAF4A")

# Merge plots:
png(filename = 'EBS_pcod/main_EBSpcod.png', width = 190, height = 160, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

