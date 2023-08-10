source('aux_fun.R')
require(ggplot2)
theme_set(theme_bw())

# Make diagram ------------------------------------------------------------
require(glue)
require(DiagrammeR)
require(htmltools)
require(DiagrammeRsvg)
require(svglite)
require(rsvg)
require(tiff)

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


# Save:
DPI = 500
WidthCM = 17
HeightCM = 8

diag1 %>% export_svg %>% charToRaw %>% 
  rsvg(width = WidthCM *(DPI/2.54), height = HeightCM *(DPI/2.54)) %>% 
  jpeg::writeJPEG("Figure-1.jpg", quality = 1)

# Now you have to modify the DPI using GIMP. Load the jpg file just created and go to
# Image > Scale Image, and change resolution (px/in) to 500

# Make data plots ---------------------------------------------------------

require(r4ss)
# GOA pollock:
mydat = readRDS('aux_data/datfile.RDS')

ebscod = r4ss::SS_output('SS_models/EBS_pcod', covar = FALSE)
goacod = r4ss::SS_output('SS_models/GOA_pcod', covar = FALSE)

jpeg(filename = 'Figure-2.jpg', width = 170, height = 150, units = 'mm', res = 500)
par(mfrow = c(2,2))
plot_data_overview(datlist = mydat, sectionCex = 0.9)
r4ss::SSplotData(replist = goacod, subplot = 2, margins = c(1.7,1,2,4.5))
title(main = 'GOA Pacific cod')
r4ss::SSplotData(replist = ebscod, subplot = 2, margins = c(1.7,1,2,4.5))
title(main = 'EBS Pacific cod')
dev.off()

# -------------------------------------------------------------------------
# GOA pollock plots -------------------------------------------------------
# -------------------------------------------------------------------------

model_names = c('Original', 'wham_ewaa', 'wham_iid', 'wham_2dar1')

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
jpeg(filename = 'GOA_pollock/Figure-3.jpg', width = 170, height = 150, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

# Plot WAA fits:
plot_waa_fit(fit = fit_b, minyr=1990, maxyr=2009)
ggsave(filename = 'GOA_pollock/summary_WAA_fit_b.jpg', width = 170, height = 140, units = 'mm', dpi = 500)
plot_waa_fit(fit = fit_c, minyr=1990, maxyr=2009)
ggsave(filename = 'GOA_pollock/summary_WAA_fit_c.jpg', width = 170, height = 140, units = 'mm', dpi = 500)
plot_waa_fit(fit = fit_b, minyr=1977, maxyr=2021, by.cohort = FALSE)
ggsave(filename = 'GOA_pollock/summary_WAA_year_fit_b.jpg', width = 170, height = 140, units = 'mm', dpi = 500)
plot_waa_fit(fit = fit_c, minyr=1977, maxyr=2021, by.cohort = FALSE)
ggsave(filename = 'GOA_pollock/summary_WAA_year_fit_c.jpg', width = 170, height = 140, units = 'mm', dpi = 500)

# Plot WAA proj:
plot_waa_proj(mods = list(proj_a, proj_b, proj_c), minyr=2020, maxyr=2024, myCols = thisCols[1:3], 
              modNames = model_names[2:4], projYear = 2021)
ggsave(filename = 'GOA_pollock/summary_WAA_proj.jpg', width = 170, height = 140, units = 'mm', dpi = 500)

# Make selectivity plot (compare ADMB and WHAM)

all_years = 1970:2021
all_ages = 1:10
all_fleets = 1:7
fleets_names = c('Fishery', 'Age 3+ Shelikof', 'NMFS BT', 'ADFG BT', 
                 'Age 1 Shelikof', 'Age 2 Shelikof', 'Summer AT')

# Read ADMB selex:
load('aux_data/admb_selex.RData')
admb_selex$type = model_names[1]

# Organize WHAM output
k = 1
tmp_df = list()
for(j in seq_along(all_fleets)) {
  for(i in seq_along(all_years)) {
    tmp_df[[k]] = data.frame(year = all_years[i], age = all_ages, 
                             selex = fit_a$rep$selAA[[j]][i,], fleet = fleets_names[j], 
                             type = model_names[2])
    k = k + 1
  }
}
for(j in seq_along(all_fleets)) {
  for(i in seq_along(all_years)) {
    tmp_df[[k]] = data.frame(year = all_years[i], age = all_ages, 
                             selex = fit_b$rep$selAA[[j]][i,], fleet = fleets_names[j], 
                             type = model_names[3])
    k = k + 1
  }
}
for(j in seq_along(all_fleets)) {
  for(i in seq_along(all_years)) {
    tmp_df[[k]] = data.frame(year = all_years[i], age = all_ages, 
                             selex = fit_c$rep$selAA[[j]][i,], fleet = fleets_names[j], 
                             type = model_names[4])
    k = k + 1
  }
}

# Merge both datasets:
plot_data = dplyr::bind_rows(tmp_df)
plot_data = rbind(admb_selex, plot_data)
plot_data$fleet = factor(plot_data$fleet, levels = fleets_names)
plot_data$age = factor(plot_data$age)

# Make plot:
ggplot(plot_data, aes(x = year, y = age)) + 
  geom_raster(aes(fill=selex)) + 
  scale_fill_viridis_c() +
  labs(x="Year", y="Age", fill = 'Selectivity') +
  scale_x_continuous(breaks = c(1970, 1980, 1990, 2000, 2010, 2020), 
                     labels = c(1970, '', 1990, '', 2010, '')) +
  theme_bw() + theme(axis.text.x=element_text(size=7.5, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11),
                     legend.position = 'top') +
  facet_grid(type ~ fleet)
ggsave(filename = 'GOA_pollock/selex_GOApollock.jpg', width = 170, height = 150, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# GOA pcod plots -------------------------------------------------------
# -------------------------------------------------------------------------

n_ages = 10
# Model names in plot:
model_names = c('SS3', 'wham')

SS_report = r4ss::SS_output(dir = 'SS_models/GOA_pcod', covar = FALSE) # from OM
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
jpeg(filename = 'GOA_pcod/Figure-5.jpg', width = 170, height = 70, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

# Plot Ecov:
jpeg(filename = 'GOA_pcod/env_index_Q.jpg', width = 85, height = 60, units = 'mm', res = 500)
plot_ecov_fit(fit_a, label = ' ', myCol = "#377EB8", yLab = 'Environmental index')
dev.off()

# Make selectivity plot (compare SS and WHAM)

all_years = 1977:2022
all_fleets = 1:5
tmp_df = list()

# Organize SS3 output:
SelMat = SS_report$sizeselex[SS_report$sizeselex$Factor == 'Lsel' & SS_report$sizeselex$Yr >= 1977 & SS_report$sizeselex$Yr <= 2022 & SS_report$sizeselex$Fleet <= 5, ]
k = 1
for(j in seq_along(all_fleets)) {
  ind_vec = numeric(length(all_years))
  for(i in seq_along(all_years)) {
    ind_yr = which(SelMat$Fleet == all_fleets[j] & SelMat$Yr == all_years[i])
    if(length(ind_yr) > 0) ind_vec[i] = ind_yr
    else ind_vec[i] = ind_vec[i-1]
    
    temp = SelMat[ind_vec[i], ]
    tmp_df[[k]] = data.frame(fleet = j, year = all_years[i], len = as.numeric(colnames(temp)[6:ncol(temp)]), 
                             selex = as.vector(as.matrix(temp[,6:ncol(temp)])), type = model_names[1])
    k = k + 1
  }
}

# Organize WHAM output
for(j in seq_along(all_fleets)) {
  for(i in seq_along(all_years)) {
    tmp_df[[k]] = data.frame(fleet = j, year = all_years[i], len = fit_a$input$data$lengths + 0.5, 
                             selex = fit_a$rep$selAL[[j]][i,], type = model_names[2])
    k = k + 1
  }
}

# Merge both datasets:
plot_data = dplyr::bind_rows(tmp_df)
plot_data$fleet = factor(plot_data$fleet, levels = 1:5, labels = c('FshTrawl', 'FshLL', 'FshPot', 'Srv', 'LLSrv'))

# Make plot:
ggplot(plot_data, aes(x = year, y = len)) + 
  geom_raster(aes(fill=selex)) + 
  scale_fill_viridis_c() +
  labs(x="Year", y="Length (cm)", fill = 'Selectivity') +
  theme_bw() + theme(axis.text.x=element_text(size=7.5, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11),
                     legend.position = 'top') +
  facet_grid(type ~ fleet)
ggsave(filename = 'GOA_pcod/selex_GOApcod.jpg', width = 170, height = 120, units = 'mm', dpi = 500)

# -------------------------------------------------------------------------
# EBS pcod plots ----------------------------------------------------------
# -------------------------------------------------------------------------

n_ages = 20
years = 1977:2022
SS_report = r4ss::SS_output(dir = 'SS_models/EBS_pcod', covar = FALSE) # from OM

# Model names in plot:
model_names = c('SS3', 'wham_ar1', 'wham_ecov')

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
jpeg(filename = 'EBS_pcod/Figure-7.jpg', width = 170, height = 150, units = 'mm', res = 500)
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

# Compare models with difference logsime for Ecov:
load('EBS_pcod/fit_a.RData')
fit_a_02 = fit_a
load('EBS_pcod/fit_a01.RData')
fit_a_01 = fit_a
load('EBS_pcod/fit_a001.RData')
fit_a_001 = fit_a
load('EBS_pcod/fit_b.RData')
fit_b_02 = fit_b
load('EBS_pcod/fit_b01.RData')
fit_b_01 = fit_b
load('EBS_pcod/fit_b001.RData')
fit_b_001 = fit_b

a1 = plot_ecov_fit(fit_a_02, label = 'obs error = 0.2', myCol = "#377EB8")
a2 = plot_ecov_fit(fit_a_01, label = 'obs error = 0.1', myCol = "#377EB8")
a3 = plot_ecov_fit(fit_a_001, label = 'obs error = 0.01', myCol = "#377EB8")
b1 = plot_ecov_fit(fit_b_02, label = 'obs error = 0.2', myCol = "#4DAF4A")
b2 = plot_ecov_fit(fit_b_01, label = 'obs error = 0.1', myCol = "#4DAF4A")
b3 = plot_ecov_fit(fit_b_001, label = 'obs error = 0.01', myCol = "#4DAF4A")

# Merge plots:
jpeg(filename = 'EBS_pcod/ecov_sigma.jpg', width = 170, height = 210, units = 'mm', res = 500)
gridExtra::grid.arrange(a1, b1, a2, b2, a3, b3, ncol = 2)
dev.off()

# Get AIC table:
all_mods = wham::compare_wham_models(mods = list(fit_a_02, fit_a_01, fit_a_001,
                                                 fit_b_02, fit_b_01, fit_b_001), 
                                     table.opts = list(calc.rho = F), do.plot = FALSE)

# Make selectivity plot (compare SS and WHAM)
all_years = 1977:2022
all_fleets = 1:2
tmp_df = list()

# Organize SS3 output:
SelMat = SS_report$sizeselex[SS_report$sizeselex$Factor == 'Lsel' & SS_report$sizeselex$Yr >= 1977 & SS_report$sizeselex$Yr <= 2022 & SS_report$sizeselex$Fleet <= 5, ]
k = 1
for(j in seq_along(all_fleets)) {
  ind_vec = numeric(length(all_years))
  for(i in seq_along(all_years)) {
    ind_yr = which(SelMat$Fleet == all_fleets[j] & SelMat$Yr == all_years[i])
    if(length(ind_yr) > 0) ind_vec[i] = ind_yr
    else ind_vec[i] = ind_vec[i-1]
    
    temp = SelMat[ind_vec[i], ]
    tmp_df[[k]] = data.frame(fleet = j, year = all_years[i], len = as.numeric(colnames(temp)[10:ncol(temp)]), 
                             selex = as.vector(as.matrix(temp[,10:ncol(temp)])), type = model_names[1])
    k = k + 1
  }
}

# Organize WHAM output
for(j in seq_along(all_fleets)) {
  for(i in seq_along(all_years)) {
    tmp_df[[k]] = data.frame(fleet = j, year = all_years[i], len = fit_a$input$data$lengths + 0.5, 
                             selex = fit_a$rep$selAL[[j]][i,], type = model_names[2])
    k = k + 1
  }
}
for(j in seq_along(all_fleets)) {
  for(i in seq_along(all_years)) {
    tmp_df[[k]] = data.frame(fleet = j, year = all_years[i], len = fit_b$input$data$lengths + 0.5, 
                             selex = fit_b$rep$selAL[[j]][i,], type = model_names[3])
    k = k + 1
  }
}

# Merge both datasets:
plot_data = dplyr::bind_rows(tmp_df)
plot_data$fleet = factor(plot_data$fleet, levels = 1:2, labels = c('Fishery', 'Survey'))

# Make plot:
ggplot(plot_data, aes(x = year, y = len)) + 
  geom_raster(aes(fill=selex)) + 
  scale_fill_viridis_c() +
  labs(x="Year", y="Length (cm)", fill = 'Selectivity') +
  theme_bw() + theme(axis.text.x=element_text(size=11, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=11),
                     plot.title=element_text(size=11),
                     legend.position = 'top') +
  facet_grid(type ~ fleet)
ggsave(filename = 'EBS_pcod/selex_EBSpcod.jpg', width = 170, height = 210, units = 'mm', dpi = 500)
