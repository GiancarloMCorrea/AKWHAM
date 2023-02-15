source('aux_fun.R')

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


# -------------------------------------------------------------------------

# Plot selex devs from SS and WHAM:
require(r4ss)

# EBS pcod:
ssmod = r4ss::SS_output('SS_models/EBS_pcod')
load('EBS_pcod/fit_a.RData')

# Get devs SS model
par3_f1 = ssmod$parameters$Value[185:230]
par6_f1 = ssmod$parameters$Value[231:276]
par1_f2 = ssmod$parameters$Value[277:317]
par3_f2 = ssmod$parameters$Value[318:358]

par(mfrow = c(2,2))
hist(par3_f1, main = 'par3_fishery')
hist(par6_f1, main = 'par6_fishery')
hist(par1_f2, main = 'par1_survey')
hist(par3_f2, main = 'par3_survey')

# Get devs WHAM model
par3_f1 = fit_a$input$par$selpars_re[1:46]
par6_f1 = fit_a$input$par$selpars_re[47:92]
par1_f2 = fit_a$input$par$selpars_re[98:138]
par3_f2 = fit_a$input$par$selpars_re[144:184]

plot1 <- hist(par3_f1, plot = FALSE)
plot2 <- hist(par6_f1, plot = FALSE)
plot3 <- hist(par1_f2, plot = FALSE)
plot4 <- hist(par3_f2, plot = FALSE)

par(mfrow = c(1,1))
plot(plot1, col = 1, xlim = c(-0.05, 0.05))
plot(plot2, col = 2, add = TRUE)

par(mfrow = c(1,1))
plot(plot3, col = 3, xlim = c(-1, 1))
plot(plot4, col = 4, add = TRUE)
