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
DAT1 [label = 'Informative data: \n Marginal length comps \n Conditional age-at-length', fillcolor = White, style = rounded]
DAT2 [label = 'Informative data: \n Observed mean \n weight-at-age', fillcolor = White, style = rounded]

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
