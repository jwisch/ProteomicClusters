library(circlize)
library('Cairo')
CairoWin() #dealing with aliasing
#BBB - blood brain barrier, N - neurodegeneration, I - inflammation, A - amyloid production, L - liver function, CVD - cardiovascular disease

df <- read.csv("./Results/Proteins_forChordDiagram.csv") #made this by hand based on results

color_map <- data.frame("Type" = c("A", "BBB", "CVD", "I", "L", "N"),
                        "Type_Full" = c("Amyloid", "BBB", "CVD", "Inflammation", "Liver", "Neurodegeneration"),
                        "Hex" = c("#FDE725FF", "#7AD151FF", "#22A884FF",
                                  "#2A788EFF", "#414487FF", "#440154FF"))

df <- merge(df, color_map, by = "Type")

# grid.col = c("Emerging vs No" = "#440154FF",
#              "AD vs Emerging" = "#238A8DFF", 
#              "AD vs No" = "#FDE725FF",
#              "BBB" = "grey13",
#              "N" = "grey25",
#              "CVD" = "grey37",
#              "I" = "grey49",
#              "A" = "grey61",
#              "L" = "grey73")
df$Classification <- as.factor(df$Classification)
levels(df$Classification) <- c("Positive vs Intermediate", "Positive vs Negative", "Intermediate vs Negative")


grid.col = structure( c(df$Hex, "#440154FF","#238A8DFF","#FDE725FF"), 
                     names = c(df$Protein, "Intermediate vs Negative", "Positive vs Intermediate","Positive vs Negative" ) )

group = structure(c(df$Protein, df$Protein), names = c(df$Type, df$Classification))

chordDiagram(df[, c("Protein", "Classification")], annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
#  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)


circos.clear()

legend(x = "bottomleft", 
       legend = c(color_map$Type_Full[c(2:3, 5, 1, 4, 6)]),  # Legend texts
       lty = rep(1, length(color_map$Type_Full[c(2:3, 5, 1, 4, 6)])),           # Line types
       col = c(color_map$Hex[c(2:3, 5, 1, 4, 6)]),           # Line colors
       lwd = 2, 
       bty = "n") 
