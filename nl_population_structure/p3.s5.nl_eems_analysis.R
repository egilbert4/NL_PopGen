# script to plot the EEMS results

library("tidyverse")
library("rEEMSplots")
library("rworldxtra")
library("rgdal")

data(countriesHigh)
countriesHigh <- fortify(countriesHigh)

coord_data <- "NL-ancestry.geog-filtered.coord"
(coord_data <- read.table(coord_data))

ggplot() +
  geom_polygon(data = countriesHigh, aes(x = long, y = lat, group = group)) +
  geom_point(data = coord_data, aes(x = V1, y = V2), colour = "red", size = 1) +
  coord_cartesian(xlim = c(-60, -52.5), ylim = c(46,52.5))

projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"

coord_merc <- sp::spTransform(SpatialPoints(as.matrix(coord_data),
                                            CRS(projection_none)), CRS(projection_mercator))
coord_merc <- coord_merc@coords

chains <- 1:10
ndemes <- 300
analysis <- "initial" # or "extension"

eems.plots(mcmcpath = paste0(analysis, "/", "NFL-ancestry.EEMS-", analysis, "-nDemes", ndemes, ".chain", chains), 
           plotpath = "eems-example",
           add.abline = TRUE, add.r.squared = TRUE, 
           add.grid = TRUE, lwd.grid = 0.5, col.grid = "grey60",
           add.demes = TRUE, col.demes = "black", min.cex.demes = 0.5, max.cex.demes = 1, 
           add.map = TRUE, lwd.map = 1, col.map = "black",
           projection.in = projection_none, projection.out = projection_mercator,
           longlat = TRUE, 
           # only use this if you want to see the original points
           #m.plot.xy = points(coord_merc, col = "red", pch = 19, cex = 0.25 ),
           out.png = TRUE, res = 600, plot.width = 7, plot.height = 8)


