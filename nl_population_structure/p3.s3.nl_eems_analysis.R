# create the coord file from the filtered plink data

# library
library("tidyverse")

# the fam data, i.e., the genotyped individuals
fam_data <- "NL-ancestry.geog-filtered.fam"
fam_data <- read_delim(file = fam_data, col_names = c("FID","IID","PID","MID","SEX","PHE"), delim = " ")

# the coord data, the per-grandparent lat/lng data
coord_data <- "NL1807.gp-birthplace_and_clusters_rescued.csv"
coord_data <- read_csv(file = coord_data)

# summarise the lat/lng data to individual
coord_data <- coord_data %>%
  group_by(FSID,IID,KCLST) %>%
  summarise(lat = mean(lat, na.rm = T), lng = mean(lng, na.rm = T),
            .groups = "drop")

# combine the coord data with fam data
coord_data <- left_join(select(fam_data, FID,IID), coord_data, by = "IID")

outfile <- "NL-ancestry.geog-filtered.EEMS.data"
write_delim(file = outfile, x = coord_data, delim = " ", col_names = F)

outfile <- "NL-ancestry.geog-filtered.coord"
write_delim(file = outfile, x = select(coord_data,lng,lat), delim = " ", col_names = F)
