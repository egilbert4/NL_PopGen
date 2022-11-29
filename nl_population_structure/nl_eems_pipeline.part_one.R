# select individuals for the NL-EEMS analysis

# library
library("tidyverse")

# the fam data, i.e., the genotyped individuals
fam_data <- "NL-ancestry.filtered-cleaned.fam"
fam_data <- read_delim(file = fam_data, col_names = c("FID","IID","PID","MID","SEX","PHE"), delim = " ")

# the coord data, the per-grandparent lat/lng data
coord_data <- "NL1807.gp-birthplace_and_clusters_rescued.csv"
coord_data <- read_csv(file = coord_data)

# summarise the lat/lng data to individual
coord_data <- coord_data %>%
  group_by(FSID,IID,KCLST) %>%
  summarise(lat = mean(lat, na.rm = T), lng = mean(lng, na.rm = T),
            .groups = "drop")

# filter fam IIDs for coord data
filter_fam_data <- fam_data %>%
  filter(IID %in% coord_data$IID) %>%
select(FID,IID)

outfile <- "NL-ancestry.geog-filtered.idfile"
write_delim(file = outfile, x = filter_fam_data, delim = " ", col_names = F)
