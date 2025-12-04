############################################################################### #
# Aim ----
#| Here is an example for Belgium data
# Content:
#| Load package
#| Load data
#| Clean data
#| Compute viral concentrations (cc)
#| Normalization
#| Smoothening
#| Aggregation at national level
#| Level of activity
############################################################################### #

# Load packages ----
# specify package location
.libPaths( "//sciensano.be/fs/1150_EPIVG_EpiInfect/15_WBE/PROJECTS/CodeLibraryR/librairies/R/4.5.2" )

# select packages
pkgs <- c("dplyr", "tidyr", "zoo", "writexl", "ggplot2")
# install packages if not yet done
install.packages(setdiff(pkgs, rownames(installed.packages())))
invisible(lapply(pkgs, FUN = library, character.only = TRUE))

# load data ----
# Belgian data are available here https://www.geo.be/catalog/details/9eec5acf-a2df-11ed-9952-186571a04de2?l=en
#| Metadata information
#| "siteName" is the name of the treatment plant
#| "collDTStart" is the date of sampling
#| "labName" is the name of the lab analysing the sample
#| "labProtocolID" is the protocol used to analyse the dample
#| "flowRate" is the flow rate measured at the inlet of the treatment plant during sampling
#| "popServ" is the population covered by the treatment plant
#| "measure" is the target measured
#| "value" is the result

# load sars-cov-2 data
df_sc <- read.csv("https://data.geo.be/ws/sciensano/wfs?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=sciensano:wastewatertreatmentplantscovid&outputFormat=csv")

# load pmmv data
df_pmmv <- read.csv("https://data.geo.be/ws/sciensano/wfs?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=sciensano:wastewatertreatmentplantspmmv&outputFormat=csv")

# join
df <- df_sc %>%
  rbind(df_pmmv)

# clean data ----
df <- df %>%
  select(siteName, collDTStart, labName, labProtocolID, flowRate, popServ, measure, value, quality)

# format date
df$date <- as.Date(df$collDTStart)

# set and subset dates
date_start <- as.Date("2024-09-01", format = "%Y-%m-%d")
date_end <- as.Date("2025-12-01", format = "%Y-%m-%d")

df <- df %>%
  filter(date > date_start & date < date_end)

# subset sars and pmmv data based on labProtocolID used betwen date_start and date_end
# display existing labProtocolID
# unique(df$labProtocolID)
df <- df %>%
  filter(labProtocolID %in%
           c("SC_COV_4.1",
             "UA_COV_4.0",
             "SC_PMMV_2.1",
             "UA_PMMV_2.0"))

# rename measures
# diplay existing measure
# unique(df$measure)
df[df$measure == "SARS-CoV-2 E gene", ]$measure <- "SARS"                                  
df[df$measure == "SARS-CoV-2 nucleocapsid gene, allele 2", ]$measure <- "SARS"             
df[df$measure == "Pepper mild mottle virus capsid protein gene region", ]$measure <- "PMMV"

# translate siteName to english 
df[df$siteName == "Bruxelles-Sud", ]$siteName <- "Brussels-South"
df[df$siteName == "Bruxelles-Nord", ]$siteName <- "Brussels-North"

# apply LOQ provided by the lab
df[df$measure == "PMMV" & df$value < 250, ]$value <- NA
df[df$measure == "SARS" & df$value < 8, ]$value <- NA

# remove outliers
df[df$quality == "Quality concerns", ]$value <- NA

# normalization ----
# compute mean of replicated analysis of each measure
df <- df %>%
  select(date, siteName, labName, flowRate, popServ, measure, value) %>%
  group_by(date, siteName, labName, flowRate, popServ, measure) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>% ungroup()

# pivot
df <- df %>%
  pivot_wider(names_from = measure, values_from = value)

# compute viral load (value_load), viral ratio (value_ratio)
df <- df %>%
  pivot_longer(cols = SARS, names_to = "measure", values_to = "value") %>%
  mutate(value_load = value*flowRate*24*1000000/popServ*100000,
         value_pmmv = value/PMMV)

# save
df_site_raw <- df

# outlier ----
# identify and exlude outliers for SARS
# tbl_outlier <- df_site_raw %>%
#   rename(SARS = value) %>%
#   filter(!is.na(SARS) & !is.na(PMMV)) %>%
#   group_by(siteName) %>%
#   mutate(sars_p95 = quantile(SARS, probs = 0.95, na.rm =T),
#          sars_p05 = quantile(SARS, probs = 0.05, na.rm =T),
#          pmmv_p95 = quantile(PMMV, probs = 0.95, na.rm =T),
#          pmmv_p05 = quantile(PMMV, probs = 0.05, na.rm =T)) %>%
#   ungroup() %>%
#   mutate(outlier = case_when(SARS > sars_p95 ~ TRUE,
#                              SARS < sars_p05 ~ TRUE,
#                              # PMMV > pmmv_p95 ~ TRUE,
#                              PMMV < pmmv_p05 ~ TRUE,
#                              .default = FALSE)) %>%
#   filter(outlier == TRUE) %>%
#   select(date, siteName, SARS, sars_p05, sars_p95, PMMV, pmmv_p05, pmmv_p95) %>%
#   arrange(siteName, desc(date))
# 
# tbl_outlier

# smoothening ----
# compute the linear extrapolation data
df <- df_site_raw %>%
  group_by(siteName) %>%
  complete(date = seq(min(date), max(date), "day")) %>%
  mutate(value_avg14d_past = na.approx(value, maxgap = 14, na.rm = FALSE),
         value_load_avg14d_past = na.approx(value_load, maxgap = 14, na.rm = FALSE),
         value_pmmv_avg14d_past = na.approx(value_pmmv, maxgap = 14, na.rm = FALSE))

# compute moving average on past 14 days
df <- df %>%
  group_by(siteName) %>%
  mutate(across(value_avg14d_past:value_pmmv_avg14d_past,
                ~ rollmean(.x, k = 14, fill = NA, na.rm = TRUE, align = "right")))

# save
df_site <- df

# national level ----
## aggregation ----
# compute weighted mean with factor being the population served by each site
df <- df_site_raw %>%
  select(date, popServ, value, value_load, value_pmmv) %>%
  mutate(siteName = "Belgium") %>%
  group_by(siteName, date) %>%
  summarise(across(value:value_pmmv, ~ weighted.mean(.x, popServ, na.rm=TRUE))) %>% ungroup()

## smoothening ----
# linear extrapolation data
df <- df %>%
  group_by(siteName) %>%
  complete(date = seq(min(date), max(date), "day")) %>%
  mutate(value_avg14d_past = na.approx(value, maxgap = 14, na.rm = FALSE),
         value_load_avg14d_past = na.approx(value_load, maxgap = 14, na.rm = FALSE),
         value_pmmv_avg14d_past = na.approx(value_pmmv, maxgap = 14, na.rm = FALSE))

# moving average on past 14 days
df <- df %>%
  group_by(siteName) %>%
  mutate(across(value_avg14d_past:value_pmmv_avg14d_past,
                ~ rollmean(.x, k = 14, fill = NA, na.rm = TRUE, align = "right")))

# save
df_nation <- df

# export data ----
# create folder if not existing
dir.create("./data", showWarnings = F)

# export as csv
write.table(df_site_raw, file = "./data/Belgium_export-site_raw.csv", sep = ";", dec = ".",
            col.names = TRUE, row.names = FALSE)

write.table(df_nation, file = "./data/Belgium_export-site.csv", sep = ";", dec = ".",
            col.names = TRUE, row.names = FALSE)

write.table(df_nation, file = "./data/Belgium_export-nation.csv", sep = ";", dec = ".",
            col.names = TRUE, row.names = FALSE)

# export as xls
write_xlsx(
  list(site_raw = df_site_raw, site = df_site, nation = df_nation),
  path = "./data/Belgium_export.xlsx"
)

# export as rds
saveRDS(list(df_site_raw, df_site, df_nation),
        file = "./data/Belgium_export.rds")

# visuals ----
# create folder if not existing
dir.create("./plot", showWarnings = F)

## viral cc ----
plot <- df_site_raw %>%
  filter(siteName %in% c("Brussels-North", "Brussels-South")) %>%
  rename(Site_name = siteName) %>%
  ggplot(aes(x = date, y = value, group = Site_name, color = Site_name)) +
  geom_point(na.rm = T) +
  geom_line(na.rm = T) +
  scale_x_date(date_labels = "%d/%m/%y (W%V)", breaks = "8 weeks") +
  labs(x = "", y = "Viral concentration (c./ml)") +
  theme_bw() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1, colour = "gray40"))

plot

# save
ggsave(file="./plot/Graph_be-bx_north_south-viral_cc.png",
       plot, width = 21, height = 12, dpi = 200)

## normalization ----
plot <- df_site_raw %>%
  filter(siteName %in% c("Brussels-North")) %>%
  mutate(value_load = value_load/2e10,
         value_pmmv = value_pmmv*1e5) %>%
  pivot_longer(cols = value:value_pmmv, names_to = "value", values_to = "result") %>%
  mutate(value = case_when(value == "value" ~ "Viral concentration (c/ml)",
                           value == "value_load" ~ "2x10E10 Viral load (c./100k inhab./24h)",
                           value == "value_pmmv" ~ "10E-5 viral ratio (c./c. PMMV)",
                           .default = value)) %>%
  ggplot(aes(x = date, y = result, group = value, color = value)) +
  geom_point(na.rm = T) +
  geom_line(na.rm = T) +
  scale_x_date(date_labels = "%d/%m/%y (W%V)", breaks = "8 weeks") +
  labs(x = "", y = "Wastewater concentration") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1, colour = "gray40"))

plot

# save
ggsave(file="./plot/Graph_be-bx_north-normalization.png",
       plot, width = 21, height = 12, dpi = 200)

## smoothening ----
plot <- df_site %>%
  filter(siteName %in% c("Brussels-North")) %>%
  ggplot() +
  geom_point(aes(x = date, y = value_pmmv), na.rm = T) +
  geom_line(aes(x = date, y = value_pmmv_avg14d_past), na.rm = T) +
  scale_x_date(date_labels = "%d/%m/%y (W%V)", breaks = "8 weeks") +
  labs(x = "", y = "Viral ratio (c./c.)") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1, colour = "gray40"))

plot

# save
ggsave(file="./plot/Graph_be-bx_north-normalization.png",
       plot, width = 21, height = 12, dpi = 200)

## aggregation ----
# aggregation brussels-north and brussels-south
# compute weighted mean with factor being the population served by each site
df1 <- df_site_raw %>%
  filter(siteName %in% c("Brussels-North", "Brussels-South")) %>%
  select(date, popServ, value, value_load, value_pmmv) %>%
  mutate(siteName = "Brussels aggregated") %>%
  group_by(siteName, date) %>%
  summarise(across(value:value_pmmv, ~ weighted.mean(.x, popServ, na.rm=TRUE))) %>% ungroup()

# smoothening
# linear extrapolation data
df1 <- df1 %>%
  group_by(siteName) %>%
  complete(date = seq(min(date), max(date), "day")) %>%
  mutate(value_avg14d_past = na.approx(value, maxgap = 14, na.rm = FALSE),
         value_load_avg14d_past = na.approx(value_load, maxgap = 14, na.rm = FALSE),
         value_pmmv_avg14d_past = na.approx(value_pmmv, maxgap = 14, na.rm = FALSE))

# moving average on past 14 days
df1 <- df1 %>%
  group_by(siteName) %>%
  mutate(across(value_avg14d_past:value_pmmv_avg14d_past,
                ~ rollmean(.x, k = 14, fill = NA, na.rm = TRUE, align = "right")))

# graph
plot <- df_site %>%
  filter(siteName %in% c("Brussels-North", "Brussels-South")) %>%
  rbind(df1) %>%
  rename(Site_name = siteName) %>%
  ggplot() +
  geom_point(aes(x = date, y = value_pmmv, group = Site_name, color = Site_name), na.rm = T) +
  geom_line(aes(x = date, y = value_pmmv_avg14d_past, group = Site_name, color = Site_name), linewidth = 1, na.rm = T) +
  scale_x_date(date_labels = "%d/%m/%y (W%V)", breaks = "8 weeks") +
  labs(x = "", y = "Viral ratio (c./c.)") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1, colour = "gray40"))

plot

# save
ggsave(file="./plot/Graph_be-bx_north_nation-viral_ratio.png",
       plot, width = 21, height = 12, dpi = 200)
