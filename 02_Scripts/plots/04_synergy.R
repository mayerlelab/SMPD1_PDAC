
library(readxl)
library(tidyverse)
library(synergyfinder)


screening_data <- readxl::read_excel("treatment summery-KPC1050-v3.xlsx", 
                                     "calculation")

screening_data <- screening_data %>%
  group_by(conc_r, conc_c) %>%
  mutate(response = median(response)) %>%
  ungroup() %>%
  unique()


res <- ReshapeData(
  data = screening_data,
  data_type = "viability",
  impute = TRUE,
  impute_method = NULL,
  noise = TRUE,
  seed = 1)

res <- CalculateSynergy(
  data = res,
  method = c("ZIP", "HSA", "Bliss", "Loewe"),
  Emin = NA,
  Emax = NA,
  correct_baseline = "non")

res1 <- CalculateSensitivity(
  data = res,
  correct_baseline = "non"
)

res1$drug_pairs

p <- Plot2DrugContour(
  data = res,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  summary_statistic = c("mean","quantile_25", "quantile_75")
) + theme_bw()  +
  geom_contour(colour = "black",
               size = 0.5) +
  scale_fill_gradientn(
    colours = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(100),
    oob = scales::squish,
    limits = c(-50, 50),
    name = "Synergy score"
  ) +
  guides(fill = guide_colourbar(
    barwidth = unit(0.3, "cm"),
    ticks.colour = "black",
    frame.colour = "black"
  ))


print(p)

ggsave("synergy_score.pdf",
       plot = p,
       width = 10, 
       height = 8,
       units = "cm",
       dpi = 300)

Plot2DrugSurface(
  data = res,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  summary_statistic = c("mean", "quantile_25", "median", "quantile_75")
) 
