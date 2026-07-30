library(tidyverse)

aesdffhjj = c(1,23,34)

filtered_ids = data.table
remove(data.table)

table(filtered_ids$celltype)

filtered_ids_rem_lq = filtered_ids |> filter(celltype != "LowQuality")

filtered_ids_rem_lq$sampleid[filtered_ids_rem_lq$sampleid == "CA1"] = "Y1"
filtered_ids_rem_lq$sampleid[filtered_ids_rem_lq$sampleid == "CA2"] = "Y2"
filtered_ids_rem_lq$sampleid[filtered_ids_rem_lq$sampleid == "CA3"] = "Y3"
filtered_ids_rem_lq$sampleid[filtered_ids_rem_lq$sampleid == "CA4"] = "Y4"

filtered_ids_rem_lq$sampleid[filtered_ids_rem_lq$sampleid == "WBI-1"] = "YWBI1"
filtered_ids_rem_lq$sampleid[filtered_ids_rem_lq$sampleid == "WBI-2"] = "YWBI2"
filtered_ids_rem_lq$sampleid[filtered_ids_rem_lq$sampleid == "WBI-3"] = "YWBI3"
filtered_ids_rem_lq$sampleid[filtered_ids_rem_lq$sampleid == "WBI-4"] = "YWBI4"

table(filtered_ids$celltype)
filtered_Lowquality = filtered_ids |> filter(celltype !="LowQuality")

filtered_lowquality$sampleid[filtered_lowquality$sampleid == "CA1"] = "Y1"
filtered_lowquality$sampleid[filtered_lowquality$sampleid == "CA2"] = "Y2"
filtered_lowquality$sampleid[filtered_lowquality$sampleid == "CA3"] = "Y3"
filtered_lowquality$sampleid[filtered_lowquality$sampleid == "CA4"] = "Y4"

filtered_lowquality$sampleid[filtered_lowquality$sampleid == "WBI-1"] = "YWBI1"
filtered_lowquality$sampleid[filtered_lowquality$sampleid == "WBI-2"] = "YWBI2"
filtered_lowquality$sampleid[filtered_lowquality$sampleid == "WBI-3"] = "YWBI3"
filtered_lowquality$sampleid[filtered_lowquality$sampleid == "WBI-4"] = "YWBI4"

filtered_lowquality_wide <- table(
  filtered_lowquality$celltype,
  filtered_lowquality$sampleid
) |>
  prop.table(margin = 2) |>
  (\(x) x * 100)() |>
  as.data.frame.matrix() |>
  rownames_to_column("celltype")

filtered_lowquality_long <- filtered_lowquality_wide |>
  pivot_longer(cols = -celltype, names_to = "sample", values_to = "abundance")


ggplot(filtered_lowquality_long,
       aes (x = sample, y = abundance, fill = celltype))+ 
  geom_bar(stat = "identity") + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) + 
  theme_classic() + 
  labs(x = "sample", y = "Cell Type Abundance (%)", fill = "Cell Type") +
  scale_fill_brewer(palette = "Spectral")



tiff("plot2", units = "in", width = 6, height = 6, res = 300)



######################


filtered_lowquality_EC = filtered_lowquality[grep("EC", filtered_lowquality$subtype), ]

filtered_lowquality_EC_filtered = filtered_lowquality_EC |> filter(subtype != "Fenestrated EC")


ggplot(filtered_lowquality_EC_filtered,
       aes (x = score, y = subtype, fill = subtype))+ 
  geom_boxplot() + 
  #scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 100, 10)) + 
  theme_classic() + 
  labs(x = "Senescent score", y = "Subtype") +
  coord_flip() + # fill = "Cell Type") +
  scale_fill_brewer(palette = "Spectral")

#########################

control = c("Y1", "Y2", "Y3", "Y4")
treatment = c("YWBI1", "YWBI2", "YWBI3", "YWBI4")

plot.data = filtered_lowquality_EC_filtered |>
  mutate(condition = case_when(
    sampleid %in% control ~ "Control",
    sampleid %in% treatment ~ "Treatment"
  ))

tail(sort(filtered_lowquality_EC_filtered$score))

threshold = 0.5 * 0.168

plot.data.filtered = plot.data |> filter(score > threshold)
plot.data.filtered = plot.data.filtered |> filter(subtype %in% c("Arterial-capillary EC", "Arterial EC", "Capillary-venous EC", "Venous EC"))

ggplot(plot.data.filtered,
       aes(x = subtype, y = score, fill = condition)) +
  #stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "% of Senescent cells", y = "Subtypes")
  #coord_flip()
  #scale_fill_brewer(palette = "Spectral")
