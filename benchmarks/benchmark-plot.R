library(tidyverse)
library(ggimage)

data <- read_csv("glymotif/benchmarks/summary-v0.14.2-vs-v0.13.1.csv")

time_data <- data |>
  select(-speedup_vs_v0.13.1) |>
  pivot_longer(
    c(v0.14.2_median_sec, v0.13.1_median_sec),
    names_to = "version",
    values_to = "median_sec",
    names_pattern = "(.*?)_median_sec"
  ) |>
  mutate(
    version = factor(version, levels = c("v0.13.1", "v0.14.2")),
    function_name = factor(
      function_name,
      levels = c(
        "have_motif", "count_motif", "match_motif",
        "have_motifs", "count_motifs", "match_motifs"
      )
    )
  )

speedup_data <- data |>
  select(function_name, speed_up = speedup_vs_v0.13.1) |>
  mutate(
    speed_up = paste0("x ", scales::label_number(accuracy = 0.01)(speed_up)),
    function_name = factor(
      function_name,
      levels = c(
        "have_motif", "count_motif", "match_motif",
        "have_motifs", "count_motifs", "match_motifs"
      )
  )) |>
  inner_join(
    time_data |> summarise(y = max(median_sec) + 10, .by = function_name),
    by = join_by(function_name)
  )

time_barplots <- ggplot(time_data) +
  geom_col(aes(function_name, median_sec, fill = version), position = "dodge") +
  scale_fill_brewer(palette = "Blues") +
  geom_text(
    data = speedup_data,
    mapping = aes(x = function_name, y = y, label = speed_up),
    size = 3.5,
  ) +
  labs(
    title = "Benchmarking",
    x = "Function Name",
    y = "Seconds",
    fill = "version"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggpreview(plot = time_barplots, width = 4.5, height = 4)
