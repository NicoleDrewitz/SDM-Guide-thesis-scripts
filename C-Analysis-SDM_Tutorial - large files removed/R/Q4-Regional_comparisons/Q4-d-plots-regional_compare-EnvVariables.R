# Created by Nicole Drewitz with Claude AI
# Created on 31 May, 2026
# paste to the end of Q4-d-regional_compare script
# create dot plots comparing mean environmental values between two climate scenarios, across various regions.

# ── Cleveland Dot Plot Matrix ──────────────────────────────────────────────────
library(ggplot2)
library(forcats)
library(ggh4x)
library(grid)        # for unit()
library(ggtext)      # for element_markdown() if you want rich text annotations

# ── A4 landscape dimensions (mm → inches, minus 20mm margins each side) ──────
# A4 landscape: 297mm × 210mm
# Minus 20mm margins: 257mm × 170mm usable
a4_w_in <- (297 - 40) / 25.4   # ≈ 10.12 inches (usable width)
a4_h_in <- (210 - 40) / 25.4   # ≈ 6.69  inches (usable height)

fig_dir <- here::here("figures")

plot_vars <- setdiff(available_vars, c( # exclude these variables
  "kg0", "Topo_northness",
  "soil_waterV33_LoamClay_0033kPa_0-5-30cm_mean_averaged",
  "soil_soc_0-5-30cm_mean_averaged",
  "soil_phh2o_0-5-30cm_mean_averaged"
))

all_data |>
  distinct(file_source) |>
  pull(file_source) |>
  walk(\(src) {
    
    group_name <- regmatches(src, regexpr("(?<=sampling-)(.+?)(?=-n)", src, perl = TRUE))
    
    regions_in_file <- all_data |>
      filter(file_source == src) |>
      distinct(region) |>
      pull(region)
    
    # ── Optional: exclude regions with "Inland" prefix ───────────────────────────
    # Comment out the next two lines to include all regions
    regions_in_file <- regions_in_file[!grepl("^Inland", regions_in_file)] # exclude regions that start with "Inland"
    
    plot_data <- final_table |>
      filter(Region %in% regions_in_file) |>
      filter(`Environmental variable` %in% var_titles[plot_vars]) |>
      mutate(
        current     = `Current climate (mean)`,
        future      = current + as.numeric(`Change with medium emissions (mean)`),
        Region      = factor(Region, levels = regions_in_file),
        facet_label = paste0(`Environmental variable`, "\n(", Unit, ")"),
        facet_label = fct_reorder(
          facet_label,
          match(`Environmental variable`, unname(var_titles))
        )
      ) |>
      select(`Environmental variable`, Unit, Region, facet_label, current, future) |>
      pivot_longer(c(current, future), names_to = "period", values_to = "value") |>
      mutate(period = factor(period, levels = c("current", "future"),
                             labels = c("Current climate", "SSP3-7.0"))) |>
      filter(!is.na(value))
    
    # ── Layout ─────────────────────────────────────────────────────────────────
    n_vars   <- length(plot_vars)
    n_cols   <- ceiling(n_vars / 2)
    n_rows   <- 2
    n_regions <- length(regions_in_file)
    
    # ── Per-variable annotation data ───────────────────────────────────────────
    # Compute range labels (current min → future max) per facet for annotation
    annotation_data <- plot_data |>
      group_by(facet_label) |>
      summarise(
        x_current = mean(value[period == "Current climate"], na.rm = TRUE),
        x_future  = mean(value[period == "SSP3-7.0"],        na.rm = TRUE),
        delta     = x_future - x_current,
        .groups   = "drop"
      ) |>
      mutate(
        # Arrow direction label
        delta_label = case_when(
          delta >  0.01 ~ paste0("+", round(delta, 1)),
          delta < -0.01 ~ as.character(round(delta, 1)),
          TRUE          ~ "≈ 0"
        ),
        # Place annotation at x = future mean, y just above top region
        y_pos = as.numeric(factor(regions_in_file[1],
                                  levels = regions_in_file)) + 0.55
      )
    
    # ── Y-scale list: labels only for left column ──────────────────────────────
    left_col_idx <- c(1, n_cols + 1)
    
    y_scales <- lapply(seq_len(n_vars), function(i) {
      if (i %in% left_col_idx) scale_y_discrete()
      else                      scale_y_discrete(labels = NULL)
    })
    
    # ── Replace the geom_segment + geom_text annotation block ────────────────────
    # Place delta arrow + label BELOW the x-axis using annotation_custom + coord,
    # or more simply: use geom_ with negative y positions outside the panel,
    # combined with coord_cartesian(clip = "off").
    
    p <- ggplot(plot_data, aes(x = value, y = Region, colour = period)) +
      
      geom_line(
        aes(group = interaction(Region, facet_label)),
        colour    = "grey70",
        linewidth = 0.4
      ) +
      
      geom_point(size = 2.5, alpha = 0.9) +
      
      # ── Δ arrow: drawn BELOW axis at y = -0.6, clipped off ───────────────────
      geom_segment(
        data = annotation_data,
        aes(
          x    = x_current, xend = x_future,
          y    = -0.55,     yend = -0.55
        ),
        colour      = "grey30",
        linewidth   = 0.5,
        arrow       = arrow(length = unit(0.12, "cm"), type = "closed"),
        inherit.aes = FALSE
      ) +
      
      # ── Δ label: sits just below the arrow ────────────────────────────────────
      geom_text(
        data = annotation_data,
        aes(
          x     = (x_current + x_future) / 2,
          y     = -1.4, # use to change vertical distance between label and arrow
          label = delta_label
        ),
        size        = 2.5,
        colour      = "grey20",
        fontface    = "bold",
        inherit.aes = FALSE
      ) +
      
      facet_wrap(~ facet_label, scales = "free_x", nrow = n_rows) +
      facetted_pos_scales(y = y_scales) +
      
      # ── Allow drawing outside panel bounds ────────────────────────────────────
      coord_cartesian(clip = "off") +
      
      scale_colour_manual(
        values = c("Current climate" = "#fdbb84", "SSP3-7.0" = "#e34a33"),
        guide  = guide_legend(
          title        = NULL,
          override.aes = list(size = 3),
          keywidth     = unit(1.2, "cm")
        )
      ) +
      
      labs(
        title    = paste0("Current climate vs SSP3-7.0 \u2014 ", group_name),
        subtitle = "Arrows and \u0394 values below each panel show mean change (current \u2192 SSP3-7.0)",
        x        = NULL, y = NULL, colour = NULL
      ) +
      
      theme_minimal(base_size = 9) +
      theme(
        strip.text         = element_text(size = 7.5, face = "bold", hjust = 0),
        strip.background   = element_rect(fill = "grey94", colour = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        panel.spacing.x    = unit(0.6, "lines"),
        panel.spacing.y    = unit(0.3, "lines"),
        legend.position    = "top",
        legend.text        = element_text(size = 8),
        legend.margin      = margin(0, 0, 2, 0),
        legend.key.height = unit(0.2, "cm"),
        plot.title         = element_text(face = "bold", size = 11),
        plot.subtitle      = element_text(size = 7.5, colour = "grey40",
                                          margin = margin(b = 4)),
        plot.caption       = element_text(size = 6.5, colour = "grey50",
                                          hjust = 0, margin = margin(t = 6)),
        axis.text.x        = element_text(size = 7),
        axis.text.y        = element_text(size = 7.5)
      )
    # ── Save at exactly A4 landscape usable area ──────────────────────────────
    ggsave(
      filename = file.path(fig_dir,
                           paste0("Q4-d-climate_change_dotplot-", group_name, ".png")),
      plot     = p,
      width    = a4_w_in,
      height   = a4_h_in,
      dpi      = 300,
      bg       = "white"
    )
    cat("Saved:", group_name, "\n")
  })

# ── Save at exactly A4 landscape usable area ──────────────────────────────
ggsave(
  filename = file.path(fig_dir,
                       paste0("Q4-d-climate_change_dotplot-", group_name, ".png")),
  plot     = p,
  width    = a4_w_in,
  height   = a4_h_in,
  dpi      = 300,
  bg       = "white"
)
cat("Save")