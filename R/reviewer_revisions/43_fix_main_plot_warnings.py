from pathlib import Path
import shutil

files = {
    "r10": Path("R/10_make_main_figure_1_pathway_remodeling.R"),
    "r11": Path("R/11_make_main_figure_2_collapse_heterogeneity.R"),
    "r12": Path("R/12_make_main_figure_3_pathology_coupling.R"),
    "r15": Path("R/15_make_main_figure_6_candidate_classification.R"),
}

backup = Path(
    "outputs/reviewer_revisions/"
    "RNA_batch_canonicalization/pre_plot_warning_cleanup_snapshot"
)
backup.mkdir(parents=True, exist_ok=True)

for p in files.values():
    shutil.copy2(p, backup / f"{p.stem}_PRE_PLOT_WARNING_FIX.R")


def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise SystemExit(
            f"{label}: expected exactly one target; found {n}."
        )
    return text.replace(old, new, 1)


# ============================================================
# Figure 1: ASCII-safe minus in axis label
# ============================================================

p = files["r10"]
s = p.read_text()

s = replace_once(
    s,
    'x = "Protein remodeling magnitude − RNA remodeling magnitude"',
    'x = "Protein remodeling magnitude - RNA remodeling magnitude"',
    "Figure 1 Unicode minus"
)

p.write_text(s)


# ============================================================
# Figure 2:
#   - shape 21 outline uses stroke, not linewidth
#   - constant significance annotations use annotate()
# ============================================================

p = files["r11"]
s = p.read_text()

s = replace_once(
    s,
    '''  geom_point(
    aes(fill = late_decline_group),
    shape = 21,
    size = 1.85,
    color = "grey35",
    linewidth = 0.18
  ) +''',
    '''  geom_point(
    aes(fill = late_decline_group),
    shape = 21,
    size = 1.85,
    color = "grey35",
    stroke = 0.18
  ) +''',
    "Figure 2 point outline"
)

s = replace_once(
    s,
    '''  geom_segment(
    aes(
      x = 1,
      xend = 2,
      y = late_effect_difference_ymax + late_effect_difference_yrange * 0.26,
      yend = late_effect_difference_ymax + late_effect_difference_yrange * 0.26
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = late_effect_difference_ymax + late_effect_difference_yrange * 0.48,
      label = paste0(late_effect_difference_stats_tbl$star, "\\n", late_effect_difference_stats_tbl$p_label)
    ),
    inherit.aes = FALSE,
    size = 2.8,
    fontface = "bold",
    lineheight = 0.95
  ) +''',
    '''  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = late_effect_difference_ymax + late_effect_difference_yrange * 0.26,
    yend = late_effect_difference_ymax + late_effect_difference_yrange * 0.26,
    color = "black",
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = 1.5,
    y = late_effect_difference_ymax + late_effect_difference_yrange * 0.48,
    label = paste0(
      late_effect_difference_stats_tbl$star,
      "\\n",
      late_effect_difference_stats_tbl$p_label
    ),
    size = 2.8,
    fontface = "bold",
    lineheight = 0.95
  ) +''',
    "Figure 2 Panel B annotation"
)

s = replace_once(
    s,
    '''  geom_segment(
    aes(
      x = 1,
      xend = 2,
      y = centrality_ymax + centrality_yrange * 0.28,
      yend = centrality_ymax + centrality_yrange * 0.28
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      x = 1.5,
      y = centrality_ymax + centrality_yrange * 0.51,
      label = paste0(centrality_stats_tbl$star, "\\n", centrality_stats_tbl$p_label)
    ),
    inherit.aes = FALSE,
    size = 2.8,
    fontface = "bold",
    lineheight = 0.95
  ) +''',
    '''  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = centrality_ymax + centrality_yrange * 0.28,
    yend = centrality_ymax + centrality_yrange * 0.28,
    color = "black",
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = 1.5,
    y = centrality_ymax + centrality_yrange * 0.51,
    label = paste0(
      centrality_stats_tbl$star,
      "\\n",
      centrality_stats_tbl$p_label
    ),
    size = 2.8,
    fontface = "bold",
    lineheight = 0.95
  ) +''',
    "Figure 2 Panel C annotation"
)

p.write_text(s)


# ============================================================
# Figure 3:
#   - geom_label border parameter: linewidth
#   - shape 21 point outline: stroke
# ============================================================

p = files["r12"]
s = p.read_text()

if s.count("label.size = 0.25") != 2:
    raise SystemExit(
        "Figure 3: expected two label.size = 0.25 occurrences; "
        f"found {s.count('label.size = 0.25')}."
    )

s = s.replace(
    "label.size = 0.25",
    "linewidth = 0.25"
)

s = replace_once(
    s,
    "label.size = 0.22",
    "linewidth = 0.22",
    "Figure 3 Panel C label border"
)

s = replace_once(
    s,
    '''    size = 2.15,
    linewidth = 0.22
  ) +''',
    '''    size = 2.15,
    stroke = 0.22
  ) +''',
    "Figure 3 point outline"
)

p.write_text(s)


# ============================================================
# Figure 6: ASCII-safe dash in panel title
# ============================================================

p = files["r15"]
s = p.read_text()

s = replace_once(
    s,
    'title = "B. Late-decline–Braak target landscape"',
    'title = "B. Late-decline-Braak target landscape"',
    "Figure 6 Unicode en dash"
)

p.write_text(s)

print("PLOT WARNING PATCH APPLIED")
print("R/10: Unicode minus -> ASCII hyphen")
print("R/11: point stroke + scalar annotations")
print("R/12: label linewidth + point stroke")
print("R/15: Unicode en dash -> ASCII hyphen")
