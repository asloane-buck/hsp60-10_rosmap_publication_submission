from pathlib import Path
import shutil

p04 = Path("R/04_build_pathway_sets_all_clients.R")
p15 = Path("R/15_make_main_figure_6_candidate_classification.R")

s04 = p04.read_text()
s15 = p15.read_text()

old04 = '''readxl::read_excel(cfg$mitocarta_xls, sheet = "A Human MitoCarta3.0")'''

new04 = '''readxl::read_excel(
      cfg$mitocarta_xls,
      sheet = "A Human MitoCarta3.0",
      col_types = "text",
      na = c("", "NA")
    )'''

if s04.count(old04) != 1:
    raise SystemExit(
        f"Expected exactly one R/04 MitoCarta read; found {s04.count(old04)}."
    )

s04 = s04.replace(old04, new04, 1)


old15_preview = '''readxl::read_excel(mitocarta_path, sheet = sh, n_max = 5)'''

new15_preview = '''readxl::read_excel(
      mitocarta_path,
      sheet = sh,
      n_max = 5,
      col_types = "text",
      na = c("", "NA")
    )'''

if s15.count(old15_preview) != 1:
    raise SystemExit(
        "Expected exactly one Figure 6 MitoCarta preview read; "
        f"found {s15.count(old15_preview)}."
    )

s15 = s15.replace(old15_preview, new15_preview, 1)


old15_full = '''raw <- readxl::read_excel(mitocarta_path, sheet = sh)'''

new15_full = '''raw <- readxl::read_excel(
    mitocarta_path,
    sheet = sh,
    col_types = "text",
    na = c("", "NA")
  )'''

if s15.count(old15_full) != 1:
    raise SystemExit(
        "Expected exactly one Figure 6 full MitoCarta read; "
        f"found {s15.count(old15_full)}."
    )

s15 = s15.replace(old15_full, new15_full, 1)


backup = Path(
    "outputs/reviewer_revisions/"
    "RNA_batch_canonicalization/pre_warning_cleanup_snapshot"
)
backup.mkdir(parents=True, exist_ok=True)

shutil.copy2(
    p04,
    backup / "04_build_pathway_sets_all_clients_PRE_MITOCARTA_TYPE_FIX.R"
)
shutil.copy2(
    p15,
    backup / "15_make_main_figure_6_candidate_classification_PRE_MITOCARTA_TYPE_FIX.R"
)

p04.write_text(s04)
p15.write_text(s15)

print("PATCH APPLIED")
print("R/04: MitoCarta read is deterministic text import")
print("R/15: MitoCarta preview/full reads are deterministic text imports")
print('Missing-value strings: "" and "NA"')
