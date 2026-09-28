import os
from pathlib import Path

import pandas as pd


# ============================================================================
# PATHS
# ============================================================================

REGIONAL_BASE_ENV = "ROSMAP_REGIONAL_BASE"
regional_base = os.environ.get(REGIONAL_BASE_ENV)

if not regional_base:
    raise RuntimeError(
        f"{REGIONAL_BASE_ENV} is not set. "
        "Set it to the Brain Region Specificity directory."
    )

BASE = Path(regional_base).expanduser().resolve()

SOURCE_FILE = (
    BASE / "Metadata/AMP-AD_DiverseCohorts_individual_metadata.csv"
)

OLD_FILE = (
    BASE
    / "REGIONAL_FIRST_PASS_outputs"
    / "AMP_individual_metadata_cleaned_for_screen.csv"
)

OUTDIR = Path(
    "outputs/reviewer_revisions/regional_braak_fix"
)

OUTDIR.mkdir(parents=True, exist_ok=True)


# ============================================================================
# HELPERS
# ============================================================================

def clean_id(series):
    return (
        series.astype("string")
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
    )


def normalize_text(series):
    return (
        series.astype("string")
        .str.strip()
        .replace({
            "missing or unknown": pd.NA,
            "nan": pd.NA,
            "NaN": pd.NA,
            "": pd.NA,
        })
    )


# ============================================================================
# LOAD DATA
# ============================================================================

print("=" * 80)
print("REGIONAL BRAAK METADATA CORRECTION")
print("=" * 80)

source = pd.read_csv(SOURCE_FILE, low_memory=False)
old = pd.read_csv(OLD_FILE, low_memory=False)

print(f"Authoritative source participants: {len(source)}")
print(f"Existing derived participants:     {len(old)}")


# ============================================================================
# PARTICIPANT IDS
# ============================================================================

if "individualID" not in source.columns:
    raise ValueError("Source metadata is missing individualID.")

if "individual_id" not in old.columns:
    raise ValueError("Derived metadata is missing individual_id.")

source["individual_id"] = clean_id(source["individualID"])
old["individual_id"] = clean_id(old["individual_id"])

if source["individual_id"].isna().any():
    raise ValueError("Source metadata contains missing individualID.")

if old["individual_id"].isna().any():
    raise ValueError("Derived metadata contains missing individual_id.")

if source["individual_id"].duplicated().any():
    dup = source.loc[
        source["individual_id"].duplicated(False),
        "individual_id"
    ]
    raise ValueError(
        "Source metadata contains duplicate participant IDs: "
        f"{dup.head(10).tolist()}"
    )

if old["individual_id"].duplicated().any():
    dup = old.loc[
        old["individual_id"].duplicated(False),
        "individual_id"
    ]
    raise ValueError(
        "Derived metadata contains duplicate participant IDs: "
        f"{dup.head(10).tolist()}"
    )


# ============================================================================
# REQUIRED SOURCE VARIABLES
# ============================================================================

for col in ["Braak", "bScore"]:
    if col not in source.columns:
        raise ValueError(
            f"Authoritative source metadata is missing required column: {col}"
        )

for col in ["braak_raw", "braak_num", "bscore_raw", "bscore_num"]:
    if col not in old.columns:
        raise ValueError(
            f"Existing derived metadata is missing required column: {col}"
        )


# ============================================================================
# MERGE AUTHORITATIVE BRAAK VARIABLES
# ============================================================================

source_subset = source[
    ["individual_id", "Braak", "bScore"]
].copy()

merged = old.merge(
    source_subset,
    on="individual_id",
    how="left",
    validate="one_to_one",
)


# ============================================================================
# VERIFY SOURCE PROVENANCE
# ============================================================================

merged["_old_braak_raw"] = normalize_text(merged["braak_raw"])
merged["_source_braak"] = normalize_text(merged["Braak"])

merged["_old_bscore_raw"] = normalize_text(merged["bscore_raw"])
merged["_source_bscore"] = normalize_text(merged["bScore"])

braak_agree = (
    merged["_old_braak_raw"].fillna("<NA>")
    == merged["_source_braak"].fillna("<NA>")
)

bscore_agree = (
    merged["_old_bscore_raw"].fillna("<NA>")
    == merged["_source_bscore"].fillna("<NA>")
)

print("\n" + "=" * 80)
print("SOURCE / DERIVED AGREEMENT")
print("=" * 80)

print(
    "Braak raw agreement:",
    f"{braak_agree.sum()} / {len(merged)}"
)

print(
    "bScore raw agreement:",
    f"{bscore_agree.sum()} / {len(merged)}"
)

if not braak_agree.all():
    bad = merged.loc[
        ~braak_agree,
        ["individual_id", "braak_raw", "Braak"]
    ]

    print("\nBraak mismatches:")
    print(bad.head(20).to_string(index=False))

    raise ValueError(
        "Existing braak_raw does not agree with authoritative source Braak."
    )

if not bscore_agree.all():
    bad = merged.loc[
        ~bscore_agree,
        ["individual_id", "bscore_raw", "bScore"]
    ]

    print("\nbScore mismatches:")
    print(bad.head(20).to_string(index=False))

    raise ValueError(
        "Existing bscore_raw does not agree with authoritative source bScore."
    )


# ============================================================================
# DEFINE CORRECT BRAAK VARIABLES
# ============================================================================

stage_map = {
    "Stage I": 1,
    "Stage II": 2,
    "Stage III": 3,
    "Stage IV": 4,
    "Stage V": 5,
    "Stage VI": 6,
}

bin3_map = {
    "Braak Stage I-II": "I-II",
    "Braak Stage III-IV": "III-IV",
    "Braak Stage V-VI": "V-VI",
}

bin3_num_map = {
    "Braak Stage I-II": 2,
    "Braak Stage III-IV": 4,
    "Braak Stage V-VI": 6,
}

merged["braak_stage_num"] = (
    merged["_source_braak"]
    .map(stage_map)
    .astype("Int64")
)

merged["braak_bin3"] = (
    merged["_source_bscore"]
    .map(bin3_map)
)

merged["braak_bin3_num"] = (
    merged["_source_bscore"]
    .map(bin3_num_map)
    .astype("Int64")
)

# Preserve old variable explicitly rather than silently overwriting it.
merged["braak_num_legacy_2v6"] = merged["braak_num"]


# ============================================================================
# VALIDATE CORRECT I-VI MAPPING
# ============================================================================

observed = (
    merged.loc[
        merged["_source_braak"].notna(),
        ["_source_braak", "braak_stage_num"]
    ]
    .drop_duplicates()
    .sort_values("braak_stage_num")
)

observed_map = {
    stage: int(value)
    for stage, value in zip(
        observed["_source_braak"],
        observed["braak_stage_num"]
    )
}

print("\n" + "=" * 80)
print("CORRECTED I-VI MAPPING")
print("=" * 80)
print(observed.to_string(index=False))

if observed_map != stage_map:
    print("\nExpected mapping:")
    print(stage_map)

    print("\nObserved mapping:")
    print(observed_map)

    raise ValueError(
        "Corrected Braak I-VI mapping failed validation."
    )

print("\nCorrected Braak I-VI mapping validation: PASSED")


# ============================================================================
# VALIDATE THREE-BIN VARIABLE
# ============================================================================

observed_bin3 = (
    merged.loc[
        merged["_source_bscore"].notna(),
        ["_source_bscore", "braak_bin3", "braak_bin3_num"]
    ]
    .drop_duplicates()
    .sort_values("braak_bin3_num")
)

observed_bin3_num_map = {
    source_value: int(value)
    for source_value, value in zip(
        observed_bin3["_source_bscore"],
        observed_bin3["braak_bin3_num"]
    )
}

if observed_bin3_num_map != bin3_num_map:
    print("\nExpected three-bin mapping:")
    print(bin3_num_map)

    print("\nObserved three-bin mapping:")
    print(observed_bin3_num_map)

    raise ValueError(
        "Three-bin Braak mapping failed validation."
    )

print("Three-bin Braak mapping validation: PASSED")


# ============================================================================
# CHECK LEGACY VARIABLE
# ============================================================================

legacy_values = sorted(
    pd.to_numeric(
        merged["braak_num_legacy_2v6"],
        errors="coerce"
    )
    .dropna()
    .unique()
    .tolist()
)

print("\nLegacy braak_num unique values:")
print(legacy_values)

if legacy_values != [2.0, 6.0]:
    raise ValueError(
        "Expected legacy braak_num to contain exactly values 2 and 6."
    )

print("Legacy 2/6 Braak provenance check: PASSED")


# ============================================================================
# EXPORT
# ============================================================================

drop_internal = [
    "_old_braak_raw",
    "_source_braak",
    "_old_bscore_raw",
    "_source_bscore",
]

final = merged.drop(columns=drop_internal)

output_file = (
    OUTDIR
    / "regional_participant_metadata_corrected_braak.csv"
)

final.to_csv(
    output_file,
    index=False
)

audit = (
    final.groupby(
        [
            "braak_raw",
            "braak_stage_num",
            "bscore_raw",
            "braak_bin3",
            "braak_bin3_num",
        ],
        dropna=False,
    )
    .size()
    .reset_index(name="n")
    .sort_values(
        "braak_stage_num",
        na_position="last",
    )
)

audit_file = OUTDIR / "braak_mapping_audit.csv"

audit.to_csv(
    audit_file,
    index=False
)


# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("CORRECTED BRAAK MAPPING SUMMARY")
print("=" * 80)

print(audit.to_string(index=False))

print("\nContinuous Braak I-VI counts:")
print(
    final["braak_stage_num"]
    .value_counts(dropna=False)
    .sort_index()
    .to_string()
)

print("\nThree-bin Braak counts:")
print(
    final["braak_bin3"]
    .value_counts(dropna=False)
    .to_string()
)

print("\nLegacy 2/6 Braak counts:")
print(
    final["braak_num_legacy_2v6"]
    .value_counts(dropna=False)
    .sort_index()
    .to_string()
)

print("\n" + "=" * 80)
print("FILES WRITTEN")
print("=" * 80)

print(output_file)
print(audit_file)

print("\nALL REGIONAL BRAAK METADATA CHECKS PASSED.")
