#!/usr/bin/env python3
"""Append Supplementary Table 13 without reserializing existing worksheets.

ST13 summarizes the Reviewer 2 complementary-direction sensitivity:
direction counts/Fisher test, abundance-matched null results,
conventional AD-versus-MCI direction counts, and the FDR-significant
shared Hsp60/10 clients with higher protein abundance in AD versus MCI.

Existing worksheet XML is preserved byte-for-byte.
"""

import argparse
import csv
import math
from pathlib import Path
import zipfile
import xml.etree.ElementTree as ET

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG = "http://schemas.openxmlformats.org/package/2006/relationships"
CT = "http://schemas.openxmlformats.org/package/2006/content-types"

ET.register_namespace("", MAIN)
ET.register_namespace("r", REL)

SHEET_NAME = "ST13_Complementary_direction"


def xmlbytes(root):
    return ET.tostring(
        root,
        encoding="utf-8",
        xml_declaration=True,
    )


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def col_letter(index):
    out = ""
    while index:
        index, rem = divmod(index - 1, 26)
        out = chr(65 + rem) + out
    return out


def numeric(value):
    if isinstance(value, bool):
        return False
    if isinstance(value, (int, float)):
        return math.isfinite(float(value))
    try:
        x = float(str(value))
        return math.isfinite(x)
    except (TypeError, ValueError):
        return False


def add_cell(row_node, row_index, col_index, value):
    if value is None or value == "":
        return

    ref = f"{col_letter(col_index)}{row_index}"

    if numeric(value):
        cell = ET.SubElement(
            row_node,
            f"{{{MAIN}}}c",
            {"r": ref},
        )
        node = ET.SubElement(
            cell,
            f"{{{MAIN}}}v",
        )
        node.text = format(float(value), ".17g")
        return

    cell = ET.SubElement(
        row_node,
        f"{{{MAIN}}}c",
        {
            "r": ref,
            "t": "inlineStr",
        },
    )
    inline = ET.SubElement(
        cell,
        f"{{{MAIN}}}is",
    )
    node = ET.SubElement(
        inline,
        f"{{{MAIN}}}t",
    )
    node.text = str(value)


def build_sheet(
    fisher,
    matched,
    conventional,
    higher,
):
    ws = ET.Element(
        f"{{{MAIN}}}worksheet"
    )

    views = ET.SubElement(
        ws,
        f"{{{MAIN}}}sheetViews",
    )
    ET.SubElement(
        views,
        f"{{{MAIN}}}sheetView",
        {"workbookViewId": "0"},
    )

    cols = ET.SubElement(
        ws,
        f"{{{MAIN}}}cols",
    )
    widths = [
        32, 24, 18, 18, 18, 18, 18, 18, 18,
        18, 18, 18, 18, 18, 18, 18,
    ]
    for i, width in enumerate(widths, 1):
        ET.SubElement(
            cols,
            f"{{{MAIN}}}col",
            {
                "min": str(i),
                "max": str(i),
                "width": str(width),
                "customWidth": "1",
            },
        )

    data = ET.SubElement(
        ws,
        f"{{{MAIN}}}sheetData",
    )

    rows = []

    rows.append([
        "Supplementary Table 13. "
        "Complementary-direction sensitivity for "
        "AD-versus-MCI protein abundance"
    ])
    rows.append([
        "Signed protein effect = AD minus MCI. "
        "Positive values indicate higher abundance in AD; "
        "negative values indicate lower abundance in AD."
    ])
    rows.append([])

    rows.append([
        "A. Directional enrichment among Hsp60/10 clients "
        "versus non-client mitochondrial proteins"
    ])
    fisher_cols = [
        "hsp_lower",
        "hsp_higher",
        "hsp_lower_fraction",
        "background_lower",
        "background_higher",
        "background_lower_fraction",
        "odds_ratio",
        "conf_low",
        "conf_high",
        "p_value_two_sided",
    ]
    rows.append(fisher_cols)
    rows.append([
        fisher[col]
        for col in fisher_cols
    ])
    rows.append([])

    rows.append([
        "B. Abundance-matched complementary-direction sensitivity"
    ])
    matched_cols = [
        "metric",
        "label",
        "primary_tail",
        "observed",
        "null_mean",
        "null_sd",
        "observed_to_null_ratio",
        "z_score",
        "primary_empirical_p_reported",
    ]
    rows.append(matched_cols)
    for source in matched:
        rows.append([
            source[col]
            for col in matched_cols
        ])
    rows.append([])

    rows.append([
        "C. Conventional shared-client AD-versus-MCI "
        "protein differential-abundance direction counts"
    ])
    conventional_cols = [
        "n_shared_clients",
        "n_protein_fdr_significant",
        "n_fdr_lower_in_ad",
        "n_fdr_higher_in_ad",
        "n_all_lower_in_ad",
        "n_all_higher_in_ad",
    ]
    rows.append(conventional_cols)
    rows.append([
        conventional[col]
        for col in conventional_cols
    ])
    rows.append([])

    rows.append([
        "D. FDR-significant shared Hsp60/10 clients "
        "with higher protein abundance in AD versus MCI"
    ])
    higher_cols = [
        "gene_symbol",
        "protein_effect",
        "protein_std_error",
        "protein_conf_low",
        "protein_conf_high",
        "protein_p_value",
        "protein_fdr",
        "rna_effect",
        "rna_p_value",
        "rna_fdr",
        "significance_class",
        "direction_concordant",
    ]
    rows.append(higher_cols)
    for source in higher:
        rows.append([
            source[col]
            for col in higher_cols
        ])

    rows.append([])
    rows.append([
        "Source",
        "outputs/reviewer_sensitivity/complementary_direction/"
        "direction_fisher_exact.csv",
    ])
    rows.append([
        "Source",
        "outputs/reviewer_sensitivity/complementary_direction/"
        "matched_null_summary.csv",
    ])
    rows.append([
        "Source",
        "outputs/reviewer_sensitivity/complementary_direction/"
        "conventional_DE_AD_vs_MCI_direction_summary.csv",
    ])
    rows.append([
        "Source",
        "outputs/supplemental_figures/tables/"
        "supplementary_figure_7_conventional_differential/"
        "Supplementary_Table_conventional_RNA_protein_clients.csv",
    ])

    for row_index, values in enumerate(rows, 1):
        row_node = ET.SubElement(
            data,
            f"{{{MAIN}}}row",
            {"r": str(row_index)},
        )
        for col_index, value in enumerate(
            values,
            1,
        ):
            add_cell(
                row_node,
                row_index,
                col_index,
                value,
            )

    return ws


def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )
    parser.add_argument(
        "--existing",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--sensitivity-dir",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--conventional-clients",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )
    args = parser.parse_args()

    if args.output.resolve() == args.existing.resolve():
        parser.error(
            "Output must differ from input."
        )
    if args.output.exists():
        parser.error(
            "Output already exists."
        )

    fisher_path = (
        args.sensitivity_dir
        / "direction_fisher_exact.csv"
    )
    matched_path = (
        args.sensitivity_dir
        / "matched_null_summary.csv"
    )
    conventional_path = (
        args.sensitivity_dir
        / "conventional_DE_AD_vs_MCI_direction_summary.csv"
    )

    for path in [
        fisher_path,
        matched_path,
        conventional_path,
        args.conventional_clients,
    ]:
        if not path.is_file():
            parser.error(
                f"Missing required ST13 source: {path}"
            )

    fisher_rows = read_csv(fisher_path)
    matched = read_csv(matched_path)
    conventional_rows = read_csv(
        conventional_path
    )
    clients = read_csv(
        args.conventional_clients
    )

    if len(fisher_rows) != 1:
        raise ValueError(
            "Expected one Fisher-summary row."
        )
    if len(conventional_rows) != 1:
        raise ValueError(
            "Expected one conventional-summary row."
        )

    fisher = fisher_rows[0]
    conventional = conventional_rows[0]

    if (
        int(float(fisher["hsp_lower"])) != 243
        or int(float(fisher["hsp_higher"])) != 63
        or int(float(fisher["background_lower"])) != 427
        or int(float(fisher["background_higher"])) != 182
    ):
        raise ValueError(
            "Directional counts do not match frozen "
            "Reviewer 2 sensitivity results."
        )

    metric_names = {
        row["metric"]
        for row in matched
    }
    expected_metrics = {
        "lower_abundance_magnitude",
        "higher_abundance_magnitude",
        "absolute_difference",
    }
    if metric_names != expected_metrics:
        raise ValueError(
            "Unexpected matched-null metrics: "
            f"{sorted(metric_names)}"
        )

    if (
        int(float(conventional["n_shared_clients"])) != 257
        or int(float(conventional["n_protein_fdr_significant"])) != 39
        or int(float(conventional["n_fdr_lower_in_ad"])) != 36
        or int(float(conventional["n_fdr_higher_in_ad"])) != 3
    ):
        raise ValueError(
            "Conventional direction counts do not match "
            "the frozen shared-client results."
        )

    higher = []
    for row in clients:
        if row.get("contrast") != "AD_vs_MCI":
            continue
        if str(row.get(
            "is_hsp60_10_client"
        )).upper() != "TRUE":
            continue
        if str(row.get(
            "protein_significant_fdr05"
        )).upper() != "TRUE":
            continue
        if float(row["protein_effect"]) <= 0:
            continue
        higher.append(row)

    higher.sort(
        key=lambda row: row["gene_symbol"]
    )

    if {
        row["gene_symbol"]
        for row in higher
    } != {"GAPDH", "GNG12", "PPIA"}:
        raise ValueError(
            "Expected exactly GAPDH, GNG12, and PPIA "
            "as FDR-significant higher-abundance "
            "shared clients."
        )

    if len(higher) != 3:
        raise ValueError(
            "Expected exactly three higher-abundance "
            "FDR-significant shared clients."
        )

    with zipfile.ZipFile(
        args.existing
    ) as archive:
        original = {
            name: archive.read(name)
            for name in archive.namelist()
        }

    wb = ET.fromstring(
        original["xl/workbook.xml"]
    )
    rels = ET.fromstring(
        original[
            "xl/_rels/workbook.xml.rels"
        ]
    )
    content = ET.fromstring(
        original["[Content_Types].xml"]
    )

    sheets = wb.find(
        f"{{{MAIN}}}sheets"
    )
    names = {
        sheet.get("name")
        for sheet in sheets
    }

    if SHEET_NAME in names:
        parser.error(
            f"{SHEET_NAME} already exists."
        )

    next_id = max(
        int(sheet.get("sheetId"))
        for sheet in sheets
    ) + 1

    used_rel_ids = {
        rel.get("Id")
        for rel in rels
    }

    rid = f"rIdTable13_{next_id}"
    while rid in used_rel_ids:
        rid += "_"

    path = (
        f"xl/worksheets/"
        f"table13_{next_id}.xml"
    )

    ET.SubElement(
        sheets,
        f"{{{MAIN}}}sheet",
        {
            "name": SHEET_NAME,
            "sheetId": str(next_id),
            f"{{{REL}}}id": rid,
        },
    )

    ET.SubElement(
        rels,
        f"{{{PKG}}}Relationship",
        {
            "Id": rid,
            "Type": REL + "/worksheet",
            "Target": (
                "worksheets/"
                + path.rsplit("/", 1)[1]
            ),
        },
    )

    ET.SubElement(
        content,
        f"{{{CT}}}Override",
        {
            "PartName": "/" + path,
            "ContentType": (
                "application/vnd.openxmlformats-"
                "officedocument.spreadsheetml."
                "worksheet+xml"
            ),
        },
    )

    result = dict(original)
    result[path] = xmlbytes(
        build_sheet(
            fisher,
            matched,
            conventional,
            higher,
        )
    )
    result["xl/workbook.xml"] = xmlbytes(
        wb
    )
    result[
        "xl/_rels/workbook.xml.rels"
    ] = xmlbytes(rels)
    result[
        "[Content_Types].xml"
    ] = xmlbytes(content)

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with zipfile.ZipFile(
        args.output,
        "w",
        zipfile.ZIP_DEFLATED,
    ) as archive:
        for name, payload in result.items():
            archive.writestr(
                name,
                payload,
            )

    changed = {
        "xl/workbook.xml",
        "xl/_rels/workbook.xml.rels",
        "[Content_Types].xml",
    }

    with zipfile.ZipFile(
        args.output
    ) as archive:
        assert all(
            archive.read(name) == payload
            for name, payload in original.items()
            if name not in changed
        )

    print("Saved:", args.output)
    print(
        "Verified: all pre-existing worksheet XML "
        "and scientific values are unchanged."
    )
    print(
        "Added:",
        SHEET_NAME,
    )
    print(
        "Higher-abundance FDR-significant clients:",
        ", ".join(
            row["gene_symbol"]
            for row in higher
        ),
    )


if __name__ == "__main__":
    main()
