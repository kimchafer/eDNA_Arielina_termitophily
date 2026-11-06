#!/usr/bin/env python3
# ==============================================================================
# PURPOSE
#   Compare MMseqs2 taxonomic assignments at 80% and 97% minimum identity.
#   Merge taxonomy by MOTU ID and retain LULU-curated per-sample abundances.
#   Produce the 21-column MOTU_comparison table in TSV and Excel formats.
#
# INPUT
#   6_mmseqs_lca4/results/taxonomy_lca4_a20_80.tsv
#   6_mmseqs_lca4/results/taxonomy_lca4_a20_97.tsv
#   5_lulu/ALL_d2/LULU_curated_MOTU_table.tsv
#
# OUTPUT
#   7_taxonomy_comparison/MOTU_comparison.tsv
#   7_taxonomy_comparison/MOTU_comparison.xlsx
#   7_taxonomy_comparison/logs/parameters.tsv
#   7_taxonomy_comparison/logs/comparison_summary.tsv
#
# SETTINGS
#   Fixed absolute paths below; no command-line arguments required.
#   Preserve the original 21-column layout and comparison rules.
#   Preserve the order of MOTUs in the LULU abundance table.
#   Match by MOTU ID; fail on missing, duplicate or unexpected IDs.
#   Require the four sample columns KPC0541-KPC0544.
#   Do not overwrite an existing output directory.
#   Empty display cells are written as em dashes (—) in TSV/XLSX.
#   Excel output contains plain cell values with default appearance.
#   No colors, borders, bold text, filters, frozen panes, table styles,
#   explicit column widths, or special row heights.
#   Python standard library only.
# ==============================================================================

import csv
import hashlib
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================

WORK_DIR = Path(
    "/storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103"
)

TAX80 = WORK_DIR / (
    "6_mmseqs_lca4/results/taxonomy_lca4_a20_80.tsv"
)

TAX97 = WORK_DIR / (
    "6_mmseqs_lca4/results/taxonomy_lca4_a20_97.tsv"
)

ABUNDANCE = WORK_DIR / (
    "5_lulu/ALL_d2/LULU_curated_MOTU_table.tsv"
)

OUT_DIR = WORK_DIR / "7_taxonomy_comparison"

# ==============================================================================
# OUTPUT TABLE SETTINGS
# ==============================================================================

HEADERS = [
    "MOTU_ID",
    "Total_reads",
    "KPC0541",
    "KPC0542",
    "KPC0543",
    "KPC0544",
    "TaxID_80",
    "Rank_80",
    "Assignment_80",
    "TaxID_97",
    "Rank_97",
    "Assignment_97",
    "Comparison",
    "Insecta_80",
    "Insecta_97",
    "Order_80",
    "Order_97",
    "Family_80",
    "Family_97",
    "Lineage_80",
    "Lineage_97",
]

NUMERIC_COLUMNS = {
    "Total_reads",
    "KPC0541",
    "KPC0542",
    "KPC0543",
    "KPC0544",
    "TaxID_80",
    "TaxID_97",
}

# ==============================================================================
# STEP 1. READ AND VALIDATE INPUT TABLES
# ==============================================================================

def read_taxonomy(path):
    data = {}

    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream, delimiter="\t")

        for lineno, row in enumerate(reader, 1):
            if len(row) != 5:
                raise ValueError(
                    f"{path}:{lineno}: expected 5 columns, "
                    f"found {len(row)}"
                )

            motu, taxid, rank, name, lineage = row

            if not motu or motu in data:
                raise ValueError(
                    f"{path}:{lineno}: empty or duplicate MOTU ID: {motu}"
                )

            try:
                taxid = int(taxid)
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{lineno}: invalid TaxID: {taxid}"
                ) from exc

            if taxid < 0:
                raise ValueError(
                    f"{path}:{lineno}: negative TaxID"
                )

            data[motu] = (taxid, rank, name, lineage)

    return data


def read_abundance(path):
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        samples = HEADERS[2:6]

        if (
            not fields
            or fields[0] != "MOTU_ID"
            or len(fields) != len(set(fields))
        ):
            raise ValueError("Invalid abundance table header")

        if fields[1:] != samples:
            raise ValueError(
                f"Expected sample columns {samples}, "
                f"found {fields[1:]}"
            )

        seen = set()
        rows = []

        for row in reader:
            motu = row["MOTU_ID"]

            if not motu or motu in seen:
                raise ValueError(
                    f"Empty or duplicate MOTU ID: {motu}"
                )

            seen.add(motu)
            counts = []

            for sample in samples:
                try:
                    value = int(row[sample])
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Invalid count: {motu}, {sample}"
                    ) from exc

                if value < 0:
                    raise ValueError(
                        f"Negative count: {motu}, {sample}"
                    )

                counts.append(value)

            if sum(counts) == 0:
                raise ValueError(
                    f"Zero-abundance MOTU: {motu}"
                )

            rows.append((motu, counts))

    if not rows:
        raise ValueError("Empty abundance table")

    return rows

# ==============================================================================
# STEP 2. PARSE TAXONOMIC LINEAGE AND COMPARE ASSIGNMENTS
# ==============================================================================

def lineage_info(lineage):
    tokens = lineage.split(";") if lineage else []

    has_insecta = any(
        token == "c_Insecta"
        for token in tokens
    )

    def first(prefix):
        return next(
            (
                token[2:]
                for token in tokens
                if token.startswith(prefix)
            ),
            "",
        )

    return (
        "Yes" if has_insecta else "No",
        first("o_"),
        first("f_"),
    )


def compare(t80, t97):
    a, b = t80[0], t97[0]

    if a == 0 and b == 0:
        return "Unclassified both"

    if a != 0 and b == 0:
        return "Lost at 97%"

    if a == 0 and b != 0:
        return "Gained at 97%"

    if a == b:
        return "Same TaxID"

    return "Changed TaxID"


def build_rows(abund, tax80, tax97):
    expected = {
        motu for motu, _ in abund
    }

    for label, tax in [
        ("80%", tax80),
        ("97%", tax97),
    ]:
        if set(tax) != expected:
            missing = sorted(expected - set(tax))
            extra = sorted(set(tax) - expected)

            raise ValueError(
                f"{label} MOTU IDs disagree: "
                f"missing={missing[:10]} "
                f"extra={extra[:10]}"
            )

    rows = []

    for motu, counts in abund:
        t80 = tax80[motu]
        t97 = tax97[motu]

        i80, o80, f80 = lineage_info(t80[3])
        i97, o97, f97 = lineage_info(t97[3])

        rows.append([
            motu,
            sum(counts),
            *counts,
            t80[0],
            t80[1],
            t80[2],
            t97[0],
            t97[1],
            t97[2],
            compare(t80, t97),
            i80,
            i97,
            o80,
            o97,
            f80,
            f97,
            t80[3],
            t97[3],
        ])

    return rows

# ==============================================================================
# STEP 3. EXPORT PLAIN TSV AND XLSX TABLES
# ==============================================================================

EM_DASH = "—"


def display_value(value):
    """Replace missing text in exported tables without changing numeric zero."""
    if value is None or (isinstance(value, str) and not value.strip()):
        return EM_DASH
    return value


def write_tsv(path, rows):
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as stream:
        writer = csv.writer(
            stream,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(HEADERS)
        writer.writerows(
            [display_value(value) for value in row] for row in rows
        )


def excel_col(n):
    result = ""

    while n:
        n, remainder = divmod(n - 1, 26)
        result = chr(65 + remainder) + result

    return result


def write_xlsx(path, rows):
    """
    Write a plain Excel worksheet using default cell appearance.

    No:
      - fill colors
      - font styles
      - borders
      - Excel tables
      - automatic filters
      - frozen panes
      - custom row heights
      - custom column widths
    """

    all_rows = [HEADERS] + rows
    xml_rows = []

    for i, values in enumerate(all_rows, 1):
        cells = []

        for j, (header, value) in enumerate(
            zip(HEADERS, values), 1
        ):
            address = f"{excel_col(j)}{i}"

            if i > 1 and header in NUMERIC_COLUMNS:
                cells.append(
                    f'<c r="{address}">'
                    f"<v>{int(value)}</v>"
                    "</c>"
                )
            else:
                cells.append(
                    f'<c r="{address}" t="inlineStr">'
                    f"<is><t>{escape(str(display_value(value)))}</t></is>"
                    "</c>"
                )

        xml_rows.append(
            f'<row r="{i}">'
            + "".join(cells)
            + "</row>"
        )

    ns = (
        "http://schemas.openxmlformats.org/"
        "spreadsheetml/2006/main"
    )

    rel_ns = (
        "http://schemas.openxmlformats.org/"
        "officeDocument/2006/relationships"
    )

    package_rel_ns = (
        "http://schemas.openxmlformats.org/"
        "package/2006/relationships"
    )

    xml_header = (
        '<?xml version="1.0" encoding="UTF-8" '
        'standalone="yes"?>'
    )

    worksheet = (
        xml_header
        + f'<worksheet xmlns="{ns}">'
        + "<sheetData>"
        + "".join(xml_rows)
        + "</sheetData>"
        + "</worksheet>"
    )

    workbook = (
        xml_header
        + f'<workbook xmlns="{ns}" '
        f'xmlns:r="{rel_ns}">'
        '<sheets>'
        '<sheet name="MOTU_comparison" '
        'sheetId="1" r:id="rId1"/>'
        '</sheets>'
        '</workbook>'
    )

    workbook_rels = (
        xml_header
        + f'<Relationships xmlns="{package_rel_ns}">'
        f'<Relationship Id="rId1" '
        f'Type="{rel_ns}/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        '</Relationships>'
    )

    package_rels = (
        xml_header
        + f'<Relationships xmlns="{package_rel_ns}">'
        f'<Relationship Id="rId1" '
        f'Type="{rel_ns}/officeDocument" '
        'Target="xl/workbook.xml"/>'
        '</Relationships>'
    )

    content_types = (
        xml_header
        + '<Types xmlns="http://schemas.openxmlformats.org/'
        'package/2006/content-types">'
        '<Default Extension="rels" '
        'ContentType="application/vnd.openxmlformats-package.'
        'relationships+xml"/>'
        '<Default Extension="xml" '
        'ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.'
        'spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.'
        'spreadsheetml.worksheet+xml"/>'
        '</Types>'
    )

    files = {
        "[Content_Types].xml": content_types,
        "_rels/.rels": package_rels,
        "xl/workbook.xml": workbook,
        "xl/_rels/workbook.xml.rels": workbook_rels,
        "xl/worksheets/sheet1.xml": worksheet,
    }

    with ZipFile(
        path,
        "w",
        compression=ZIP_DEFLATED,
    ) as bundle:
        for filename, content in files.items():
            bundle.writestr(filename, content)

# ==============================================================================
# STEP 4. RECORD PARAMETERS AND VALIDATE OUTPUTS
# ==============================================================================

def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def write_log_tsv(path, pairs):
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(["metric", "value"])
        writer.writerows(pairs)


def main():
    print(
        "[INFO] Start: "
        f"{datetime.now().astimezone().isoformat(timespec='seconds')}"
    )
    print(f"[INFO] Work directory: {WORK_DIR}")
    print("[INFO] Taxonomy conditions: 80% and 97%")

    # --------------------------------------------------------------------------
    # STEP 1. Validate and read all input files.
    # --------------------------------------------------------------------------

    print("[STEP 1] Validate input taxonomy and LULU abundance tables")

    for path in (TAX80, TAX97, ABUNDANCE):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(
                f"Missing or empty input: {path}"
            )

    if OUT_DIR.exists():
        raise FileExistsError(
            f"Output directory already exists: {OUT_DIR}"
        )

    abundance = read_abundance(ABUNDANCE)
    tax80 = read_taxonomy(TAX80)
    tax97 = read_taxonomy(TAX97)

    # --------------------------------------------------------------------------
    # STEP 2. Join taxonomy by MOTU ID.
    # --------------------------------------------------------------------------

    print("[STEP 2] Merge taxonomy by MOTU ID and reconstruct 21 columns")

    rows = build_rows(
        abundance,
        tax80,
        tax97,
    )

    n_motus = len(rows)
    n_reads = sum(row[1] for row in rows)

    changes = Counter(
        row[12] for row in rows
    )

    insecta80 = sum(
        row[13] == "Yes" for row in rows
    )

    insecta97 = sum(
        row[14] == "Yes" for row in rows
    )

    print(
        f"[QC PASS] Input IDs agree: "
        f"{n_motus} MOTUs, {n_reads:,} reads"
    )

    # --------------------------------------------------------------------------
    # STEP 3. Export comparison tables and metadata.
    # --------------------------------------------------------------------------

    print("[STEP 3] Write comparison outputs and QC logs")

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=False,
    )

    LOG_DIR = OUT_DIR / "logs"
    LOG_DIR.mkdir()

    tsv_path = OUT_DIR / "MOTU_comparison.tsv"
    xlsx_path = OUT_DIR / "MOTU_comparison.xlsx"

    write_tsv(tsv_path, rows)
    write_xlsx(xlsx_path, rows)

    write_log_tsv(
        LOG_DIR / "parameters.tsv",
        [
            ("work_dir", str(WORK_DIR)),
            ("taxonomy_80_path", str(TAX80)),
            ("taxonomy_97_path", str(TAX97)),
            ("abundance_path", str(ABUNDANCE)),
            ("taxonomy_80_sha256", sha256(TAX80)),
            ("taxonomy_97_sha256", sha256(TAX97)),
            ("abundance_sha256", sha256(ABUNDANCE)),
            (
                "script_sha256",
                sha256(Path(__file__).resolve()),
            ),
            ("python_version", sys.version.split()[0]),
            ("comparison_key", "MOTU_ID"),
        ],
    )

    write_log_tsv(
        LOG_DIR / "comparison_summary.tsv",
        [
            ("total_MOTUs", n_motus),
            ("total_reads", n_reads),
            ("Insecta_80_MOTUs", insecta80),
            ("Insecta_97_MOTUs", insecta97),
            *[
                (status, changes[status])
                for status in (
                    "Same TaxID",
                    "Lost at 97%",
                    "Gained at 97%",
                    "Changed TaxID",
                    "Unclassified both",
                )
            ],
        ],
    )

    # --------------------------------------------------------------------------
    # STEP 4. Validate output files and report summary.
    # --------------------------------------------------------------------------

    print("[STEP 4] Validate output files and report summary")

    for path in (
        tsv_path,
        xlsx_path,
        LOG_DIR / "parameters.tsv",
        LOG_DIR / "comparison_summary.tsv",
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing or empty output: {path}"
            )

    if sum(changes.values()) != n_motus:
        raise RuntimeError(
            "Comparison status accounting mismatch"
        )

    print(
        f"[QC PASS] Total MOTUs: {n_motus}; "
        f"reads: {n_reads:,}"
    )

    print(
        f"[QC] Insecta: "
        f"80% = {insecta80}; "
        f"97% = {insecta97}"
    )

    for status in (
        "Same TaxID",
        "Lost at 97%",
        "Gained at 97%",
        "Changed TaxID",
        "Unclassified both",
    ):
        print(f"[QC] {status}: {changes[status]}")

    print(f"[OUTPUT] {tsv_path}")
    print(f"[OUTPUT] {xlsx_path}")

    print(
        "[DONE] "
        f"{datetime.now().astimezone().isoformat(timespec='seconds')}"
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)
