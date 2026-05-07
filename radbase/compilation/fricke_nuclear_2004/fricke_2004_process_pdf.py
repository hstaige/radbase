#!/usr/bin/env python3
"""
Extract Fricke et al. muonic atom data tables from a PDF into a pandas DataFrame.

The target tables are the element-wise muonic-atom tables labeled Z.1.3.1
(where Z is the atomic number), e.g. "20.1.3.1" for calcium.

Dependencies:
    pip install pymupdf pandas openpyxl

Usage:
    python extract_fricke_muonic.py Fricke2004.pdf -o fricke_muonic.csv
    python extract_fricke_muonic.py Fricke2004.pdf -o fricke_muonic.xlsx

Notes:
- The script keeps values as strings so uncertainty notation such as 4.4628(11;18)
  and unicode minus signs are not lost. Additional *_value columns can be added by
  applying pandas.to_numeric to cleaned strings if desired.
- The layout parser is coordinate-based. It handles the two main table families in
  the compilation: ordinary transition-energy tables and deformed-nucleus Barrett
  radius parameter tables.
"""

import argparse
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Literal, Optional, Sequence

import fitz  # PyMuPDF
import pandas as pd

DASHES = str.maketrans({
    "\u2013": "-",  # en dash
    "\u2014": "-",  # em dash
    "\u2212": "-",  # minus sign
    "\u00ad": "",   # soft hyphen
})

# A conservative list of atomic symbols for sanity checking element headings.
SYMBOLS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
}

acount=0


@dataclass
class Word:
    x0: float
    y0: float
    x1: float
    y1: float
    text: str


@dataclass
class Line:
    page_index: int
    y0: float
    y1: float
    words: list[Word]
    text: str
    abs_index: int = -1


@dataclass
class Section:
    z: int
    label: str
    page_start: int
    abs_start: int
    abs_end: int
    title: str
    symbol: Optional[str] = None
    element: Optional[str] = None


@dataclass
class Layout:
    kind: Literal['energy', 'deformed']
    columns: list[str]
    x_positions: dict[str, float]
    boundaries: dict[str, tuple[float, float]]
    header_abs_index: int
    title: str


def norm_text(s: str) -> str:
    return re.sub(r"\s+", " ", s.translate(DASHES)).strip()


def is_number_like(s: str) -> bool:
    """True for values such as 1255.849(63), -2.853, 0.1, 3.1(2;4)."""
    s = norm_text(s).replace(" ", "")
    return bool(re.fullmatch(r"[-+]?\d+(?:\.\d+)?(?:\([^)]*\))?", s))


def get_lines(page: fitz.Page, page_index: int, y_tol: float = 2.0) -> list[Line]:
    """Extract page words and group them into horizontal text lines."""
    raw = page.get_text("words")  # x0, y0, x1, y1, word, block, line, word_no
    words = [Word(x0, y0, x1, y1, t, ) for x0, y0, x1, y1, t, *_ in raw]
    words.sort(key=lambda w: ((w.y0 + w.y1) / 2, w.x0))

    groups: list[list[Word]] = []
    group_centers: list[float] = []
    for w in words:
        yc = (w.y0 + w.y1) / 2
        if groups and abs(group_centers[-1] - yc) <= y_tol:
            groups[-1].append(w)
            group_centers[-1] = sum((x.y0 + x.y1) / 2 for x in groups[-1]) / len(groups[-1])
        else:
            groups.append([w])
            group_centers.append(yc)

    lines: list[Line] = []
    for grp in groups:
        grp.sort(key=lambda w: w.x0)
        text = norm_text(" ".join(w.text for w in grp))
        if not text:
            continue
        lines.append(Line(
            page_index=page_index,
            y0=min(w.y0 for w in grp),
            y1=max(w.y1 for w in grp),
            words=grp,
            text=text,
        ))
    return lines


def all_lines(doc: fitz.Document) -> list[Line]:
    lines: list[Line] = []
    for i, page in enumerate(doc):
        for line in get_lines(page, i):
            line.abs_index = len(lines)
            lines.append(line)
    return lines


def page_element_from_lines(page_lines: Sequence[Line]) -> tuple[Optional[int], Optional[str], Optional[str]]:
    """Find headings such as '20-Ca Calcium' on a page."""
    for line in page_lines[:20]:
        m = re.search(r"\b(?P<z>\d{1,3})-(?P<sym>[A-Z][a-z]?)\s+(?P<name>[A-Z][A-Za-z-]+)\b", line.text)
        if m and m.group("sym") in SYMBOLS:
            return int(m.group("z")), m.group("sym"), m.group("name")
    return None, None, None


def find_muonic_sections(lines: Sequence[Line]) -> list[Section]:
    """Find section headers that start with Z.1.3.1."""
    starts: list[tuple[int, int, str]] = []
    for line in lines:
        m = re.match(r"^\s*(?P<label>(?P<z>\d{1,3})\.1\.3\.1)\b\s*(?P<title>.*)$", line.text)
        if m:
            starts.append((line.abs_index, int(m.group("z")), norm_text(line.text)))

    sections: list[Section] = []
    # Collect lines by page
    by_page: DefaultDict[int, list[Line]] = defaultdict(list)
    for line in lines:
        by_page[line.page_index].append(line)

    for n, (start_abs, z, title) in enumerate(starts):

        # These sections have highly irregular structures due to the lack of data. Easier to exclude them.
        if z in (65, 67, 71):
            continue

        label = f"{z}.1.3.1"
        # End at the first next heading in the same element that leaves subsection .1,
        # or at the next Z.1.3.1 start, whichever comes first. This prevents pulling in
        # the .1.3.2 differences table.
        next_start_abs = starts[n + 1][0] if n + 1 < len(starts) else len(lines)
        end_abs = min(len(lines), next_start_abs)
        stop_re = re.compile(rf"^\s*{z}\.(?:1\.3\.[2-9]|1\.[4-9]|[23])\b")
        for j in range(start_abs + 1, end_abs):
            if stop_re.match(lines[j].text):
                end_abs = j
                break

        page_start = lines[start_abs].page_index
        _, sym, name = page_element_from_lines(by_page[page_start])
        sections.append(Section(z=z, label=label, page_start=page_start, abs_start=start_abs, abs_end=end_abs, title=title, symbol=sym, element=name))
    return sections


def infer_layout(header_line: Line, header_window: Sequence[Line], section: Section) -> Layout:
    """Infer table column x-positions from the header line and nearby unit lines."""
    header_text = norm_text(" ".join(line.text for line in header_window))
    kind = "deformed" if ("R0" in header_text or "β" in header_text or "deformed" in section.title) and "Eexp" not in header_text else "energy"

    x: dict[str, float] = {}
    words: list[Word] = [word for line in header_window for word in line.words]
    possible_columns = ("A", "Isotope", "Eexp", "Etheo", "Npol", "c", "R0", "a", "β", "β2", "<", "α", "k", "R", "C", "Cz", "Ref")
    aliases = {"Isotope": "A", "β": "beta2", "β2": "beta2", "<": "r2_model", "α": "alpha", "C": "Cz", "R": "Rkalpha"}

    # Use the actual header line first; the nearby lines help find unit-only columns.
    for w in words:
        t = norm_text(w.text).strip(".,")
        if t in possible_columns:
            x.setdefault(aliases.get(t, t), w.x0)

    if kind == "energy":
        # Decide whether an Etheo and c column exists in this particular table.
        if "Etheo" in x:
            columns = ["A", "Eexp", "Etheo", "Npol", "c", "r2_model", "alpha", "k", "Cz", "Rkalpha", "Ref"]
            template = {"A": 85, "Eexp": 106, "Etheo": 162, "Npol": 205, "c": 234, "r2_model": 278, "alpha": 320,
                        "k": 352, "Cz": 385, "Rkalpha": 434, "Ref": 487}
        else:
            columns = ["A", "Eexp", "Npol", "c", "r2_model", "alpha", "k", "Cz", "Rkalpha", "Ref"]
            template = {"A": 85, "Eexp": 127, "Npol": 187, "r2_model": 217, "alpha": 273, "k": 308, "Cz": 351,
                        "Rkalpha": 406, "Ref": 464}
            # Some light-nucleus tables do not have c; remove it if it was not found.
            if "c" not in x:
                columns.remove("c")
    else:
        columns = ["A", "Npol", "R0", "a", "beta2", "r2_model", "alpha", "k", "Cz", "Rkalpha", "Ref"]
        template = {"A": 85, "Npol": 110, "R0": 157, "a": 199, "beta2": 233, "r2_model": 268, "alpha": 313, "k": 346, "Cz": 378, "Rkalpha": 429, "Ref": 479}

    # Boundaries are halfway between sorted column starts.
    ordered = sorted([(col, x[col]) for col in columns], key=lambda kv: kv[1])
    boundaries: dict[str, tuple[float, float]] = {}
    for idx, (col, xpos) in enumerate(ordered):
        left = -10_000 if idx == 0 else (ordered[idx - 1][1] + xpos) / 2
        right = 10_000 if idx == len(ordered) - 1 else (xpos + ordered[idx + 1][1]) / 2
        boundaries[col] = (left, right)

    return Layout(kind=kind, columns=columns, x_positions=x, boundaries=boundaries, header_abs_index=header_line.abs_index, title=section.title)


def line_to_cells(line: Line, layout: Layout) -> dict[str, str]:
    # Assign words to the corresponding columns based on layout boundaries.
    cells: dict[str, list[str]] = {c: [] for c in layout.columns}
    for w in sorted(line.words, key=lambda z: z.x0):  # Want first text left-to-right that fits in column.
        for col, (left, right) in layout.boundaries.items():
            if left <= w.x0 < right:
                cells[col].append(w.text)
                break
    return {col: norm_text(" ".join(parts)) for col, parts in cells.items()}


def find_table_headers(segment_lines: Sequence[Line]) -> list[int]:
    """Return indexes of table header rows in segment_lines."""
    idxs: list[int] = []
    for i, line in enumerate(segment_lines):
        t = line.text
        # Main data-table header: must contain A or Isotope and Ref.; avoid K X-ray A-A' tables.
        if re.search(r"\b(A|Isotope)\b", t) and "Ref" in t and "A - A" not in t and "A – A" not in t:
            if any(key in t for key in ("Eexp", "Npol", "R0", "Cz", "c [fm]", "R µ", "R k")):
                idxs.append(i)
    return idxs


def looks_like_stop_line(line: Line) -> bool:
    # Checks if we hit the end of the page or a footnote.
    t = line.text
    if t.startswith("Landolt-Börnstein") or t.startswith("New Series"):  # End of page
        return True
    if re.match(r"^\s*\d+\s*\)", t):  # Footnote, e.g. 1)
        return True
    return False


def is_dataish(cells: dict[str, str], layout: Layout) -> bool:
    if cells.get("A") and re.fullmatch(r"\d{1,3}|nat[A-Za-z]+", cells["A"]):
        return True
    if layout.kind == "energy":
        return bool(cells.get("Eexp") or cells.get("Etheo") or cells.get("Npol") or cells.get("c") or cells.get("r2_model") or cells.get("Rkalpha") or cells.get("alpha"))
    else:
        return bool(cells.get("Npol") or cells.get("R0") or cells.get("a") or cells.get("beta2") or cells.get("r2_model") or cells.get("Rkalpha") or cells.get("alpha"))


def has_real_data_value(cells: dict[str, str]) -> bool:
    # Keep closing parentheses on new lines, happens for Z=81
    if cells.get("Eexp", "") and cells['Eexp'].strip() == ')':
        return True
    if has_a(cells):
        return True

    def cell_starts_with_number(s: str) -> bool:
        s = norm_text(s)
        return bool(re.match(r"^[-+]?\d+(?:\.\d+)?(?:\([^)]*\))?\*?$", s) or re.match(r"^[-+]?\d+(?:\.\d+)?", s))

    data_cols = ["Eexp", "Etheo", "Npol", "c", "R0", "a", "beta2", "r2_model", "alpha", "k", "Cz", "Rkalpha"]
    return any(cell_starts_with_number(cells.get(col, "")) for col in data_cols)


def transition_name(section_title: str, index: int, kind: str) -> str:
    t = norm_text(section_title)
    if kind == "deformed":
        if index == 1:
            return "2p1/2-1s1/2 c.g./component"
        if index == 2:
            return "2p3/2-1s1/2 c.g./component"
        return f"muonic component {index}"
    if "2p1/2" in t and "2p3/2" in t:
        return "2p1/2-1s1/2" if index == 1 else "2p3/2-1s1/2" if index == 2 else f"transition {index}"
    if "2p3/2" in t:
        return "2p3/2-1s1/2"
    if "2p1/2" in t:
        return "2p1/2-1s1/2"
    if "2p" in t and "1s" in t:
        return "2p-1s"
    return "muonic atom"


def has_unclosed_parenthesis(s: str) -> bool:
    return s.count("(") > s.count(")")


def has_a(cells: dict[str, str]) -> bool:
    return bool(cells.get("A") and re.fullmatch(r"\d{1,3}|nat[A-Za-z]+", cells.get("A", "")))


def is_pure_parameter_row(row: dict[str, str], layout: Layout) -> bool:
    if layout.kind == "energy":
        return bool(row.get("c") or row.get("r2_model")) and not any(row.get(c) for c in ["Eexp", "Etheo", "Npol", "alpha", "k", "Cz", "Rkalpha", "Ref"])
    return bool(row.get("R0") or row.get("a") or row.get("beta2") or row.get("r2_model")) and not any(row.get(c) for c in ["Npol", "alpha", "k", "Cz", "Rkalpha", "Ref"])


def build_energy_transition_rows(cluster: list[dict[str, str]]) -> list[dict[str, str]]:
    """Merge split Eexp fragments, then return physical transition rows.

    In a few heavy-element tables the value printed in Eexp is split across
    three PDF text lines, for example '5726.140(670' / main row / ')'.
    This function attaches those fragments to the row that contains Etheo/Rkalpha.
    """
    transition_rows: list[dict[str, str]] = []
    prefix_eexp: list[str] = []
    current_main: Optional[dict[str, str]] = None
    non_eexp_cols = ["Etheo", "Npol", "alpha", "k", "Cz", "Rkalpha", "Ref"]

    for row in cluster:
        has_only_eexp = bool(row.get("Eexp")) and not any(row.get(c) for c in non_eexp_cols + ["c", "r2_model"])
        if has_only_eexp:
            if current_main is not None and has_unclosed_parenthesis(current_main.get("Eexp", "")):
                current_main["Eexp"] = norm_text(current_main.get("Eexp", "") + row.get("Eexp", ""))
            else:
                prefix_eexp.append(row["Eexp"])
            continue

        is_main = bool(any(row.get(c) for c in non_eexp_cols) or row.get("Eexp")) and not is_pure_parameter_row(row, Layout("energy", [], {}, {}, -1, ""))
        if is_main:
            if prefix_eexp:
                parts = prefix_eexp + ([row.get("Eexp", "")] if row.get("Eexp") else [])
                row["Eexp"] = norm_text("".join(parts))
                prefix_eexp.clear()
            transition_rows.append(row)
            current_main = row
    return transition_rows


def finalize_cluster(cluster: list[dict[str, str]], section: Section, layout: Layout, page_index: int) -> list[dict[str, str | int]]:
    if not cluster:
        return []
    A = next((row.get("A", "") for row in cluster if row.get("A")), "")
    if not A or not re.fullmatch(r"\d{1,3}|nat[A-Za-z]+", A):
        return []

    # Values that apply to all transition rows in the isotope cluster.
    shared: dict[str, str] = {}
    shared_cols = ["c", "r2_model"] if layout.kind == "energy" else ["R0", "a", "beta2", "r2_model"]
    for col in shared_cols:
        for row in cluster:
            if row.get(col):
                shared[col] = row[col]
                break

    first_ref = next((row.get("Ref", "") for row in cluster if row.get("Ref")), "")

    out: list[dict[str, str | int]] = []
    transition_rows: list[dict[str, str]] = []
    if layout.kind == "energy":
        transition_rows = build_energy_transition_rows(cluster)
    else:
        for row in cluster:
            if row.get("Rkalpha") or row.get("alpha") or row.get("Npol"):
                # Avoid treating pure parameter rows as transition rows.
                if is_pure_parameter_row(row, layout):
                    continue
                transition_rows.append(row)

    for idx, row in enumerate(transition_rows, start=1):
        rec: dict[str, str | int] = {
            "Z": section.z,
            "element_symbol": section.symbol or "",
            "element_name": section.element or "",
            "section_label": section.label,
            "section_title": section.title,
            "table_kind": layout.kind,
            "page": page_index + 1,
            "A": A,
            "transition_index": idx,
            "transition": transition_name(section.title, idx, layout.kind),
        }
        for col in layout.columns:
            if col == "A":
                continue
            value = row.get(col, "") or shared.get(col, "")
            if col == "Ref" and not value:
                value = first_ref
            rec[col] = value
        # Ensure all common columns exist, regardless of table family.
        for col in ["Eexp", "Etheo", "Npol", "c", "R0", "a", "beta2", "r2_model", "alpha", "k", "Cz", "Rkalpha", "Ref"]:
            rec.setdefault(col, "")
        out.append(rec)
    return out


def parse_table_after_header(segment_lines: Sequence[Line], header_i: int, section: Section) -> list[dict[str, str]]:
    # Header windows are the header row plus a few unit/subscript rows below it.
    header_line = segment_lines[header_i]
    header_window = segment_lines[header_i: min(len(segment_lines), header_i + 5)]
    layout = infer_layout(header_line, header_window, section)

    # Convert text lines to candidate data rows first. This makes it easier to handle
    # heavy-element cases where Eexp starts one PDF line above the line containing A.
    candidates: list[tuple[Line, dict[str, str]]] = []
    for line in segment_lines[header_i + 1:]:
        if looks_like_stop_line(line):
            break
        # Ignore repeated header/unit lines if a table continues on another page.
        if re.search(r"\b(A|Isotope)\b", line.text) and "Ref" in line.text and any(key in line.text for key in ("Eexp", "Npol", "R0")):
            continue
        # Matches [ ] to find units
        if re.fullmatch(r"(?:\[[^]]+\]\s*)+", line.text):
            continue
        cells = line_to_cells(line, layout)
        # Exclude rows where obviously no data is present.
        if not is_dataish(cells, layout):
            continue
        # Exclude rows that don't have numbers or closing parenthesis for at least one column.
        if not has_real_data_value(cells):
            continue
        joined = " ".join(cells.values())
        if "[keV]" in joined or "[fm]" in joined or "[1/fm]" in joined or "[10-" in joined:
            continue
        candidates.append((line, cells))

    def only_open_eexp(cells: dict[str, str]) -> bool:
        # Check if this line is just Eexp with unbalanced parenthesis (e.g. Z=81, '5.906(670' ).
        if not (layout.kind == "energy" and cells.get("Eexp") and has_unclosed_parenthesis(cells["Eexp"])):
            return False
        return not any(cells.get(c) for c in ["Etheo", "Npol", "c", "r2_model", "alpha", "k", "Cz", "Rkalpha", "Ref"])

    rows: list[dict[str, str]] = []
    cluster: list[dict[str, str]] = []
    pending_before_A: list[dict[str, str]] = []
    page_index = header_line.page_index

    # Data can appear across multiple lines, so we collect data into clusters. Once we think a new cluster is starting,
    # we consolidate our current cluster and start collecting the new one. The pending_before_A syntax helps catch odd
    # elements that appear across multiple lines.
    for idx, (line, cells) in enumerate(candidates):
        if has_a(cells):  # New isotope is starting, finish current cluster if any and then start new cluster.
            if cluster:
                rows += finalize_cluster(cluster, section, layout, page_index)
            cluster = pending_before_A + [cells]
            pending_before_A = []
            page_index = line.page_index
            continue

        if cluster:
            # Specifically to handle Z=81, where energies appear on multiple lines and may appear before the new A.
            if only_open_eexp(cells) and idx + 1 < len(candidates) and has_a(candidates[idx + 1][1]):
                rows += finalize_cluster(cluster, section, layout, page_index)
                cluster = []
                pending_before_A = [cells]
                page_index = line.page_index
            else:
                cluster.append(cells)
        else:
            pending_before_A.append(cells)

    if cluster:
        rows += finalize_cluster(cluster, section, layout, page_index)
    return rows


def extract_muonic_tables(pdf_path: str | Path) -> pd.DataFrame:
    doc = fitz.open(str(pdf_path))
    lines = all_lines(doc)
    sections = find_muonic_sections(lines)

    all_records: list[dict[str, str]] = []
    for section in sections:
        segment = lines[section.abs_start:section.abs_end]
        header_idxs = find_table_headers(segment)
        if not header_idxs:  # Z=77 data contains only rms radii.
            continue
        for header_i in header_idxs:  # Z=48 spans multiple pages, so we parse after each header.
            all_records += parse_table_after_header(segment, header_i, section)

    df = pd.DataFrame(all_records)
    if df.empty:
        return df

    # Useful ordering and de-duplication for occasional duplicate header detections.
    columns = [
        "Z", "element_symbol", "element_name", "section_label", "section_title", "table_kind", "page", "A",
        "transition_index", "transition", "Eexp", "Etheo", "Npol", "c", "R0", "a", "beta2", "r2_model",
        "alpha", "k", "Cz", "Rkalpha", "Ref",
    ]
    df = df[[c for c in columns if c in df.columns] + [c for c in df.columns if c not in columns]]
    df = df.sort_values(["Z", "A", "transition_index"], kind="stable").reset_index(drop=True)
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Fricke Z.1.3.1 muonic atom tables from a PDF.")
    parser.add_argument("pdf", help="Input Fricke PDF")
    parser.add_argument("-o", "--output", default="fricke_muonic_tables.csv", help="CSV/XLSX output path")
    args = parser.parse_args()

    df = extract_muonic_tables(args.pdf)

    out = Path(args.output)
    if out.suffix.lower() in {".xlsx", ".xlsm"}:
        df.to_excel(out, index=False)
    else:
        df.to_csv(out, index=False)
    print(f"Extracted {len(df):,} transition rows from {df['section_label'].nunique() if not df.empty else 0} Z.1.3.1 tables -> {out}")


if __name__ == "__main__":
    main()
