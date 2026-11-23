"""Serialize complete standard MARCS atmosphere records as .mod or .mod.gz."""

import argparse
from dataclasses import dataclass
import gzip
import math
from pathlib import Path


# The column order follows the public MARCS readmarcs.f reader. The abundance
# block is fixed-width Fortran 10F7.2. Depth rows are read list-directed.
STRUCTURE = (
    " k lgTauR  lgTau5    Depth     T        Pe          Pg         Prad       Pturb",
    (".2f", ".4f", ".3E", ".1f", ".4E", ".4E", ".4E", ".4E"),
)
THERMODYNAMICS = (
    " k lgTauR    KappaRoss   Density   Mu      Vconv   Fconv/F      RHOX",
    (".2f", ".4E", ".4E", ".3f", ".3E", ".5f", ".6E"),
)
PARTIAL_PRESSURES = (
    " k lgPgas   H I    H-     H2     H2+    H2O    OH     CH     CO     CN     C2",
    " k    N2     O2     NO     NH     TiO   C2H2    HCN    C2H    HS     SiH    C3H",
    " k    C3     CS     SiC   SiC2    NS     SiN    SiO    SO     S2     SiS   Other",
)
TABLES = (STRUCTURE, THERMODYNAMICS) + tuple(
    (header, (".3f",) + (".2f",) * 10 if i == 0 else (".2f",) * 11)
    for i, header in enumerate(PARTIAL_PRESSURES)
)
LEADING_LINES = ("Model structure", None, "Assorted logarithmic partial pressures", None, None)


@dataclass(frozen=True)
class MarcsDocument:
    """Complete MARCS model.

    The header holds the eleven fixed-format records from the identifier
    through the mass fractions. Depth arrays are ordered outermost first
    and omit the row index. See the module's table headers for column order.
    Provide all metadata and thermodynamic/chemical quantities for the same
    model when constructing a document.
    """

    header: tuple
    abundances: tuple
    structure: tuple
    thermodynamics: tuple
    partial_pressures: tuple


def _open(path, mode):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, mode, encoding="ascii", newline="")


def _read_table(lines, count, columns):
    header = next(lines, "").split()
    if not header or header[0] != "k":
        raise ValueError("Missing MARCS depth-table header")
    rows = []
    for index in range(1, count + 1):
        fields = next(lines, "").split()
        if len(fields) != columns + 1 or fields[0] != str(index):
            raise ValueError(f"Invalid MARCS depth row {index}")
        rows.append(tuple(map(float, fields[1:])))
    return tuple(rows)


def read_document(path):
    """Read a complete standard MARCS model suitable for re-export."""
    with _open(Path(path), "rt") as stream:
        lines = iter(stream.read().splitlines())

    header = tuple(next(lines, "") for _ in range(11))
    if not next(lines, "").startswith("Logarithmic chemical number abundances"):
        raise ValueError("Missing MARCS chemical abundance block")
    abundances = tuple(float(x) for _ in range(10) for x in next(lines, "").split())
    count_record = next(lines, "")
    if "Number of depth points" not in count_record:
        raise ValueError("Missing MARCS depth count")
    count = int(count_record.split()[0])

    blocks = []
    for title, (_, formats) in zip(LEADING_LINES, TABLES):
        if title and next(lines, "").strip() != title:
            raise ValueError(f"Missing MARCS section {title}")
        blocks.append(_read_table(lines, count, len(formats)))
    if any(line.strip() for line in lines):
        raise ValueError("Unexpected content after the MARCS partial-pressure tables")
    return MarcsDocument(header, abundances, blocks[0], blocks[1], tuple(blocks[2:]))


def write_document(model, path):
    """Write a complete MARCS model to a new plain-text or gzip file."""
    if len(model.header) != 11 or any(
        not line or "\n" in line or "\r" in line for line in model.header
    ):
        raise ValueError("Provide eleven complete MARCS header records")
    if len(model.abundances) != 92 or model.abundances[0] != 12 or not all(
        math.isfinite(value) and len(f"{value:7.2f}") == 7 for value in model.abundances
    ):
        raise ValueError("Provide 92 fixed-width logarithmic abundances with hydrogen set to 12")

    count = len(model.structure)
    blocks = (model.structure, model.thermodynamics, *model.partial_pressures)
    if len(blocks) != len(TABLES) or not 2 <= count <= 9999:
        raise ValueError("Provide both structure tables, three partial-pressure tables, and 2–9999 layers")
    records = list(model.header)
    records.append("Logarithmic chemical number abundances, H always 12.00")
    records.extend(
        "".join(f"{value:7.2f}" for value in model.abundances[i:i + 10])
        for i in range(0, 92, 10)
    )
    records.append(f"{count:4d} Number of depth points")
    for title, (header, formats), rows in zip(LEADING_LINES, TABLES, blocks):
        if len(rows) != count:
            raise ValueError("MARCS tables must all share the same depth count")
        if title:
            records.append(title)
        records.append(header)
        for index, row in enumerate(rows, start=1):
            if len(row) != len(formats) or not all(math.isfinite(value) for value in row):
                raise ValueError(f"Invalid {header.strip()} row {index}")
            fields = (format(value, spec) for value, spec in zip(row, formats))
            records.append(f"{index:4d} " + " ".join(f"{field:>11}" for field in fields))

    output = "\n".join(records) + "\n"
    output.encode("ascii")
    with _open(Path(path), "xt") as stream:
        stream.write(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Complete MARCS .mod or .mod.gz model")
    parser.add_argument("output", help="New .mod or .mod.gz file")
    args = parser.parse_args()
    try:
        write_document(read_document(args.input), args.output)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
