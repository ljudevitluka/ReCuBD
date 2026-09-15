"""Readers for the three taxon-specific input files, plus table I/O helpers."""
import re

import pandas as pd

INT_COLS = {"length", "n_cds"}
BOOL_COLS = {"included", "refseq", "is_phage", "from_host_taxon", "placeholder_taxon",
             "def_complete", "def_cds_complete"}


def read_config_tsv(path, required):
    """Read a tab-separated config file, ignoring '#' comment lines."""
    rows = []
    header = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = [f.strip() for f in fields]
                missing = [c for c in required if c not in header]
                if missing:
                    raise ValueError(path + ": missing required column(s) " + ", ".join(missing))
                continue
            fields += [""] * (len(header) - len(fields))
            rows.append({h: fields[i].strip() for i, h in enumerate(header)})
    if not rows:
        raise ValueError(path + ": no data rows found")
    return rows


def load_search_terms(path):
    """-> {'query': [(term, field)], 'name': [term]}"""
    rows = read_config_tsv(path, ["term", "use"])
    query, name = [], []
    for r in rows:
        use = (r.get("use") or "both").lower()
        if use in ("query", "both"):
            query.append((r["term"], r.get("field") or "All Fields"))
        if use in ("name", "both"):
            name.append(r["term"])
    if not query:
        raise ValueError(path + ": no terms with use=query or use=both")
    return {"query": query, "name": name}


def load_host_accept(path):
    """-> {'clades': [...], 'vernacular': [...]}"""
    rows = read_config_tsv(path, ["kind", "value"])
    out = {"clades": [], "vernacular": []}
    for r in rows:
        kind = r["kind"].lower()
        if kind == "clade":
            out["clades"].append(r["value"])
        elif kind == "vernacular":
            out["vernacular"].append(r["value"])
        else:
            raise ValueError(path + ": unknown kind '" + r["kind"] + "' (expected clade or vernacular)")
    if not out["clades"] and not out["vernacular"]:
        raise ValueError(path + ": no clade or vernacular rows")
    return out


def load_tier2_groups(path):
    """-> {'environmental': set(), 'phage': set()}"""
    rows = read_config_tsv(path, ["group", "kind"])
    out = {"environmental": set(), "phage": set()}
    for r in rows:
        kind = r["kind"].lower()
        if kind not in out:
            raise ValueError(path + ": unknown kind '" + r["kind"] + "' (expected environmental or phage)")
        out[kind].add(r["group"])
    return out


def term_regex(terms):
    """Word-boundary, case-insensitive regex over a list of terms."""
    if not terms:
        return re.compile(r"(?!x)x")
    parts = [re.escape(t.strip()) for t in terms if t.strip()]
    return re.compile(r"\b(?:" + "|".join(parts) + r")", re.I)


def build_entrez_query(search_terms, virus_taxid_root):
    clauses = []
    for term, field in search_terms["query"]:
        clauses.append('"' + term + '"[' + field + ']' if field else '"' + term + '"')
    return "txid" + str(virus_taxid_root) + "[Organism:exp] AND (" + " OR ".join(clauses) + ")"


def write_table(df, path):
    df.to_csv(path, sep="\t", index=False)


def read_table(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[])
    for c in df.columns:
        if c in INT_COLS:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
        elif c in BOOL_COLS:
            df[c] = df[c].astype(str).str.lower().isin(["true", "1", "yes"])
    return df
