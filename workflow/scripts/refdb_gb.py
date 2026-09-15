"""Parse GenBank XML (GBSet) into flat records."""
import glob
import os
import xml.etree.ElementTree as ET

COLUMNS = ["accession", "version", "definition", "organism", "lineage", "length", "moltype",
           "topology", "division", "keywords", "create_date", "update_date", "db_source",
           "taxid", "host", "isolation_source", "country", "collection_date", "strain",
           "segment", "note", "n_cds", "products", "pubmed"]


def _quals(feature):
    out = {}
    for q in feature.findall("./GBFeature_quals/GBQualifier"):
        name = q.findtext("GBQualifier_name")
        if name:
            out.setdefault(name, []).append(q.findtext("GBQualifier_value") or "")
    return out


def parse_gbseq(seq):
    feats = seq.findall("./GBSeq_feature-table/GBFeature")
    src, cds = {}, []
    for f in feats:
        key = f.findtext("GBFeature_key")
        if key == "source" and not src:
            src = _quals(f)
        elif key == "CDS":
            cds.append(_quals(f))
    taxid = ""
    for xref in src.get("db_xref", []):
        if xref.startswith("taxon:"):
            taxid = xref.split(":", 1)[1]
    return {
        "accession": seq.findtext("GBSeq_primary-accession") or "",
        "version": seq.findtext("GBSeq_accession-version") or "",
        "definition": seq.findtext("GBSeq_definition") or "",
        "organism": seq.findtext("GBSeq_organism") or "",
        "lineage": seq.findtext("GBSeq_taxonomy") or "",
        "length": int(seq.findtext("GBSeq_length") or 0),
        "moltype": seq.findtext("GBSeq_moltype") or "",
        "topology": seq.findtext("GBSeq_topology") or "",
        "division": seq.findtext("GBSeq_division") or "",
        "keywords": ";".join((e.text or "") for e in seq.findall("./GBSeq_keywords/GBKeyword")),
        "create_date": seq.findtext("GBSeq_create-date") or "",
        "update_date": seq.findtext("GBSeq_update-date") or "",
        "db_source": seq.findtext("GBSeq_source-db") or "",
        "taxid": taxid,
        "host": ";".join(src.get("host", [])),
        "isolation_source": ";".join(src.get("isolation_source", [])),
        "country": ";".join(src.get("geo_loc_name", src.get("country", []))),
        "collection_date": ";".join(src.get("collection_date", [])),
        "strain": ";".join(src.get("strain", []) + src.get("isolate", [])),
        "segment": ";".join(src.get("segment", [])),
        "note": ";".join(src.get("note", [])),
        "n_cds": len(cds),
        "products": " | ".join(sorted({(c.get("product", [""])[0] or "")[:60] for c in cds})),
        "pubmed": ";".join(sorted({e.text for e in seq.findall(".//GBReference_pubmed") if e.text})),
    }


def parse_dir(xmldir, pattern="*.xml"):
    records = []
    for path in sorted(glob.glob(os.path.join(xmldir, pattern))):
        root = ET.parse(path).getroot()
        for seq in root.findall("GBSeq"):
            records.append(parse_gbseq(seq))
    return records


def lineage_rank(lineage, suffix):
    for token in str(lineage).split("; "):
        if token.endswith(suffix):
            return token
    return ""


def add_viral_taxonomy(df):
    df["virus_class"] = df.lineage.apply(lambda s: lineage_rank(s, "viricetes"))
    df["virus_order"] = df.lineage.apply(lambda s: lineage_rank(s, "virales"))
    df["virus_family"] = df.lineage.apply(lambda s: lineage_rank(s, "viridae") or "unclassified")
    df["lineage_terminal"] = df.lineage.str.split("; ").str[-1]
    return df
