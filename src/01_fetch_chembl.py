"""Step 1 - pull raw bioactivity from ChEMBL for InhA and M. tuberculosis whole-cell.

Two targets are pulled:
  CHEMBL1849  Enoyl-[acyl-carrier-protein] reductase [NADH], M. tuberculosis  (the enzyme)
  CHEMBL360   Mycobacterium tuberculosis                                      (whole-cell)

The enzyme target is exhaustible by straightforward pagination. The whole-cell organism
target carries >200k activities and the activity endpoint offers a minimum-potency filter
but no usable offset past a few thousand rows, so the potency axis is partitioned into
pChEMBL bands and each band paginated separately. Records are de-duplicated on
`activity_id` when the bands are merged.

Output: data/chembl_raw.json with keys `meta`, `inha`, `mtb`.

Provenance note: the run behind the committed dataset went through a ChEMBL MCP connector
rather than direct HTTPS, because the analysis sandbox could not reach www.ebi.ac.uk. The
records and field names are identical either way - the connector is a transport, not a
transform - and this script is the portable version that anyone can run.
"""
import json, os, time
import requests

BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "chembl_raw.json")
TARGETS = {"inha": ("CHEMBL1849", "Enoyl-[acyl-carrier-protein] reductase [NADH] / M. tuberculosis"),
           "mtb": ("CHEMBL360", "Mycobacterium tuberculosis (whole-cell)")}
TYPES = ["IC50", "Ki", "Kd", "MIC", "MIC90", "EC50"]
BANDS = [(3, 5), (5, 6), (6, 7), (7, 8), (8, 20)]      # pChEMBL bands
FIELDS = ["activity_id", "molecule_chembl_id", "canonical_smiles", "standard_type",
          "standard_value", "standard_units", "standard_relation", "pchembl_value",
          "assay_chembl_id", "assay_description", "assay_type", "data_validity_comment",
          "document_year", "document_chembl_id", "target_chembl_id", "potential_duplicate",
          "molecule_pref_name"]


def page(params, cap=20000, pause=0.15):
    """Follow ChEMBL's `page_meta.next` cursor until exhausted or `cap` rows."""
    rows, url = [], BASE
    while url and len(rows) < cap:
        r = requests.get(url, params=params if url == BASE else None, timeout=60)
        r.raise_for_status()
        d = r.json()
        rows += d["activities"]
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
        params = None
        time.sleep(pause)
    return rows


def main():
    out = {"meta": {"pulled_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "source": "ChEMBL activity endpoint"}}
    for key, (tid, desc) in TARGETS.items():
        seen = {}
        for t in TYPES:
            for lo, hi in BANDS:
                p = dict(target_chembl_id=tid, standard_type=t, limit=1000,
                         pchembl_value__gte=lo, pchembl_value__lt=hi)
                for a in page(p):
                    seen[a["activity_id"]] = {k: a.get(k) for k in FIELDS}
        out[key] = list(seen.values())
        out["meta"][f"{key}_target"] = f"{tid} / {desc}"
        out["meta"][f"n_{key}"] = len(out[key])
        print(f"{key} ({tid}): {len(out[key])} unique activities")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
