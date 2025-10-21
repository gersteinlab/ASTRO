#!/usr/bin/env python3
import csv
import requests
from io import StringIO

INFILE = "data/GSE151334set.humanmat.txt"
OUTFILE = "cache/GSE151334set.humanmat_with_srr.txt"

def srx_to_srr(srx: str):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "sra", "term": srx, "retmode": "json"}
    r = requests.get(url, params=params, timeout=30).json()
    idlist = r.get("esearchresult", {}).get("idlist", [])
    if not idlist:
        return []
    uid = idlist[0]
    url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params2 = {"db": "sra", "id": uid, "rettype": "runinfo", "retmode": "text"}
    r2 = requests.get(url2, params=params2, timeout=60)
    reader = csv.DictReader(StringIO(r2.text))
    return [row["Run"] for row in reader if row["Run"].startswith("SRR")]

with open(INFILE, newline="") as fin, open(OUTFILE, "w", newline="") as fout:
    reader = csv.DictReader(fin, delimiter="\t")
    fieldnames = reader.fieldnames + ["SRR Accession"]
    writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for row in reader:
        srx = row.get("SRA Accession", "").strip()
        if srx:
            srrs = srx_to_srr(srx)
            row["SRR Accession"] = ",".join(srrs) if srrs else "NA"
        else:
            row["SRR Accession"] = "NA"
        writer.writerow(row)

print('end')
