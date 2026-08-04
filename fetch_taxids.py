import sys
import time
import requests

def fetch_taxid(accession):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "nucleotide",
        "term": accession,
        "retmode": "json"
    }
    r = requests.get(url, params=params)
    r.raise_for_status()
    data = r.json()
    idlist = data.get("esearchresult", {}).get("idlist", [])
    if not idlist:
        return None

    # Now get taxid for first id
    uid = idlist[0]

    url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params2 = {
        "db": "nucleotide",
        "id": uid,
        "retmode": "json"
    }
    r2 = requests.get(url2, params=params2)
    r2.raise_for_status()
    data2 = r2.json()
    docsum = data2.get("result", {}).get(uid, {})
    taxid = docsum.get("taxid", None)
    return taxid

def main(input_file, output_file):
    with open(input_file) as f, open(output_file, "w") as out:
        for line in f:
            acc = line.strip()
            if not acc:
                continue
            try:
                taxid = fetch_taxid(acc)
                if taxid is None:
                    print(f"{acc}\tNA")
                else:
                    print(f"{acc}\t{taxid}")
                    out.write(f"{acc}\t{taxid}\n")
            except Exception as e:
                print(f"Error fetching {acc}: {e}", file=sys.stderr)
            time.sleep(0.34)  # NCBI limit ~3 requests per second

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fetch_taxids.py <input_accessions.txt> <output_taxid_map.tsv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
