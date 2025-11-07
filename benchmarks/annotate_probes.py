import csv, json, urllib.request, urllib.parse, time
with open(r"C:\Users\felix\downloads\test\probes.csv") as f:
    probes = [row["probe_id"] for row in csv.DictReader(f)]
mapping = {}
for i in range(0, len(probes), 200):
    batch = probes[i:i+200]
    body = urllib.parse.urlencode({"q": ",".join(batch), "scopes": "reporter",
                                    "species": "human", "fields": "symbol",
                                    "size": 200}).encode()
    req = urllib.request.Request("https://mygene.info/v3/query", data=body,
              headers={"Content-Type": "application/x-www-form-urlencoded"})
    with urllib.request.urlopen(req, timeout=30) as r:
        for hit in json.loads(r.read()):
            if not hit.get("notfound") and hit.get("symbol"):
                mapping[hit["query"]] = hit["symbol"]
    time.sleep(0.3)
with open(r"C:\Users\felix\downloads\test\probe_symbols.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["probe_id", "symbol"])
    for p in probes:
        w.writerow([p, mapping.get(p, "")])
print(f"Mapped {len(mapping)}/{len(probes)} probes")
