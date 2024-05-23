import requests
import time
from tqdm import tqdm
import subprocess
import pandas as pd

#############################
# Read in domain descriptions

organism = "Homo sapiens"
organism_domains = {}
with open("data/cath-domain-description-file.txt") as f:

    domain = None
    pdb_id = None
    source = None
    for line in f:
        # skip comments
        if line.startswith("#"):
            continue

        if line.startswith("DOMAIN"):
            splits = line.split(" ")
            splits = list(filter(None, splits))  # delete empty strings
            domain = splits[1]
            print(domain)
            pdb_id = domain[0:4]

        if line.startswith("SOURCE"):
            splits = line.split(" ")
            splits = list(filter(None, splits))  # delete empty strings
            source = splits[1].split(".")[0]
            print(source)

        if domain is not None and source is not None:
            if source == organism:
                organism_domains[pdb_id] = organism

            domain = None
            pdb_id = None
            source = None

            while not line.startswith("//"):
                line = f.readline()
                print(line)

            break


print(organism_domains)



#############################
# Read in domains list

domains = {}
with open("data/cath-domain-list.txt") as f:
    i = 0
    for line in f:
        # skip comments
        if line.startswith("#"):
            continue

        splits = line.split(" ")

        # delete empty strings
        splits = list(filter(None, splits))

        # extract information
        pdb_id = splits[0][0:4]
        chain = splits[0][4:]
        superfamily = f"{splits[1]}.{splits[2]}.{splits[3]}.{splits[4]}"

        # collect
        if pdb_id in domains.keys():
            domains[pdb_id]["chains"][chain] = superfamily
        else:
            i += 1
            domains[pdb_id] = {"chains": {chain: superfamily}, "mapped_ids": [], "failed": False,
                               "organism": [], "gene_name": []}

        if i == 1026:
            break

#############################


def chunks_size_n(data, n):
    for i in range(0, len(data), n):
        yield data[i:i + n]


PDB_mapping_url = "https://rest.uniprot.org/idmapping"

# Post a request
ids = list(domains.keys())
n = 5
r = None

results = {"results": [], "failedIds": [], "obsoleteCount": 0}
for ids in tqdm(chunks_size_n(ids, n)):
    r = requests.post(PDB_mapping_url + "/run", params={"ids": ids, "from": "PDB", "to": "UniProtKB"})
    jobId = r.json()["jobId"]

    # Fetch result
    unfinished = True
    while unfinished:
        r = subprocess.run(["curl", "-i", PDB_mapping_url + f"/status/{jobId}"], text=True, capture_output=True, check=True)
        time.sleep(0)
        if r.stdout[-10:-2] == "FINISHED":
            unfinished = False

    r = requests.get(PDB_mapping_url + f"/status/{jobId}")
    result = r.json()
    for i in range(len(result["results"])):
        organism = result["results"][i]["to"]["organism"]["scientificName"]

        gene_name = "-"
        if "genes" in result["results"][i]["to"].keys():
            if "geneName" in result["results"][i]["to"]["genes"][0].keys():
                gene_name = result["results"][i]["to"]["genes"][0]["geneName"]["value"]

        # fill gene name and organism in domains
        domains[result["results"][i]["from"]]["organism"].append(organism)
        domains[result["results"][i]["from"]]["gene_name"].append(gene_name)

    r = requests.get(PDB_mapping_url + f"/results/{jobId}")
    result = r.json()

    results["results"].extend(result["results"])
    results["obsoleteCount"] += result["obsoleteCount"]
    if "failedIds" in result.keys():
        results["failedIds"].extend(result["failedIds"])

# More resulting IDs since some are mapped to multiple IDs

# Process results
for result in results["results"]:
    domains[result["from"]]["mapped_ids"].append(result["to"])

for failedID in results["failedIds"]:
    domains[failedID]["failed"] = True

df = pd.DataFrame.from_dict(domains, orient="index")
df.index.name = "PDB_ID"
df.to_csv(f"output/domains_{n}_chunks", sep="\t")

