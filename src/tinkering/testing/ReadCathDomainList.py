
domains = {}
with open("data/cath-domain-list.txt") as f:
    for line in f:
        # skip comments
        if line.startswith("#"):
            continue

        splits = line.split(" ")

        # delete empty strings
        splits = list(filter(None, splits))

        # extract information
        domain = splits[0]
        superfamily = f"{splits[1]}.{splits[2]}.{splits[3]}.{splits[4]}"

        # collect
        domains[domain] = superfamily


i = 500478 / 100000
print(i)