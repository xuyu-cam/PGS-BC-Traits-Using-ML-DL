def get_trait_abbrs(trait_abbr_file):
    #get the full trait list
    trait_abbr =[]
    with open(trait_abbr_file) as f:
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            trait_abbr.append(line[0])
    return trait_abbr