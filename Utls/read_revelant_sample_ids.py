def read_related_sampleIDs(file):
    # read the sample ids in INTERVAL that are related to samples in UKB
    sampleIDs = []
    with open(file) as f:
        for line in f:
            sampleIDs.append(int(line.strip()))
    return sampleIDs