

def get_trait_vector(pheno_name,pheno_file_path):
    # function that gets the adjusted trait value data via a vector
    col_name_index_map = {}
    pheno_vec = []
    with open(pheno_file_path) as f:
        col_name_line = f.readline().strip().split()
        for i in range(len(col_name_line)):
            col_name_index_map[col_name_line[i]]=i
        target_index = col_name_index_map[pheno_name+ '_gwas_normalised']
        for line in f:
            line = line.strip().split()
            pheno_vec.append(line[target_index])
    return pheno_vec  #still a string list


