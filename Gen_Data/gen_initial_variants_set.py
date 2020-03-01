

import os,sys

def get_top_P_val_variants(trait,cond_variants_files_path,top_1M_variants_file_path,N):
    variants = []
    with open(cond_variants_files_path+trait+'_condsig') as f:
        for line in f:
            line = line.strip()
            variants.append(line)
    with open(top_1M_variants_file_path + trait + '_gwas_1M_SNPs_full.out') as f:
        for line in f:
            line = line.strip().split('\t')
            rsid = line[-1]
            if rsid not in variants:
                variants.append(rsid)
            if len(variants) >= N:
                break
    return variants

def write_variants_file(sel_variants,temp_variants_files_path):
    var_file_path = temp_variants_files_path + trait + '_condsig'
    with open(var_file_path,'w') as f:
        for var in sel_variants:
            f.write(var+'\n')
    return var_file_path


if __name__ == "__main__":

    # Generate xdata ydata numpy matrix for all the traits and the top N p-val SNPs excluding those conditional analysis SNPs
    trait_start_index = int(sys.argv[1])
    trait_end_index = int(sys.argv[2])
    map_file = "/home/yx322/tests/TRAIT_MAP.tsv"
    cond_variants_files_path = '/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_variants_com/'
    top_1M_variants_file_path = '/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/500k_blood_traits_GWAS/top_p_val_gwas_1M/'
    plink = "/home/yx322/plink_2.0/plink"

    raw_variants_files_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_and_top_p_variants_raw/"
    if os.path.isdir(raw_variants_files_path) == False:
           os.mkdir(raw_variants_files_path)
    N = 7500
    traits = []
    with open(map_file) as f:
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            traits.append(line[0])

    for i in range(trait_start_index,trait_end_index):

        trait = traits[i]
        print("Processing top p val variants - trait {} ".format(trait))
        top_variants = get_top_P_val_variants(trait,cond_variants_files_path,top_1M_variants_file_path,N)
        trait_variants_file = write_variants_file(top_variants,raw_variants_files_path)



