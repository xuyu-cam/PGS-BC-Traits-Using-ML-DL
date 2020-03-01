import glob,ntpath,datetime


def compare_variants_files(variants_file,compare_variants_files_path, trait):

    compare_file = compare_variants_files_path + '/' + trait + '_interval_ukb_common_MAF_0.1'
    vars = []
    vars_comp= []
    with open(variants_file) as f1:
        for line in f1:
           vars.append(line.strip())

    with open(compare_file) as f2:
        for line in f2:
            line = line.strip().split()
            line = [line[0]] + line[3:6]
            id = '_'.join(line)
            vars_comp.append(id)

    print("variants not in vars_comp but in the new vars set:")
    for var in vars:
        if var not in vars_comp:
            print(var)

    print("variants not in the new variants set but in old set:")
    for var in vars_comp:
        if var not in vars:
            print(var)




if __name__ == "__main__":

    variants_files_path = '/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_variants_com_maf_0.01'
    compare_variants_files_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/interval_blood_traits_genetics/condsig_com_UKB_INT__MAF_0.1/variants_full_info"


    for filepath in glob.glob(variants_files_path + '/' + '*'):
        filename = ntpath.basename(filepath)
        trait = filename.strip().replace('_condsig','')
        print("\n{} - Processing trait {}".format(datetime.datetime.now(),trait))
        variants_file = filepath

        # gen_variants_files_with_dosage_files(trait, dosage_path, out_path)
        # gen_variants_files_with_rsid_gwas(trait, variants_file, gwas_path, out_path)
        # gen_variants_files_with_rsid_gwas_sub(trait, sub_gwas_path, out_path)

        compare_variants_files(variants_file,compare_variants_files_path, trait)