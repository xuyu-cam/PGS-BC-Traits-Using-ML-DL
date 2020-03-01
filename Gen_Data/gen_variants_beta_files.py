
import glob,ntpath,datetime,os


def read_variants_file(variants_file):
    vars = []
    with open(variants_file) as f:
        for line in f:
            line = line.strip()
            vars.append(line)
    return vars

def gen_variants_files_with_rsid_gwas(trait,variants_file,gwas_path,out_path):
    rsids = read_variants_file(variants_file)
    out = out_path + '/' + trait + '_beta'
    fw = open(out, 'w')
    gwas_file = gwas_path + '/' + trait  +'_gwas_normalised_imputed_full_panel.out'
    print("{}-Processing file {}".format(datetime.datetime.now(),gwas_file))
    ids_beta = {}
    with open(gwas_file) as f:
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            cur_rsid  = '_'.join([line[1],line[2],line[4],line[5]])
            ids_beta[cur_rsid] = line[10]


    for cur_rsid in rsids:
        beta = ids_beta[cur_rsid]
        line = [cur_rsid,beta]
        fw.write(" ".join(line) + "\n")
    fw.close()


if __name__ == "__main__":
    variants_files_path = '/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_and_top_p_variants_com_maf_0.01'
    out_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_and_top_p_variants_com_maf_0.01_beta"
    gwas_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/500k_blood_traits_GWAS"

    if os.path.isdir(out_path) == False:
           os.mkdir(out_path)

   # sub_gwas_path = "/scratch/yx322/ukb_blood_traits_genetics/igv_condsig_variants_grs"

    for filepath in glob.glob(variants_files_path + '/' + '*'):
        filename = ntpath.basename(filepath)
        trait = filename.strip().replace('_condsig','')
        if trait in ['mchc']:
            print("{} - Processing trait {}".format(datetime.datetime.now(),trait))
            variants_file = filepath

         #   dosage_path = "/scratch/yx322/ukb_blood_traits_genetics/condsig_MAF_0.1/genetics_dosage/" + trait

            # gen_variants_files_with_dosage_files(trait, dosage_path, out_path)
            gen_variants_files_with_rsid_gwas(trait, variants_file, gwas_path, out_path)
            # gen_variants_files_with_rsid_gwas_sub(trait, sub_gwas_path, out_path)

         #   compare_variants_files(out_path, trait)
