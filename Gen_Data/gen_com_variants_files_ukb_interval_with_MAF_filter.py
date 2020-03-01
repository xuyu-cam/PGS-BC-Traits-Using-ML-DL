
import pandas as pd
import datetime
import glob
import ntpath
import os

def get_name_abbr_map():
    trait_abbr_map ={}
    with open("E:\\Blood Cell Trait Data\\TRAIT_MAP.tsv") as f:
        f.readline()
        for line in f:
            line=line.strip().split('\t')
            trait_abbr_map[line[1]] = line[0]
    return trait_abbr_map

def write_con_variants_files():
    trait_abbr_map = get_name_abbr_map()
    df = pd.read_excel('E:\\Blood Cell Trait Data\\ukbb500k_condsig.xlsx', sheet_name='BCX_final_output')
    trait_variants = {}
    for trait in trait_abbr_map:
        trait_variants[trait_abbr_map[trait]] = []

    for i in range(len(df['Associated Blood Index'])):
        trait_abb = trait_abbr_map[df['Associated Blood Index'][i]]
        rs_id = df['rsID (where available)'][i]
        new_rsid = '_'.join([str(df['Chr (GRCh37)'][i]),str(df['BP (GRCh37)'][i]),str(df['REF (GRC37)'][i]),str(df['ALT (GRC37)'][i])])

        if df['Chr (GRCh37)'][i]!='X' and df['Chr (GRCh37)'][i]!='XY' and df['Minor Allele Frequency'][i]>0.01:
            trait_variants[trait_abb].append(new_rsid)

    path = "E:\\Blood Cell Trait Data\\condsig_variants_new_rsids_maf_0.01\\"

    for trait in trait_variants:
        variants = trait_variants[trait]
        file_name = path+trait+'_condsig'
        with open(file_name,'w') as f:
            for v in variants:
                print(v)
                f.write(v+'\n')

    return trait_variants

def read_variants_files_by_chr(variant_file):
    # read variants of a trait by chromosome
    # e.g. variants = [[],[],....,[]]; variants[i] stores the variants of chromosome i+1 for the given trait
    variants =  [[] for x in range(22)]
    with open(variant_file) as f:
        for line in f:
            chr = int(line.strip().split('_')[0])
            variants[chr-1].append(line.strip())
    return variants

def read_variants_by_order_in_file(variant_file):
    # read variants of a trait into a list
    vars = []
    with open(variant_file) as f:
        for line in f:
            vars.append(line.strip())
    return vars


def gen_com_variants_ukb_interval_by_chr(trait_variants,trait_variants_order,ukb_bim_path,interval_bim_path):
    # Only keep variants available on both UKB and Interval
    # ukb_bim_path to locate bim files of UKB cohort
    # interval_bim_path to locate bim files of Interval cohort

    trait_com_vars = {}
    trait_com_vars_order = {}
    for trait in trait_variants:
        trait_com_vars[trait] = []
        trait_com_vars_order[trait] = []

    for i in range(1,23):
        ukb_bim_file = ukb_bim_path + 'ukb_imp_chr' + str(i) + '.bim'                        #name of the bim files - changable
        interval_bim_file = interval_bim_path + 'interval_impute_chr' + str(i) + '.bim'      #name of the bim files - changable
        ukb_vars = []
        interval_vars = []
        print("{} Reading variants of Chr{} of UKB and Interval Cohorts".format(datetime.datetime.now(),i))
        print("UKB bim file: {}".format(ukb_bim_file))
        print("Interval bim file: {}".format(interval_bim_file))
        with open(ukb_bim_file) as f1,open(interval_bim_file) as f2:
            for line in f1:
                var = line.strip().split('\t')[1]
                ukb_vars.append(var)
            for line in f2:
                var = line.strip().split('\t')[1]
                interval_vars.append(var)

        print("{}-Comparing common variants of chr{} in both UKB and Interval for all the traits".format(datetime.datetime.now(), i))

        for trait in trait_variants:
            variants = trait_variants[trait]
            print("{} Processing trait: {}".format(datetime.datetime.now(),trait))
            for var in variants[i-1]:
                if var in ukb_vars and var in interval_vars:
                    trait_com_vars[trait].append(var)

    for trait in trait_variants_order:
        for var in trait_variants_order[trait]:
            if var in trait_com_vars[trait]:
                trait_com_vars_order[trait].append(var)

    return trait_com_vars_order


def filter_variants_with_maf(trait_variants,ukb_bim_path,maf_th):
    # filtering out variants with the MAF threshold set at maf_th
    # Using ukb_bim_path to access the MAF of variants in UKB (.afreq files - PLINK2)

    for i in range(1,23):
        # chrs 1-22
        frq_file = ukb_bim_path + 'ukb_imp_chr' + str(i) + ".afreq"   # name of the variants frequency files
        print("UKB MAF file: {}".format(frq_file))
        ukb_vars_freq = {}
        with open(frq_file) as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                ukb_vars_freq[line[1]] = line[4]

        for trait in trait_variants:
            new_vars = list(trait_variants[trait][i-1])

            print("{} Processing trait: {}".format(datetime.datetime.now(), trait))
            for var in trait_variants[trait][i-1]:
                frq1 = float(ukb_vars_freq[var])
                frq2 = 1-frq1
                if (frq1 < maf_th or frq2 <maf_th):
                    new_vars.remove(var)
                 #   print(var, frq1, frq2)
            trait_variants[trait][i-1] = new_vars
    return trait_variants



def gen_com_variants_ukb_interval_with_MAF_filter(raw_variants_file_path,ukb_bim_path,interval_bim_path,write_path,maf_th,max_n):
    # for each trait, this generates a variant file that stores variants that have MAF > maf_th and avaialble on both UKB and Interval cohorts under the folder write_path;
    # Under raw_variants_file_path folder, each file stores a set of initially selected variants for a trait with file name "trait_condsig"

    if os.path.isdir(write_path) == False:
           os.mkdir(write_path)

    trait_variants = {}
    trait_variants_order = {}
    for filepath in glob.glob(raw_variants_file_path + '*'):
        filename = ntpath.basename(filepath)
        trait = filename.strip().replace('_condsig','')
        print("Reading variants of trait: {} with file {}".format(trait,filepath))
        trait_variants_file = filepath
        variants = read_variants_files_by_chr(trait_variants_file)
        variants_order = read_variants_by_order_in_file(trait_variants_file)

        trait_variants[trait] = variants
        trait_variants_order[trait] = variants_order

    print("------ \n Finished reading all the listed variants of all the traits \n------")

    trait_variants = filter_variants_with_maf(trait_variants, ukb_bim_path, maf_th)

    print("------ \n Finished MAF variants filtering of all the traits \n------")

    trait_com_vars = gen_com_variants_ukb_interval_by_chr(trait_variants, trait_variants_order,ukb_bim_path, interval_bim_path)

    print("{}- Writing selected variants by trait".format(datetime.datetime.now()))
    for trait in trait_com_vars:
        write_file = write_path + trait + '_condsig'
        variants = trait_com_vars[trait]

        if len(variants)>max_n:
            variants = variants[:max_n]
        with open(write_file,'w') as f:
            for each in variants:
                f.write(each + '\n')


if  __name__ == "__main__":
    raw_variants_file_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_and_top_p_variants_raw/"
    ukb_bim_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_genetics_new_rsids/"
    interval_bim_path =  "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/interval_genetics_new_rsids/"
    write_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_and_top_p_variants_com_maf_0.01/"
    maf_th = 0.01
    max_n = 5000
    gen_com_variants_ukb_interval_with_MAF_filter(raw_variants_file_path, ukb_bim_path, interval_bim_path, write_path, maf_th,max_n)