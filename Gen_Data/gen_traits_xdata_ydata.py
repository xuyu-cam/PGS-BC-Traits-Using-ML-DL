

import numpy as np
import os
import glob
import subprocess
import datetime

from collections import defaultdict
from itertools import izip
import ntpath
import pandas as pd
import shutil
import sys


#Generating plink genotype files of a given set of variants of a trait
def extract_variants_bed_files(trait,plink,cohort_bed_files_path,cohort_file_name,trait_variants_file,temp_bed_files,log_path):
    print(" {} Extrating variants of trait {} to {} using Plink bed files at {}".format(datetime.datetime.now(),trait,temp_bed_files,cohort_bed_files_path))
    with open(log_path +'/' +trait + '_write_bed_log.txt', 'w') as f:
        for i in range(1,23):
            out = temp_bed_files + '/' + trait +'_chr'+str(i)
            ukb_bed_files_chr = cohort_bed_files_path +'/'+ cohort_file_name +str(i)
            subprocess.call([plink, '--bfile', ukb_bed_files_chr,'--extract',trait_variants_file,'--make-bed', '--out', out],stdout=f)
    print(" {} Finished variants extraction".format(datetime.datetime.now()))

def get_sample_ids(ids):
    new_ids = []
    for each in ids:
        id = each.strip().split('_')[0]
        new_ids.append(id)
    return new_ids

def convert_to_dosage_one_chr(trait,bfile_path,outfile_path,plink_path,log_path,chr_num):

  bfile = bfile_path
  out = outfile_path
  plink = plink_path

  with open(log_path + '/' + trait + '_write_dosage_log_chr'+str(chr_num), 'w') as f:
      subprocess.call([plink, '--bfile', bfile,'--export', 'A-transpose', '--out', out],stdout=f)
      subprocess.call([plink, '--bfile', bfile, '--freq', '--out', out],stdout=f  )

  diter = iter(x.split() for x in open(out + ".traw"))
  fiter = iter(x.split() for x in open(out + ".afreq"))
  biter = iter(x.split() for x in open(bfile + ".bim"))

  buff = defaultdict(list)

  diter.next()
  fiter.next()

  for (dcols, fcols, bcols) in izip(diter, fiter, biter):
    nline = [fcols[1]]

    MAF = 1- float(fcols[4])
    nline = nline + [str(MAF*2) if e == "NA" else e for e in dcols[6:]]

    buff[fcols[0]].append(" ".join(nline) + "\n")

  for curCHR, lines in buff.iteritems():
    with open(out +'_chr' +curCHR + ".txt", "w") as ofile:
      ofile.writelines(lines)

  os.remove(out + ".afreq")
  os.remove(out + ".traw")
  os.remove(out + ".log")
  if os.path.isfile(out + ".nosex"):
      os.remove(out + ".nosex")

#generating dosage files by using a given set of plink bed files
def gen_dosage_files(trait,plink,temp_bed_files,temp_dosage_files,xdata_path,log_path):
    print(" {} Generating xdata (samples-snps) csv file of trait {} at: {}".format(datetime.datetime.now(),trait,xdata_path))
    df = pd.DataFrame()
    with open(log_path + '/'+trait + '_write_xdata_log.txt', 'w') as f: # log file
        for i in range(1, 23): # for each chromosome
            plink_path = plink
            bfile_path = temp_bed_files + "/" + trait + "_chr" + str(i)
            outfile_path = temp_dosage_files + "/" + trait

            if os.path.isfile(bfile_path + '.bed'):
                bfile = bfile_path
                out = outfile_path
                plink = plink_path

                with open(log_path + '/' + trait + '_write_dosage_log_chr' + str(i), 'w') as f:
                    subprocess.call([plink, '--bfile', bfile, '--export', 'A-transpose', '--out', out], stdout=f)
                    subprocess.call([plink, '--bfile', bfile, '--freq', '--out', out], stdout=f)

                diter = iter(x.split() for x in open(out + ".traw"))
                fiter = iter(x.split() for x in open(out + ".afreq"))

                header = diter.next()
                fiter.next()
                sample_ids = get_sample_ids(header[6:])
                if i ==1:
                    df["sample_ids"] = sample_ids

                for (dcols, fcols) in izip(diter, fiter):
                    rsid = fcols[1]
                    MAF = 1 - float(fcols[4])
                    nline = []
                    nline =  nline + [str(MAF * 2) if e == "NA" else e for e in dcols[6:]]
                    df[rsid]= nline

                os.remove(out + ".afreq")
                os.remove(out + ".traw")
                os.remove(out + ".log")
                if os.path.isfile(out + ".nosex"):
                    os.remove(out + ".nosex")

            else:
                pass
    save_file_name =  xdata_path + "/" + trait + '_xdata' + ".csv.gz"
    # df.to_pickle(save_file_name)
    df.to_csv(save_file_name,index=False,compression='gzip')

    print(" {} Finished xdata csv file generation for trait {} \n-------".format(datetime.datetime.now(),trait))


#Generate the X and Y matrixes, given a variants list for a trait
def gen_xdata_ydata_per_trait(trait,plink,cohort_bed_files_path,cohort_file_name,trait_variants_file,temp_bed_files,temp_dosage_files,xdata_path,log_path):

    #generating the bed files of the given set of variants for a trait
    extract_variants_bed_files(trait, plink, cohort_bed_files_path, cohort_file_name, trait_variants_file,temp_bed_files, log_path)

    #generating dosage files (easy to be read using dataframe) of a given set of plink bed files
    gen_dosage_files(trait, plink, temp_bed_files,temp_dosage_files, xdata_path, log_path)

if __name__ == "__main__":

    # input for which trait to be calculated
    trait_start_index = int(sys.argv[1])
    trait_end_index = trait_start_index +1
    # trait_end_index = int(sys.argv[2])

    #files stores a list of traits being studied
    map_file = "/TRAIT_MAP.tsv"

    #The folder path stores files having the variants of each trait
    variants_files_path = "/variants/condsig_variants_com"

    #Plink path
    plink = "/plink_2.0/plink2"

    #Folder stores the Plink genotype files of the whole UKB
    ukb_bed_files = "/ukb_genetics_new_rsids"

    #Bed file prefix
    bed_file_name = 'ukb_imp_chr'

    #A folder saves teporarily generated files
    temp_path = "/temp/"

    #A folder saves log files - an empty folder
    log_path = "/temp/temp_logs"

    #Folder saves the extracted genotype files
    xdata_path = "/condsig_variants_com"

    if os.path.isdir(xdata_path) == False:
        os.mkdir(xdata_path)

    traits = []
    #read the trait list
    with open(map_file) as f:
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            traits.append(line[0])

    # loop for the traits from index trait_start_index, to index trait_end_index
    for i in range(trait_start_index,trait_end_index):
        #Current trait
        trait = traits[i]

        #the file stores all the variants of a trait to be extracted - each row a variant
        trait_variants_file = variants_files_path + '/' + trait + '_condsig'
        print(" ------- \n {} \n Processing trait: {} \n variant file: {} \n ------".format(datetime.datetime.now(),trait, trait_variants_file))

        #Generated temporary bed files path and dosage files path for the current trait
        temp_bed_files = temp_path + trait + "_temp_beds"
        temp_dosage_files = temp_path + trait + "_temp_dosage"

        if os.path.isdir(temp_bed_files) == False:
            os.mkdir(temp_bed_files)
        else:
            shutil.rmtree(temp_bed_files)
            os.mkdir(temp_bed_files)

        if os.path.isdir(temp_dosage_files) == False:
            os.mkdir(temp_dosage_files)
        else:
            shutil.rmtree(temp_dosage_files)
            os.mkdir(temp_dosage_files)

        #Extracting funtion for a trait
        gen_xdata_ydata_per_trait(trait, plink, ukb_bed_files, bed_file_name, trait_variants_file, temp_bed_files,temp_dosage_files,xdata_path, log_path)

        shutil.rmtree(temp_bed_files)
        shutil.rmtree(temp_dosage_files)