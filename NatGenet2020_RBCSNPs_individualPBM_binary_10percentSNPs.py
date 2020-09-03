#!/usr/bin/env python

### Defining reverse complement function

def reverse_complement(DNA):
    reverse = DNA.upper()[::-1]
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return(''.join([complement[base] for base in reverse]))

### Defining function that returns 8 base pair sliding window from input of 15 bp sequence, major allele, minor allele

def eightbp_major_minor(seq_15, major, minor):
    major_15 = seq_15[0:7] + major + seq_15[8:]
    minor_15 = seq_15[0:7] + minor + seq_15[8:]
    list_of_major_minor_pairs = []
    for i in range(0, 8):
        list_of_major_minor_pairs.append([major_15[i:(i+8)].upper(), minor_15[i:(i+8)].upper()])
    return(list_of_major_minor_pairs)

### Getting all common SNPs in dbSNP, with indels removed
### Filtering out any SNPs in protein coding regions (NSF, NSM, NSN, SYN)

import re

file_all_rsIDs = open("SNPs_caf.recode.vcf", "r")

count = 0
count_NSF = 0
count_NSM = 0
count_NSN = 0
count_SYN = 0
count_nsSNP = 0
hash_rsIDs_common_SNPs = {}

for readLine in file_all_rsIDs:
    count += 1
    if (count >= 58):
        list_of_rsID_info = readLine.split()
        match_NSF = re.search("NSF=1", list_of_rsID_info[7])
        if match_NSF:
            count_NSF += 1
            continue
        match_NSM = re.search("NSM=1", list_of_rsID_info[7])
        if match_NSM:
            count_NSM += 1
            continue
        match_NSN = re.search("NSN=1", list_of_rsID_info[7])
        if match_NSN:
            count_NSN += 1
            continue
        match_SYN = re.search("SYN=1", list_of_rsID_info[7])
        if match_SYN:
            count_SYN += 1
            continue
        count_nsSNP += 1
        hash_rsIDs_common_SNPs[readLine.split()[2]] = (list_of_rsID_info[0], list_of_rsID_info[1],
                                                       list_of_rsID_info[2], list_of_rsID_info[3],
                                                       list_of_rsID_info[4], list_of_rsID_info[7])

file_all_rsIDs.close()
print(count)
print(count_NSF)
print(count_NSM)
print(count_NSN)
print(count_SYN)
print(count_nsSNP)

hash_rsIDs_10percent_SNPs = {}
count_ref_allele = 0
count_1_alt_allele = 0
count_2_alt_allele = 0
count_3_alt_allele = 0

for key in hash_rsIDs_common_SNPs:
    count_ref_allele = 0
    if len(hash_rsIDs_common_SNPs[key][4]) == 1:
        count_1_alt_allele += 1
        match_test = re.search("CAF=(.*),(.*)$", hash_rsIDs_common_SNPs[key][5])
        if (float(match_test.group(1)) > 0.1) & (float(match_test.group(2)) > 0.1):
            hash_rsIDs_10percent_SNPs[key] = hash_rsIDs_common_SNPs[key]
        continue
    if len(hash_rsIDs_common_SNPs[key][4]) == 3:
        count_2_alt_allele += 1
        match_test = re.search("CAF=(.*),(.*),(.*)$", hash_rsIDs_common_SNPs[key][5])
        count_check_10percent = 0
        if match_test.group(1) != ".":
            if (float(match_test.group(1)) > 0.1):
                count_check_10percent += 1
        if match_test.group(2) != ".":
            if (float(match_test.group(2)) > 0.1):
                count_check_10percent += 1
        if match_test.group(3) != ".":
            if (float(match_test.group(3)) > 0.1):
                count_check_10percent += 1
        if count_check_10percent > 1:
            hash_rsIDs_10percent_SNPs[key] = hash_rsIDs_common_SNPs[key]
        continue
    if len(hash_rsIDs_common_SNPs[key][4]) == 5:
        count_3_alt_allele += 1
        match_test = re.search("CAF=(.*),(.*),(.*),(.*)$", hash_rsIDs_common_SNPs[key][5])
        count_check_10percent = 0
        if match_test.group(1) != ".":
            if (float(match_test.group(1)) > 0.1):
                count_check_10percent += 1
        if match_test.group(2) != ".":
            if (float(match_test.group(2)) > 0.1):
                count_check_10percent += 1
        if match_test.group(3) != ".":
            if (float(match_test.group(3)) > 0.1):
                count_check_10percent += 1
        if match_test.group(4) != ".":
            if (float(match_test.group(4)) > 0.1):
                count_check_10percent += 1
        if count_check_10percent > 1:
            hash_rsIDs_10percent_SNPs[key] = hash_rsIDs_common_SNPs[key]
        continue

print(len(hash_rsIDs_common_SNPs))
hash_rsIDs_common_SNPs = {}
print(len(hash_rsIDs_10percent_SNPs))
print(count_ref_allele)
print(count_1_alt_allele)
print(count_2_alt_allele)
print(count_3_alt_allele)

### Getting individual variant alleles out as single bases in tuples, so that sequence substitution can be performed

count_ref_allele = 0
count_1_alt_allele = 0
count_2_alt_allele = 0
count_2_alt_allele_multiple_3 = 0
count_3_alt_allele = 0
count_3_alt_allele_multiple_3 = 0
count_3_alt_allele_multiple_4 = 0
hash_individual_rsIDs_10percent_SNPs = {}

for key, value in hash_rsIDs_10percent_SNPs.items():
    count_ref_allele += 1
    if len(hash_rsIDs_10percent_SNPs[key][4]) == 1:
        hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4]] = (value[0], value[1], value[2], value[3], value[4])
        count_1_alt_allele += 1
        continue
    if len(hash_rsIDs_10percent_SNPs[key][4]) == 3:
        count_2_alt_allele += 1
        match_test = re.search("CAF=(.*),(.*),(.*)$", hash_rsIDs_10percent_SNPs[key][5])
        count_check_10percent = 0
        if (match_test.group(1) != ".") & (match_test.group(2) != ".") :
            if (float(match_test.group(1)) > 0.1) & (float(match_test.group(2)) > 0.1):
                hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][0]] = (value[0], value[1], 
                                                                                 value[2], value[3], value[4][0])
                count_check_10percent += 1
        if (match_test.group(1) != ".") & (match_test.group(3) != ".") :
            if (float(match_test.group(1)) > 0.1) & (float(match_test.group(3)) > 0.1):
                hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][2]] = (value[0], value[1], 
                                                                                 value[2], value[3], value[4][2])
                count_check_10percent += 1
        if count_check_10percent == 2:
            count_2_alt_allele_multiple_3 += 1
            print(hash_rsIDs_10percent_SNPs[key])
            print(hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][0]])
            print(hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][2]])
        continue
    if len(hash_rsIDs_10percent_SNPs[key][4]) == 5:
        count_3_alt_allele += 1
        match_test = re.search("CAF=(.*),(.*),(.*),(.*)$", hash_rsIDs_10percent_SNPs[key][5])
        count_check_10percent = 0
        if (match_test.group(1) != ".") & (match_test.group(2) != ".") :
            if (float(match_test.group(1)) > 0.1) & (float(match_test.group(2)) > 0.1):
                hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][0]] = (value[0], value[1],
                                                                                 value[2], value[3], value[4][0])
                count_check_10percent += 1
        if (match_test.group(1) != ".") & (match_test.group(3) != ".") :
            if (float(match_test.group(1)) > 0.1) & (float(match_test.group(3)) > 0.1):
                hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][2]] = (value[0], value[1], 
                                                                                 value[2], value[3], value[4][2])
                count_check_10percent += 1
        if (match_test.group(1) != ".") & (match_test.group(4) != ".") :
            if (float(match_test.group(1)) > 0.1) & (float(match_test.group(4)) > 0.1):
                hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][4]] = (value[0], value[1], 
                                                                                 value[2], value[3], value[4][4])
                count_check_10percent += 1
        if count_check_10percent == 2:
            count_3_alt_allele_multiple_3 += 1
            print(hash_rsIDs_10percent_SNPs[key])
            if str(key + "_" + value[4][0]) in hash_individual_rsIDs_10percent_SNPs:
                print(hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][0]])
            if str(key + "_" + value[4][2]) in hash_individual_rsIDs_10percent_SNPs:
                print(hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][2]])
            if str(key + "_" + value[4][4]) in hash_individual_rsIDs_10percent_SNPs:
                print(hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][4]])
            continue
        if count_check_10percent == 3:
            count_3_alt_allele_multiple_4 += 1
            print(hash_rsIDs_10percent_SNPs[key])
            print(hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][0]])
            print(hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][2]])
            print(hash_individual_rsIDs_10percent_SNPs[key + "_" + value[4][4]])
        continue
    
print(len(hash_individual_rsIDs_10percent_SNPs))
print(count_ref_allele)
print(count_1_alt_allele)
print(count_2_alt_allele)
print(count_3_alt_allele)

print(count_2_alt_allele_multiple_3)
print(count_3_alt_allele_multiple_3)
print(count_3_alt_allele_multiple_4)

### Using Biopython SeqIO module to obtain 15-bp window with each SNP at the centre - 
### for getting sliding window for modules

testcount = 0
biallele_count = 0
multi_2_allele_count = 0
multi_3_allele_count = 0
hash_rsID_8mer_major_minor = {}

from Bio import SeqIO

list_chr = []

for i in range(1, 23):
	list_chr.append(i)
list_chr.append("X")
list_chr.append("Y")

for str_chr in list_chr:
    str_file_name_chr = "chr" + str(str_chr) + ".fa"
    for seq_record in SeqIO.parse(str_file_name_chr, "fasta"):
        for k, v in hash_individual_rsIDs_10percent_SNPs.items():
            if v[0] == str(str_chr):
                chrpos = int(v[1])
                seq_15_current = str(seq_record.seq[(chrpos-8):(chrpos+7)]).upper()
                if len(v[4]) == 1:
                    biallele_count += 1
                    hash_rsID_8mer_major_minor[k] = eightbp_major_minor(seq_15_current, v[3], v[4])
                else:
                	print(hash_individual_rsIDs_10percent_SNPs[k])

### Creating dictionary versions of TF PBM files of interest

import os

list_PBM_files = os.listdir("Individual_PBMs_18092018/")
hash_all_PBMs_with_Escores = {}

for i in range(0, len(list_PBM_files)):
    file = list_PBM_files[i]
    m = re.match("(.*)_Escore\.txt", file)
    str_PBM_name = (m.group(1))
    hash_all_PBMs_with_Escores[str_PBM_name] = {}
    PBM_file = open("Individual_PBMs_18092018/" + file, "r")
    for readLine in PBM_file:
        key = readLine.rstrip("\n").split("\t")[0]
        value = readLine.rstrip("\n").split("\t")[1]
        key_rc = reverse_complement(key)
        hash_all_PBMs_with_Escores[str_PBM_name][key] = value
        hash_all_PBMs_with_Escores[str_PBM_name][key_rc] = value
    PBM_file.close()

### Getting average GATA PBM dataset to represent GATA3, GATA4, GATA5 and GATA6 datasets

hash_all_PBMs_with_Escores["GATA_average"] = {}
file_GATA_average = open("GATA_average_Escore.txt", "w")

for i in hash_all_PBMs_with_Escores["Gata3_1024_contig8mers"]:
    Escore_Gata3 = float(hash_all_PBMs_with_Escores["Gata3_1024_contig8mers"][i])
    Escore_Gata4 = float(hash_all_PBMs_with_Escores["Gata4_8mers_11111111"][i])
    Escore_Gata5 = float(hash_all_PBMs_with_Escores["Gata5_3768_contig8mers"][i])
    Escore_Gata6 = float(hash_all_PBMs_with_Escores["Gata6_3769_contig8mers"][i])
    average_Escore_Gata = (Escore_Gata3 + Escore_Gata4 + Escore_Gata5 + Escore_Gata6)/4
    hash_all_PBMs_with_Escores["GATA_average"][i] = average_Escore_Gata
    file_GATA_average.write(i + "\t" + str(average_Escore_Gata) + "\n")
file_GATA_average.close()

### Adapting individual PBMs to scan common SNPs with >10% allele frequency - what types of binding changes do we see? 
### Used criteria of E-score > 0.35 for one allele and E-score < 0.3 for other allele

hash_PBM_statistics = {}
file_to_output = open("file_to_output_Sept2018.txt", "w")
	
for i in range(0, len(list_PBM_files)):
    file = list_PBM_files[i]
    m = re.match("(.*)_Escore\.txt", file)
    str_PBM_name = (m.group(1))

    count_loss_binding = 0
    count_surprise_no_loss_binding = 0
    count_gain_binding = 0
    count_surprise_no_gain_binding = 0
    list_eightmer_change = []

    for k, v in hash_rsID_8mer_major_minor.items():
    	eight8mer_change = 0
    	for i in range(0, 8):
    		allele1_key = v[i][0]
    		allele2_key = v[i][1]
    		if allele1_key not in hash_all_PBMs_with_Escores[str_PBM_name]:
    			print(allele1_key)
    			continue
    		if allele2_key not in hash_all_PBMs_with_Escores[str_PBM_name]:
    			print(allele2_key)
    			continue
    		if (float(hash_all_PBMs_with_Escores[str_PBM_name][allele1_key]) > 0.35) & (float(hash_all_PBMs_with_Escores[str_PBM_name][allele2_key]) < 0.3):
    			count_loss_binding += 1
    			eight8mer_change += 1
    		if (float(hash_all_PBMs_with_Escores[str_PBM_name][allele2_key]) > 0.35) & (float(hash_all_PBMs_with_Escores[str_PBM_name][allele1_key]) < 0.3):
    			count_gain_binding += 1
    			eight8mer_change += 1
    	list_eightmer_change.append(eight8mer_change)
	
    testing_write = open(str_PBM_name + ".txt", "w")
    for item in list_eightmer_change:
        testing_write.write("%s\n" % item)
    testing_write.close()
    hash_PBM_statistics[str_PBM_name] = (str_PBM_name, len(hash_rsID_8mer_major_minor)*8, count_loss_binding, count_gain_binding, sum(list_eightmer_change))
    print(hash_PBM_statistics[str_PBM_name])
    file_to_output.write(str_PBM_name + "\t" + str(len(hash_rsID_8mer_major_minor)*8)  + "\t" + str(count_loss_binding)  + "\t" + str(count_gain_binding))

file_to_output.close()