# -*- coding: utf-8 -*-

import requests
import json
import re
import time

api_rootURL = 'https://api.ncbi.nlm.nih.gov/variation/v0/'
headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.90 Safari/537.36'}

class SNPData:
    def __init__(self, rsid, continuance=False):
        self.rsid = rsid
        
        if isinstance(self.rsid, int):
            self.rsid = self.rsid
        else:
            m_rsid = re.fullmatch('(rs)?(?P<rsid>[1-9]\d*)', self.rsid)
            if m_rsid:
                self.rsid = int(m_rsid.group('rsid'))
        if continuance == True:
            self.run_snp()
                
    def run_snp(self):
        url = api_rootURL + 'beta/refsnp/' + str(self.rsid)
        self.req = requests.get(url, headers=headers)
        self.rs = json.loads(self.req.text)
        try:
            data = self.rs['primary_snapshot_data']
        except:
            data = self.rs
        
        # allele_annot = []
        # for annot in data['allele_annotations']:
        #     for clininfo in annot['clinical']:
        #         allele_annot.append(clininfo['clinical_significances'])
        
        return data
    
    def get_disease_data(self, rs_obj, genome_alleles):
        if "error" in rs_obj:
            rs_error = {}
            rs_error["rsid"] = "rs" + str(self.rsid)
            rs_error["error"] = "Rsid not found" 
            return rs_error
        
        rsid_fulldict = {}
        rsid_fulldict["rsid"] = "rs" + str(self.rsid)
        try:
            allele_annotations = rs_obj["allele_annotations"]
        except:
            rs_error = {}
            rs_error["rsid"] = "rs" + str(self.rsid)
            rs_error["error"] = "Rsid not found" 
            return rs_error
        
        clinical_data_list = []
        
        #Iterates through every variation of the SNP
        for annot in allele_annotations:
            
            #Checks to see if the variation has any clinical significance
            #If it doesn't --> discard
            if not annot["clinical"]:
                continue
            clinical_data = {}
            
            #if it does --> get the severity and the disease
            #Gets the nucleotide of the mutation when it finds one
            if annot["frequency"] != []:
                mutated_nucleotide = annot["frequency"][0]["observation"]["inserted_sequence"]
            else:
                continue
                
            if mutated_nucleotide not in genome_alleles:
                continue
            #Gets the disease name and the severity(clinical_significance)
            clin_sig_list = []
            disease_name_list = []
            for clininfo in annot['clinical']:
                disease_names = clininfo["disease_names"]
                clinical_significances = clininfo['clinical_significances']
                clin_sig_list.append(clinical_significances[0])
                disease_name_list.append(disease_names[0])
            
            #adds mutated nucleotide, list of diseases, and list of clinical signficances to clinical data dictionary
            clinical_data["mutated_nucleotide"] = mutated_nucleotide
            clinical_data["disease_names"] = disease_name_list
            clinical_data["clinical_significances"] = clin_sig_list
            #print("didnt continue")
            
            clinical_data_list.append(clinical_data)
            
        #adds all clinical data to large dictionary
        rsid_fulldict["genome_alleles"] = genome_alleles
        rsid_fulldict["clinical_data"] = clinical_data_list
        return rsid_fulldict
        
    def remove_empty(self, res_obj):
        entire_obj = []
        for i in range(len(res_obj)):
            new_res_obj = []
            if "error" in res_obj[i]:
                file2= open("errors.txt", "a")
                file2.write(f"error: {res_obj[i]['rsid']}\n")
            else:
                for j in range(len(res_obj[i]["clinical_data"])):
                    if not res_obj[i]["clinical_data"][j]:
                        continue
                    else:
                        if res_obj[i] in new_res_obj:
                            continue
                        # elif "not specified" in res_obj[i]["clinical_data"][j]["disease_names"] or "not provided" in res_obj[i]["clinical_data"][j]["disease_names"]:
                        #     continue
                        else:
                            new_res_obj.append(res_obj[i])
                            entire_obj.append(res_obj[i])
                        
        return entire_obj
    
    
    #Outputs data line by line into snp_output.txt (rsid, patient allelels, study allele, disease names)
    def format_disease_data(self, rsid_list, output_file, firstopen=False):
        #If it is the first time to be opened, erase existing data and write over it
        if firstopen:
            file = open(output_file, "w")
            file.write("rsid;patient_alleles;study_allele;disease_names;Pathogenicity\n")
            #If its not, append the data
        else:
            file = open(output_file, "a")
            
        #Sort through every item in the rs_obj returned
        for item in rsid_list:
            rsid = item["rsid"]
            genome_alleles = item["genome_alleles"]
            i = 0
            
            #Search for alleles in each rs_obj
            for alleles in item["clinical_data"]:
                j = 0
                #Search through all the disease names associated with each allele
                for diseases in item["clinical_data"][i]["disease_names"]:
                    #For every disease name in the allele data, get the allele nucleotide, disease name, and pathogenicity
                    mutated_nucleotide = item["clinical_data"][i]["mutated_nucleotide"]
                    disease_names = item["clinical_data"][i]["disease_names"][j]
                    clinical_significance = item["clinical_data"][i]["clinical_significances"]
                    #Check the pathogenicity of the allele and assign it a variable
                    if "pathogenic" in clinical_significance:
                        pathogenic = True
                    else:
                        pathogenic = False
                    if "benign" in clinical_significance or "likely-benign" in clinical_significance:
                        benign = True
                    else:
                        benign = False
                    #Write the result into a variable
                    if benign and pathogenic:
                        pathogenicity = "Pathogenic or Benign"
                    elif benign and not pathogenic:
                        pathogenicity = "Benign"
                    elif not benign and pathogenic:
                        pathogenicity = "Pathogenic"
                    else:
                        pathogenicity = "Unknown"
                    
                    #If it is the first disease name and allele, print the rsid, allele, disease name, and pathogenicity
                    if i == 0 and j == 0:
                        first_line = f"{rsid};{genome_alleles};{mutated_nucleotide};{disease_names};{pathogenicity}\n"
                        file.write(first_line)
                    else:
                        #If it is the first allele, print the allele, disease names, and pathogenicity
                        if j == 0:
                            line = f";;{mutated_nucleotide};{disease_names};{pathogenicity}\n"
                        #If it is the second disease name of the allele, only print the disease name
                        else:
                            line = f";;;{disease_names}\n"
                        file.write(line)
                    j+=1
                i+=1
        file.close()
        
        
        