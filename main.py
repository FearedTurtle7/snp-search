import json
from SNPDataGet import *
import time
import threading
import os

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

def perform_snp_search(input_file, output_file, step_amount, cores, resume):
    start_time = time.time()
    genome_file = input_file
    
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    print(f"Step amount: {step_amount}")
    print(f"Cores: {cores}")
    
    if os.path.exists(output_file) and resume:
        originalOutput = open(output_file, "r")
        outputLines = originalOutput.readlines()
        list_of_outputs = []
        for outputLine in outputLines:
            s = outputLine.replace("\n", "").split(";")
            if s[0].replace("rs", "").isdecimal():
                list_of_outputs.append(s[0])
        last_output = list_of_outputs[-1]
        new_output = False
        print("Resuming with " + last_output)
    else:
        last_output = 0
        new_output = True
    
    additionalLines = 0
    #Splitting genome file by line
    genome = open(genome_file, "r")
    Lines = genome.readlines()
    #list_of_rsids = [["rs328", "CG"],["rs334", "AT"]]
    list_of_rsids = []
    for line in Lines:
        s = line.replace("\n", "").split("\t")
        if s[0].replace("rs", "").isdecimal():
            if not new_output:
                if last_output == s[0]:
                    additionalLines = len(list_of_rsids)
                    list_of_rsids.clear()
                    continue
            if len(s) == 5:
                genome_alleles = s[3] + s[4]
            elif len(s) == 4:
                genome_alleles = s[3]
            rsid_name = s[0]
            rsid_allele_list = []
            rsid_allele_list.append(rsid_name)
            rsid_allele_list.append(genome_alleles)
            list_of_rsids.append(rsid_allele_list)
    

    
    
    print("starting program")
    
    #introduces list_of_rsids in steps
    for rsid_index in range(len(list_of_rsids)):
        #Engages the 
        if rsid_index % step_amount == 0:
            step_time = time.time()
            
            list_of_returns = []
            def execute_code(rsid_number, genome_alleles):
                newvar = SNPData(rsid_number, True)
                rs_obj = newvar.run_snp()
                returns = newvar.get_disease_data(rs_obj, genome_alleles)
                list_of_returns.append(returns)
            
            
            
            x = 0
            y = 0
            templist_totallen = 0
            temp_list = []
            thread_time = time.time()
            for rsid in list_of_rsids[rsid_index:(rsid_index + step_amount)]:
                temp_list.append(rsid)
                if x % cores == 0 or x == step_amount - 1:
                    y+=1
                    rsid_list = []
                    allele_list = []
                    
                    for value in temp_list:
                        rsid_list.append(value[0])
                        allele_list.append(value[1])
                    threads = [threading.Thread(target=execute_code, args= (rsid_tuple[0],rsid_tuple[1])) for rsid_tuple in temp_list]
                    templist_totallen = templist_totallen + len(temp_list)
                    temp_list = []
                    thread_time = time.time()
                    for thread in threads:
                        thread.start()
                    for thread in threads:
                        thread.join()
                x+=1
            
            newvar = SNPData("")
            
            list_of_returns = newvar.remove_empty(list_of_returns)
        
            if list_of_rsids[rsid_index] == list_of_rsids[0] and new_output:
                end_seconds = (time.time()-step_time) * (len(list_of_rsids)/step_amount)
                m, s = divmod(end_seconds, 60)
                h, m = divmod(m, 60)
                d, h = divmod(h, 24)
                print(f"[0]: Expected to take {int(round(d, 0))}d {int(round(h, 0))}h {int(round(m, 0))}m {round(s, 2)}s with an average time of {round((time.time() - start_time)/step_amount, 2)}s per SNP")
                newvar.format_disease_data(list_of_returns, output_file, True)
            elif list_of_rsids[rsid_index] == list_of_rsids[0]:
                end_seconds = (time.time()-step_time) * (len(list_of_rsids)/step_amount)
                m, s = divmod(end_seconds, 60)
                h, m = divmod(m, 60)
                d, h = divmod(h, 24)
                print(f"[0]: Expected to take {int(round(d, 0))}d {int(round(h, 0))}h {int(round(m, 0))}m {round(s, 2)}s with an average time of {round((time.time() - start_time)/step_amount, 2)}s per SNP")
                newvar.format_disease_data(list_of_returns, output_file)
            else:
                newvar.format_disease_data(list_of_returns, output_file)
            
            end_seconds = (time.time()-step_time) * (len(list_of_rsids)/step_amount)
            m, s = divmod(end_seconds, 60)
            h, m = divmod(m, 60)
            d, h = divmod(h, 24)
            os.system('cls' if os.name == 'nt' else 'clear')
            print("------------------------------------------------------")
            print(f"Input file: {input_file}")
            print(f"Output file: {output_file}")
            print(f"Step amount: {step_amount}")
            print(f"Cores: {cores}")
            print("------------------------------------------------------")
            print("\n\n")
            print(f"Progress:                 {rsid_index + step_amount + additionalLines}/{len(list_of_rsids) + additionalLines}")
            print(f"Percentage Progress:      {round((rsid_index + step_amount + additionalLines)/(len(list_of_rsids) + additionalLines)*100, 2)}%")
            print(f"Last Average SNP rate:    {round((time.time() - step_time)/step_amount, 2)}s per SNP")
            print(f"Estimated remaining time: {int(round(d, 0))}d {int(round(h, 0))}h {int(round(m, 0))}m {round(s, 0)}s")
            print("\n\n")
            printProgressBar((rsid_index + step_amount + additionalLines), len(list_of_rsids) + additionalLines, length=50)
    
    end_seconds = time.time()-start_time
    m, s = divmod(end_seconds, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    print(f"[END]: {len(list_of_rsids)} SNPs looked through")
    print(f"[END]: Finished program in {int(round(h, 0))}d {int(round(h, 0))}h {int(round(m, 0))}m {round(s, 2)}s with an average time of {round((time.time() - start_time)/len(list_of_rsids), 2)}s per SNP")