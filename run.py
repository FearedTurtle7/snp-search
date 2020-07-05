import sys, getopt
from main import *

def main(argv):
    inputfile = ''
    outputfile = ''
    step_amount = 0
    cores = 0
    resume = False
    
    try:
        opts, args = getopt.getopt(argv, "hi:o:s:c:r")
        print(opts)
        print(args)
    except getopt.GetoptError:
        print("run.py -i <inputfile> -o <outputfile> -s <checkpoint length> -c <cores>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("run.py -i <inputfile> -o <outputfile> -s <step amount> -c <cores>")
            sys.exit()
        elif opt == "-i":
            inputfile = arg
        elif opt == "-o":
            outputfile = arg
        elif opt == "-c":
            cores= int(arg)
        elif opt == "-s":
            step_amount = int(arg)
        elif opt == "-r":
            resume = True
        
    
    perform_snp_search(inputfile, outputfile, step_amount, cores, resume)
        
        
if __name__ == "__main__":
    main(sys.argv[1:])