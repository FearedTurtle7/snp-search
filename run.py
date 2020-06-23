import sys, getopt
from main import *

def main(argv):
    inputfile = ''
    outputfile = ''
    step_amount = 0
    cores = 0
    
    try:
        opts, args = getopt.getopt(argv, "hi:o:s:c:")
        print(opts)
    except getopt.GetoptError:
        print("run.py -i <inputfile> -o <outputfile> -s <checkpoint length> (optional) -c <cores> (optional)")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("run.py -i <inputfile> -o <outputfile> -s <step amount> (optional) -c <cores> (optional)")
            sys.exit()
        elif opt == "-i":
            inputfile = arg
        elif opt == "-o":
            outputfile = arg
        elif opt == "-s":
            cores= int(arg)
        elif opt == "-c":
            step_amount = int(arg)
        
    
    perform_snp_search(inputfile, outputfile, step_amount, cores)
        
        
if __name__ == "__main__":
    main(sys.argv[1:])