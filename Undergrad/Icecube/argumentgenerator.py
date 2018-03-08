#!/usr/bin/env python

import glob
import argparse

desc = "Generates a list of arguments for extractor.py that can be passed " + \
       "to condor so jobs can be run"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-d", "--directory", type = str, 
                    help = 'specify directory to generate arguments for')
parser.add_argument("-o", "--outfile", type = str, default = "arguments",
                    help = "name of the oufile")
parser.add_argument("-s", "--savedir", type=str, help="name of save directory")
parser.add_argument("-y", "--yaxis", type=str, choices=["exp","rawexp","true"],
                    help = "scale of y axis (exp, rawexp, or true)")
parser.add_argument("-n", "--number", type = int, default = 0,
                    help = "number of files to generate args for, default all")
args = parser.parse_args() 

if __name__ == "__main__":
    
    #get list of files to generate arguments for
    files = glob.glob("%s/ic*.root" % args.directory)
    files.sort()

    #set number of arguments to make
    if args.number:
        filenums = range(args.number)
    else:
        filenums = range(len(files))

    outfiles = []
    arguments = []

    #open file, write user specified arguments
    argfile = open("/home/bdigiorgio/%s.txt" % args.outfile, "w")
    for i in filenums:

        #generate filename, write arguments
        outfiles.append("%sextracted-%s" % (args.savedir, files[i][-15:-5]))
        arguments.append("-d %s -f %d -o %s -y %s" 
                     % (args.directory, filenums[i], outfiles[i], args.yaxis))
        argfile.write(arguments[i])
        argfile.write("\n")

    argfile.close()
