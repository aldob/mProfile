#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import sys
from re import sub
from multiprocessing import Pool
try:
    from itertools import izip as zip
except ImportError:
    pass



def argypargy():
    parser = ArgumentParser(description='TransloCapture -1 input_read1.fastq -2 input_read2.fastq -o output.csv')
    req_args = parser.add_argument_group('Required arguments')
    add_args = parser.add_argument_group('Additional arguments')
    req_args.add_argument("--input", "-i", help="Input fastq file for SR sequencing")
    req_args.add_argument("--read1", "-1", help="Fastq read 1 from PE sequencing")
    req_args.add_argument("--read2", "-2", help="Fastq read 2 from PE sequencing")
    req_args.add_argument("--output", "-o", help="Output file to write to, format is csv")
    add_args.add_argument("--control", "-c", help="The fastq you want to normalise to (e.g. untreated).\nIf unspecified, will not normalise.")
    add_args.add_argument("--control1", "-c1", help="Read 1 of the fastq you want to normalise to (e.g. untreated).\nIf unspecified, will not normalise.")
    add_args.add_argument("--control2", "-c2", help="Read 2 of the fastq you want to normalise to (e.g. untreated).\nIf unspecified, will not normalise.")
    req_args.add_argument("--primers", "-p", help="A 3 column .csv file of the name, foward primer sequence and reverse primer sequence (reverse complement) for each site to be analysed.")
    add_args.add_argument("--preproc", "-pp", help="If specified, --input (-i) and --control (-c) must be already quantified TransloCapture matrices.\nOutput will be a new matrix that is the differential of input-control.", action='store_true')
    add_args.add_argument("--fastqout", "-fo", help="Fastq file to write translocated sequences to.\n If unspecified, will not write")
    add_args.add_argument("--fastqout1", "-fo1", help="Fastq file to write read1 of translocated sequences to.\n If unspecified, will not write.")
    add_args.add_argument("--fastqout2", "-fo2", help="Fastq file to write read2 of translocated sequences to.\n If unspecified, will not write.")
    add_args.add_argument("--quiet", "-q", help="Removes all messages.", action='store_true')
    args = parser.parse_args()
    if len(sys.argv)==1: # If no arguments are given, print help information.
        parser.print_help()
        sys.exit()
    if args.input is None and args.read1 is None: # Input file is required
        print("\nTransloCapture ERROR: Please provide an input file with --input (-i) or with --read1 and --read2 (-1 -2) for PE seqeuencing.\n")
        sys.exit()
    if args.output is None: # Output file is required
        print("\nTransloCapture ERROR: Please provide an output file with --output (-o).\n")
        sys.exit()
    if args.input is not None and args.read1 is not None: # Don't crossover the SR and PE options
        print("\nTransloCapture ERROR: --input (-i) is for single-read sequencing and --read1 --read2 (-1 -2) are for PE seqeuencing. They cannot be used together.\n")
        sys.exit()
    if args.control is not None and args.control1 is not None: # Don't crossover the SR and PE options, control edition
        print("\nTransloCapture ERROR: --control (-c) is for single-read sequencing and --control1 --control2 (-c1 -c2) are for PE seqeuencing. They cannot be used together.\n")
        sys.exit()                
    if args.read1 is not None and args.read2 is None or args.read1 is None and args.read2 is not None: # Need both reads specified for PE
        print("\nTransloCapture ERROR: If --read1 (-1) or --read2 (-2) are specified you must also supply the other. For single read sequencing use --input (-i) instead.\n")
        sys.exit()
    if args.control1 is not None and args.control2 is None or args.control1 is None and args.control2 is not None: # Need both reads specified for PE, control edition
        print("\nTransloCapture ERROR: If --control1 (-c1) or --control2 (-c2) are specified you must also supply the other. For single read sequencing use --control (-c) instead.\n")
        sys.exit()
    if args.read1 is not None and args.control is not None or args.input is not None and args.control1 is not None: # Don't crossover the SR and PE options, mix match addition
        print("\nTransloCapture ERROR: --read1 (-1)/--control1 (-c1) must be used alongside each other, not alongside --control (-c)/--input (-i).\n")
        sys.exit()
    if args.fastqout1 is not None and args.read1 is None or args.fastqout is not None and args.input is None: # Don't crossover the SR and PE options, fastq output edition
        print("\nTransloCapture ERROR: To write translocated reads, for single-read sequencing use --input (-i) and --fastqout (-fo), for PE sequencinig use --read1 --read2 (-1 -2) and --fastqout1 --fastqout2 (-fo1 -fo2). Don't mix and match.\n")
        sys.exit()
    if args.fastqout1 is not None and args.fastqout2 is None or args.fastqout1 is None and args.fastqout2 is not None: # Need both reads specified for PE, fastq output edition
        print("\nTransloCapture ERROR: If --fastqout1 (-fo1) or --fastqout2 (-fo2) are specified you must also supply the other. For single read sequencing use --input (-i) instead.\n")
        sys.exit()
    if args.preproc == False and args.input is not None and args.input.endswith(".csv"): # .csv input suggests they want --preproc
        print("\nDetected translocation matrix.csv input instead of fastq, activating --preproc (-pp).\n")
        args.preproc = True
    if args.preproc == True and args.control is None: # Preproc needs a control
        print("\nTransloCapture ERROR: --preproc (-pp) also needs --control (-c) to calculate a differential to the --input (-i) sample.\n")
        sys.exit()
    if args.preproc == True and args.read1 is not None: # Need to use SR options for preproc
        print("\nTransloCapture ERROR: --read1/2 (-1/2) and --control1/2 (-c1/2) are for paired fastq files.\nPlease use --input (-i) and --control (-c) with --preproc (-pp).\n")
        sys.exit()        
    if args.preproc == True and args.fastqout is not None or args.preproc == True and args.fastqout1 is not None: # Need to use SR options for preproc
        print("\nTransloCapture ERROR: --preproc (-pp) cannot be used alongside --fastqout (-fo) or --fastqout1/2 (-fo1/2) because no fastq is being analysed with preproc.\n")
        sys.exit()        
    if args.primers is None and args.preproc is None: # Need primer sequences unless using preproc
        print("\nTransloCapture ERROR: --primers (-p) is needed to identify translocated sequences in fastq files.\n")
        sys.exit()
    return(args)
def numsafe(anum):
    try:
        float(anum)
        return(True)
    except ValueError:
        return(False)
def TransloCapture(argues):
    fastq=argues[0]
    fastq1=argues[1]
    fastq2=argues[2]
    site_file=argues[3]
    out_fastq=argues[4]
    out_fastq1=argues[5]
    out_fastq2=argues[6]
    # First, generate lists of the primer targets and their sequences to identify each site in the fastqs
    primer_names = list()
    fw_primer_list = list()
    rv_primer_list = list()
    with open(site_file) as sites:
        for site in sites:
            primer_names.append(site.split(",")[0])
            fw_primer_list.append(site.split(",")[1].upper())
            rv_primer_list.append(sub("\n|\r", "", site.split(",")[2]).upper())
    # Make an empty dict and fill it with all the possible crossover events 
    samp_dict = {}
    for donor in primer_names:
        for acceptor in primer_names:
            if donor != acceptor:
                samp_dict[donor+"-"+acceptor] = 0 
            else:
                samp_dict[donor+"-"+acceptor] = "NA"
    # Loop over each read and identify which primers generated it
    # First identify the forward primer used, then identify the reverse
    # If it is a crossover event, increase the value of that event in the dict by 1
    # If it is a canonical target then increase then increase the readcoutn for that target as this is then used for normalisation
    count = 0
    readcounts = {key:0 for key in primer_names}
    lines1=list()
    lines2=list()
    if out_fastq is not None:
        output = open(out_fastq, "w")
    elif out_fastq1 is not None:
        output1 = open(out_fastq1, "w")
        output2 = open(out_fastq2, "w")
    if fastq1 is None:
        with open(fastq) as samp_fq:
            for rd in samp_fq:
                count+=1
                lines1.append(rd.strip("\n"))
                if count%4==0:
                    set1=lines1
                    read=set1[1]
                    count=0
                    lines1=list()
                    foundcheck=False
                    for fw, rv, donor in zip(fw_primer_list, rv_primer_list, primer_names):
                        if fw in read[0:len(fw)+1]:
                            if rv in read[-(len(rv)+1):0]:
                                readcounts[donor] += 1
                            else:
                                for all_rv, acceptor in zip(rv_primer_list, primer_names):
                                    if rv != all_rv and not foundcheck:
                                        if all_rv in read[-(len(all_rv)+1):0]: 
                                            readcounts[donor] += 1
                                            if str(donor+"-"+acceptor) in samp_dict:
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if out_fastq is not None:
                                                    output.write("\n".join(set1)+"\n")
                                for all_fw, acceptor in zip(fw_primer_list, primer_names):
                                    if all_fw in read[-(len(all_rv)+1):0] and not foundcheck:
                                        readcounts[donor] += 1
                                        if str(donor+"-"+acceptor) in samp_dict:
                                            foundcheck=True
                                            samp_dict[str(donor+"-"+acceptor)] += 1
                                            if out_fastq is not None:
                                                output.write("\n".join(set1)+"\n")
                        elif rv in read[0:len(fw)+1] and not foundcheck:
                            if fw in read[-(len(fw)+1):0]:
                                readcounts[donor] += 1
                            else:
                                for all_fw, acceptor in zip(fw_primer_list, primer_names):
                                    if fw != all_fw and not foundcheck:
                                        if all_fw in read[-(len(all_fw)+1):0]: 
                                            readcounts[donor] += 1
                                            if str(donor+"-"+acceptor) in samp_dict:
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if out_fastq is not None:
                                                    output.write("\n".join(set1)+"\n")
                                for all_rv, acceptor in zip(rv_primer_list, primer_names):
                                    if all_rv in read[-(len(all_fw)+1):0] and not foundcheck:
                                        readcounts[donor] += 1
                                        if str(donor+"-"+acceptor) in samp_dict:
                                            foundcheck=True
                                            samp_dict[str(donor+"-"+acceptor)] += 1
                                            if out_fastq is not None:
                                                output.write("\n".join(set1)+"\n")
    else:
        with open(fastq1) as fq1, open(fastq2) as fq2:
            dubcount=0
            for rd1, rd2 in zip(fq1, fq2):
                count+=1
                lines1.append(rd1.strip("\n"))
                lines2.append(rd2.strip("\n"))
                if count%4==0:
                    set1=lines1
                    set2=lines2
                    read1=set1[1]
                    read2=set2[1]
                    count=0
                    lines1=list()
                    lines2=list()
                    foundcheck=False
                    for fw, rv, donor in zip(fw_primer_list, rv_primer_list, primer_names):
                        if fw in read1[0:len(fw)+1]:
                            if rv in read2[0:len(rv)+1]:
                                readcounts[donor] += 1
                            else:
                                for all_rv, acceptor in zip(rv_primer_list, primer_names):
                                    if rv != all_rv and not foundcheck:
                                        if all_rv in read2[0:len(all_rv)+1]:
                                            readcounts[donor] += 1
                                            if str(donor+"-"+acceptor) in samp_dict:
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if out_fastq1 is not None:
                                                    output1.write("\n".join(set1)+"\n")
                                                    output2.write("\n".join(set2)+"\n")
                                for all_fw, acceptor in zip(fw_primer_list, primer_names):
                                    if all_fw in read2[0:len(all_fw)+1] and not foundcheck:
                                        readcounts[donor] += 1
                                        if str(donor+"-"+acceptor) in samp_dict:
                                            foundcheck=True
                                            samp_dict[str(donor+"-"+acceptor)] += 1
                                            if out_fastq1 is not None:
                                                output1.write("\n".join(set1)+"\n")
                                                output2.write("\n".join(set2)+"\n")

                        elif rv in read1[0:len(fw)+1] and not foundcheck:
                            if fw in read2[0:len(fw)+1]:
                                readcounts[donor] += 1
                            else:
                                for all_fw, acceptor in zip(fw_primer_list, primer_names):
                                    if fw != all_fw and not foundcheck:
                                        if all_fw in read2[0:len(all_fw)+1]: 
                                            readcounts[donor] += 1
                                            if str(donor+"-"+acceptor) in samp_dict:
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if out_fastq1 is not None:
                                                    output1.write("\n".join(set1)+"\n")
                                                    output2.write("\n".join(set2)+"\n")
                                for all_rv, acceptor in zip(rv_primer_list, primer_names):
                                    if all_rv in read2[0:len(all_rv)+1] and not foundcheck:
                                        readcounts[donor] += 1
                                        if str(donor+"-"+acceptor) in samp_dict:
                                            foundcheck=True
                                            samp_dict[str(donor+"-"+acceptor)] += 1
                                            if out_fastq1 is not None:
                                                output1.write("\n".join(set1)+"\n")
                                                output2.write("\n".join(set2)+"\n")

    if out_fastq is not None:
        output.close()
    elif out_fastq1 is not None:
        output1.close()
        output2.close()
    # Normalise all counts to readcount of the canonical donor target
    samp_dict_norm = {}
    for key,val in samp_dict.items():
        if val != "NA":
            if readcounts[key.split("-")[0]] > 0:
                samp_dict_norm[key]=float((val/(readcounts[key.split("-")[0]]))*100)
            else:
                samp_dict_norm[key]=0
        else:
            samp_dict_norm[key]="NA"
    return(samp_dict_norm)
def dict_diff(ctrl_dict, treat_dict):
    diff_dict = {}
    for (key1,val1), (key2,val2) in zip(sorted(ctrl_dict.items()), sorted(treat_dict.items())):
        if val1 != "NA":
            diff_dict[key2] = float(val2)-float(val1)
        elif val1 == "NA":
            diff_dict[key2] = val2
    return(diff_dict)
def translomap_write(tc_dict="", tc_output="", names=""):
    with open(tc_output, 'w') as outputfile:
        outputfile.write(str(","+','.join(names)+"\n"))
        for acceptor in names:
            outputfile.write(str(acceptor+","))
            for donor in names:
                outputfile.write(str(tc_dict[str(donor+"-"+acceptor)])+",")
            outputfile.write("\n")
def main(args=argypargy()):
    if args.preproc == False:
        with open(args.primers) as sites:
            primer_names = [site.split(",")[0] for site in sites]
        if args.control is not None or args.control1 is not None:
            if args.quiet == False:
                print("\nIdentifying translocated sequences in treated and control.\n")
            p=Pool(2)
            both_dicts=p.map(TransloCapture, [[args.control, args.control1, args.control2, args.primers, args.fastqout, None, None], [args.input, args.read1, args.read2, args.primers, args.fastqout, args.fastqout1, args.fastqout2]])
            p.close()
            if args.quiet == False:
                print("\nQuantifying differential and writing output file.\n")
            diff_dict = dict_diff(both_dicts[0], both_dicts[1])
            translomap_write(tc_dict=diff_dict, tc_output=args.output, names=primer_names)
        else:
            if args.quiet == False:
                print("\nIdentifying translocated sequences.\n")
            treat_dict = TransloCapture([args.input, args.read1, args.read2, args.primers, args.fastqout, args.fastqout1, args.fastqout2])
            translomap_write(tc_dict=treat_dict, tc_output=args.output, names=primer_names)
    elif args.preproc == True:
        if args.quiet == False:
            print("\nQuantifying differential and writing output file.\n")
        with open(args.control) as ctrl, open(args.input) as treat, open(args.output, "w") as outputfile:
            for line1, line2 in zip(ctrl, treat):
                for val1, val2 in zip(line1.split(","), line2.split(",")):
                    if numsafe(val1):
                        outputfile.write(str(float(val2)-float(val1))+",")
                    else:
                        outputfile.write(val1.strip("\n")+",")
                outputfile.write("\n")

