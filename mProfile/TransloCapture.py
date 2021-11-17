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
    add_args.add_argument("--translocated", "-t", help="Fastq file to write translocated sequences to.\n If unspecified, will not write")
    add_args.add_argument("--translocated1", "-t1", help="Fastq file to write read1 of translocated sequences to.\n If unspecified, will not write.")
    add_args.add_argument("--translocated2", "-t2", help="Fastq file to write read2 of translocated sequences to.\n If unspecified, will not write.")
    add_args.add_argument("--fastqout", "-fq", help="Fastq file to write non-translocated sequences to.\n If unspecified, will not write")
    add_args.add_argument("--fastqout1", "-fq1", help="Fastq file to write read1 of non-translocated sequences to.\n If unspecified, will not write.")
    add_args.add_argument("--fastqout2", "-fq2", help="Fastq file to write read2 of non-translocated sequences to.\n If unspecified, will not write.")
    add_args.add_argument("--sensitivity", "-s", help="Pads window at start of reads for identify the primer used for amplification.\nLarger numbers increase detection, but reduce specificity. default=2, max=10.", default=2)
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
    if args.translocated1 is not None and args.read1 is None or args.translocated is not None and args.input is None: # Don't crossover the SR and PE options, fastq output edition
        print("\nTransloCapture ERROR: To write translocated reads, for single-read sequencing use --input (-i) and --translocated (-t), for PE sequencinig use --read1 --read2 (-1 -2) and --translocated1 --translocated2 (-t1 -t2). Don't mix and match.\n")
        sys.exit()
    if args.translocated1 is not None and args.translocated2 is None or args.translocated1 is None and args.translocated2 is not None: # Need both reads specified for PE, fastq output edition
        print("\nTransloCapture ERROR: If --translocated1 (-t1) or --translocated2 (-t2) are specified you must also supply the other. For single read sequencing use --input (-i) instead.\n")
        sys.exit()
    if args.fastqout1 is not None and args.read1 is None or args.fastqout is not None and args.input is None: # Don't crossover the SR and PE options, fastq output edition
        print("\nTransloCapture ERROR: To write reads, for single-read sequencing use --input (-i) and --fastqout (-fq), for PE sequencinig use --read1 --read2 (-1 -2) and --fastqout1 --fastqout2 (-fq1 -fq2). Don't mix and match.\n")
        sys.exit()
    if args.fastqout1 is not None and args.fastqout2 is None or args.fastqout1 is None and args.fastqout2 is not None: # Need both reads specified for PE, fastq output edition
        print("\nTransloCapture ERROR: If --fastqout1 (-fq1) or --fastqout2 (-fq2) are specified you must also supply the other. For single read sequencing use --input (-i) instead.\n")
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
    if args.preproc == True and args.translocated is not None or args.preproc == True and args.translocated1 is not None: # Need to use SR options for preproc
        print("\nTransloCapture ERROR: --preproc (-pp) cannot be used alongside --translocated (-t) or --translocated1/2 (-t1/2) because no fastq is being analysed with preproc.\n")
        sys.exit()
    if args.preproc == True and args.fastqout is not None or args.preproc == True and args.fastqout2 is not None: # Need to use SR options for preproc
        print("\nTransloCapture ERROR: --preproc (-pp) cannot be used alongside --fastqout (-fq) or --fastqout1/2 (-fq1/2) because no fastq is being analysed with preproc.\n")
        sys.exit()  
    if args.primers is None and args.preproc is None: # Need primer sequences unless using preproc
        print("\nTransloCapture ERROR: --primers (-p) is needed to identify translocated sequences in fastq files.\n")
        sys.exit()
    try:
        args.sensitivity = int(args.sensitivity)
    except ValueError:
        print("\nTransloCapture ERROR: --sensitivity (-s) must be a number from 0-10.\n")
        sys.exit()
    if args.sensitivity < 0: # Sensitivity below zero error
        print("\nTransloCapture ERROR: --sensitivity (-s) must be a number from 0-10.\n")
        sys.exit()
    elif 10 >= args.sensitivity > 5: # Low specificity warning for sensitivity > 5
        print("\nTransloCapture WARNING: --sensitivity (-s) greater than 5 greatly reduces specificity of the analysis.\n")
    elif args.sensitivity > 10: # Sensitivity limit warning
        print("\nTransloCapture WARNING: --sensitivity (-s) must not exceed 10, TransloCapture will run with a sensitivity of 10.\n")
        args.sensitivity = 10
    return(args)
def numsafe(anum):
    try:
        float(anum)
        return(True)
    except ValueError:
        return(False)
def rev_comp(seq):
    newseq = seq.replace("A", "t")
    newseq = newseq.replace("T", "a")
    newseq = newseq.replace("G", "c")
    newseq = newseq.replace("C", "g")
    return(newseq.upper()[::-1])
def TransloCapture(argues):
    fastq=argues[0]
    fastq1=argues[1]
    fastq2=argues[2]
    site_file=argues[3]
    translocated=argues[4]
    translocated1=argues[5]
    translocated2=argues[6]
    regular=argues[7]
    regular1=argues[8]
    regular2=argues[9]
    sens=argues[10]
    # First, generate lists of the primer targets and their sequences to identify each site in the fastqs
    primer_names = list()
    fw_primer_list = list()
    rv_primer_list = list()
    fw_lens = list()
    rv_lens = list()
    with open(site_file) as sites:
        for site in sites:
            primer_names.append(site.split(",")[0])
            fw_primer_list.append(site.split(",")[1].upper()[1:])
            rv_primer_list.append(sub("\n|\r", "", site.split(",")[2]).upper()[1:])
            fw_lens.append(len(site.split(",")[1])+sens)
            rv_lens.append(len(sub("\n|\r", "", site.split(",")[2]))+sens)
    # Make an empty dict and fill it with all the possible crossover events 
    samp_dict = {}
    for donor in primer_names:
        for acceptor in primer_names:
            samp_dict[donor+"-"+acceptor] = 0 
    # Loop over each read and identify which primers generated it
    # First identify the forward primer used, then identify the reverse
    # If it is a crossover event, increase the value of that event in the dict by 1
    # If it is a canonical target then increase then increase the readcoutn for that target as this is then used for normalisation
    count = 0
    readcounts = {key:0 for key in primer_names}
    lines1=list()
    lines2=list()
    if translocated is not None:
        output = open(translocated, "w")
    elif translocated1 is not None:
        output1 = open(translocated1, "w")
        output2 = open(translocated2, "w")
    if regular is not None:
        regout = open(regular, "w")
    elif regular1 is not None:
        regout1 = open(regular1, "w")
        regout2 = open(regular2, "w")
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
                    for rv, fw, rv_l, fw_l, donor in zip(rv_primer_list, fw_primer_list, rv_lens, fw_lens, primer_names):
                        if foundcheck:
                            break
                        else:
                            if fw in read[:fw_l]:
                                if rev_comp(rv) in read[-(rv_l):]:
                                    readcounts[donor] += 1
                                    break
                                else:
                                    for all_rv, all_fw, arv_l, afw_l, acceptor in zip(rv_primer_list, fw_primer_list, rv_lens, fw_lens, primer_names):
                                        if donor != acceptor:
                                            if rev_comp(all_rv) in read[-(arv_l):]: 
                                                readcounts[donor] += 1
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if translocated is not None:
                                                    set1[0] = set1[0].strip("\n") + " " + donor + "_fw:" + acceptor + "_rv"
                                                    output.write("\n".join(set1)+"\n")
                                                break
                                            elif rev_comp(all_fw) in read[-(afw_l):]:
                                                readcounts[donor] += 1
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if translocated is not None:
                                                    set1[0] = set1[0].strip("\n") + " " + donor + "_fw:" + acceptor + "_fw"
                                                    output.write("\n".join(set1)+"\n")
                                                break

                            elif rv in read[:rv_l]:
                                if rev_comp(fw) in read[-(fw_l):]:
                                    readcounts[donor] += 1
                                    break
                                else:
                                    for all_rv, all_fw, arv_l, afw_l, acceptor in zip(rv_primer_list, fw_primer_list, rv_lens, fw_lens, primer_names):
                                        if donor != acceptor:
                                            if rev_comp(all_fw) in read[-(afw_l):]: 
                                                readcounts[donor] += 1
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if translocated is not None:
                                                    set1[0] = set1[0].strip("\n") + " " + donor + "_rv:" + acceptor + "_rv"
                                                    output.write("\n".join(set1)+"\n")
                                                break
                                            elif rev_comp(all_rv) in read[-(arv_l):]: 
                                                readcounts[donor] += 1
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if translocated is not None:
                                                    set1[0] = set1[0].strip("\n") + " " + donor + "_rv:" + acceptor + "_fw"
                                                    output.write("\n".join(set1)+"\n")
                                                break
                    if not foundcheck and regular is not None:
                        regout.write("\n".join(set1)+"\n") 

    else:
        with open(fastq1) as fq1, open(fastq2) as fq2:
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
                    for rv, fw, rv_l, fw_l, donor in zip(rv_primer_list, fw_primer_list, rv_lens, fw_lens, primer_names):
                        if foundcheck:
                            break
                        else:
                            if fw in read1[:fw_l]:
                                if rv in read2[:rv_l]:
                                    readcounts[donor] += 1                                    
                                    break
                                else:
                                    for all_rv, all_fw, arv_l, afw_l, acceptor in zip(rv_primer_list, fw_primer_list, rv_lens, fw_lens, primer_names):
                                        if donor != acceptor:
                                            if all_rv in read2[:arv_l]:
                                                readcounts[donor] += 1
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if translocated1 is not None:
                                                    set1[0] = set1[0].strip("\n") + " " + donor + "_fw:" + acceptor + "_rv"
                                                    set2[0] = set2[0].strip("\n") + " " + donor + "_fw:" + acceptor + "_rv"
                                                    output1.write("\n".join(set1)+"\n")
                                                    output2.write("\n".join(set2)+"\n")
                                                break
                                            elif all_fw in read2[:afw_l]:
                                                readcounts[donor] += 1
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if translocated1 is not None:
                                                    set1[0] = set1[0].strip("\n") + " " + donor + "_fw:" + acceptor + "_fw"
                                                    set2[0] = set2[0].strip("\n") + " " + donor + "_fw:" + acceptor + "_fw"
                                                    output1.write("\n".join(set1)+"\n")
                                                    output2.write("\n".join(set2)+"\n")
                                                break                                                

                            elif rv in read1[:rv_l]:
                                if fw in read2[:fw_l]:
                                    readcounts[donor] += 1                                      
                                    break
                                else:
                                    for all_rv, all_fw, arv_l, afw_l, acceptor in zip(rv_primer_list, fw_primer_list, rv_lens, fw_lens, primer_names):
                                        if donor != acceptor:
                                            if all_rv in read2[:arv_l]:
                                                readcounts[donor] += 1
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if translocated1 is not None:
                                                    set1[0] = set1[0].strip("\n") + " " + donor + "_rv:" + acceptor + "_rv"
                                                    set2[0] = set2[0].strip("\n") + " " + donor + "_rv:" + acceptor + "_rv"
                                                    output1.write("\n".join(set1)+"\n")
                                                    output2.write("\n".join(set2)+"\n")
                                                break
                                            elif all_fw in read2[:afw_l]:
                                                readcounts[donor] += 1
                                                foundcheck=True
                                                samp_dict[str(donor+"-"+acceptor)] += 1
                                                if translocated1 is not None:
                                                    set1[0] = set1[0].strip("\n") + " " + donor + "_rv:" + acceptor + "_fw"
                                                    set2[0] = set2[0].strip("\n") + " " + donor + "_rv:" + acceptor + "_fw"
                                                    output1.write("\n".join(set1)+"\n")
                                                    output2.write("\n".join(set2)+"\n")
                                                break
                    if not foundcheck and regular1 is not None:
                        regout1.write("\n".join(set1)+"\n")
                        regout2.write("\n".join(set2)+"\n")                                            

    if translocated is not None:
        output.close()
    elif translocated1 is not None:
        output1.close()
        output2.close()
    if regular is not None:
        regout.close()
    elif regular1 is not None:
        regout1.close()
        regout2.close()
    # Normalise all counts to readcount of the canonical donor target
    samp_dict_norm = {}
    for donor, acceptor in zip(primer_names, primer_names):
        if donor != acceptor:
            if readcounts[donor] + readcounts[acceptor] > 0:
                fw_rate = samp_dict[donor+"-"+acceptor]
                rv_rate = samp_dict[acceptor+"-"+donor]
                samp_dict_norm[donor+"-"+acceptor] = (fw_rate/(readcounts[donor] + readcounts[acceptor]))*100
                samp_dict_norm[acceptor+"-"+donor] = (rv_rate/(readcounts[donor] + readcounts[acceptor]))*100
            else:
                samp_dict_norm[donor+"-"+acceptor] = 0
                samp_dict_norm[acceptor+"-"+donor] = 0
        else:
            samp_dict_norm[donor+"-"+acceptor] = "NA"
            samp_dict_norm[acceptor+"-"+donor] = "NA"
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
            outputfile.write(','.join([str(tc_dict[str(donor+"-"+acceptor)]) for donor in names])+"\n")
def main(args=argypargy()):
    if args.preproc == False:
        with open(args.primers) as sites:
            primer_names = [site.split(",")[0] for site in sites]
        if args.control is not None or args.control1 is not None:
            if args.quiet == False:
                print("\nIdentifying translocated sequences in treated and control.\n")
            p=Pool(2)
            both_dicts=p.map(TransloCapture, [[args.control, args.control1, args.control2, args.primers, args.translocated, None, None, args.fastqout, None, None, args.sensitivity], [args.input, args.read1, args.read2, args.primers, args.translocated, args.translocated1, args.translocated2, args.fastqout, args.fastqout1, args.fastqout2, args.sensitivity]])
            p.close()
            if args.quiet == False:
                print("\nQuantifying differential and writing output file.\n")
            diff_dict = dict_diff(both_dicts[0], both_dicts[1])
            translomap_write(tc_dict=diff_dict, tc_output=args.output, names=primer_names)
        else:
            if args.quiet == False:
                print("\nIdentifying translocated sequences.\n")
            treat_dict = TransloCapture([args.input, args.read1, args.read2, args.primers, args.translocated, args.translocated1, args.translocated2, args.fastqout, args.fastqout1, args.fastqout2, args.sensitivity])
            translomap_write(tc_dict=treat_dict, tc_output=args.output, names=primer_names)
    elif args.preproc == True:
        if args.quiet == False:
            print("\nQuantifying differential and writing output file.\n")
        with open(args.control) as ctrl, open(args.input) as treat, open(args.output, "w") as outputfile:
            outputfile.write(",")
            for line1, line2 in zip(ctrl, treat):
                lengt=len(line1.split(","))
                count=0
                for val1, val2 in zip(line1.split(","), line2.split(",")):
                    count+=1
                    if val1 or val2:
                        if count == lengt:
                            extra=""
                        else:
                            extra=","
                        if numsafe(val1):
                            outputfile.write(str(float(val2)-float(val1))+extra)
                        else:
                            outputfile.write(val1.strip("\n")+extra)
                outputfile.write("\n")


