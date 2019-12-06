#!/usr/bin/env python
from __future__ import division
import argparse
import sys
import os
import re
try:
    from itertools import izip as zip
except ImportError:
    pass

def argypargy():
    parser = argparse.ArgumentParser(description="TransloCapture -1 <input_read1> -2 <input_read2> -o <output>")
    group1 = parser.add_mutually_exclusive_group(required=True)
    group2 = parser.add_mutually_exclusive_group()
    group1.add_argument("--input", "-i", help="The fastq file you wish to process")
    group1.add_argument("--read1", "-1", help="Read 1 from PE sequencing")
    parser.add_argument("--read2", "-2", help="Read 2 from PE sequencing")
    parser.add_argument("command", help="Read 2 from PE sequencing")
    parser.add_argument("--output", "-o", help="The file you want the mprofile file to be written to", required=True)
    group2.add_argument("--control", "-c", help="The fastq you want to normalise to (e.g. untreated).\nIf unspecified, will not normalise.")
    group2.add_argument("--control1", "-c1", help="Read 1 of the fastq you want to normalise to (e.g. untreated).\nIf unspecified, will not normalise.")
    parser.add_argument("--control2", "-c2", help="Read 2 of the fastq you want to normalise to (e.g. untreated).\nIf unspecified, will not normalise.")
    parser.add_argument("--primers", "-p", help="A .csv file containing the name and foward and reverse primer (reverse complement) sequences for each site to be analysed on each line.")
    parser.add_argument("--preproc", "-pp", help="A .csv file containing the name and foward and reverse primer (reverse complement) sequences for each site to be analysed on each line.", action='store')
    args = parser.parse_args()
    if args.preproc is None:
        if args.input.endswith(".csv"):
            print("Detected translocation matrix input instead of fastq, activating --preproc.")
            args.preproc = True
    if args.preproc is not None:
        if args.control is None and args.control1 is None:
            print("To calculate a differential you must supply a control input to normalise to.")
    return(args)
def numsafe(anum):
    try:
        float(anum)
    except ValueError:
        return(False)
    else:
        return(True)
def TransloCapture(samp_fq="NULL", fq1="NULL", fq2="NULL", site_file="NULL"):
    # First, generate lists of the primer targets and their sequences to identify each site in the fastqs
    primer_names = list()
    fw_primer_list = list()
    rv_primer_list = list()
    with open(site_file) as sites:
        for site in sites:
            primer_names.append(site.split(",")[0])
            fw_primer_list.append(site.split(",")[1])
            rv_primer_list.append(site.split(",")[2])
    # Make an empty dict and fill it with all the possible crossover events 
    samp_dict = {}
    for donor in primer_names:
        for acceptor in primer_names:
            if donor == acceptor:
                samp_dict[donor+"-"+acceptor] = "NA"
            elif donor != acceptor:
                samp_dict[donor+"-"+acceptor] = 0 
    # Loop over each read and identify which primers generated it
    # First identify the forward primer used, then identify the reverse
    # If it is a crossover event, increase the value of that event in the dict by 1
    # If it is a canonical target then increase then increase the readcoutn for that target as this is then used for normalisation
    count = 2
    readcounts = {key:0 for key in primer_names}
    if fq1=="NULL":
        for read in samp_fq:
            count +=1
            if count%4==0:
                for fw, rv, donor in zip(fw_primer_list, rv_primer_list, primer_names):
                    if fw in read[0:len(fw)+5]:
                        if rv in read[-(len(rv)+5):0]:
                            readcounts[donor] += 1
                        elif rv not in read[-(len(rv)+5):0]:
                            for all_rv, acceptor in zip(rv_primer_list, primer_names):
                                if rv != all_rv:
                                    if all_rv in read[-(len(all_rv)+5):0]: 
                                        readcounts[donor] += 1
                                        if str(donor+"-"+acceptor) in samp_dict:
                                            samp_dict[str(donor+"-"+acceptor)] += 1
                    elif rv in read[0:len(fw)+5]:
                        if fw in read[-(len(fw)+5):0]:
                            readcounts[donor] += 1
                        elif fw not in read[-(len(fw)+5):0]:
                            for all_fw, acceptor in zip(fw_primer_list, primer_names):
                                if fw != all_fw:
                                    if all_fw in read[-(len(all_fw)+5):0]: 
                                        readcounts[donor] += 1
                                        if str(donor+"-"+acceptor) in samp_dict:
                                            samp_dict[str(donor+"-"+acceptor)] += 1
    if fq1!="NULL":     
        for read1, read2 in zip(fq1, fq2):
            count +=1
            if count%4==0:
                for fw, rv, donor in zip(fw_primer_list, rv_primer_list, primer_names):
                    if fw in read1[0:len(fw)+5]:
                        if rv in read2[0:len(rv)+5]:
                            readcounts[donor] += 1
                        elif rv not in read2[0:len(fw)+5]:
                            for all_rv, acceptor in zip(rv_primer_list, primer_names):
                                if rv != all_rv:
                                    if all_rv in read2[0:len(all_rv)+5]: 
                                        readcounts[donor] += 1
                                        if str(donor+"-"+acceptor) in samp_dict:
                                            samp_dict[str(donor+"-"+acceptor)] += 1
                    elif rv in read1[0:len(fw)+5]:
                        if fw in read2[0:len(fw)+5]:
                            readcounts[donor] += 1
                        elif fw not in read2[0:len(fw)+5]:
                            for all_fw, acceptor in zip(fw_primer_list, primer_names):
                                if fw != all_fw:
                                    if all_fw in read2[0:len(all_fw)+5]: 
                                        readcounts[donor] += 1
                                        if str(donor+"-"+acceptor) in samp_dict:
                                            samp_dict[str(donor+"-"+acceptor)] += 1
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
def main_process(args):
    if args.preproc is None:
        with open(args.primers) as sites:
            primer_names = [site.split(",")[0] for site in sites]
        print("Quantifying translocation events")
        if args.control is not None or args.control1 is not None:
            if args.read1 is not None:
                with open(args.control1) as ctrl1, open(args.control2) as ctrl2, open(args.read1) as treat1, open(args.read2) as treat2:
                    ctrl_dict = TransloCapture(fq1=ctrl1, fq2=ctrl2, site_file=args.primers)
                    treat_dict = TransloCapture(fq1=treat1, fq2=treat2, site_file=args.primers)
                diff_dict = dict_diff(ctrl_dict, treat_dict)
                translomap_write(tc_dict=diff_dict, tc_output=args.output, names=primer_names)
            elif args.input is not None:
                with open(args.control) as ctrl, open(args.input) as treat:
                    ctrl_dict = TransloCapture(samp_fq=ctrl, site_file=args.primers)
                    treat_dict = TransloCapture(samp_fq=treat, site_file=args.primers)
                diff_dict = dict_diff(ctrl_dict, treat_dict)
                translomap_write(tc_dict=diff_dict, tc_output=args.output, names=primer_names)
        elif args.control is None and args.control1 is None:
            if args.read1 is not None:
                with open(args.read1) as treat1, open(args.read2) as treat2:
                    treat_dict = TransloCapture(fq1=treat1, fq2=treat2, site_file=args.primers)
                translomap_write(tc_dict=treat_dict, tc_output=args.output, names=primer_names)
            elif args.input is not None:
                with open(args.input) as treat:
                    treat_dict = TransloCapture(samp_fq=treat, site_file=args.primers)
                translomap_write(tc_dict=treat_dict, tc_output=args.output, names=primer_names)
    elif args.preproc is not None:
        print("Quantifying translocation events")
        with open(args.control) as ctrl, open(args.input) as treat, open(args.output, "w") as outputfile:
            for line1, line2 in zip(ctrl, treat):
                for val1, val2 in zip(line1.split(","), line2.split(",")):
                    if numsafe(val1):
                        outputfile.write(str(float(val2)-float(val1))+",")
                    else:
                        outputfile.write(val1.strip("\n")+",")
                outputfile.write("\n")

