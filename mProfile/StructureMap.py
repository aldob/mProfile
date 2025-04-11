
#!/usr/bin/env python
# ToDo: Use pysam.header to retrieve sorting information, .to_dict()?

from __future__ import division
from argparse import ArgumentParser
import sys
import pysam as ps

def argypargy():
    parser = ArgumentParser(description="StructureMap -i input.sam -o output.sprofile")
    req_args = parser.add_argument_group('Required arguments')
    add_args = parser.add_argument_group('Additional arguments')
    req_args.add_argument("--input", "-i", help="Input alignment file of SAM/BAM/CRAM format MUST be paired end and sorted")
    req_args.add_argument("--output", "-o", help="Output structure profile, referred to as .sprofile files, is a tab delimited table.")
    add_args.add_argument("--insert_max", "-im", help="Maximum size of input DNA fragments, over which deletions are annotated, default=2000", nargs='?', default=2000)
    add_args.add_argument("--fragment_sizes", "-fs", help="Optional output file for the fragment size of correctly mapped reads")
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if args.input is None:
        print("\nINPUT ERROR: Please provide an input file with --input (-i).\n")
        sys.exit()
    if args.output is None:
        print("\nOUTPUT ERROR: Please provide an output file with --output (-o).\n")
        sys.exit()
    if not args.input.lower().endswith(("sam", "cram", "bam")):
        print("\nINPUT ERROR: Please provide an input file with the correct file extension (.sam, .cram, .bam).\n")
        sys.exit()

    try:
        args.insert_max=int(args.insert_max)
    except ValueError:
        print("\n--insert_max is not a valid number, running as default (2000).\n")
        args.large = 2000

    return(args)

args = argypargy()



def alignprocess(input=args.input, output=args.output, insertmax=args.insert_max, sizes=args.fragment_sizes):
    # BAM file flag interpretation dict
    # BAM file must be sorted; the following categories rely on the first read of each pair that is found in the file having smaller coordinate, as this is the only one processed.
    flags = {
        "unmap":[73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137, 77, 141], # 0x4 and/or 0x10, segment and/or next segment unmapped. no other categories have 0x4 or 0x10
        "map":[99, 163], # 0x20 or 0x40 (not both), seq or next seq reverse complemented; 0x2, both segments properly aligned; reported where coordinate of read (163 or 99) is less its reverse complement (83 or 147), i.e. expected orientation. reported insert size is the distance between the far ends of each read, i.e. the actual fragment size.
        "large_insert": [161, 97], # 0x20 or 0x40 (not both), seq or next seq reverse complemented; no 0x2, not properly aligned; reported where coordinate of read (161 or 97) is less its reverse complement (81 or 145), i.e. expected orientation. reported insert size is the distance between the far ends of each read, i.e. the actual fragment size.
        "diverging":[81, 145], # 0x20 or 0x40 (not both), seq or next seq reverse complemented; no 0x2, not properly aligned; reported where coordinate of reverse complement read (81 or 145) is less than other read (161 or 97), i.e. opposite to expected orientation. reported insert size is the gap between the reads, including neither of the actual read lengths.
        "costrand":[67, 131, 115, 179, 65, 129, 113, 177], # 0x20 and 0x40, or neither 0x20 or 0x40, both seq and next seq reverse complemented or not; can include 0x2, each segment properly aligned
        "other_map": [147, 83] # mapped and proprly paired; but reported where coordinate of reverse read (83 or 147) is less than other read (163 or 99). Should not usually see these as properly paired and bam is ordered (occurs when a read starts from 'within the fragemnt' that is seen from the read pair). reported insert size is the overlap between the two reads, i.e. only the part which was sequenced in both directions.
        }

    # BAM file flag strandedness interpretation 
    flagstrands = {
        99:["+", "-"], 147:["-", "+"], 83:["-", "+"], 163:["+", "-"],
        67:["+", "+"], 131:["+", "+"], 115:["-", "-"], 179:["-", "-"],
        81:["-", "+"], 161:["+", "-"], 97:["+", "-"], 145:["-", "+"],
        65:["+", "+"], 129:["+", "+"], 113:["-", "-"], 177:["-", "-"]
        }

    if sizes:
        size_file = open(args.fragment_sizes, 'w')

    past = set()
    outlist = list() # Reads to be written to output file
    rc = 0 # Counts number of aligned reads to normalise to
    error = 0

    with ps.AlignmentFile(input) as seqfil:
        for line in seqfil:
            if line.query_name not in past:
                past.add(line.query_name)
                flag = line.flag

                if flag > 255:  # added to skip supp alignment, not passed filters, duplicates, secondary alignments
                    continue

                elif flag in flags['unmap']:
                    continue

                rc+=1
                chr_x = line.reference_name
                chr_y = line.next_reference_name 
                coord_x = str(line.reference_start)
                coord_y = str(line.next_reference_start)
                try:
                    strands = flagstrands[flag]
                except KeyError as e:
                    raise KeyError(f"Flag {flag} not found in flagstrands dictionary on read number {rc}, {line.query_name}.\nError: {e}")
                insert = line.template_length

                if flag in flags['map']:
                    if sizes:
                        temp = size_file.write(str(abs(line.template_length)) +"\n")
                    if insert > insertmax: # Is a large deletion
                        outlist.append("\t".join([chr_x, chr_y, coord_x, coord_y, strands[0], strands[1], str(insert), "Deletion", str(flag)]) +"\n")

                else:
                    if chr_y != chr_x: 
                        outlist.append("\t".join([chr_x, chr_y, coord_x, coord_y, strands[0], strands[1], "0", "Inter-chromosomal translocation", str(flag)]) +"\n")
                        continue
                
                    elif flag in flags['large_insert']: # opposite strands, relative position and orientation of reads as expected, but too far apart to be properly mapped 
                        if sizes:
                            temp = size_file.write(str(abs(line.template_length)) +"\n")
                        if insert > insertmax: 
                            outlist.append("\t".join([chr_x, chr_y, coord_x, coord_y, strands[0], strands[1], str(insert), "Deletion", str(flag)]) +"\n")
                        else:
                            outlist.append("\t".join([chr_x, chr_y, coord_x, coord_y, strands[0], strands[1], str(insert), "Large insert", str(flag)]) +"\n")

                    elif flag in flags['costrand']: # Same strand, could be classed as inversion or intrachromosomal translocation
                        outlist.append("\t".join([chr_x, chr_y,coord_x,coord_y, strands[0], strands[1], str(insert), "Same strand", str(flag)]) +"\n")

                    elif flag in flags['diverging']: # opposite strands, relative orientation as expected but position of reverse read is less than forward read
                        outlist.append("\t".join([chr_x, chr_y,coord_x,coord_y, strands[0], strands[1], str(insert), "Diverging", str(flag)]) +"\n")

                    elif flag in flags['other_map']: # opposite strands, relative orientation as expected but position of reverse read is less than forward read
                        outlist.append("\t".join([chr_x, chr_y, coord_x, coord_y, strands[0], strands[1], str(insert), "Other", str(flag)]) +"\n")
                        
                    else: # I don't know, probably don't exist
                        outlist.append("\t".join([chr_x, chr_y,coord_x,coord_y, strands[0], strands[1], str(insert), "Unknown", str(flag)]) +"\n")
    if sizes:
        size_file.close()

    with open(args.output, 'w') as outfil:
        temp = outfil.write("Chr FW\tChr RV\tCoord FW\tCoord RV\tStrand FW\tStrand RV\tInsert bp\tAbberation\tFlag, reads="+str(rc)+"\n")
        for line in outlist:
            temp = outfil.write(line)


def main(args=argypargy()):
    alignprocess()

    