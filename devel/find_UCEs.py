#!/usr/bin/python

import sys,warnings
import argparse

VERSION="1.0"
DESCRIPTION='''
Find UCE candidates in Mummer output
version: v.%s


###############
# Disclaimer: #

Tested in Python 2.7.12
Requires Biopython

The script is currently in the beta phase. Any feedback is much appreciated.

##############

''' %VERSION


def reformat_mummer(mummer_out_file):
    
    outlist = []
    temp_list = []
    
    for line in open(mummer_out_file, 'r'):
        line=line.strip()
        if line.startswith(">"):
            temp_list = [line.replace("> ","")]
        else:
            temp_list.append('\t'.join(line.split()))
            if len(temp_list) == 3:
                outlist.append("\t".join(temp_list))
                temp_list = [temp_list[0]]
            
    return outlist
            
            
def output_mummer_tabular(outfile, result_list):
    
    fh = open(outfile,'w')
    for res in result_list:
        fh.write(res+'\n')
        

def get_ctg_lengths(fasta):
    
    from Bio import SeqIO
    lengths = {}

    for seq in SeqIO.parse(open(fasta,'r'),'fasta'):
        lengths[seq.id]=len(seq.seq)
        
    return lengths


def evaluate_dist_from_end(dist_from_edge, pos, mode, total_length=0):
    
    import sys
    
    ok = True
    
#    dist_from_edge = int(dist_from_edge)
#    pos=int(p)
    if mode == 'start':
        if not pos >= (dist_from_edge): #if match is not far enough from the start of the sequence
            ok = False
        
    elif mode == 'end':
     
        if not total_length > 0:
            sys.exit('need to know the length of the contig')
            
        if not (total_length-pos) >= dist_from_edge:
            ok = False
        
    if ok:
        return True
    else:
        return False

def reformat_match(ref_id, ref_start_pos, baselist):
    
    #ref_id is the unique id of the reference contig
    #ref_start_pos is the coordinate of the start position on the reference
    #baselist is a list containing [query_id, query_start_pos, match_length, sequence]
    #returns a dictionary of format:
    #{'ref':{'id':'xyz', 'start': '123', 'end': '234', 'or':'+', 'seq':'AGCT'}, 
    #'que':{'id':'xyz', 'start': '123', 'end': '234', 'or':'-', 'seq':'AGCT'}}
    
    temp_dict = {'ref':{}, 'que':{}}
    
    temp_dict['ref']['id'] = ref_id
    temp_dict['ref']['start'] = ref_start_pos
    temp_dict['ref']['end'] = ref_start_pos + int(baselist[2])
    temp_dict['ref']['seq'] = baselist[3].upper()
    temp_dict['ref']['or'] = '+'
    temp_dict['que']['id'] = baselist[0]
    temp_dict['que']['start'] = int(baselist[1])
    temp_dict['que']['end'] = int(baselist[1]) + int(baselist[2])
    temp_dict['que']['seq'] = baselist[3].upper()
    temp_dict['que']['or'] = '+'
                
    if baselist[0].endswith('Reverse'):
        temp_dict['que']['id'] = temp_dict['que']['id'].replace('_Reverse','')
        temp_dict['que']['or'] = '-'

    return temp_dict


def write_binned_match_details(outfile, master_dict):
    
    order = ['id','start','end','or','seq']
        
    ## Write out coordinates of candidates
    outfh = open(outfile,'w')
    for group in master_dict.keys():
        for uce in master_dict[group]:
            linestart = "%s\t%s\t" %(uce,group)
            if isinstance(master_dict[group][uce], list):
#                print uce,len(binned[group][uce])
                for rep in master_dict[group][uce]:
                    temp = []
                    for target in ['ref','que']:
                        for key in order:
                            temp.append(str(rep[target][key]))
                    line=linestart+"\t".join(temp)
                    outfh.write(line+'\n')
            else:
                temp = []
                for target in ['ref','que']:
                    for key in order:
                        temp.append(str(master_dict[group][uce][target][key]))
                line=linestart+"\t".join(temp)
                outfh.write(line+'\n')
    
    outfh.close()


def write_UCE_coordinates(outfile, master_dict):
    
    order = ['id','start','end','or','seq']
        
    ## Write out coordinates of candidates
    outfh = open(outfile,'w')
    for group in master_dict.keys():
        for uce in master_dict[group]:
            line = "%s\t%s\t" %(uce,group)
            temp = []
            for target in ['ref','que']:
                for key in order:
                    temp.append(str(master_dict[group][uce][target][key]))
            line+="\t".join(temp)
            outfh.write(line+'\n')
    
    outfh.close()

def extract_UCEs(prefix, candidates, target='ref', extension=0, fasta='', minlength=0):
    """
    Write UCE's to fasta file optionally including and extension of x bp up and downstream of the UCE
    """

    import sys
    from Bio import SeqIO

    if extension and not fasta:
        sys.exit("\nYou want to extract UCE's with and extension - need a fasta file for that\n")
    
    if not extension:
        outfh = open(prefix+'.fasta','w')
        for group in candidates.keys():
            for uce in candidates[group]:
                if not len(candidates[group][uce]['ref']['seq']) >= minlength:
                        continue
                header_suffix = '%s|%s|%s' %(uce, candidates[group][uce][target]['id'], candidates[group][uce][target]['start'])
                outfh.write(">%s|%s\n%s\n" %(group, header_suffix, candidates[group][uce][target]['seq']))
    
        outfh.close()
    
    else:
        #identify ref contigs for which I need to extract sequences
        ctgs = {}

        for group in candidates.keys():
            for uce in candidates[group]:
                
                ctg_id = candidates[group][uce][target]['id']
                if not ctg_id in ctgs:
                    ctgs[ctg_id] = []
                ctgs[ctg_id].append(uce+' '+group)
        
        
        outfh = open(prefix+'.fasta','w')
        for seq in SeqIO.parse(open(fasta,'r'),'fasta'):
            if seq.id in ctgs:
        #        print seq.id
                for uce_with_group in ctgs[seq.id]:
                    (uce,group) = uce_with_group.split(" ")
                    if not len(candidates[group][uce]['ref']['seq']) >= minlength:
                        continue
                        
                    startadd = ''
                    endadd = ''
                    
        #            print uce,group
                    start = candidates[group][uce][target]['start']-extension-1
                    end = candidates[group][uce][target]['end']+extension-1
                    if start < 0:
#                        print group,uce,candidates[group][uce][target]
                        startadd = 'N' * (start*-1)
                        start = 0
                
                    if end > len(seq.seq):
    #                    print uce,candidates[group][uce][target]
    #                    print end,len(seq.seq)
                        endadd = 'N' * (end-len(seq.seq))
                    
                    if candidates[group][uce][target]['or'] == '+':
                        uc_seq = str(seq.seq[start:end])
                    else:
                        uc_seq = str(seq.seq.reverse_complement()[start:end])
        #            print "\t%s,%s,%s" %(uce, candidates['merged'][uce]['que']['start'], candidates['merged'][uce]['que']['end'])
                    
                    header_suffix = '%s|%s|%s' %(uce, candidates[group][uce][target]['id'], candidates[group][uce][target]['start'])
                    outfh.write(">%s|%s\n%s\n" %(group, header_suffix, candidates[group][uce][target]['seq']))
                    outfh.write(">%s|%s\n%s%s%s\n" %(group, header_suffix, startadd, uc_seq, endadd))
    
        outfh.close()



parser = argparse.ArgumentParser(description=DESCRIPTION, prog='find_UCEs.py', 
	formatter_class=argparse.RawDescriptionHelpFormatter, 
	epilog='''   examples: 
	# Minimum command
	./find_UCEs.py --mummer_raw mummer.out --ref_fasta ref.fa --que_fasta que.fa 
	
	# apply some reasonable filters
	./find_UCEs.py --mummer_raw mummer.out --ref_fasta ref.fa --que_fasta que.fa --min_length 100 --max_merge_distance 40 --unique_distance 1000 --dist_from_end 10 --max_merge_distance_diff 5 --extend 500 --output_length_dists
	
	''')

parser.add_argument("--mummer_raw", help="Mummer output file.", metavar="<FILE>", action="store")
parser.add_argument("--ref_fasta", help="Reference assembly in fasta format.", metavar="<FILE>", action="store")
parser.add_argument("--que_fasta", help="Query assembly in fasta format.", metavar="<FILE>", action="store")
parser.add_argument("--min_length", help="Minimum length of candidates to be retained. Default=0, i.e. no filter", metavar="<INT>", default=0, type=int, action="store")
parser.add_argument("--max_merge_distance", help="Maximum allowed distance between Mummer matches on reference sequence to be merged into a single candidate. default=0, i.e. don't merge.", metavar="<INT>", default=0, type=int, action="store")
parser.add_argument("--unique_distance", help="Minimum distance between two adjacent Mummer matches on reference sequence to be consider the matches as separate candidates. default=1000", metavar="<INT>", default=1000, type=int, action="store")
parser.add_argument("--dist_from_end", help="Minimum distance of Mummer match from end of Reference contig/scaffold to be considered as candidate. default=10", metavar="<INT>", default=10, type=int, action="store")
parser.add_argument("--max_merge_distance_diff", help="Maximum difference of distances between adjacent matches to be merged on the reference and query. default=0.", metavar="<INT>", default=0, type=int, action="store")
parser.add_argument("--extend", help="Extract candidates with additional up/downstream sequence of this length into separate file. Default=0, i.e. don't extract extended sequences. Default=0", metavar="<INT>", default=0, type=int, action="store")
parser.add_argument("--output_length_dists", help="Output to screen summaries of the candidate length distributions.", action="store_true", default=False)

args = parser.parse_args()

if len(sys.argv) < 2:	#if the script is called without any arguments display the usage
    print
    print "%s\n" %DESCRIPTION
    parser.print_usage()
    print
    sys.exit(1)


if not args.mummer_raw:
	sys.exit("\n\nCan't work without a Mummer output file\n\n")

if not args.ref_fasta:
	sys.exit("\n\nCan't work without a Reference fasta file\n\n")

if not args.que_fasta:
	sys.exit("\n\nCan't work without a Query fasta file\n\n")


### MAIN

mummer_result = args.mummer_raw
#Tab delimited file
#6 columns
#col 1: query contig id
#col 2: REf contig id
#col 3: ref start pos
#col 4: query start pos
#col 5: length of match
#col 6: nucleotide sequence
#ref = subject in mummer, i.e. the fasta that was specified first in the mummer command
#query = query in mummer, i.e. the fasta that was specified second in the mummer command
#example mummer command: mummer -maxmatch -l 120 -F ref.fasta query.fas
#The input file was created from mummer output with a custom script

#Input files - reference and query genome assembly
quefasta = args.que_fasta
reffasta = args.ref_fasta

max_merge_distance = args.max_merge_distance #maximum distance between 2 unique perfect matches from mummer to be merged to a single bait
unique_distance = args.unique_distance #minimum distance between 2 candidates (if distance is shorter only consider the longer perfect match)
dist_from_edge = args.dist_from_end
max_diff = args.max_merge_distance_diff #maximum gap size difference between ref and query in a case of adjacent candidates 

extension=args.extend
minlength=args.min_length

###########################################

mummerlist = reformat_mummer(mummer_result)

mummer_result_prefix = mummer_result.split("/")[-1]
print mummer_result_prefix+'.tsv'
output_mummer_tabular(mummer_result_prefix+'.tsv', mummerlist)


mummer_fh = open(mummer_result_prefix+'.tsv','r')
#mummer_fh = open('../reformating_mummer_results/test.out.tsv','r')

per_ref_ctg = {}
total = 0

candidates = {'singlematch':{}, 'unique':{}, 'merged':{}}
binned = {'position':{}, 'too-close':{}, 'overlap-filter':{}, 'skip': {}, 'irregular-merge':{}, 'merge-dist-diff':{}}


que_prefix = ".".join(quefasta.split("/")[-1].split(".")[:-1])
ref_prefix = ".".join(reffasta.split("/")[-1].split(".")[:-1])
global_prefix = mummer_result_prefix

lengths = get_ctg_lengths(reffasta)

#read in reformatted mummer output
#The loop produces a dictionary with the following structure:
#key = reference contig id
#value = dictionary with key = start position on reference and value is list with [queryid, querystart, matchlength, sequence]

for line in mummer_fh:
    llist = line.strip().split("\t")
    total+=1
#    print llist
    if not llist[1] in per_ref_ctg:
        per_ref_ctg[llist[1]] = {}
    if not int(llist[2]) in per_ref_ctg[llist[1]]:
        per_ref_ctg[llist[1]][int(llist[2])] = []
    per_ref_ctg[llist[1]][int(llist[2])].append([llist[0],llist[3],llist[4],llist[5]])

print "Total number of MUMMER matches: %i" %total
#print per_ref_ctg    
    
#Filter 'multiple query hits with identical start position on reference'
#searches for multiple matches with the same startposition on the same reference contig, i.e. more than one matches with the query start at the 
#same position in the reference

filtercount=0
for ref_ctg_id in sorted(per_ref_ctg):
#    print ref_ctg_id
    for ref_ctg_pos in sorted(per_ref_ctg[ref_ctg_id]):
        replicates = len(per_ref_ctg[ref_ctg_id][ref_ctg_pos])
        if replicates > 1:
#            print "\nReference startpoint '%s' on reference ctg id '%s' is hit by %s queries -> filter" %(ref_ctg_pos,ref_ctg_id,replicates)

            matchid = ref_ctg_id+"|%s" %ref_ctg_pos
            binned['position'][matchid] = []
            for rep in per_ref_ctg[ref_ctg_id][ref_ctg_pos]:

                binned['position'][matchid].append(reformat_match(ref_id=ref_ctg_id, ref_start_pos=ref_ctg_pos, baselist=rep))
                
                filtercount+=1
                
            del per_ref_ctg[ref_ctg_id][ref_ctg_pos]
            if len(per_ref_ctg[ref_ctg_id]) == 0:
                del per_ref_ctg[ref_ctg_id]

                
print "Filter 'multiple query hits with identical start position on reference' removed %i Mummer hits" %filtercount
total-=filtercount
print "Leaves us with %i hits to process" %total


##extract uniques
#print "Number of contigs in dictionary before unique matches were removed: %s" %len(per_ref_ctg)
#identify the matches 
for ref_ctg_id in sorted(per_ref_ctg):
#    print ref_ctg_id,len(per_ref_ctg[ref_ctg_id]),sorted(per_ref_ctg[ref_ctg_id])
    if len(per_ref_ctg[ref_ctg_id]) == 1:
        ref_ctg_pos = per_ref_ctg[ref_ctg_id].keys()[0]
        temp = per_ref_ctg[ref_ctg_id][ref_ctg_pos][0]

        
        #check if distance from ctg start is ok
#        if not ref_ctg_pos >= (dist_from_edge): #if match is not far enough from the start of the sequence
        start = evaluate_dist_from_end(dist_from_edge=dist_from_edge, pos=ref_ctg_pos, mode='start')
        end = evaluate_dist_from_end(dist_from_edge=dist_from_edge, pos=ref_ctg_pos, mode='end', total_length=lengths[ref_ctg_id])
        if not start or not end:
#        if not evaluate_dist_from_end(dist_from_edge=dist_from_edge, pos=ref_ctg_pos, mode='start'):
            temp.append('too-close')
#            print "too close to edge - %s" %ref_ctg_id
            del per_ref_ctg[ref_ctg_id]
            continue
            
        else: #move match to candidate status
            
            matchid = ref_ctg_id+"|%s" %ref_ctg_pos
#            print "sending to candidates -> %s" %matchid

            candidates['singlematch'][matchid] = reformat_match(ref_id=ref_ctg_id, ref_start_pos=ref_ctg_pos, baselist=temp)
    
#            print "unique - %s" %ref_ctg_id
                    
#        print "%s\n" %candidates['singlematch'][matchid]

        #remove match from original dictionary
        del per_ref_ctg[ref_ctg_id]

#    print "more than one - %s" %ref_ctg_id


#print "Number of contigs in dictionary after unique matches were removed: %s" %len(per_ref_ctg)
################################################

to_merge = {}


#Identify matches in close proximity
#At this stage the per_ref_ctg dictionary contains exclusively reference contigs that had more than one match with the query, but none 
#that started at the same position (see filter 'multiple query hits with identical start position on reference' above)



for ref_ctg_id in sorted(per_ref_ctg):
    pos_list = sorted(per_ref_ctg[ref_ctg_id])
#    print pos_list
#    print ref_ctg_id,len(per_ref_ctg[ref_ctg_id]),sorted(per_ref_ctg[ref_ctg_id])
#    print "\n### %s" %ref_ctg_id #,len(per_ref_ctg[ref_ctg_id]),sorted(per_ref_ctg[ref_ctg_id])

    for i in range(len(pos_list)-1):
#        print per_ref_ctg[ref_ctg_id][pos_list[i]][0]
        if len(per_ref_ctg[ref_ctg_id][pos_list[i]][0]) == 4:
            if pos_list[i] >= (dist_from_edge):
                per_ref_ctg[ref_ctg_id][pos_list[i]][0].append('unique') #only if it's the first match on a contig ][pos_list[i]][0][4]
            else:
                per_ref_ctg[ref_ctg_id][pos_list[i]][0].append('too-close')
                
#        print "\n[%i]: compare '%i' vs. %i" %(i, pos_list[i], pos_list[i+1])

        #calculate the distance between the current match [i] and the next match [i+1] on the reference contig = distance
        current_length = int(per_ref_ctg[ref_ctg_id][pos_list[i]][0][2])
        distance = pos_list[i+1] - (pos_list[i]+current_length)
#        print "%s: %s %i %s %s" %(ref_ctg_id,pos_list[i],current_length, pos_list[i+1], distance)

        #if the distance between two adjacent matches is larger than a defined value (unique_distance) the match will be considered as unique
        if distance >= unique_distance:
#            print "Unique distance"
            per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('unique')
#            if per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] == 'unique':
#                print pos_list[i],per_ref_ctg[ref_ctg_id][pos_list[i]][0]
#                continue
        
        elif distance < 0:
#            print "overlapping match -> exclude"
            per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] = 'overlap-filter'
            per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('overlap-filter')
        
        elif distance > 0:
            if distance <= max_merge_distance:
#            print "Matches within merge distance"
                if per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] == 'unique':
                    per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] = 'merge-'+ref_ctg_id+'|'+str(pos_list[i])
                    if not ref_ctg_id+'|'+str(pos_list[i]) in to_merge:
                        to_merge[ref_ctg_id+'|'+str(pos_list[i])] = []
                    to_merge[ref_ctg_id+'|'+str(pos_list[i])].append(pos_list[i])
                    
                    per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('merge-'+ref_ctg_id+'|'+str(pos_list[i]))
                    to_merge[ref_ctg_id+'|'+str(pos_list[i])].append(pos_list[i+1])
                elif per_ref_ctg[ref_ctg_id][pos_list[i]][0][4].startswith('merge'): #merge label found -> append same merge label to next 
                    per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append(per_ref_ctg[ref_ctg_id][pos_list[i]][0][4])
                    to_merge[per_ref_ctg[ref_ctg_id][pos_list[i]][0][4].replace('merge-','')].append(pos_list[i+1])
                elif per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] == 'overlap-filter':
                    per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('overlap-filter')
                elif per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] == 'too-close':
                    per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('too-close')
                elif per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] == 'skip':
#                    print "### SKIP CASE"
#                    print pos_list[i+1],per_ref_ctg[ref_ctg_id][pos_list[i+1]],pos_list[i-1],per_ref_ctg[ref_ctg_id][pos_list[i-1]]
#                    print "distance: %s vs. %s" %(pos_list[i+1] - (pos_list[i-1]+int(per_ref_ctg[ref_ctg_id][pos_list[i-1]][0][2])), unique_distance)
                    if (pos_list[i+1] - (pos_list[i-1]+int(per_ref_ctg[ref_ctg_id][pos_list[i-1]][0][2]))) >= unique_distance:
                        per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('unique')
#                        print "i+1 will be unique"
                    else:
                        per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('skip')
#                        print "i+1 will be skip"
                        
            else:
                if per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] == 'unique' or per_ref_ctg[ref_ctg_id][pos_list[i]][0][4].startswith('merge'):
                    per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('skip')
                elif per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] == 'skip':
#                    print "### SKIP CASE 2"
#                    print pos_list[i+1],per_ref_ctg[ref_ctg_id][pos_list[i+1]],pos_list[i-1],per_ref_ctg[ref_ctg_id][pos_list[i-1]]
#                    print "distance: %s vs. %s" %(pos_list[i+1] - (pos_list[i-1]+int(per_ref_ctg[ref_ctg_id][pos_list[i-1]][0][2])), unique_distance)
                    if (pos_list[i+1] - (pos_list[i-1]+int(per_ref_ctg[ref_ctg_id][pos_list[i-1]][0][2]))) >= unique_distance:
                        per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('unique')
#                        print "i+1 will be unique 2"
                    else:
                        per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('skip')
#                        print "i+1 will be skip 2"
                else:
                    per_ref_ctg[ref_ctg_id][pos_list[i+1]][0].append('unique')
                                    
                
#        print pos_list[i],per_ref_ctg[ref_ctg_id][pos_list[i]][0]

        if per_ref_ctg[ref_ctg_id][pos_list[i]][0][4] == 'unique':
            matchid = ref_ctg_id+"|%s" %pos_list[i]
            temp = per_ref_ctg[ref_ctg_id][pos_list[i]][0]
#            print "sending to candidates -> %s" %matchid
        
            candidates['unique'][matchid] = reformat_match(ref_id=ref_ctg_id, ref_start_pos=pos_list[i], baselist=temp)
            
#            print "%s\n" %candidates['unique'][matchid]
    
    #output the decision for the last match on the contig
#    print "[%i]: %i - %s" %(len(pos_list), pos_list[-1], per_ref_ctg[ref_ctg_id][pos_list[-1]][0]) 
    if per_ref_ctg[ref_ctg_id][pos_list[-1]][0][4] == 'unique':
        matchid = ref_ctg_id+"|%s" %pos_list[-1]
        temp = per_ref_ctg[ref_ctg_id][pos_list[-1]][0]
#        print "sending to candidates -> %s" %matchid

        candidates['unique'][matchid] = reformat_match(ref_id=ref_ctg_id, ref_start_pos=pos_list[-1], baselist=temp)            
#        print "%s\n" %candidates['unique'][matchid]


print "\nMERGING\n"

#print to_merge

for mergers in sorted(to_merge):
    ok = True
#    print "\n"+mergers,to_merge[mergers]
    ctg_id = "|".join(mergers.split("|")[:-1])
#    print ctg_id
    temp = []
    query_ids = []

    for start in to_merge[mergers]:
#        print start,per_ref_ctg[ctg_id][start][0]
        temp.append(per_ref_ctg[ctg_id][start][0])
        query_ids.append(per_ref_ctg[ctg_id][start][0][0])

    if len(list(set(query_ids))) == 1: #this is true if all query matches are on the same contig
#        print "single queryID"
        pass
    else:
#        print "non unique query ctg"
#        print "\n"+mergers,to_merge[mergers]
        binned['irregular-merge'][mergers] = []
        for start2 in to_merge[mergers]:
#            print ctg_id,start2,per_ref_ctg[ctg_id][start2][0]
            binned['irregular-merge'][mergers].append(reformat_match(ref_id=ctg_id, ref_start_pos=start2, baselist=per_ref_ctg[ctg_id][start2][0]))

        continue
    
    for i in range(len(temp)-1):
#        print temp[i]
        distance_ref = to_merge[mergers][i+1] - (to_merge[mergers][i] + int(temp[i][2]))
#        print "REFERENCE DISTANCE: %i" %distance_ref
        distance_que = int(temp[i+1][1]) - (int(temp[i][1]) + int(temp[i][2]))
#        print "QUERY DISTANCE: %i" %distance_que

#        diff = distance_que-distance_ref
	#calcuate absolute difference between distance of two adjacent matches on the reference and the query
	diff = abs(distance_que-distance_ref)
#        print "DIFF: %i\n" %diff
        
#        if diff >= min_diff and diff <= max_diff:
	if diff <= max_diff:
            pass
        else:
#            print "Diff conflict - diff = %i" %diff
	    binned['merge-dist-diff'][mergers] = []
            for start2 in to_merge[mergers]:
            	binned['merge-dist-diff'][mergers].append(reformat_match(ref_id=ctg_id, ref_start_pos=start2, baselist=per_ref_ctg[ctg_id][start2][0]))
            ok = False
            break
    
    if ok:
#        print "GO MERGE"
        pass
    else:
        continue
    
    ##Merge
    candidates['merged'][mergers] = {'ref':{}, 'que':{}}
#    print candidates['merged'][mergers]
    
    candidates['merged'][mergers]['ref']['id'] = ctg_id
    candidates['merged'][mergers]['ref']['start'] = to_merge[mergers][0]
    candidates['merged'][mergers]['ref']['end'] = to_merge[mergers][-1] + int(temp[-1][2])
    candidates['merged'][mergers]['ref']['seq'] = ''
    candidates['merged'][mergers]['ref']['or'] = '+'
    candidates['merged'][mergers]['que']['id'] = query_ids[0]
    candidates['merged'][mergers]['que']['start'] = int(temp[0][1])
    candidates['merged'][mergers]['que']['end'] = int(temp[-1][1]) + int(temp[-1][2])
    candidates['merged'][mergers]['que']['seq'] = ''
    candidates['merged'][mergers]['que']['or'] = '+'
    
    if query_ids[0].endswith('Reverse'):
        candidates['merged'][mergers]['que']['or'] = '-'
        candidates['merged'][mergers]['que']['id'] = query_ids[0].replace('_Reverse','')

#    print candidates['merged'][mergers]
        

    
## Extract merged sequences from fasta


#identify ref contigs for which I need to extract sequences
merge_ref_ctgs = {}
for merged in candidates['merged']:
    ref_ctg = "|".join(merged.split("|")[:-1])
    if not ref_ctg in merge_ref_ctgs:
        merge_ref_ctgs[ref_ctg] = []
    merge_ref_ctgs[ref_ctg].append(merged)
    
#for ctg in sorted(merge_ref_ctgs):
#    print ctg,merge_ref_ctgs[ctg]
    
from Bio import SeqIO

for seq in SeqIO.parse(open(reffasta,'r'),'fasta'):
    if seq.id in merge_ref_ctgs:
#        print seq.id
        for uce in merge_ref_ctgs[seq.id]:
            if candidates['merged'][uce]['ref']['or'] == '+':
                uc_seq = str(seq.seq[candidates['merged'][uce]['ref']['start']-1:candidates['merged'][uce]['ref']['end']-1])
            else:
                uc_seq = str(seq.seq.reverse_complement()[candidates['merged'][uce]['ref']['start']-1:candidates['merged'][uce]['ref']['end']-1])
#            print "\t%s,%s,%s" %(uce, candidates['merged'][uce]['ref']['refstart'], candidates['merged'][uce]['ref']['refend'])
            candidates['merged'][uce]['ref']['seq'] = uc_seq
#            print candidates['merged'][uce]['ref']

            
#Identify query contigs for which I need to extract sequences
merge_que_ctgs = {}
for merged in candidates['merged']:
    ctg = candidates['merged'][merged]['que']['id']
    
    if not ctg in merge_que_ctgs:
        merge_que_ctgs[ctg] = []
    merge_que_ctgs[ctg].append(merged)

    
#for ctg in sorted(merge_que_ctgs):
#    print ctg,merge_que_ctgs[ctg]


for seq in SeqIO.parse(open(quefasta,'r'),'fasta'):
    if seq.id in merge_que_ctgs:
#        print seq.id
        for uce in merge_que_ctgs[seq.id]:
            if candidates['merged'][uce]['que']['or'] == '+':
                uc_seq = str(seq.seq[candidates['merged'][uce]['que']['start']-1:candidates['merged'][uce]['que']['end']-1])
            else:
                uc_seq = str(seq.seq.reverse_complement()[candidates['merged'][uce]['que']['start']-1:candidates['merged'][uce]['que']['end']-1])
#            print "\t%s,%s,%s" %(uce, candidates['merged'][uce]['que']['start'], candidates['merged'][uce]['que']['end'])
            candidates['merged'][uce]['que']['seq'] = uc_seq

#            print "%s\n%s" %(candidates['merged'][uce]['ref']['seq'],candidates['merged'][uce]['que']['seq'])
            
            #test
#            ctg_id = "|".join(uce.split("|")[:-1])
#            for start in to_merge[uce]:
#                print per_ref_ctg[ctg_id][start][0][3]
            

#move invalid matches to bin
for bingroup in ['too-close','overlap-filter','skip']:
#    print bingroup
    for ref_ctg_id in sorted(per_ref_ctg):
    
        for start in per_ref_ctg[ref_ctg_id]:
#            print bingroup,ref_ctg_id,start
#            print per_ref_ctg[ref_ctg_id][start][0]
            temp = per_ref_ctg[ref_ctg_id][start][0]
            
            if bingroup in temp:
                matchid = ref_ctg_id+"|%s" %start
                
                binned[bingroup][matchid] = reformat_match(ref_id=ref_ctg_id, ref_start_pos=start, baselist=temp)


#output length distribution per candidate group
if args.output_length_dists:
    print "\nLength distributions for candidate groups:"
    for group in candidates.keys():
        print "\n%s\t%s" %(group,len(candidates[group]))
        length_dist = {}
        for uce in candidates[group]:
            length = candidates[group][uce]['ref']['end'] - candidates[group][uce]['ref']['start']
            if not length in length_dist:
                length_dist[length] = 0
            length_dist[length]+=1
        for l in sorted(length_dist):
            print l,length_dist[l]

#Summarize results
print "\n### SUMMARY ###"
               
#output summary for binned matches
print "\nCandidates dropped:"
for bingroup in binned:
    print "%s\t%s" %(bingroup,len(binned[bingroup]))
    
##summarize candidates
total = 0
print "\nCandidates retained:"
for group in candidates.keys():
    print "%s\t%s" %(group,len(candidates[group]))
    total += len(candidates[group])
print "_________________\nGrand TOTAL: %s\n=================" %total
        

print "\n#### PRODUCING OUTPUT FILES ####"
print "\nWRITING GLOBAL OUTPUTS with prefix: '%s'" %global_prefix

write_binned_match_details(global_prefix+'.bin.reason.tsv', master_dict=binned)

write_UCE_coordinates(global_prefix+'.candidate_coordinates.tsv', master_dict=candidates)

name_body = '.candidate_UCEs'
if minlength:
    name_body += '-minlength'+str(minlength)
        
#Extract only UCE's from reference with minimum length of 100
print "WRITING: %s%s.fasta" %(ref_prefix,name_body)
extract_UCEs(prefix = ref_prefix+name_body, candidates=candidates, target='ref', minlength=minlength)
#Extract only UCE's from query with minimum length of 100
print "WRITING: %s%s.fasta" %(que_prefix,name_body)
extract_UCEs(prefix = que_prefix+name_body, candidates=candidates, target='que', minlength=minlength)

if extension:
    name_body += '-extended'+str(extension)
    #Extract UCE's from reference with minimum length of 100 with 500 bp extension
    print "WRITING: %s%s.fasta" %(ref_prefix,name_body)
    extract_UCEs(prefix = ref_prefix+name_body, candidates=candidates, target='ref', minlength=minlength, extension=extension, fasta=reffasta)
    #Extract UCE's from query with minimum length of 100 with 500 bp extension
    print "WRITING: %s%s.fasta" %(que_prefix,name_body)
    extract_UCEs(prefix = que_prefix+name_body, candidates=candidates, target='que', minlength=minlength, extension=extension, fasta=quefasta)
    


