# import sys
# import os
import vcf
import pandas as pd
from Bio import SeqIO
import csv
import subprocess
# import time
from intervaltree import Interval, IntervalTree
import argparse

# input files
parser = argparse.ArgumentParser()
parser.add_argument("--vcf_file", "-v", type=str, required=True)
parser.add_argument("--out", "-o", type=str, required=True)
parser.add_argument("--amp_size", "-a", type=int, required=True)
args = parser.parse_args()
out_file = args.out
a_size = args.amp_size


def design_amplicons(snp_vcf):
    """Returns a list of possible amplicons coordinates"""
    repeat_bed = "/path/to/low_comlexity_regions.bed"
    pls = 100  # primer landing site = length in bases to fit a primer (i.e 80-100 bases before first target SNP)
    vcf_reader = vcf.Reader(open(snp_vcf, 'r'))
    snp_list = []
    for record in vcf_reader:
        record = (record.CHROM, record.POS)  # extract snp chrom and Position info from vcf
        snp_list.append(record)  # save snp info to list

    snp_dict = dict()  # create empty snp dictionary

    # convert snp_list a dict and save it to snp_dict
    for chrm, pos in snp_list: snp_dict.setdefault(chrm, []).append(pos)

    # save represented chrms in the vcf to a list
    chrms = []
    for chrm in snp_dict.keys(): chrms.append(chrm)

    # make dataframe for the space (in bp) between adjacent snps
    spaced_dfs = [] # list of dataframes, one for each chrm in the vcf
    # make dataframe for each chrm using the snp positions in the vcf
    for chrm, pos in snp_dict.items():
        chrm_df = pd.DataFrame.from_dict((list(filter(None, pos))))
        chrm_df.columns = ['snp_pos']
        chrm_df.loc[-1] = [0]  # make 0 the first snp position such that the first snp is captured
        chrm_df.index = chrm_df.index + 1
        chrm_df = chrm_df.sort_index().reset_index(drop=True)  # reset index such that 0 is at index = 0
        low_idx = 0 # lower index
        high_idx = 1 # high index
        snp1_list = []
        snp2_list = []
        snp_count = [-999] # number of snps captured
        last_idx = len(chrm_df) - 1 # index of the last row in the df
        last_diffindex = -1
        snp1 = chrm_df.loc[low_idx, 'snp_pos']
        snp2 = chrm_df.loc[high_idx, 'snp_pos']
        diff = snp2 - snp1
        while high_idx < last_idx:
            if diff >= pls:
                if last_diffindex == -1:
                    snp_count_left = -999
                else:
                    snp_count_left = (low_idx - last_diffindex + 1)
                    snp_count.append(snp_count_left)
                last_diffindex = high_idx
                snp1_list.append(snp1)
                snp2_list.append(snp2)
            low_idx = low_idx + 1
            high_idx = high_idx + 1
            snp1 = chrm_df.loc[low_idx, 'snp_pos']
            snp2 = chrm_df.loc[high_idx, 'snp_pos']
            diff = snp2 - snp1
        spaced_df = pd.DataFrame({'SNP1': snp1_list, 'SNP2': snp2_list, 'snp_count': snp_count})
        spaced_dfs.append(spaced_df)

    # use order of chrms to specify which intervals to use in the tree
    # make an interval tree for each chromosome, then query [start_primer1:end_primer2]
    amplicondfs = []
    for i, spaced_df in enumerate(spaced_dfs):
        chrom = chrms[i]
        chrm_intervals = []
        with open(repeat_bed) as rptfile:
            reader = csv.DictReader(rptfile, delimiter='\t')
            for line in reader:
                if line['chrm'] == chrom:
                    begin = int(line['begin'])
                    end = int(line['end'])
                    interval = (begin, end)
                    chrm_intervals.append(interval)
            rptfile.close()
        tree = IntervalTree.from_tuples(chrm_intervals)

        index_list = list(spaced_df.index)
        amplicondf = pd.DataFrame(columns=['primer1_start', 'primer1_end', 'target_snp',
                                           'primer2_start', 'primer2_end', 'gap', 'num_snps'])
        # loop through all first primer positions
        for primer1_index in index_list:
            # end of possible placements for primer 1
            end_primer1 = (spaced_df.loc[primer1_index, 'SNP2'] - 1)
            # possible start of primer1 is
            start_primer1 = (end_primer1 - pls)
            a_start = int(start_primer1)
            primer2_index = primer1_index + 1
            a_snps = 0  # amplicon snps
            # loop through primer 2 positions
            while primer2_index < index_list[-1]:
                a_snps = a_snps + spaced_df.loc[primer2_index, 'snp_count']
                start_primer2 = (spaced_df.loc[primer2_index, 'SNP1'] + 1)
                # end_primer2 = (spaced_df.loc[primer2_index, 'snp2'] - 1)
                end_primer2 = (start_primer1 + a_size)
                if start_primer2 - end_primer1 > a_size - (2 * pls):
                    # if the start of primer 2 region is too far, skip to the end of this loop
                    primer2_index = index_list[-1]
                else:
                    # otherwise, store this as  a possible amplicon
                    a_stop = int(end_primer2)
                    if len(sorted(tree[a_start:a_stop])) < 1:
                        target = spaced_df.loc[primer1_index, 'SNP2']
                        gap = start_primer2 - (end_primer1 + 1)
                        amplicon_info = {'primer1_start': start_primer1, 'primer1_end': end_primer1,
                                         'target_snp': target, 'primer2_start': start_primer2,
                                         'primer2_end': end_primer2, 'gap': gap, 'num_snps': a_snps}
                        amplicondf = amplicondf.append(amplicon_info, ignore_index=True)
                primer2_index = primer2_index + 1
        rptfile.close()
        amplicondfs.append(amplicondf)

    # save each chrm's amplicon regions to a file
    for chrm, amplicondf in zip(chrms, amplicondfs):
        amplicondf.to_csv(chrm+"_"+out_file+"_"+str(a_size)+"_amplicons.tsv",
                          index=False, header=True, sep='\t', chunksize=1000000)
    return chrms


# ----------------designing primers and checking their specificity-----------------------------------------------------
def design_primers(chrms):
    """Returns a list of specific primer sets for amplicon coordinates"""
    fastafile = "/path/to/reference.fa"
    # primer parameters
    primer_opt = 19
    primer_min = 17
    primer_max = 21
    primer_ns = 1
    product_size = "%i-%i" % ((a_size - 50), a_size)
    min_tm = 58
    max_tm = 62
    maxdiff_tm = 5
    # blast parameters
    # '-evalue', evalue,
    # evalue = str(30000)
    fmt = str(6)  # no header
    wordsize = str(7)
    vcf_name = args.vcf_file
    db = vcf_name[:8] + str(a_size) + "genomeDB" # name database with a prefix of the vcf file input

    # output lists
    primer_id = []
    primer_chrm = []
    primer_seq = []
    primer_start = []
    primer_end = []
    primer_length = []
    primer_tm = []
    primer_gc = []
    primer_selfcomp = []
    primer_selfendcomp = []
    pair_prdt = []

    # format DB for the BLAST runs
    make_db = ['makeblastdb',
               '-in', fastafile,
               '-dbtype', 'nucl',
               '-out', db,
               '-parse_seqids']
    subprocess.call(make_db)

    # make primers for the amplicons, run the primers against the fastafile, pick only the specific primers
    for chrm in chrms:
        # parse fastafile for each chrom's seq
        for record in SeqIO.parse(fastafile, "fasta"):
            if record.id == chrm:
                chrm_seq = record.seq
                # print(record.id)
                # open corresponding file with amplicon
                with open(chrm+"_"+out_file+"_"+str(a_size)+"_amplicons.tsv") as ampfile:
                    reader = csv.DictReader(ampfile, delimiter='\t')
                    # each Amplicon is a line in the file
                    for k, amplicon in enumerate(reader):
                        amplicon_id = chrm + '_A' + str(k)
                        start = int(amplicon['primer1_start']) - 1
                        end = int(amplicon['primer2_end'])
                        amplicon_seq = chrm_seq[start:end]
                        target = str(int(amplicon['target_snp']) - start + 1) + ',' + amplicon['gap']
                        # target is relative to the ampliconSeq, includes the first snp and gap
                        # which includes the remaining gaps
                        # define contents of boulderfile as a dict
                        boulder = {
                            "SEQUENCE_ID": amplicon_id,
                            "SEQUENCE_TEMPLATE": amplicon_seq,
                            "SEQUENCE_TARGET": target,
                            "PRIMER_TASK": "generic",
                            "PRIMER_OPT_SIZE": primer_opt,
                            "PRIMER_MIN_SIZE": primer_min,
                            "PRIMER_MAX_SIZE": primer_max,
                            "PRIMER_MAX_NS_ACCEPTED": primer_ns,
                            "PRIMER_PRODUCT_SIZE_RANGE": product_size,
                            "PRIMER_MIN_TM": min_tm,
                            "PRIMER_MAX_TM": max_tm,
                            "PRIMER_PAIR_MAX_DIFF_TM": maxdiff_tm,
                            "P3_FILE_FLAG": 1,
                            "SEQUENCE_INTERNAL_EXCLUDED_REGION": target,
                            "PRIMER_EXPLAIN_FLAG": 1,
                        }
                        # each amplicon has a separate boulder file
                        boulderfile = amplicon_id + ".boulderio"
                        # write the boulder contents to the amplicon's file
                        with open(boulderfile, "w") as boulderf:
                            for param in boulder.keys():
                                # add "=" btn parameter and value, then "newline" after each parameter
                                boulderf.write(str(param + "=" + (str(boulder[param]) + "\n")))
                            boulderf.write("=" + "\n")
                        # Run primer3_core via commandline (outside script)
                        # print("running Primer3")
                        primer3out = subprocess.check_output(["primer3_core", boulderfile])
                        primeroutdict = {}
                        primerinfolist = primer3out.decode().strip().split("\n")  # decode the binary checkoutput format
                        for line in primerinfolist:
                            entry = line.strip().split("=")
                            # entry[0] is the key, entry[1] has the output values
                            primeroutdict[entry[0]] = entry[1]
                        # sometimes zero pairs are found
                        if primeroutdict.get('PRIMER_PAIR_NUM_RETURNED') != '0':
                            # print(ampliconID, "has primers")
                            primer_file = open(amplicon_id + "_primers.fa", "w")  # check file the primers were saved
                            amp_primers_ids = []
                            pairs = int(primeroutdict.get('PRIMER_PAIR_NUM_RETURNED'))
                            for num in range(0, pairs):
                                left_head = amplicon_id + "_left_" + str(num)
                                right_head = amplicon_id + "_right_" + str(num)
                                amp_primers_ids.extend(
                                    (left_head, right_head))  # save only the the headers of the fasta
                                left_seq = primeroutdict.get("PRIMER_LEFT_" + str(num) + "_SEQUENCE")
                                right_seq = primeroutdict.get("PRIMER_RIGHT_" + str(num) + "_SEQUENCE")
                                # write fasta for the primer sequences to a file
                                primer_file.write('>{0}\n{1}\n>{2}\n{3}\n'.format(left_head, left_seq, right_head,
                                                                                  right_seq))
                            primer_file.close()
                            # run blast with query as the amplicon primers that were saved in the .fa
                            query_seq = amplicon_id + "_primers.fa"  # primers make the query sequence
                            blast_output = amplicon_id + "_blastout"  # specifies blast output file
                            query_line = ['blastn',
                                          '-query', query_seq,
                                          '-word_size', wordsize,
                                          '-out', blast_output,  # save BLAST output to a file
                                          '-outfmt', fmt,
                                          '-db', db]
                            subprocess.call(query_line)
                            hits_df = pd.read_csv(blast_output, delimiter='\t')
                            hits_df.columns = ["query", "subject", "identity", "alignment_length", "mismatches",
                                               "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue",
                                               "bit_score"]
                            # group the primers by their last digit (SM_V7_1_A1_left_0', 'SM_V7_1_A1_right_0',
                            # 'SM_V7_1_A1_left_1', 'SM_V7_1_A1_right_1)
                            d = {}  # make temporary dict
                            for primer in amp_primers_ids:
                                d.setdefault(primer[-1], []).append(primer)
                            # make a list of lists each with a specific pair e.g [['SM_V7_1_A1_left_0',
                            # 'SM_V7_1_A1_right_0'], ['SM_V7_1_A1_left_1', 'SM_V7_1_A1_right_1']]
                            primer_pairs = [d[n] for n in sorted(d, reverse=False)]

                            # pick details for each primer, retreive the blast-hits and sort the specific ones
                            for i, pair in enumerate(primer_pairs):
                                pair_prdt_key = "PRIMER_PAIR_" + str(i) + "_PRODUCT_SIZE"
                                # left primer details
                                left_primer = pair[0]
                                left_df = hits_df[(hits_df['query'] == left_primer) & (hits_df["identity"] == 100.0)]
                                left_seq_key = "PRIMER_LEFT_" + str(i) + "_SEQUENCE"
                                left_pos_key = "PRIMER_LEFT_" + str(i)
                                left_tm_key = "PRIMER_LEFT_" + str(i) + "_TM"
                                left_gc_key = "PRIMER_LEFT_" + str(i) + "_GC_PERCENT"
                                left_self_key = "PRIMER_LEFT_" + str(i) + "_SELF_ANY_TH"
                                left_selfend_key = "PRIMER_LEFT_" + str(i) + "_SELF_END_TH"
                                # subset hits df for all hits whose aligment length is close to primer length by <= 1
                                left_df = left_df[left_df['alignment_length'] >=
                                                  (len(primeroutdict.get(left_seq_key)) - 1)]

                                # right primer details
                                right_primer = pair[1]
                                right_df = hits_df[(hits_df['query'] == right_primer) & (hits_df["identity"] == 100.0)]
                                right_seq_key = "PRIMER_RIGHT_" + str(i) + "_SEQUENCE"
                                right_pos_key = "PRIMER_RIGHT_" + str(i)
                                right_tm_key = "PRIMER_RIGHT_" + str(i) + "_TM"
                                right_gc_key = "PRIMER_RIGHT_" + str(i) + "_GC_PERCENT"
                                right_self_key = "PRIMER_RIGHT_" + str(i) + "_SELF_ANY_TH"
                                right_selfend_key = "PRIMER_RIGHT_" + str(i) + "_SELF_END_TH"
                                # subset hits df for all hits whose aligment length is close to primer length by <= 1
                                right_df = right_df[right_df['alignment_length'] >=
                                                    (len(primeroutdict.get(right_seq_key)) - 1)]

                                if len(left_df) == 1 and len(right_df) == 1:
                                    left_query = left_df['query'].iloc[0]
                                    right_query = right_df['query'].iloc[0]
                                    left_subject = left_df['subject'].iloc[0]
                                    right_subject = right_df['subject'].iloc[0]
                                    if left_subject in left_query and right_subject in right_query:
                                        primer_id.extend([left_primer, right_primer])
                                        primer_chrm.extend(["_".join(left_primer.split("_")[:3]),
                                                            "_".join(right_primer.split("_")[:3])])
                                        primer_seq.extend([str(primeroutdict.get(left_seq_key)),
                                                           str(primeroutdict.get(right_seq_key))])
                                        left_primer_start = (start + int(primeroutdict.get(left_pos_key).split(",")[0]))
                                        right_primer_start = (start +
                                                              int(primeroutdict.get(right_pos_key).split(",")[0]))
                                        primer_start.extend([left_primer_start, right_primer_start])
                                        left_primer_length = primeroutdict.get(left_pos_key).split(",")[1]
                                        right_primer_length = primeroutdict.get(right_pos_key).split(",")[1]
                                        primer_length.extend([left_primer_length, right_primer_length])
                                        primer_end.extend([left_primer_start + int(left_primer_length),
                                                           right_primer_start + int(right_primer_length)])
                                        primer_tm.extend([primeroutdict.get(left_tm_key),
                                                          primeroutdict.get(right_tm_key)])
                                        primer_gc.extend([primeroutdict.get(left_gc_key),
                                                          primeroutdict.get(right_gc_key)])
                                        primer_selfcomp.extend([primeroutdict.get(left_self_key),
                                                                primeroutdict.get(right_self_key)])
                                        primer_selfendcomp.extend([primeroutdict.get(left_selfend_key),
                                                                   primeroutdict.get(right_selfend_key)])
                                        pair_prdt.extend([primeroutdict.get(pair_prdt_key),
                                                          primeroutdict.get(pair_prdt_key)])
                        # print(amplicon_id, ": done")
                        subprocess.call(["rm", amplicon_id + ".for"])
                        subprocess.call(["rm", amplicon_id + ".rev"])
                        subprocess.call(["rm", amplicon_id + ".boulderio"])
                        subprocess.call(["rm", amplicon_id + "_primers.fa"])
                        subprocess.call(["rm", amplicon_id + "_blastout"])
                ampfile.close()
    subprocess.call(["rm", db + "*"])
    specific_df = pd.DataFrame(list(zip(primer_chrm, primer_id, primer_start, primer_end, primer_seq, primer_length,
                                        primer_tm, primer_gc, primer_selfcomp, primer_selfendcomp, pair_prdt)),
                               columns=['Chr', 'ID', 'Start', 'End', 'Primer', "Length", "TM", "%GC", "SelfComp",
                                        "Self3'Comp", "PrdtSize"])
    specific_df.to_csv(out_file+"_"+str(a_size)+"_primers.tsv", index=False, header=True, sep='\t', chunksize=1000000)


def main():
    vcf_file = args.vcf_file
    chrm_list = design_amplicons(vcf_file)
    design_primers(chrm_list)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('Interrupted')
