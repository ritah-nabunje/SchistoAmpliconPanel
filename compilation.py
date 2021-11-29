import vcf
import csv
from intervaltree import Interval, IntervalTree
import pandas as pd
import itertools
from itertools import groupby
import subprocess
from subprocess import PIPE
from Bio import SeqIO
import math
import datetime
import glob, os


# ----------------design possible genome-wide amplicons----------------------
def design_genome_wide_amplicons(snp_vcf, a_size, exclude_bed, out_file, pls=60):
    """Returns a list of possible amplicons coordinates"""
    if a_size > 1000:
        flex = 500  # flexibility of amplicon size e.g +/- flex
    else:
        flex = 50
    # global a_stop, a_start
    # repeat_bed = "/nfs/users/nfs_r/rn8/schisto_amplicon_panel/all_problematic_regions012021.bed"
    # repeat_bed = "/Users/rnabunje/Projects/SchistoAmplicons/output/2kb032021/low_complexity_regions_032021.bed"

    # pls = 100  # primer landing site = length in bases to fit a primer (i.e 60-100 bases before first target SNP)
    vcf_reader = vcf.Reader(open(snp_vcf, 'r'))  # read vcf for the snp info
    snp_list = []  # make list to which snp records will be added
    # save each snp record to the list
    for record in vcf_reader:
        snp_info = (record.CHROM, record.POS)  # extract snp chrom and Position info from vcf
        snp_list.append(snp_info)

    snp_dict = dict()  # create empty snp dictionary
    # convert snp_list to a dict and save it the newly created empty snp_dict
    for chrom, pos in snp_list: snp_dict.setdefault(chrom, []).append(pos)
    # save represented chromosomes in the vcf to a list
    chroms = []
    for chrom in snp_dict.keys(): chroms.append(chrom)

    # make dataframe for the space (in bp) between adjacent snps
    spaced_dfs = []  # list of dataframes, one for each chromosome in the vcf
    # make dataframe for each chromosome using the snp positions in the vcf
    for chrom, pos in snp_dict.items():  # dict keys are the chromosomes
        chrm_df = pd.DataFrame.from_dict((list(filter(None, pos))))  # make a df of one column of snp positions
        chrm_df.columns = ['snp_pos']  # name the column as snp_pos for easy reference
        chrm_df.loc[-1] = [0]  # make 0 the first snp position such that the first snp is captured
        chrm_df.index = chrm_df.index + 1  # indices increase by one having place 0 at the beginning
        chrm_df_sorted = chrm_df.sort_index().reset_index(drop=True)  # reset index such that position 0 is at index 0
        low_idx = 0  # lower index
        high_idx = 1  # high index
        snp1_list = []  # make list for the first snp position
        snp2_list = []  # list for the second snp position
        snp_count = [-999]  # number of snps captured btn two points
        last_diffindex = -1  # the index at which the difference was captured
        snp1 = chrm_df_sorted.loc[low_idx, 'snp_pos']  # first snp at the lower index = 0
        snp2 = chrm_df_sorted.loc[high_idx, 'snp_pos']  # second snp is at index 1
        diff = snp2 - snp1  # make difference to record space btn the two snps
        while high_idx < len(chrm_df_sorted) - 1:  # repeat until the last index
            if diff >= pls:
                if last_diffindex == -1:  # this skips the first iteration
                    snp_count_left = -999
                else:
                    snp_count_left = (low_idx - last_diffindex + 1)  # count number of snps
                    snp_count.append(snp_count_left)  # save count of snps btn the two points
                last_diffindex = high_idx  # keep this so as to update the number of snps covered
                snp1_list.append(snp1)  # make record of the snps if the space btn them is enough to fit a primer
                snp2_list.append(snp2)
            low_idx = low_idx + 1
            high_idx = high_idx + 1
            snp1 = chrm_df_sorted.loc[low_idx, 'snp_pos']
            snp2 = chrm_df_sorted.loc[high_idx, 'snp_pos']
            diff = snp2 - snp1
        spaced_df = pd.DataFrame({'SNP1': snp1_list, 'SNP2': snp2_list, 'snp_count': snp_count})
        spaced_dfs.append(spaced_df)  # save the df to the list of spaced_df's

    amplicondfs = []  # list to store df's of amplicons
    # use order of chromosomes to specify which intervals to use to make an interval tree for repeat regions
    for i, spaced_df in enumerate(spaced_dfs):
        chrom = chroms[i]  # for each separate chrom
        chrm_intervals = []
        with open(exclude_bed) as exfile:
            reader = csv.DictReader(exfile, delimiter='\t', fieldnames=['chr', "begin", "end"])
            for line in reader:
                if line['chr'] == chrom:
                    begin = int(line['begin'])
                    end = int(line['end'])
                    interval = (begin, end)
                    chrm_intervals.append(interval)
            exfile.close()
        # make an interval tree for each chromosome, for which a queries [start_primer1:end_primer2] will be made
        tree = IntervalTree.from_tuples(chrm_intervals)
        index_list = list(spaced_df.index)  # list of indices in a chromosome's spaced_df
        # define an amplicon's details as col names
        amplicondf = pd.DataFrame(columns=['Chr', 'primer1_start', 'primer1_end', 'first_snp',
                                           'primer2_start', 'primer2_end', 'gap', 'num_snps', 'length'])
        # loop through all first primer positions
        for primer1_index in index_list:
            first_snp = spaced_df.loc[primer1_index, 'SNP2']
            # end of possible placements for primer 1
            end_primer1 = first_snp - 1  # just before the first snp
            primer2_index = primer1_index + 1  # possible placement of primer2 is on the next row
            a_snps = 0  # amplicon snps
            # loop through primer 2 positions
            while primer2_index < index_list[-1]:  # until the primer 2 index is equal to the last index
                a_snps = a_snps + spaced_df.loc[primer2_index, 'snp_count']  # get snps covered btn primer 1 and 2
                start_primer2 = (spaced_df.loc[primer2_index, 'SNP1'] + 1)  # possible start of primer2
                if start_primer2 - end_primer1 > a_size - (2 * pls):
                    # if the start of primer 2 region is too far, skip to the end of this loop
                    primer2_index = index_list[-1]
                else:
                    length_of_snp_region = start_primer2 - end_primer1
                    max_start_primer1 = max(end_primer1 - (a_size - length_of_snp_region - pls),
                                            spaced_df.loc[primer1_index, 'SNP1'])
                    max_end_primer2 = min(start_primer2 + (a_size - length_of_snp_region - pls),
                                          spaced_df.loc[primer2_index, "SNP2"])
                    amplicon_limit = round((a_size - length_of_snp_region - (2 * pls)) / 2)
                    possible_start_primer1 = end_primer1 - pls - amplicon_limit
                    possible_end_primer2 = start_primer2 + pls + amplicon_limit
                    start_primer1 = 0
                    end_primer2 = 0
                    # check the validity of the primer regions
                    if max_start_primer1 > possible_start_primer1 and max_end_primer2 < possible_end_primer2:
                        primer2_index = index_list[-1]
                    elif max_start_primer1 < possible_start_primer1 and max_end_primer2 < possible_end_primer2:
                        extra = possible_end_primer2 - max_end_primer2
                        start_primer1 = max((possible_start_primer1 - extra), max_start_primer1)
                        end_primer2 = min(possible_end_primer2, max_end_primer2)


                    elif max_start_primer1 > possible_start_primer1 and max_end_primer2 > possible_end_primer2:
                        extra = max_start_primer1 - possible_start_primer1
                        start_primer1 = max(possible_start_primer1, max_start_primer1)
                        end_primer2 = min((possible_end_primer2 + extra), max_end_primer2)
                        # a_start = int(start_primer1)
                        # a_stop = int(end_primer2)
                        # amplicon = end_primer2 - start_primer1

                    else:
                        start_primer1 = max(possible_start_primer1, max_start_primer1)
                        end_primer2 = min(possible_end_primer2, max_end_primer2)
                        # a_start = int(start_primer1)
                        # a_stop = int(end_primer2)
                        # amplicon = end_primer2 - start_primer1

                    a_start = int(start_primer1)
                    a_stop = int(end_primer2)
                    amplicon = end_primer2 - start_primer1
                    # check tree for overlap with low complexity regions if amplicon is within range
                    if amplicon in range((a_size - flex), (a_size + flex)) and len(sorted(tree[a_start:a_stop])) < 1:
                        # store this as  a possible amplicon
                        gap = length_of_snp_region
                        length = a_stop - a_start
                        amplicon_info = {'Chr': chrom, 'primer1_start': int(start_primer1),
                                         'primer1_end': int(end_primer1), 'first_snp': int(first_snp),
                                         'primer2_start': int(start_primer2),
                                         'primer2_end': int(end_primer2), 'gap': int(gap),
                                         'num_snps': int(a_snps), "length": int(length)}
                        amplicondf = amplicondf.append(amplicon_info, ignore_index=True)

                primer2_index = index_list[-1]
        amplicondf["index"] = amplicondf.index
        amplicondf["Amp_ID"] = amplicondf['Chr'] + '_A' + amplicondf['index'].astype(str)
        amplicondf.pop("index")
        amplicondfs.append(amplicondf)

    # save the amplicon regions to a  file
    amplicons = pd.concat(amplicondfs)
    out_name = out_file + "_" + str(a_size) + "bp_amplicons_" + datetime.datetime.now().strftime(
        '%m_%d_%Y_%H_%M') + ".tsv"
    amplicons.to_csv(out_name, index=False, header=True, sep='\t', chunksize=1000000)
    return out_name


# ----------------designing target specific primers -----------------------------------------------------
def design_primers(amplicons_file, a_size, ref_seq, out_file):
    """Returns a list of specific primer sets for amplicon coordinates"""
    if a_size > 500:
        flex = 500  # flexibility of amplicon
    else:
        flex = 50
    # fastafile = "/lustre/scratch118/infgen/team133/rn8/panel/schistosoma_mansoni.PRJEA36577.WBPS15.genomic.fa"
    # fastafile = "/Users/rnabunje/Projects/SchistoAmplicons/data/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa"
    # primer parameters
    primer_opt = 19
    primer_min = 17
    primer_max = 21
    primer_ns = 1
    product_size = "%i-%i" % ((a_size - flex), a_size)
    min_tm = 59.5
    max_tm = 60.5
    maxdiff_tm = 1
    # blast parameters
    # '-evalue', evalue,
    # evalue = str(30000)
    # fmt = str(6)  # no header
    fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send score bitscore"
    wordsize = str(7)
    db = out_file + "refDB"  # name BLAST database

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
    snps_covered = []
    amp_id = []
    amp_start = []
    amp_end = []

    # format DB for the BLAST runs
    make_db = ['makeblastdb',
               '-in', ref_seq,
               '-dbtype', 'nucl',
               '-out', db,
               '-parse_seqids']
    subprocess.call(make_db)
    print(f'makeblastdb -n {ref_seq} -dbtype nucl -out {db} -parse_seqids')
    # make primers for the amplicons, run the primers against the fastafile, pick only the specific primers
    # parse the list of files provided

    with open(amplicons_file) as ampfile:  # open file with amplicon coordinates
        reader = csv.DictReader(ampfile, delimiter='\t')
        # each Amplicon is a line in the file
        for k, amp in enumerate(reader):
            # parse fastafile for the corresponding chrom seq
            for record in SeqIO.parse(ref_seq, "fasta"):
                if record.id == amp['Chr']:
                    chrm_seq = record.seq
            amplicon_id = amp['Amp_ID']
            start = int(amp['primer1_start']) - 1
            end = int(amp['primer2_end'])
            amplicon_seq = chrm_seq[start:end]
            target = str(int(amp['first_snp']) - start + 1) + ',' + amp['gap']
            snps = amp['num_snps']
            cor_amplicon_start = amp['primer1_start']
            cor_amplicon_end = amp['primer2_end']

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
            primerinfolist = primer3out.decode().strip().split("\n")
            # decode the binary checkoutput format
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
                print(f'query line: blastn -query {query_seq} -word_size {wordsize} -out {blast_output} -outfmt {fmt} -db {db}')
                hits_df = pd.read_csv(blast_output, delimiter='\t')
                hits_df.columns = ["query", "subject", "identity", "alignment_length", "mismatches",
                                   "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue",
                                   "bit_score"]
                # group the primers by their last digit (SM_V7_1_A1_left_0', 'SM_V7_1_A1_right_0',
                # 'SM_V7_1_A1_left_1', 'SM_V7_1_A1_right_1)
                d = {}  # make temporary dict
                for a_primer in amp_primers_ids:
                    d.setdefault(a_primer[-1], []).append(a_primer)
                # make a list of lists each with a specific pair e.g [['SM_V7_1_A1_left_0',
                # 'SM_V7_1_A1_right_0'], ['SM_V7_1_A1_left_1', 'SM_V7_1_A1_right_1']]
                pairs = [d[prm] for prm in sorted(d, reverse=False)]

                # pick details for each primer, retrieve the blast-hits and sort the specific ones
                for i, pair in enumerate(pairs):
                    pair_prdt_key = "PRIMER_PAIR_" + str(i) + "_PRODUCT_SIZE"
                    # left primer details
                    left_primer = pair[0]
                    # left primer hits are those with 100% identity alignment to the reference
                    left_df = hits_df[(hits_df['query'] == left_primer) & (hits_df['identity'] >= 90.0)]
                    left_seq_key = "PRIMER_LEFT_" + str(i) + "_SEQUENCE"
                    left_pos_key = "PRIMER_LEFT_" + str(i)
                    left_tm_key = "PRIMER_LEFT_" + str(i) + "_TM"
                    left_gc_key = "PRIMER_LEFT_" + str(i) + "_GC_PERCENT"
                    left_self_key = "PRIMER_LEFT_" + str(i) + "_SELF_ANY_TH"
                    left_selfend_key = "PRIMER_LEFT_" + str(i) + "_SELF_END_TH"
                    # subset hits df for all hits whose alignment length is at least 90% primer length
                    left_df = left_df[left_df['alignment_length'] >= (len(primeroutdict.get(left_seq_key)) * 0.9)]
                    print(f'left primer hits: {left_df}')
                    # pick those whose alignment goes till the end of the primer sequence
                    # left_df = left_df[left_df['q.end'] >= (len(primeroutdict.get(left_seq_key)) -5)]
                    # print(f'left primer hits: {left_df}')
                    # get those with mismatches less than 4 (assuming they'll hydridise)
                    # left_df = left_df[left_df['mismatches'] < 3]
                    # print(f'left primer hits: {left_df}')
                    # right primer details
                    right_primer = pair[1]
                    right_df = hits_df[(hits_df['query'] == right_primer) & (hits_df['identity'] >= 90.0)]
                    right_seq_key = "PRIMER_RIGHT_" + str(i) + "_SEQUENCE"
                    right_pos_key = "PRIMER_RIGHT_" + str(i)
                    right_tm_key = "PRIMER_RIGHT_" + str(i) + "_TM"
                    right_gc_key = "PRIMER_RIGHT_" + str(i) + "_GC_PERCENT"
                    right_self_key = "PRIMER_RIGHT_" + str(i) + "_SELF_ANY_TH"
                    right_selfend_key = "PRIMER_RIGHT_" + str(i) + "_SELF_END_TH"
                    # subset hits df for all hits whose alignment length is at least 90% primer length
                    right_df = right_df[right_df['alignment_length'] >= (len(primeroutdict.get(right_seq_key)) * 0.9)]
                    print(f'right primer hits: {right_df}')
                    # pick those whose alignment goes till the end of the primer sequence
                    # right_df = right_df[right_df['q.end'] >= (len(primeroutdict.get(right_seq_key)) -5)]
                    # print(f'right primer hits: {right_df}')
                    # get those with mismatches less than 4 (assuming they'll hydridise)
                    # right_df = right_df[right_df['mismatches'] < 3]
                    # print(f'right primer hits: {right_df}')
                    # specific primer = one of the primers have one hit with 90% identity for 90% of its length
                    if 0 < len(left_df) < 2 or 0 < len(right_df) < 2:
                        # get the blast hits query and subject detail
                        left_query = left_df['query'].iloc[0]
                        right_query = right_df['query'].iloc[0]
                        left_subject = left_df['subject'].iloc[0]
                        right_subject = right_df['subject'].iloc[0]
                        # if the subject e.g "SM_V7_1" in query = ampliconid e.g "SM_V7_1_Amp20"
                        if left_subject in left_query and right_subject in right_query:
                            # save the primers and the details to corresponding lists
                            primer_id.extend([left_primer, right_primer])
                            primer_chrm.extend(["_".join(left_primer.split("_")[:3]),
                                                "_".join(right_primer.split("_")[:3])])
                            primer_seq.extend([str(primeroutdict.get(left_seq_key)),
                                               str(primeroutdict.get(right_seq_key))])
                            left_primer_start = (start + int(primeroutdict.get(left_pos_key).split(",")[0]) + 1)
                            right_primer_start = (start + int(primeroutdict.get(right_pos_key).split(",")[0]) + 1)
                            primer_start.extend([left_primer_start, right_primer_start])
                            left_primer_length = primeroutdict.get(left_pos_key).split(",")[1]
                            right_primer_length = primeroutdict.get(right_pos_key).split(",")[1]
                            primer_length.extend([left_primer_length, right_primer_length])
                            primer_end.extend([left_primer_start + int(left_primer_length) - 1,
                                               right_primer_start - int(right_primer_length) + 1])
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
                            snps_covered.extend([snps, snps])
                            amp_id.extend([amplicon_id, amplicon_id])
                            amp_start.extend([cor_amplicon_start, cor_amplicon_start])
                            amp_end.extend([cor_amplicon_end, cor_amplicon_end])

            # remove the intermediate files
            subprocess.call(["rm", amplicon_id + ".for"])
            subprocess.call(["rm", amplicon_id + ".rev"])
            subprocess.call(["rm", amplicon_id + ".boulderio"])
            # subprocess.call(["rm", amplicon_id + "_primers.fa"])
            subprocess.call(["rm", amplicon_id + "_blastout"])
        ampfile.close()
    # remove all the database files created
    for dbfile in glob.glob(f'{db}.*'):
        os.remove(dbfile)
    # save the lists as a dataframe
    specific_df = pd.DataFrame(list(zip(primer_chrm, primer_id, primer_start, primer_end, primer_seq, primer_length,
                                        primer_tm, primer_gc, primer_selfcomp, primer_selfendcomp, pair_prdt,
                                        snps_covered, amp_id, amp_start, amp_end)),
                               columns=['Chr', 'ID', 'Start', 'End', 'Primer', "Length", "TM", "%GC", "SelfComp",
                                        "Self3'Comp", "PrdtSize", "SNPs_covered", "Amplicon_ID", "Amplicon_start",
                                        "Amplicon_end"])
    # write dataframe to a file
    out_name2 = out_file + "_"+str(a_size) + "bp_primers"+ datetime.datetime.now().strftime('%m_%d_%Y_%H_%M') + ".tsv"
    specific_df.to_csv(out_name2, index=False, header=True, sep='\t', chunksize=1000000)
    return out_name2


# ------------------------------spacing genome-wide-primers -----------------------
# for genome-wide amplicons only
def space_amplicons(primers_file, ref_seq, total_amplicons, out_file):
    # genome = "/lustre/scratch118/infgen/team133/rn8/panel/schistosoma_mansoni.PRJEA36577.WBPS15.genomic.fa"
    # genome = "/Users/rnabunje/Projects/SchistoAmplicons/data/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa"
    # genome = "/Users/rnabunje/Projects/SchistoAmplicons/data/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa"
    # get chroms in the ref
    primers = pd.read_csv(primers_file, sep='\t')  # read the primer file
    chroms = list(set(primers['Chr'].tolist()))
    # chrms = ["SM_V7_1", "SM_V7_ZW", "SM_V7_3", "SM_V7_2", "SM_V7_4", "SM_V7_5", "SM_V7_6","SM_V7_7"]
    chrom_lengths = []
    for chrom in chroms:
        chrom_length = [len(record) for record in SeqIO.parse(ref_seq, "fasta") if record.id == chrom]
        chrom_lengths.extend(chrom_length)
    genome_size = sum(chrom_lengths)
    chromdict = {chroms[c]: chrom_lengths[c] for c in range(len(chroms))}  # dictionary "chrm" : length
    pool_chrm_primers_dfs = []  # empty list of dataframes for the pool primers
    for chrom in chromdict.keys():
        spaced_amp_positions = [] # list for the a chrom's spaced amplicon positions
        chrom_primers = primers[primers['Chr'] == chrom]
        chrom_length = chromdict.get(chrom)  # length of chrm
        # number of amplicons from the chrm
        chrom_amplicons = math.floor((total_amplicons * (chrom_length / genome_size)) + 0.5)
        chrom_spaces = chrom_amplicons + 1  # spaces btn adjacent amplicons
        spacing_size = math.floor((chrom_length / chrom_spaces) + 0.5)  # length between points
        chrom_points = []  # where amplicons should be
        for idx in range(0, chrom_amplicons):
            point = spacing_size * (idx + 1)
            chrom_points.append(point)
        # get_midpoint for an amplicon
        chrom_primers["midpoint"] = chrom_primers[['Amplicon_start', 'Amplicon_end']].mean(axis=1)
        # priority goes to 3+ snp amplicons
        best_primers = chrom_primers[chrom_primers["SNPs_covered"] >= 3]
        others = chrom_primers[chrom_primers["SNPs_covered"] < 3]
        amplicon_positions = list(set(best_primers['midpoint'].tolist()))
        other_positions = list(set(others['midpoint'].tolist()))
        for point in chrom_points:
            if len(amplicon_positions) != 0:
                closest_amplicon_position = min(amplicon_positions, key=lambda x: abs(x - point))
                spaced_amp_positions.append(closest_amplicon_position)
                amplicon_positions.remove(closest_amplicon_position)
            else:
                closest_amplicon_position = min(other_positions, key=lambda x: abs(x - point))
                spaced_amp_positions.append(closest_amplicon_position)
                other_positions.remove(closest_amplicon_position)
        # compile the spaced amplicons
        spaced_amp_df = chrom_primers[chrom_primers['midpoint'].isin(spaced_amp_positions)]
        pool_chrm_primers_dfs.append(spaced_amp_df)
    spaced_primers_df = pd.concat(pool_chrm_primers_dfs)
    out_name3 = out_file + "_" + "_spaced_primers.tsv" + datetime.datetime.now().strftime('%m_%d_%Y_%H_%M')
    spaced_primers_df.to_csv(out_name3, index=False, header=True, sep='\t', chunksize=1000000)
    return out_name3


# ---------------------------------------primer compatibility check-----------------------------------------------
# score primers by their compatibility with the others
def score_primer_compatibility(primers_file, out_file, dg_threshold=-6000):
    # dg_threshold in kcal/mol more negative values suggest stronger binding
    primers = pd.read_csv(primers_file, sep='\t')
    primers["score"] = 0  # pre-set the dimer score to 0
    primer_seqs = list(primers["Primer"])  # make a list of primers
    primer_ids = list(primers["ID"])  # list of IDs

    # make combinations for primers in the list without repetition
    for s1, s2 in itertools.combinations_with_replacement(primer_ids, 2):
        # check if the primers are for the same amplicon
        if "_".join(s1.split("_")[:5]) == "_".join(s2.split("_")[:5]):
            dg = -1
        else:
            seq1 = primer_seqs[primer_ids.index(s1)]
            seq2 = primer_seqs[primer_ids.index(s2)]
            dimer_check = subprocess.check_output(["ntthal", "-s1", seq1, "-s2", seq2])
            info = dimer_check.decode().strip().split("\t")
            dg = int(float(info[3].strip().split("=")[1]))
        if dg <= dg_threshold:
            primers.loc[primers['ID'] == s1, 'score'] += 999
            primers.loc[primers['ID'] == s2, 'score'] += 999
    primers.reset_index(inplace=True)
    out_name4 = out_file + "_compatibility_scored_primers_" + datetime.datetime.now().strftime(
        '%m_%d_%Y_%H_%M') + ".tsv"
    primers.to_csv(out_name4, index=False, header=True, sep='\t', chunksize=1000000)
    return out_name4


# ----------------------secondary structure formation between primers and Illumina adapters---------------------------
def score_adapter_compatibility(primers_file, forward_ad, reverse_ad, dg_threshold, out_file):
    # forward_ad = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
    # reverse_ad = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
    primers = pd.read_csv(primers_file, sep='\t')
    primer_seqs = list(primers["Primer"])  # make a list of primers
    primer_ids = list(primers["ID"])  # list of IDs
    for i, id_ in enumerate(primer_ids):
        if "left" in id_:
            adapter = forward_ad
        else:
            adapter = reverse_ad
        primer = primer_seqs[i]
        adapter_check = subprocess.check_output(["ntthal", "-s1", primer, "-s2", adapter])
        info = adapter_check.decode().strip().split("\t")
        dg = int(float(info[3].strip().split("=")[1]))
        if dg <= dg_threshold:
            primers.loc[primers['ID'] == id_, 'score'] += 999
    out_name5 = out_file + "_adapter_scored_primers_" + datetime.datetime.now().strftime('%m_%d_%Y_%H_%M') + ".tsv"
    primers.to_csv(out_name5, index=False, header=True, sep='\t', chunksize=1000000)
    return out_name5


# --------------------pick best scoring primer pairs--------------------
def best_scoring_pairs(scored_primers_file, out_file):
    primers = pd.read_csv(scored_primers_file, sep='\t')
    # add column with total score
    primers['total_score'] = primers['SelfComp'] + primers["Self3'Comp"] + primers['score']
    # make list of final amplicons
    best_amplicon_pairs = []
    pool_amplicons = list(set(primers['Amplicon_ID'].tolist()))
    for amplicon in pool_amplicons:
        amplicon_primersdf = primers[primers['Amplicon_ID'] == amplicon]
        amplicon_primers = [x for x in amplicon_primersdf["ID"].to_list()]
        temp = sorted(amplicon_primers, key=lambda x: x[-1])  # sort primers IDs by their last digit in a temp
        # make list of lists grouped by the last digit
        primer_pairs = [list(elt) for i, elt in groupby(temp, lambda x: x[-1])]
        # make list of sums of scores for each pair of the amplicon
        pair_scores = [sum([(primers.loc[primers['ID'] == pair[0], 'total_score'].iloc[0]),
                            (primers.loc[primers['ID'] == pair[1], 'total_score'].iloc[0])]) for pair in primer_pairs]
        best_pair = primer_pairs[pair_scores.index(min(pair_scores))]  # pair with min scores
        best_amplicon_pairs.extend(best_pair)
    final_primers = primers[primers['ID'].isin(best_amplicon_pairs)]
    final_primers = final_primers.drop(['index', 'score', 'total_score', 'Amplicon_ID', 'midpoint'], axis=1)
    out_name6 = out_file + "_best_scoring_pairs_" + datetime.datetime.now().strftime('%m_%d_%Y_%H_%M') + ".tsv"
    final_primers.to_csv(out_name6, index=False, header=True, sep='\t', chunksize=1000000)
    return out_name6


# ---------------------------------add TrueSeq 5' tail sequence-----------------------
def add_tails(best_primers_file, leftprimer_tail, rightprimer_tail, out_file):
    # leftprimer_tail = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
    # rightprimer_tail = "TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
    primers = pd.read_csv(best_primers_file, sep='\t')
    primers["TrueSeq_tail"] = "None"
    primers.loc[primers['ID'].str.contains("left"), "TrueSeq_tail"] = leftprimer_tail
    primers.loc[primers['ID'].str.contains("right"), "TrueSeq_tail"] = rightprimer_tail
    primers["tailedPrimer"] = primers["TrueSeq_tail"] + primers["Primer"]
    # output
    out_name7 = out_file + "_tailed_primers_" + datetime.datetime.now().strftime('%m_%d_%Y_%H_%M') + ".tsv"
    primers.to_csv(out_name7, index=False, header=True, sep='\t', chunksize=1000000)
