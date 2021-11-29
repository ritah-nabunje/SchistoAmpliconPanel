import argparse
import compilation

# --------------------------------------------Arguments-------------------------------------------------------------
# Initiate the parser
parser = argparse.ArgumentParser(description='This scripts designs panels for SNP-based Targeted Amplicon Sequencing'
                                             ' of Schistosoma mansoni Populations')
# Add input arguments
parser.add_argument("--vcf_file", "-vcf", type=str, help='file in vcf format with the variants/snps to be targeted.',
                    required=False)
parser.add_argument("--exclude_regions", "-exclude", nargs='?', type=str, help="file in BED format with positions of "
                                                                               "the low complexity regions of the "
                                                                               "genome.",
                    required=False)
parser.add_argument("--ref_genome", "-ref", type=str, help="genome's reference sequence in FASTA format.",
                    required=False)
parser.add_argument("--primer_region", "-plr", type=int, default=60,  help="primer landing region.",
                    required=False)
parser.add_argument("--output_prefix", "-out", type=str, help='prefix for all output files', required=True)
parser.add_argument("--amplicon_length", "-alen", type=int, help='chosen length of amplicon in base pairs.',
                    required=True)
parser.add_argument("--total_amplicons", "-amps", type=int, help="number of genome-wide amplicons.",
                    required=False)
parser.add_argument("--dg_threshold", "-dg", type=int, default=-6000, help="maximum accepted delta-g value.",
                    required=False)
parser.add_argument('--left_adapter', '-ad1', type=str, help='adapter for left primers. To be used for scoring of '
                                                             'secondary structures formed between the left primers '
                                                             'and adapter sequences',
                    required=False)
parser.add_argument('--right_adapter', '-ad2', type=str, help='adapter for right primers. To be used for scoring of '
                                                              'secondary structures formed between the right primers '
                                                              'and adapter sequences', required=False)
parser.add_argument('--left_tail', '-tail1', type=str,
                    help='left TrueSeq tail sequence. To be added at the 5- end of the left primer.',
                    required=False)
parser.add_argument('--right_tail', '-tail2', type=str, nargs='?',
                    help='right TrueSeq tail sequence. To be added at the 3- end of the left primer.',
                    required=False)
parser.add_argument('--mode', '-m', type=str, choices=('genome_wide', 'locus_specific'),
                    help='design a genome-wide panel or design a panel for a given locus',
                    required=True)
parser.add_argument('--flag', '-f', type=str, choices=("run_all", 'design_amplicons', "design_primers",
                                                                  "space_amplicons", "score_dimers", "adapter_comp",
                                                                  "best_scoring", "add_tails"),
                    help='run all the steps or run a specific step at a time',
                    required=True)
parser.add_argument('--input', '-in', type=str, help='input file for specific steps', required=False)

# Read arguments from the command line
args = parser.parse_args()
# ----------------------------------------------------------------------------------------------


def main():
    # Check arguments
    vcf = args.vcf_file
    exclude = args.exclude_regions
    ref = args.ref_genome
    plr = args.primer_region
    out = args.output_prefix
    alen = args.amplicon_length
    amps = args.total_amplicons
    dg = args.dg_threshold
    ad1 = args.left_adapter
    ad2 = args.right_adapter
    tail1 = args.left_tail
    tail2 = args.right_tail
    file = args.input
    if args.mode == "genome_wide":
        if args.flag == "design_amplicons":
            return compilation.design_genome_wide_amplicons(snp_vcf=vcf, a_size=alen, exclude_bed=exclude,
                                                            out_file=out, pls=plr)
        elif args.flag == "design_primers":
            return compilation.design_primers(amplicons_file=file, a_size=alen, ref_seq=ref, out_file=out)
        elif args.flag == "space_amplicons":
            return compilation.space_amplicons(primers_file=file, ref_seq=ref, total_amplicons=amps, out_file=out)
        elif args.flag == "score_primer_compatibility":
            return compilation.score_primer_compatibility(primers_file=file, dg_threshold=dg, out_file=out)
        elif args.flag == "score_adapter_compatibility":
            return compilation.score_adapter_compatibility(primers_file=file, forward_ad=ad1,
                                                           reverse_ad=ad2, dg_threshold=dg, out_file=out)
        elif args.flag == "best_scoring_pairs":
            return compilation.best_scoring_pairs(scored_primers_file=file, out_file=out)
        elif args.flag == "add_tails":
            return compilation.add_tails(best_primers_file=file, leftprimer_tail=tail1, rightprimer_tail=tail2,
                                         out_file=out)
        else:
            amplicons = compilation.design_genome_wide_amplicons(snp_vcf=vcf, a_size=alen, exclude_bed=exclude,
                                                                 out_file=out, pls=plr)
            primers = compilation.design_primers(amplicons_file=amplicons, a_size=alen, ref_seq=ref, out_file=out)
            spaced_primers = compilation.space_amplicons(primers_file=primers, ref_seq=ref, total_amplicons=amps,
                                                         out_file=out)
            dimer_checked = compilation.score_primer_compatibility(primers_file=spaced_primers, dg_threshold=dg,
                                                                   out_file=out)
            adapter_compatible = compilation.score_adapter_compatibility(primers_file=dimer_checked, forward_ad=ad1,
                                                                         reverse_ad=ad2, dg_threshold=dg, out_file=out)
            best_primers = compilation.best_scoring_pairs(scored_primers_file=adapter_compatible, out_file=out)
            tailed_primers = compilation.add_tails(best_primers_file=best_primers, leftprimer_tail=tail1,
                                                   rightprimer_tail=tail2, out_file=out)
            return tailed_primers
    else:
        print("locus_specific mode not included yet")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print(' Oops! Panel design process Interrupted.')
