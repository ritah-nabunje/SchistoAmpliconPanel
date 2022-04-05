=================== SchistoAmpliPan ===================

SchistoAmpliPan is an automated amplicon panel designing module that optimally selects evenly spaced potential amplicons in highly polymorphic regions of the genome. It has been tested and validation will be done using Schistosoma mansoni.

## Content of the repository
SchistoAmpliPan - Folder containing the python scripts of the module

## Required third party software
In order for the module to work correctly, it requires that Primer3, and NCBI-BLAST softwares are installed and can be accessed by the main directory. 

## Usage

The program can be invoked with the -h option to get help and more information.
```
./run_ampliPan.py -h
```

```
usage: run_ampliPan.py [-h] [--vcf_file VCF_FILE] [--exclude_regions [EXCLUDE_REGIONS]] [--ref_genome REF_GENOME]
                       [--primer_region PRIMER_REGION] --output_prefix OUTPUT_PREFIX --amplicon_length AMPLICON_LENGTH
                       [--total_amplicons TOTAL_AMPLICONS] [--dg_threshold DG_THRESHOLD] [--left_adapter LEFT_ADAPTER]
                       [--right_adapter RIGHT_ADAPTER] [--left_tail LEFT_TAIL] [--right_tail [RIGHT_TAIL]] --mode
                       {genome_wide,locus_specific} --flag
                       {run_all,design_amplicons,design_primers,space_amplicons,score_dimers,adapter_comp,best_scoring,add_tails}
                       [--input INPUT]

This scripts designs panels for SNP-based Targeted Amplicon Sequencing of Schistosoma mansoni Populations

optional arguments:
  -h, --help            show this help message and exit
  --vcf_file VCF_FILE, -vcf VCF_FILE
                        file in vcf format with the variants/snps to be targeted.
  --exclude_regions [EXCLUDE_REGIONS], -exclude [EXCLUDE_REGIONS]
                        file in BED format with positions of the low complexity regions of the genome.
  --ref_genome REF_GENOME, -ref REF_GENOME
                        genome's reference sequence in FASTA format.
  --primer_region PRIMER_REGION, -plr PRIMER_REGION
                        primer landing region.
  --output_prefix OUTPUT_PREFIX, -out OUTPUT_PREFIX
                        prefix for all output files
  --amplicon_length AMPLICON_LENGTH, -alen AMPLICON_LENGTH
                        chosen length of amplicon in base pairs.
  --total_amplicons TOTAL_AMPLICONS, -amps TOTAL_AMPLICONS
                        number of genome-wide amplicons.
  --dg_threshold DG_THRESHOLD, -dg DG_THRESHOLD
                        maximum accepted delta-g value.
  --left_adapter LEFT_ADAPTER, -ad1 LEFT_ADAPTER
                        adapter for left primers. To be used for scoring of secondary structures formed between the left primers
                        and adapter sequences
  --right_adapter RIGHT_ADAPTER, -ad2 RIGHT_ADAPTER
                        adapter for right primers. To be used for scoring of secondary structures formed between the right
                        primers and adapter sequences
  --left_tail LEFT_TAIL, -tail1 LEFT_TAIL
                        left TrueSeq tail sequence. To be added at the 5- end of the left primer.
  --right_tail [RIGHT_TAIL], -tail2 [RIGHT_TAIL]
                        right TrueSeq tail sequence. To be added at the 3- end of the left primer.
  --mode {genome_wide,locus_specific}, -m {genome_wide,locus_specific}
                        design a genome-wide panel or design a panel for a given locus
  --flag {run_all,design_amplicons,design_primers,space_amplicons,score_dimers,adapter_comp,best_scoring,add_tails}, -f {run_all,design_amplicons,design_primers,space_amplicons,score_dimers,adapter_comp,best_scoring,add_tails}
                        run all the steps or run a specific step at a time
  --input INPUT, -in INPUT
                        input file for specific steps
```

#### Example to design potential genome-wide amplicons of 250 nt
```
./run_ampliPan.py -m genome_wide -f design_amplicons -vcf /path/to/vcf -exclude /path/to/low_complexity_regions/bed -out output_prefix -alen 250
```

#### Example to design primers for 250 nt amplicon coordinates 
```
./run_ampliPan.py -m genome_wide -f design_primers -ref /path/to/ref.fa -in amplicon_coordinates -out output_prefix -alen 250
```

#### Example to combine primer files generated
```
./concatenate_primer_files.py -in ch1_primers.tsv ch2_primers.tsv ch3_primers.tsv ch4_primers.tsv ch5_primers.tsv ch6_primers1.tsv -out all_primers
```

#### Example to obtain 160 spaced potential amplicons across the genome
```
./run_ampliPan.py -m genome_wide -f space_amplicons -in primers_file.tsv -amps 160 -ref /path/to/ref.fa -out output_prefix
```


