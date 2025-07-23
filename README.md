# Saturation-Mutagenesis-Pipeline

This workflow is used for analysis of directed evolution of saturation mutagenesis libraries, as seen in “Evolution of a Functionally Intact but Antigenically Distinct DENV Fusion Loop” (eLife, 2023).

Required:
•	Linux terminal
•	Python
•	Perl
•	R

•	Anylyze sequencing libraries using CAMseqv4.pl
  •	Usage: perl CAMseqv4.pl file.fastq flank5’ flank3’
  •	Optional input - R can be added as fourth input on command line to reverse complement sequences before analyzing. 
  •	Note: flanking sequences must be in same direction as sequences as represented in original fastq file. It is recommended to use 5-7 nucleotides on either side of the saturation mutagenesis site.
  •	CAMseqv4.pl forces a specific number of nucleotides between flanking regions. To allow other lengths, use CAMseqv3.pl (Tse et al. PNAS 2017).

•	Prepare output for plotting using merge.py
  •	usage: merge.py selected_aminoacids.txt unselected_aminoacids.txt outfile.txt
•	If interested in calculating amino acid distance from wildtype sequence, use merge_AA.py
  •	Edit file to include wildtype amino acid sequence
  •	usage: merge.py selected_aminoacids.txt unselected_aminoacids.txt outfile.txt

•	Plot data in R using bubbplotplot_enriched.r and piechart.r
  •	Code will need to be modified with input file, and for any changes to plotting options.
•	If interested in visualizing amino acid distance from wildtype sequence, use bubbleplot_enriched_AA.r
  • Code will need to be modified with input file, and for any changes to plotting options. Additionally, fuzzy percentages to call binning of amino acid differences must be adjusted.
