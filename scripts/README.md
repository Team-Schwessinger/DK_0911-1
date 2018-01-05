# genome v01 = original
# genome v03 = manual reassignment of pwoh as pwh
• While writing 'DK_0911_v04_filter_transposable_elements', an incomplete reassignment error (pwoh to htgs) made earlier was discovered.
• DK_0911_v03 was found to have an error regarding reassignment in the gff3 file. While the seqid was changed (pcontig_xxx -> hcontig_xxx_xxx), the ID tag was not (pcontig_xxx.x instead of hcontig_xxx_xxx.x). Furthermore, some fasta files were generated that had a locus tag instead of a header. These were then mapped from locus tag to ID using the incompletely reassigned gff3 file.
• This error was manually fixed on 04/1/18, and the fasta-generating notebook ('DK_0911_v03_generate_fa_files_from_gff3') was also changed so that it would use ID instead of locus tags as soon as they are made (rather than changed down-stream).
• The original 'DK_0911_gff3_file_reassignment' was updated so that proper reassignment occurred in the gff3 files.
• The 'new' (bug-fixed) genome generated with the updated gff3 file is genome_v031.
• Since some lengthy transposable-element-related dataframes were already generated with this incorrect labelling, we 'hack-fixed' by filtering based on the incorrect labels in 'DK_0911_v04_filter_transposable_elements' (same as what was inputted into the transposable-elements comparison programs), and then fix the labels after.
# genome v03.1 = manual reassignment of pwoh as pwh - bugfix
# genome v04 = manual reassignment of pwoh as pwh + filtered out transposable elements



N.B.
• v03 and v031 .fasta files are identical except for headers (v03: locus_tag/ID mix, v031: ID).
• v03 and v031 .gff3 files differ in their attributes column. When using v03 gff3 files to rename v03 headers, incorrect analysis may ensue.

Script running order:
    1) DK_0911_contig_analysis
        -Inputs: v01 (original) .fasta files
        -Programs: NUCmer & MUMmer (alignment of nucleotide sequences)
        -Purpose: determine which pwohs should be reassigned as htgs by seeing how well they belong to other primary contigs.
    2) DK_0911_contig_lengths
        -Inputs: v031 (post-reassignment) .fasta files
        -Programs: N/A
        -Purpose: calculate genomic information pertaining to number and size of contigs. PWH_SIZE is required later for DK_0911_assemblytics_analysis
    3) DK_0911_pwoh_reassignment
        -Inputs: v01 (pre-reassignment) .fasta files, reassignment pairs
        -Programs: N/A
        -Purpose: reassign pwoh as haplotigs in fasta files
    4) DK_0911_gff3_file_reassignment.ipynb
        -Inputs: v01 (original) .gff3 files, reassignment pairs
        -Programs: N/A
        -Purpose: reassign pwoh as haplotigs in .gff3 files
    5) DK_0911_nucmer_and_assemblytics_mapping
        -Inputs: v03 (post-reassignment) .fasta files
        -Programs: NUCmer (alignment) & assemblytics (detection of variants)
        -Purpose: extract meaningful data pertaining to variants in the genome to be further analysed in DK_0911_assemblytics_analysis
        -Notes: Added bug-catching function that checks if mapping folders contain all the analytical files that they should, and with this added file_counting.ipynb that can be imported by other notebooks.
    6) DK_0911_assemblytics_analysis
        -Inputs: Outputs from DK_0911_nucmer_and_assemblytics_mapping
        -Programs: N/A
        -Purpose: Graphical representation of genomic variants.
        -Notes: Added function that automatically labels bar charts at the correct position above the bars (rather than hand-picking values by trial-and-error), included x-axis scientific notation flexibility, relocated titles, and streamlined font-sizing options.
    7) DK_0911_defining_alleles
        -Inputs: assymbletics outputs, v031 gff3 files (post-reassignment & fixed) & protein and gene .fasta files
        -Programs: BLAST
        -Purpose: Identifies alleles and classifies them into different types based on BLAST results, on both the protein and gene level.
    8) DK_0911_
    
    
    
    
    
    
    
    
    