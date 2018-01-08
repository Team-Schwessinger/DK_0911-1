# genome v01 = original
# genome v03 = manual reassignment of pwoh as pwh
• While writing 'DK_0911_v04_filter_transposable_elements', an incomplete reassignment error (pwoh to htgs) made earlier was discovered.
• DK_0911_v03 was found to have an error regarding reassignment in the gff3 file. While the seqid was changed (pcontig_xxx -> hcontig_xxx_xxx), the ID tag was not (pcontig_xxx.x instead of hcontig_xxx_xxx.x). Furthermore, some fasta files were generated that had a locus tag instead of a header. These were then mapped from locus tag to ID using the incompletely reassigned gff3 file.
• This error was manually fixed on 04/1/18, and the fasta-generating notebook ('DK_0911_v03_generate_fasta_files_from_gff3') was also changed so that it would use ID instead of locus tags as soon as they are made (rather than changed down-stream).
• The original 'DK_0911_gff3_file_reassignment' was updated so that proper reassignment occurred in the gff3 files.
• The 'new' (bug-fixed) genome generated with the updated gff3 file is genome_v031.
• Since some lengthy transposable-element-related dataframes were already generated with this incorrect labelling, we 'hack-fixed' by filtering based on the incorrect labels in 'DK_0911_v04_filter_transposable_elements' (same as what was inputted into the transposable-elements comparison programs), and then fix the labels after.
# genome v03.1 = manual reassignment of pwoh as pwh - bugfix
• Bug-fixed gff3 files
# genome v03.2 = manual reassignment of pwoh as pwh - bugfix
• Discovered bug in fasta files, where identifier for reassigned contigs was 'hcontig_033_100 pcontig_100' or 'hcontig_074_103 reverse complement'. This was changed.
• Furthermore, one contig was manually reassigned as its reverse complement, but this manual reverse complementing would not make the contig compatible with the gff3 file. Thus, this was changed back.
• Bug-fixed fasta files
• from DK_0911_generate_fasta_files_from_gff3, fasta files for both proteins and genes had the same identifier (evm.model.xxx) which is incompatible with DK_0911_defining_alleles. Changed DK_0911_dictionaries to name genes evm.TU.xxx and protein evm.model.xxx.
# genome v04 = manual reassignment of pwoh as pwh + filtered out transposable elements

Script running order:
    1) DK_0911_contig_analysis
        -Inputs: v01 (original) .fasta files
        -Programs: NUCmer & MUMmer (alignment of nucleotide sequences)
        -Purpose: determine which pwohs should be reassigned as htgs by seeing how well they belong to other primary contigs.
    2) DK_0911_fasta_file_reassignment
        -Inputs: original .fasta files, reassignment pairs
        -Programs: N/A
        -Purpose: reassign pwoh as haplotigs in fasta files
    3) DK_0911_gff3_file_reassignment
        -Inputs: original .gff3 files, reassignment pairs
        -Programs: N/A
        -Purpose: reassign pwoh as haplotigs in .gff3 files
    4) DK_0911_nucmer_assemblytics_mapping
        -Inputs: .fasta files
        -Programs: NUCmer (alignment) & assemblytics (detection of variants)
        -Purpose: extract meaningful data pertaining to variants in the genome to be further analysed in DK_0911_assemblytics_analysis
        -Notes: Added bug-catching function that checks if mapping folders contain all the analytical files that they should, and with this added file_counting.ipynb that can be imported by other notebooks.
    5) DK_0911_assemblytics_analysis
        -Inputs: Outputs from DK_0911_nucmer_and_assemblytics_mapping
        -Programs: N/A
        -Purpose: Graphical representation of genomic variants.
        -Notes: Added function that automatically labels bar charts at the correct position above the bars (rather than hand-picking values by trial-and-error), included x-axis scientific notation flexibility, relocated titles, and streamlined font-sizing options.
    5)i) DK_0911_contig_lengths
        -Inputs: .fasta files
        -Programs: N/A
        -Purpose: calculate genomic information pertaining to number and size of contigs. PWH_SIZE is required for DK_0911_assemblytics_analysis
    6) DK_0911_generate_fasta_files_from_gff3
        -Inputs: .fasta files, gff3
        -Programs: N/A
        -Purpose: Creates .gene, .cds, .protein fasta files to be used in DK_0911_defining_alleles and DK_0911_post_allele_analysis
    6)i) DK_0911_dictionaries
        -Inputs: gff3 file
        -Programs: N/A
        -Purpose: Creates a dictionary that that maps 'locus_tag' to 'id' from the 'attributes' column of gff3 files. 
    7) DK_0911_defining_alleles
        -Inputs: assemblytics outputs, .gff3 files (post-reassignment & fixed) & protein and gene .fasta files
        -Programs: BLAST & proteinortho
        -Purpose: Identifies alleles and classifies them into different types based on BLAST results, on both the protein and gene level.
    7)i) DK_0911_proteinortho
        -Inputs: 
    7)ii) DK_0911_defining_alleles_no_proteinortho
        -Inputs: assemblytics outputs, .gff3 files (post-reassignment & fixed) & protein and gene .fasta files
        -Programs: BLAST
        -Purpose: filtering alleles on both BLAST & proteinortho can be too aggressive; this notebook only filters using BLAST results
    8) DK_0911_post_allele_analysis
        -Inputs: output from DK_0911_defining_alleles & ph-protein/gene/cds files from DK_0911_generate_fasta_files_from_gff3.
        -Programs: N/A
        -Purpose: generate and save a DataFrame containing dN/dS information (number of nonsynonymous substitutions per non-synonymous site to the number of synonymous substitutions per synonymous site), as well as Hamming & Levenshtein distances (measures of % identity).
    9) DK_0911_filter_transposable_elements
        -Inputs: BLAST & transposonPSI DataFrames, .gff3 files.










