# DK_0911 - Genome Analysis

These scripts were written during the 2017-2018 Summer break while on a Summer Research Scholarship at the Australian National University.

Supervisor: **Dr Benjamin Schwessinger**

Lab Group: Rathjen Lab, Research School of Biology

## Scripts - in semi running order
1. `DK_0911_contig_analysis`
   - Inputs: v01 (original) .*fasta* files
   - Programs: **NUCmer** & **MUMmer** (alignment of nucleotide sequences)
   - Purpose: Determine which pwohs should be reassigned as htgs by:
     1. aligning them individually against all primary contigs and seeing how well they align with each primary contigs, and
     2. for each pwoh that aligns well with primary contigs, check whether the alignment occurs over an area already occupied by another haplotig.
     - If (A) the pwoh aligns well with another primary contig and (B) the alignment area is not already occupied by another haplotig, the pwohis very likely a haplotig.
2. `DK_0911_fasta_file_reassignment`
   - Inputs: original .*fasta* files, reassignment pairs
   - Programs: N/A
   - Purpose: reassign pwoh as haplotigs in *fasta* files
3. `DK_0911_gff3_file_reassignment`
   - Inputs: original .*gff3* files, reassignment pairs
   - Programs: N/A
   - Purpose: reassign pwoh as haplotigs in .*gff3* files
4. `DK_0911_nucmer_assemblytics_mapping`
   - Inputs: .*fasta* files
   - Programs: **NUCmer** (alignment) & **Assemblytics** (detection of variants)
   - Purpose: extract meaningful data pertaining to variants in the genome to be further analysed in `DK_0911_assemblytics_analysis`
   - Notes: added bug-catching function that checks if mapping folders contain all the analytical files that they should, using a notebook `file_counting.ipynb` that can also be imported by other notebooks.
 `DK_0911_Pst104E_assemblytics_plot`
   - This notebook generates the combined assemblytics analysis as shown in Figure 1. 
5. `DK_0911_assemblytics_analysis`
   - Inputs: outputs from `DK_0911_nucmer_and_assemblytics_mapping`
   - Programs: N/A
   - Purpose: graphical representation of genomic variants.
   - Notes: added function that automatically labels bar charts at the correct position above the bars (rather than hand-picking values by trial-and-error), included x-axis scientific notation flexibility, relocated titles, and streamlined font-sizing options.
   - 5(i) `DK_0911_contig_lengths`
     - Inputs: .*fasta* files
     - Programs: N/A
     - Purpose: calculate genomic information pertaining to number and size of contigs. PWH_SIZE is required for `DK_0911_assemblytics_analysis`.
6. `DK_0911_generate_fasta_files_from_gff3`
   - Inputs: .*fasta* files, *gff3*
   - Programs: **bedtools**, **biocode** (*write_fasta_from_gff.py*)
   - Purpose: creates *gene*, *cds*, *protein* *.fasta* files (for primary contigs and haplotigs, independently and combined) to be used in `DK_0911_defining_alleles_v02` and `DK_0911_post_allele_analysis`.
   - 6(i) `DK_0911_dictionaries`
     - Inputs: *gff3* file
     - Programs: N/A
     - Purpose: creates a dictionary that that maps **locus_tag** to **id** from the **attributes** column of *gff3* files.
7. `DK_0911_defining_alleles_v02`
   - Inputs: **Assemblytics** outputs, .*gff3* files (post-reassignment & fixed) & protein and gene .*fasta* files
   - Programs: **BLAST** & **proteinortho**
   - Purpose: identifies alleles and classifies them into different types based on **BLAST** results, on the protein level. Writes out detailed DataFrames of: (1) all alleles identified by **proteinortho**, (2) best **BLAST** (p on h) alleles out of those not identified by **proteinortho**, and (3) best reciprocal **BLAST** (h on p) alleles based on haplotig proteins that were not classified as alleles in either steps (1) or (2).
   - 7(i) `DK_0911_proteinortho`
     - Inputs: .*gff* files, .*faa* files
     - Programs: **proteinortho**
     - Purpose: **proteinortho** implements a blast-based approach to determine sets of (co-)orthologous proteins or nucleic acid sequences that generalises the reciprocal best alignment heuristic
     - Notes: uses `file_counting.ipynb` to determine whether or not **proteinortho** has been previously run.
   - 7(ii) `DK_0911_defining_alleles_no_proteinortho`
     - OBSOLETE: DK_0911_defining_alleles now uses **proteinortho** by default and any filtering can be done in DK_0911_post_allele_analysis
     - Inputs: **Assemblytics** outputs, .*gff3* files (post-reassignment & fixed) & protein and gene .*fasta* files
     - Programs: **BLAST**
     - Purpose: filtering alleles on both **BLAST** & **proteinortho** can be too aggressive; this notebook only filters using **BLAST** results
8. `DK_0911_post_allele_analysis_v02`
   - Inputs: output from `DK_0911_defining_alleles_v02` & ph-protein/gene/cds files from `DK_0911_generate_fasta_files_from_gff3`.
   - Programs: **MUSCLE**, **PAML**
   - Purpose: generate and save a DataFrame containing dN/dS information (number of nonsynonymous substitutions per non-synonymous site to the number of synonymous substitutions per synonymous site), as well as Hamming & Levenshtein distances (measures of % identity). Also provides visualisations of some of this data.
   - Notes: uses `file_counting.ipynb` to determine whether or not **PAML** has been previously run.
9. `DK_0911_filter_transposable_elements`
   - Inputs: **BLAST** & **transposonPSI** DataFrames, .*gff3* files.
   - Programs: N/A
   - Purpose: generate new *gff3* files without transposable elements. After this, re-run code from step 6. `DK_0911_generate_fasta_files_from_gff3` onwards.
10. `DK0911_TE_filtering_and_summary_p_contigs_fix.ipynb, DK0911_TE_filtering_and_summary_h_contigs_fix.ipynb, Pst_104E_v14_TE_filtering_and_summary_p_contigs_fix.ipynb, Pst_104E_v13_TE_filtering_and_summary_h_contigs_fix.ipynb`
   - TO-DO description
11. `DK0911_TE_variation_analysis.ipynb`
   - TO-DO description
   
12. `Scripts in protein_annotation`
   - TO-DO
   
13. `DK0911_vs_Pst104E_gene_pair_analysis`
   - This notebooks describes the comparative gene pair analysis. The inputs are the outputs of `DK_0911_defining_alleles_v02`, `DK_0911_post_allele_analysis_v02` `Pst_104E_defining_alleles_v02_RT`, and `Pst_104E_defining_alleles_v02_RT`.
   These previous notebooks describe the assignment of alleles and the caculation of different pairwise measures including Levensthein distances and dN/dS ratios.

14. `DK0911_Pst104E_comperative_coverage_analysis`
   - This notebook plots Figure 2 representing the relative fully phased, homozygous collapsed, and hemizygous regions. It takes the input of `DK0911_SRM_cov_DK0911_on_DK0911` and `Pst_104E_SRM_cov_Pst_104E_on_Pst_104E`
   
15. `DK0911v04_Pst104E_presence_absence`
   - This notebook explores the presence absence polymorphisms between teh two genomes based on reciprocal short read mapping and whole genome alignments. The variable regions are then explored to see the genes contained within these low coverage or no coverage regions. These are the regions unique to one genome vs. the others. It looks at genes and at TEs and preforms permutation tests for both. This takes the input of the following other notebooks `DK0911_SRM_cov_DK0911_on_DK0911`, `Pst_104E_SRM_cov_Pst_104E_on_Pst_104E`, `DK0911_SRM_cov_Pst_104E_on_DK0911`, `Pst_104E_SRM_cov_Pst_104E_on_Pst_104E`, and `Mummer_DK0911_Pst_104E`.
   
16. `DK0911_vs_Pst104E_reciprocal_gene_pair_analysis`
   - This notebook explores the reciprocal gene pair analysis between the two isolates and caclualtes a set of different statistics for each gene pair including Levensthein distances and dN/dS ratios. It generates Figure 4. Takes input from `Pst104E_vs_DK0911_synteny` and `DK0911_vs_Pst104E_synteny`.

17. `DK0911v04_Pst104E_presence_absence_orthology` 
   - This notebook explores the presence absence polymorphisms between the two isolates based on Orthofinder define gene families (Orthogroups). It does several statistical tests for enrichment of certain gene groups in the presence absence gene. It also pulls in the synteny analysis by MCScanX and identifies genes within synteny blocks at different gap sizes and the maximal synteny block for genes. This notebooks generates Figure 5.Takes input from `Pst104E_vs_DK0911_synteny` and `DK0911_vs_Pst104E_synteny`.
   
18. `DK0911_v04_effector_analysis`
   - This notebook generates Supplemental Figure 6 with identical analysis published in mBio https://mBio.asm.org/content/9/1/e02275-17.
   
19. `DK0911_tables_and_subfigures`
   - This notebook is a general utility notebook that helps generating some figures.

## Genome Versions
### genome_v01
* Original genome.

### genome_v03
* Manual reassignment of pwoh (primray contigs without haplotigs) as haplotigs.

### genome_v03.1
* While writing `DK_0911_v04_filter_transposable_elements`, an incomplete reassignment error (pwoh to htgs) made earlier was discovered.
* DK_0911_v03 was found to have an error regarding reassignment in the *gff3* file. While the **seqid** was changed (*pcontig_xxx* -> *hcontig_xxx_xxx*), the **ID** was not (*pcontig_xxx.x* instead of *hcontig_xxx_xxx.x*). This caused downstream errors: some *fasta* files were generated that had a **locus_tag** instead of a header, and were then mapped from **locus_tag** to **ID** using the incompletely reassigned *gff3* file.
* This error was manually fixed on 04/1/18, and the *fasta*-generating notebook (`DK_0911_v03_generate_fasta_files_from_gff3`) was also changed so that it would use **ID** instead of **locus_tag** as soon as they are made (rather than changed further downstream).
* The original `DK_0911_gff3_file_reassignment` was updated so that proper reassignment occurred in the *gff3* files.
* The 'new' (bug-fixed) genome generated with the updated *gff3* files is genome DK_0911_v031.
* Since some lengthy transposable-element-related dataframes were already generated with this incorrect labelling, we 'hack-fixed' by filtering based on the incorrect labels in `DK_0911_v04_filter_transposable_elements` (same as what was inputted into the transposable-elements comparison programs), and then fix the labels after. This is hack-fix is commented in detail in `DK_0911_v04_filter_transposable_elements`.

### genome_v03.2
* Discovered bug in *fasta* files, where identifier for manually reassigned contigs was *hcontig_033_100 pcontig_100* or *hcontig_074_103 reverse complement*. `DK_0911_fasta_file_reassignment.ipynb` was changed to generate *fasta* files with the appropriate identifier (e.g. simply *hcontig_033_100* or *hcontig_074_103*, without the superfluous descriptor).
* Furthermore, one contig was manually reassigned as its reverse complement, but this manual reverse complementing would not make the contig compatible with the *gff3* file. Thus, feature was removed in `DK_0911_fasta_file_reassignment.ipynb`.
* From `DK_0911_generate_fasta_files_from_gff3`, *fasta* files for both proteins and genes had the same identifier (*evm.model.xxx...*) which is incompatible with `DK_0911_defining_alleles`. Changed `DK_0911_dictionaries` to map **locus_tag** to **ID** in the form of *evm.TU.xxx...* (genes) and *evm.model.xxx...* (proteins).

### genome_v04
* Genome without transposable elements (*gff3* files have their transposable elements removed in `DK_0911_filter_transposable_elements`).
* Manual reassignment of pwoh (primary contigs without haplotigs) as haplotigs from previous genome versions (> v01) remain.






