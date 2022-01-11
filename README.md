# nf for axolotl visium processing with kallisto
The original codes were from Tomas and Ashley from Barbara group ETH and highly appreciate that they were willing to share.


## Documentation of nf visium initially by Tomas.

IMPORTANT:
in the kallisto_pipeline.nf file, there is a command for running spaceranger (you can Ctrl+F "spaceranger" to find it).
you will need to edit the path for the mock fastq reads to include where you're saving those that I'm sending.
For example the current path there is
"/links/groups/treutlein/USERS/tomasgomes/gene_refs/human/refdata-gex-GRCh38-2020-A/mock_fastq",
this will need to be changed to the folder now containing the "mock_fastq" folder I'm sending.

Example of a command to quantify Visium data with the Axolotl transcriptome:

	nextflow run kallisto_pipeline.nf
		--transcriptome /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial.fa
		--transindex AmexT_v47_artificial.kalid
		--t2g /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial_genenames_t2g.txt
		--white /links/groups/treutlein/USERS/tomasgomes/gene_refs/other/visium-v1_whitelist_kallisto.txt --samplename "A1_limb"
		--outdir /links/groups/treutlein/USERS/tomasgomes/data/axolotl/
		--protocol visiumv1
		--reads "/links/groups/treutlein/DATA/sequencing/20200821_P1288_ASHLEY_VISIUM_axolotl_visium_control_11dpa/raw/A1_Animal1_Control/*.fastq.gz"
		--images "V19S23-109"
		--imagear "A1"
		--imagef "/links/groups/treutlein/DATA/sequencing/20200821_P1288_ASHLEY_VISIUM_axolotl_visium_control_11dpa/image/A1_large_image1.jpg"
		--imageal "/links/groups/treutlein/DATA/sequencing/20200821_P1288_ASHLEY_VISIUM_axolotl_visium_control_11dpa/alignment_files/V19S23-109-A1.json"


Command breakdown:

nextflow run kallisto_pipeline.nf

	# this is the fasta transcriptome of axolotl. in this case it also includes some sequences for eGFP etc, although that might not be needed.
	# the fasta headers just have the isoform name (e.g. ">AMEX60DD301000001.2")
	# I originally got this from Sergej, but can also be obtained from the axolotl website
	--transcriptome /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial.fa

	# this is the indexed transcriptome for kallisto. if it doesn't exist, it'll be generated from the fasta file, and saved in the same folder
	# if it does exist, the fasta file still has to be specified, and they have to be in the same folder
	--transindex AmexT_v47_artificial.kalid

	# this is a table with the isoform name in one column and the gene in another column. it's required to summarise the kallisto results by gene
	# can use either the gene name or gene ID, however you prefer to summarise the data
	# example row: AMEX60DD301000001.1	ZNF568
	--t2g /links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/cDNA_transcripts/AmexT_v47_artificial_genenames_t2g.txt

	# this is the whitelist for the visium barcodes (present in spaceranger, also needed for kallisto)
	--white /links/groups/treutlein/USERS/tomasgomes/gene_refs/other/visium-v1_whitelist_kallisto.txt

	# this is the name given to the sample
	--samplename "A1_limb"

	# the directory where to output the results
	--outdir /links/groups/treutlein/USERS/tomasgomes/data/axolotl/

	# the protocol being used (since the pipeline also works for other types of data)
	--protocol visiumv1

	# the path for the reads. this has to be passed like in this example (i.e. all fastq files in the same folder,
	# containing either R1 or R2 in the name, just passed with the "*" to automatically list all of them)
	# example file names: A1_Animal1_Control_S1_L001_R1_001.fastq.gz and A1_Animal1_Control_S1_L001_R2_001.fastq.gz
	--reads "/links/groups/treutlein/DATA/sequencing/20200821_P1288_ASHLEY_VISIUM_axolotl_visium_control_11dpa/raw/A1_Animal1_Control/*.fastq.gz"

	# these next parameters are all for spaceranger
	# Visium slide serial number, for example 'V10J25-015'
	--images "V19S23-109"

	# Visium area identifier, for example 'A1'
	--imagear "A1"

	# path to the H&E brightfield image file
	--imagef "/links/groups/treutlein/DATA/sequencing/20200821_P1288_ASHLEY_VISIUM_axolotl_visium_control_11dpa/image/A1_large_image1.jpg"

	# Alignment file produced by the manual Loupe alignment step
	# this is something that Ashley did (manually aligning the images and the spots in Loupe), if you need help you can ask her
	--imageal "/links/groups/treutlein/DATA/sequencing/20200821_P1288_ASHLEY_VISIUM_axolotl_visium_control_11dpa/alignment_files/V19S23-109-A1.json"
