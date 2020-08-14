PMID: 31746739

Species: Drosophila

Tissue: Brain

Paper title: Single cell transcriptome atlas of the Drosophila larval brain

Paper link: https://pubmed.ncbi.nlm.nih.gov/31746739/

GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134722

Note:

	Supplementary providing two conditions in this paper, see raw data folder.

		StarvationCondition_finalaggr: 17493 genes, 4645 cells 

		NormalCondition_finalaggr: 17493 genes, 4708 cells

	Sudhir contacted author getting two metadata folders, also contain normalized matrix, only using metadata for processing, see raw data folder.

		Condition based data processing

			Using NormalMetadata for NormalCondition_finalaggr, getting 4407 cells (the number is consistent with paper) 

			StarvationMetadata only contains information about merged condition, so using StarvationMetadata match StarvationCondition_finalaggr, getting 4347 cells

		Merged condition

			Using StarvationMetadata to match StarvationCondition_finalaggr and NormalCondition_finalaggr, getting 4347 and 4349 cells, 8696 in total (which is consistent with paper)


Result location: /n/groups/flyrnai/Yue/scseq_data_formatting/processing_data/31746739_Drosophila_Brain/results/


Note: reprocess the results by using the same metadata file (17 clusters)


