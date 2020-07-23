PMID: 32162708

Species: Drosophila

Tissue: Blood

Paper title: Temporal specificity and heterogeneity of Drosophila immune cells

Paper link: https://www-embopress-org.ezp-prod1.hul.harvard.edu/doi/full/10.15252/embj.2020104486

ArrayExpress link: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8698/

This paper contains two samples: Wild type, Wasp infestation.

| Source Name | Characteristics[stimulus]                 | Derived Array Data File          | Comment [Derived ArrayExpress FTP file]                      | Derived Array Data  File | Comment [Derived ArrayExpress FTP file]                      |
| ----------- | ----------------------------------------- | -------------------------------- | ------------------------------------------------------------ | ------------------------ | ------------------------------------------------------------ |
| NI dataset  | none                                      | raw_feature_bc_matrix_RRCZ22.tgz | ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8698/E-MTAB-8698.processed.3.zip | NI_cell_cluster_ID.txt   | ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8698/E-MTAB-8698.processed.5.zip |
| WI dataset  | infested by the wasp Leptopilina boulardi | raw_feature_bc_matrix_RRCZ23.tgz | ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8698/E-MTAB-8698.processed.3.zip | WI_cell_cluster_ID.txt   | ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8698/E-MTAB-8698.processed.5.zip |



Process them seperately. Add a combine version.

The difference between `E-MTAB-8698_Wildtype.Rmd` and `E-MTAB-8698_Waspinfestation.Rmd` is the input sample name, other parts are identical.

