# Comprehensive Discovery and Characterization of Driver Genes and Mutations in Human Cancers

__TCGA Driver and Gene Essentiality Working Group__

__Supplement: Exploration of Clinical Actionability in TCGA mc3 with PHIAL and the TARGET database__

Brendan Reardon, Nathanael Moore, Eliezer Van Allen

VanAllen Lab, Dana-Farber Cancer Institute, Broad Institute of MIT & Harvard

Last updated: 31 July 2017

## Introduction
Precision heuristics for the alteration landscape (PHIAL)<sup>1</sup> was developed as a heuristic clinical interpretation algorithm to identify clinically actionable or biologically relevant alterations in patient tumor sequencing data. This is accomplished by leveraging curated databases relevant to cancer, such as TARGET and COSMIC. A database of tumor alterations relevant for genomics-driven therapy (TARGET) was co-developed with PHIAL to synthesize the current landscape of precision medicine research, through a strategy of manual curation of literature and the use of expert opinions on the therapeutic and prognostic implications of alterations for patients with cancer. Although PHIAL was developed to study each patient individually, we have applied it to the TCGA mc3 cohort in parallel to gain insight into clinically relevant somatic alterations across TCGA mc3 and individual tissue types. 

## Data Collection
Single nucleotide variants (SNVs) and insertions and deletions (InDels) were provided to identify clinically relevant somatic alterations within the TCGA mc3 cohort across tissue types. Copy number variations (CNVs) called from GISTIC<sup>2</sup> were provided to identify clinically relevant somatic copy number alterations. 

## Sample Processing and Methods
[PHIAL v1.2.0](https://github.com/broadinstitute/phial/tree/v1.2.0) (currently private) was run utilizing [TARGET v1.4.2] (https://github.com/vanallenlab/2017-tcga-mc3_phial/blob/master/reference/target/TARGET_db_v1.4.2_05312017.txt) and COSMIC v79<sup>3</sup> on TCGA mc3 SNV and InDels. Code was modified to map the column names of inputs to required fields for PHIAL v1.2.0, which is dependent on Oncotator for annotations and format. 

As a result of input data not being annotated with Oncotator, UniProt region information is not available and thus score adjustments based on variants involving a kinase domain are unavailable for this evaluation. 

We discuss variants that matched gene and feature type to a clinical assertion in TARGET annotated as “PUTATIVELY ACTIONABLE” (PHIAL Score bins: “Potential TARGET Actionability”, “Investigate TARGET Actionability”) and variants that matched only genes to clinical assertions annotated as “BIOLOGICALLY RELEVANT” (PHIAL Score bin: “High Priority”). 

Thresholded copy number variations (CNVs) called from GISTIC were received for TCGA mc3 data (syn5049520). A manual query of TARGET database was performed for copy number alterations since PHIAL is not currently compatible with GISTIC outputs. Copy number calls were subsetted for genes present in TARGET v1.4 and further restricted to consider only amplifications and deletions, as defined by the GISTIC value +2, -2, for the purpose of this analysis.

This analysis will treat SNVs and InDels entirely separately from the copy number alterations due to the number of samples considered for each data type. SNV and InDels previously received filtering after pathology review and removal of hypermutated samples, multiple samples from the same individual, and other sequencing artifacts. 8,577 samples are considered for SNV and InDels and ~10,700 samples are considered for copy number alterations.

## Cohort Description 
TCGA mc3 contains samples derived from 33 tissue types. A total of 9,079 samples produced 1,457,702 SNV and InDel variant calls, and 10,713 samples contained 97,798,661 copy number events. The intersection of samples yields 8775 samples. 

## Cohort Results
A total of 9,039 putatively actionable SNVs and InDels were discovered across TCGA mc3, accounting for 0.62 % percent of all variant calls, as well as 13,953 amplifications and 7,011 deletions. The most commonly mutated genes across samples were PIK3CA (11.5% of samples, 1176 variants), BRAF (6.96%, 677 variants), ARID1A (6.78%, 640 variants), KRAS (6.1%, 628 variants), and IDH1 (5.45%, 496 variants). Genes most subject to copy number alterations were CDKN2A/B (16.68% / 16.10%, 1429 / 1379 samples), MYC (10.21%, 874), PIK3CA (7.61%, 651), CCND1 (7.36%, 629 ), and EGFR (6.28%, 537). 

Taken together, we observe that 55.6% of SNV & InDel samples and 48.1% of CNV samples had at least one putatively actionable variant, with a mean number of 0.94 SNVs/InDels and 1.09 CNVs per sample. If we expand our scope to also include biologically relevant alterations, 76.2% of SNV/InDel samples and 57.8% of CNV samples had at least 1 TARGET related alteration, with a mean number of 1.81 SNVs/InDels and 1.96 CNVs per sample. 
