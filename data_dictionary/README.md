# Data Dictionary for the Childhood Cancer and Structural Birth Defects datasets
This data dictionary contains descriptions about data sets that have been added to the Data Distillery graph (DDKG), which can be found [here](https://github.com/nih-cfde/data-distillery/blob/main/DataDistillery29August2025/DD_29August2025_data_dictionary.md). We have 4 types of data that we've added on top of the DDKG: 
- `Germline variants`  
- `Somatic variants`  
- `Molecular Targets Program (MTP) Tumor Expression data`  
- `single-cell data` (coming soon)  




## Germline and Somatic Variants Datasets

### Germline and Somatic Variant datasets summary

| Kids First Name | Kids First Code | dbGAP ID | Domain | Germline Variants | Somatic Variants |
|-----------------|-----------------|----------|--------|------------------|------------------ |
| Comprehensive Genomic Profiling to Improve Prediction of Clinical Outcome for Children with T-cell Acute Lymphoblastic Leukemia | KF-TALL | [phs002276](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002276.v4.p1) | Cancer | 5478 | 13250 |
| Kids First: Genetics at the Intersection of Childhood Cancer and Birth Defects | KF-GNINT | [phs001846](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001846.v1.p1) | Cancer & Birth Defects | 4323 | 0 |
| National Heart, Lung, and Blood Institute (NHLBI) Bench to Bassinet Program: The Gabriella Miller Kids First Pediatric Research Program of the Pediatric Cardiac Genetics Consortium (PCGC) | KF-CHD | [phs001138](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001138.v4.p2) | Birth Defects | 2084 | 0 |
| Discovering the Genetic Basis of Human Neuroblastoma: A Gabriella Miller Kids First Pediatric Research Program (Kids First) Project | KF-NBL | [phs001436](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001436.v1.p1) | Cancer | 4232 | 1591 |
| Kids First: Whole Exome, Genome, and RNA Sequencing in Recessive Structural Brain Defects in Children | KF-RSBD | [phs002590](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002590.v2.p1) | Birth Defects | 827 | 0 |
| Kids First: Germline and Somatic Variants in Myeloid Malignancies in Children | KF-MMC | [phs002187](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002187.v1.p1) | Cancer |137 | 0 |
----------------------------------------------------------

### Germline Variants Summary 

|   |   |
|---|---|
| Dataset SAB(s) | KFGLCHD, KFGLNBL, KFGLTALL, KFGLGNINT, KFGLRSBD, KFGLMMC |
| Description |    Germline variants from 6 WES/WGS cohorts have been ingested into the graph. The cohorts come from the Gabriella Miller Kids First Data Portal and additional details can be found in the table above, along with the total number of germline vs somatic variants that have been ingested. Only variants that were scored as HIGH from Ensembl's Variant Effect Predictor algorthim were included.  |
| Purpose |     Including germline variants in a graph database allows you to represent complex relationships between patients, variants, genes, and biological pathways in a way that supports flexible, multi-hop querying. By integrating variant data from pediatric cancer and congenital defect cohorts, the graph structure enables discovery of shared genetic architecture across disorders. It makes it possible to identify variants that connect to common genes, pathways, or functional annotations that may underlie both cancer susceptibility and developmental abnormalities. The graph format also supports advanced algorithms—such as link prediction, community detection, and GNN-based inference—that can highlight pleiotropic or convergent risk loci. Ultimately, this structure accelerates the identification of germline variants that may be causal or contributory to both childhood cancers and structural birth defects.        |
| Schema Organization | All of the germline datasets follow the exact same simple schema. Variants, with HGVSG identifiers, are connected to Cohort nodes, which define what cohort they belong to, Gene nodes, and Transcript nodes. Variants are also attached to their corresponding Population nodes which identifies which demographic the variant is found in. A schema diagram showing the nodes and edges of the germline schema, along with the Cypher query used to generate the figure can be found below. | 
| Website | https://kidsfirstdrc.org (dbGap study pages can be found in the dbGAP ID column of the Germline and Somatic Variant datasets summary table.) |
---

### Schema 
![](https://github.com/U24-CC-BD-Taylor-2024/Childhood_Cancer_Birth_Defects_Taylor_U24/blob/main/data_dictionary/images/Screenshot%202025-11-02%20at%204.58.37%20PM.png)

### Cypher Query
```cypher
with 'KFGLCHD' as sab
match (var:Code {SAB:'HGVSG'})-[:HAS_CODE]-(var_cui:Concept)-[r:cohort_has_variant {SAB:sab}]-(c2:Concept)-[:HAS_CODE]-(co:Code)
match  (var_cui)-[:gene_has_variant {SAB:sab}]-(c3:Concept)-[:HAS_CODE]-(gene:Code) 
match (var_cui)-[:has_population_frequency {SAB:sab}]-(c4:Concept)-[:HAS_CODE]-(co4:Code)
match (var_cui)-[:transcript_has_variant_impact_vep {SAB:sab}]-(c5:Concept)-[:HAS_CODE]-(co5:Code)
match (var_cui)-[:transcript_has_variant_impact_autogvp {SAB:sab}]-(c6:Concept)-[:HAS_CODE]-(co6:Code)
match (var_cui)-[:transcript_has_variant_impact_polyphen {SAB:sab}]-(c7:Concept)-[:HAS_CODE]-(co7:Code)
match (var_cui)-[:transcript_has_variant_impact_sift {SAB:sab}]-(c8:Concept)-[:HAS_CODE]-(co8:Code)
//match (var_cui)-[s]-(c9:Concept)-[:HAS_CODE]-(co9:Code {SAB: 'HSCLO'})
return * LIMIT 1
```

#### Node Counts
| Entity               | Node SAB   |   KFGLCHD |   KFGLNBL |   KFGLTALL |   KFGLGNINT |   KFGLRSBD |   KFGLMMC |
|:---------------------|:-----------|----------:|----------:|-----------:|------------:|-----------:|----------:|
| Variant              | HGVSG      |      2084 |      4232 |       5478 |        4323 |        827 |       137 |
| Cohort               | COHORT     |         1 |         1 |          1 |           1 |          1 |         1 |
| Gene                 | ENSEMBL    |      1359 |      1900 |       2563 |        2166 |        641 |       123 |
| Transcript           | ENSEMBL    |      1359 |      1899 |       2563 |        2168 |        641 |       123 |
| Population           | POPULATION |         7 |         7 |          7 |           7 |          7 |         7 |
| Chromosomal Location | HSCLO      |      2084 |      4232 |       5478 |        4323 |        827 |       137 |
---------------

#### Edge Counts
| Subject SAB   | Predicate                              | Object SAB   |   KFGLCHD |   KFGLNBL |   KFGLTALL |   KFGLGNINT |   KFGLRSBD |   KFGLMMC |
|:--------------|:---------------------------------------|:-------------|----------:|----------:|-----------:|------------:|-----------:|----------:|
| COHORT        | cohort_has_variant                     | HGVSG        |      2084 |      4232 |       5478 |        4323 |        827 |       137 |
| ENSEMBL       | gene_has_variant                       | HGVSG        |     13813 |     27052 |      35716 |       28906 |       5380 |       921 |
| POPULATION    | has_population_frequency               | HGVSG        |      8638 |      6706 |      18207 |       17689 |       2247 |       539 |
| HGVSG         | transcript_has_variant_impact_autogvp  | ENSEMBL      |      2084 |      4235 |       5478 |        4323 |        827 |       137 |
| HGVSG         | transcript_has_variant_impact_polyphen | ENSEMBL      |       677 |       584 |       1458 |        1439 |        305 |        55 |
| HGVSG         | transcript_has_variant_impact_sift     | ENSEMBL      |       672 |       578 |       1446 |        1430 |        301 |        55 |
| HGVSG         | transcript_has_variant_impact_vep      | ENSEMBL      |      2084 |      4231 |       5478 |        4323 |        827 |       137 |
----------------------------




## Somatic Datasets

|   |   |
|---|---|
| Dataset SAB(s) | KFSOMNBL, KFSOMTALL, KFSOMCBTN |
| Description |    Somatic variants from 3 WES/WGS cohorts have been ingested into the graph. Like the germline cohorts, the somatic variant cohorts come from the Gabriella Miller Kids First Data Portal. Only variants that were scored as HIGH from Ensembl's Variant Effect Predictor algorthim were included in the ingestion.  |
| Purpose |    Including somatic variants in the graph database allows for capturing tumor-specific alterations and connecting them to the germline landscape, patients, genes, pathways, and clinical features. By integrating somatic data from the three cancer cohorts, we can uncover interactions between inherited risk and acquired mutations that may jointly drive pediatric tumor development. The graph structure supports identifying recurrent somatic events and shared mutational signatures across cancers and structural birth defects. It also facilitates multi-hop analyses linking somatic variants to functional pathways, therapeutic targets, or co-occurring germline variants that might reveal combined genetic mechanisms. Overall, incorporating somatic variants strengthens the ability to discover biologically meaningful relationships that differentiate or unify the cancer cohorts.        |
| Schema Organization | Similar to the germline organiza, all of the somatic datasets follow the same schema. Variants, (again with HGVSG identifiers, are connected to Cohort nodes, which define what cohort they belong to, Gene nodes, and Transcript nodes. The somatic variants are also attached to their respective Protein nodes. A schema diagram showing the nodes and edges of the somatic schema, along with the Cypher query used to generate the figure can be found below. | 
| Website | https://kidsfirstdrc.org (dbGap study pages can be found in the dbGAP ID column of the Germline and Somatic Variant datasets summary table.) |
------------------------


### Schema 

![](https://github.com/U24-CC-BD-Taylor-2024/Childhood_Cancer_Birth_Defects_Taylor_U24/blob/main/data_dictionary/images/Screenshot%202025-11-06%20at%201.06.08%20PM.png)

### Cypher Query
```cypher
match (hgvsg:Code {SAB:'HGVSG'})-[:HAS_CODE]-(n:Concept)-[r {SAB:"KFSOMNBL"}]-(m:Concept)-[:HAS_CODE]-(hsclo:Code {SAB:'HSCLO'})  
match (n)-[r2:belongs_to_cohort]-(o:Concept)-[:HAS_CODE]-(cohort:Code {SAB:"KFCOHORT"})
match (o)-[r5]-(s:Concept)-[:HAS_CODE]-(kfstudy:Code {SAB:"KFSTUDY"})
match (n)-[r3:related_to_gene]-(p:Concept)-[:HAS_CODE]-(gene:Code {SAB:"ENSEMBL"})
match (n)-[r4:has_protein]-(q:Concept)-[:HAS_CODE]-(protein:Code {SAB:"ENSEMBL"})
return * LIMIT 1
```

#### Node Counts
| Entity               | Node SAB   |   KFSOMNBL |   KFSOMTALL |   KFSOMCBTN |
|:---------------------|:-----------|-----------:|------------:|------------:|
| Variant              | HGVSG      |       1591 |       13250 |        8655 |
| Gene                 | ENSEMBL    |       1508 |        7300 |        6028 |
| Cohort               | KFCOHORT   |          1 |           1 |           1 |
| Study                | KFSTUDY    |          1 |           1 |           1 |
| Chromosomal Location | HSCLO      |       1591 |       13250 |        8655 |

#### Edge Counts

| Subject SAB   | Predicate         | Object SAB   |   KFSOMNBL |   KFSOMTALL |   KFSOMCBTN |
|:--------------|:------------------|:-------------|-----------:|------------:|------------:|
| HGVSG         | related_to_gene   | ENSEMBL      |       1619 |       13253 |        8657 |
| HGVSG         | belongs_to_cohort | KFCOHORT     |          1 |           1 |           1 |
| HGVSG         | has_location      | HSCLO        |       2103 |       14477 |       10871 |
| KFCOHORT      | study_has_cohort  | KFSTUDY      |          1 |           1 |           1 |


--------------------------------------------------------
## MTP Expression data

### Description

`...`

### Schema
![](https://github.com/U24-CC-BD-Taylor-2024/Childhood_Cancer_Birth_Defects_Taylor_U24/blob/main/preprocessing/MTP_expression/Screenshot%202025-07-14%20at%2011.11.14%20AM.png)
### Cypher Queries
```
# Query to generate schema figure
MATCH (co1:Code {SAB:'MTPEXP'})-[:CODE]-(c1:Concept)-[:gene]-(c2:Concept)-[:CODE]-(co2:Code {SAB:'ENSEMBL'})-[]-(t1:Term) 
WHERE co1.CodeID CONTAINS 'Neuroblastoma'
MATCH (c1)-[:disease]-(c3:Concept)-[:CODE]-(co3:Code {SAB: 'MONDO'})-[:PT]-(t2:Term)
MATCH (c1)-[:tumor_expression]-(c4:Concept)-[:CODE]-(co4:Code {SAB:'EXPBINS'})
RETURN * limit 1
```

#### Node Counts
| Entity          | Node SAB   |   MTPEXP |
|:----------------|:-----------|---------:|
| MTP Expression  | MTPEXP     |  1064061 |
| Gene            | ENSEMBL    |    50635 |
| Disease         | MONDO      |       84 |
| Expression Bins | EXPBINS    |       18 |
-------------------------------------------

#### Edge Counts

| Subject SAB   | Predicate         | Object SAB   |   MTPEXP |
|:--------------|:------------------|:-------------|---------:|
| MTPEXP        | related_to_gene   | ENSEMBL      |  1064061 |
| MTPEXP        | belongs_to_cohort | KFCOHORT     |  3424928 |
| MTPEXP        | has_location      | HSCLO        |   579887 |
-------------
### Data Preprocessing scripts
[MTPEXP preprocessing script](https://github.com/U24-CC-BD-Taylor-2024/Childhood_Cancer_Birth_Defects_Taylor_U24/blob/main/preprocessing/MTP_expression/gene_counts_rsem_expected_count_collapsed_deseq.ipynb)


## Data from large scale Kids Firsts cohorts (processed through AutoGVP)

Congenital Heart Defects; KF-CHD (697 probands)  
Neuroblastoma; KF-NBL (460 probands)  
Recessive Structural Birth Defects; KF-SBD  (193 probands)  
Genetics at the Intersection of Childhood Cancer and Birth Defects; KF-GNINT (1,279 probands)  
MMC – (38 probands)  
TALL – (1,310 probands)  





