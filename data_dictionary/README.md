# Data Dictionary for the Childhood Cancer and Structural Birth Defects U24 project
This data dictionary contains descriptions about data sets that have been added to the Data Distillery graph (DDKG), which can be found [here](https://github.com/nih-cfde/data-distillery/blob/main/DataDistillery29August2025/DD_29August2025_data_dictionary.md).




## Germline Datasets

### Description
Germline variants from six WES/WGS cohorts have been ingested in the graph. 
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



## Somatic Datasets

### Description
Somatic variants from three cohorts have been ingested into the graph, including KF-NBL, KF-TALL and CBTN.

### Schema 
`...`
### Cypher Query
```cypher
match (hgvsg:Code {SAB:'HGVSG'})-[:HAS_CODE]-(n:Concept)-[r {SAB:"KFSOMNBL"}]-(m:Concept)-[:HAS_CODE]-(hsclo:Code {SAB:'HSCLO'})  
match (n)-[r2:belongs_to_cohort]-(o:Concept)-[:HAS_CODE]-(cohort:Code {SAB:"KFCOHORT"})
match (o)-[r5]-(s:Concept)-[:HAS_CODE]-(kfstudy:Code {SAB:"KFSTUDY"})
match (n)-[r3:related_to_gene]-(p:Concept)-[:HAS_CODE]-(gene:Code {SAB:"ENSEMBL"})
match (n)-[r4:has_protein]-(q:Concept)-[:HAS_CODE]-(protein:Code {SAB:"ENSEMBL"})
return * LIMIT 1
```

#### KF-NBL (Neuroblastoma) Node Counts
| SAB |  Count | 
| :------- | :------: | 
| Cell 1A  | Cell 2A  | 
| Cell 1B  | Cell 2B  | 

#### KF-NBL (Neuroblastoma) Node Counts

| Subject SAB |Predicate |Object SAB | Count |
|----------|----------|----------|----------|
| Row 1 Col 1 | Row 1 Col 2 | Row 1 Col 3 | Row 1 Col 4 |
| Row 2 Col 1 | Row 2 Col 2 | Row 2 Col 3 | Row 2 Col 4 |
| Row 2 Col 1 | Row 2 Col 2 | Row 2 Col 3 | Row 2 Col 4 |
| Row 2 Col 1 | Row 2 Col 2 | Row 2 Col 3 | Row 2 Col 4 |

### KF-TALL (T-cell Acute Lymphoblastic Leukemia)

### CBTN ()



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
### Data Preprocessing scripts
[MTPEXP preprocessing script](https://github.com/U24-CC-BD-Taylor-2024/Childhood_Cancer_Birth_Defects_Taylor_U24/blob/main/preprocessing/MTP_expression/gene_counts_rsem_expected_count_collapsed_deseq.ipynb)


## Data from large scale Kids Firsts cohorts (processed through AutoGVP)

Congenital Heart Defects; KF-CHD (697 probands)  
Neuroblastoma; KF-NBL (460 probands)  
Recessive Structural Birth Defects; KF-SBD  (193 probands)  
Genetics at the Intersection of Childhood Cancer and Birth Defects; KF-GNINT (1,279 probands)  
MMC – (38 probands)  
TALL – (1,310 probands)  





