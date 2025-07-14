# Data Dictionary for the Childhood Cancer and Structural Birth Defects U24 project
This data dictionary contains descriptions about data sets that have been added to the Data Distillery graph. Descriptions about data sets in the Data Distillery graph can be found in the Data Distillerys data dictionary


## MTP Expression data

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
[MTPEXP script](preprocessing/MTP_expression/gene_counts_rsem_expected_count_collapsed_deseq.ipynb)
## -------------------------------
## MTP Variant data

## Data from large scale Kids Firsts cohorts (processed through AutoGVP)

Congenital Heart Defects; KF-CHD (697 probands)  
Neuroblastoma; KF-NBL (460 probands)  
Recessive Structural Birth Defects; KF-SBD  (193 probands)  
Genetics at the Intersection of Childhood Cancer and Birth Defects; KF-GNINT (1,279 probands)  
MMC – (38 probands)  
TALL – (1,310 probands)  

## Somatic Datasets
