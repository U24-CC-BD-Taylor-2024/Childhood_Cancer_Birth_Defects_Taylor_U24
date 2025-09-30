# Submission of MTP variants to the Allele Registry

This module outlines steps to obtain HGVSg IDs for all MTP variants, and submit IDs to the Allele Registry

# Run Rscript to obtain HGVSg IDs

```bash
cd preprocessing/allele-regsitry
Rscript --vanilla get-mtp-variant-HGVSg-ids.R
```

# Create Allele Registry account

- Navigate to [ClinGen Allele Registry website](https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/landing) and create a count by clicking "Register"
- Create username and password to be used in subsequent steps

# Submit HGVSg IDs to Allele Registry via API

Specifications for submitting to Allele Registry can be found [here](https://ldh.clinicalgenome.org/redmine/projects/clingen-allele-registry/wiki/ClinGen_Allele_Registry_Performance_and_Guidelines)

```bash
bash request_with_payload.sh "http://reg.test.genome.network/alleles?file=hgvs" GMKF-variants.tsv <USERNAME> <PASSWORD> > post_test_GMKF_variants_out.json

bash request_with_payload.sh "http://reg.test.genome.network/alleles?file=hgvs" PBTA-variants.tsv <USERNAME> <PASSWORD> > post_test_PBTA_variants_out.json

bash request_with_payload.sh "http://reg.test.genome.network/alleles?file=hgvs" TARGET-variants.tsv <USERNAME> <PASSWORD> > post_test_TARGET_variants_out.json
```