# VarChAMP Allele Collection

This README.md documents the detailed information about the allele collection of VarChAMP project we have generated up to December 31st, 2025. This is the complete allele collection we are planning to publish and share with the scientific community in our 1% VarChAMP manuscript. 

The file allele collection file to be used for all the analyses are located at:
```varchamp_seq_confirm_alleles_clinvar_gnomad_seq_final.tsv.gz``` and ```VarChampSeqConfirmationResult_final.tsv.gz```. 

```varchamp_seq_confirm_alleles_clinvar_gnomad_seq_final.tsv.gz``` includes the complete nucleotide and amino acid sequences of all the variants we have generated in our library, while ```VarChampSeqConfirmationResult_final.tsv.gz``` contains the allele sequence confirmation results we used for conducting the quality control of these variants.

## Email Notes from Tong 

### Dec. 23, 2025

```
Luke just found out that the SPDIs I generated for our variants is off by one base. This is due to that SPDI is 0 based while genomic positions from clinivar are 1 based. I’ve fixed the SPDI column in the files attached.
```

### Dec. 10, 2025

```
I've provided missing Ensembl gene IDs and Ensembl protein IDs in file attached. 
Please note that the Ensembl protein IDs provided are the best matched Ensemble proteins for the ORFs and not necessarily with 100% identity. You’ll still see the mismatches. 
```