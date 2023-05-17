## Assessing assembly quality
Annotate metazoan-conserved genes from BUSCO 5.2.2 database with Metaeuk gene predictor (Manni et al., 2021):
```
busco -i assemblies/ -o buscoout -l metazoa_odb10 --update-data -m genome --cpu 30
```
Calculate assembly metrics with Quast 5.0.2 (Gurevich et al., 2013):
```
python3 ~/softwares/quast-5.0.2/quast.py assembly.fa -e -s -o assembly.fa_quast_report
```
