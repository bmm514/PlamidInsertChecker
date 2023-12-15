# PlamidInsertChecker
Code for SWbioDTP Data Science and Machine Learning module

This is an installable python package that can be used to interrogate DNA sequences for restriction enzyme cutting sites. 

## Installation

```pip install plasmidin```

## Useage
1a. Import the python library
```from plasmidin import RSInserter```

1. Input fasta file or bio python Seq() for plasmid and insert
```
plasmid_seq = 'path/to/plasmid_fasta.fa'
plasmid_linear = False
insert_seq = 'path/to/insert_fasta.fa'
insert_linear = True
remove_ambiguous_enzymes = True

rsinserter = RSInserter(plasmid_seq, plasmid_linear, insert_seq, insert_linear)
```
2. Run analysis on these to find restriction enzyme cut sites
show the enzymes of interest
3. Select from the list appropriate restriction enzymes to cut
```
plasmid_cut_enzymes = ('EcoRI', 'EcoRI')
insert_cut_enzymes = ('A', 'B')
```
4. Integrate insert sequence into plasmid sequence 
```
rsinserter.integrate_seq(plasmid_cut_enzymes, insert_cut_enzymes)
```
5. Analyse the output restriction sites and full sequence

6. Create plasmid maps and tables

## Examples

See example.py for code to run