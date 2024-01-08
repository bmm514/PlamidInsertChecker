# PlamidInsertChecker
Code for SWbioDTP Data Science and Machine Learning module

This is an installable python package that can be used to interrogate DNA sequences for restriction enzyme cutting sites. It can produce: 
    i) a table of cut sites with useful information on the number of cut sites, the cut site locations, the enzyme cut overhang, and the commercial suppliers
    ii) a figure showing the cut locations on a plasmid or linear section of DNA

## Installation

Create a virtual environment using your favourite tool i.e.
```virtualenv env_plasmidin```

Clone the PlasmidInsertChecker repository
```git clone git@github.com:bmm514/PlamidInsertChecker.git```

The install from PyPI using the command:
```pip install plasmidin```

## Basic useage
1a. Import the python library
From in python:
```
from plasmidin.plasmidin import RSInserter
from plasmidin.plasmid_diagrams import PlasmidDrawer

1. Input fasta file or bio python Seq() for plasmid and insert
```
plasmid_seq = 'path/to/plasmid_fasta.fa'
plasmid_linear = False
insert_seq = 'path/to/insert_fasta.fa'
insert_linear = True
remove_ambiguous_enzymes = True
```
2. Run analysis on these to find restriction enzyme cut sites
show the enzymes of interest
```
rsinserter = RSInserter(plasmid_seq, plasmid_linear, insert_seq, insert_linear, remove_ambiguous_enzymes)
rsinserter_XbaI_BamHI.shared_single_enzymes
```

3. Select from the list appropriate restriction enzymes to cut
```
plasmid_cut_enzymes = ('EnzymeA', 'EnzymeB')
insert_cut_enzymes = ('EnzymeA', 'EnzymeB')
```

4. Integrate insert sequence into plasmid sequence 
```
rsinserter.integrate_seq(plasmid_cut_enzymes, insert_cut_enzymes)
```
5. Analyse the output restriction sites and save a table of the restriction sites present
integrated_table = '/path/to/restriction_enzymes.csv'
rsinserter.integrated_rsfinder.save_enzyme_table(integrated_table, delimiter = ',')

6. Create plasmid map for the integrated sequence
```
input_seq = rsinserter.integrated_rsfinder.input_seq
feature_info = rsinserter.integrated_rsfinder.feature_info
integrated_figure = '/path/to/integrated_restriction_map.pdf'

plasmid_drawer = PlasmidDrawer(input_seq, 'IntegratedSeq', feature_info)
plasmid_drawer.draw_gd_diagram(integrated_figure, 'circular', {'pagesize' : 'A4', 'circle_core' : 0.5, 'track_size' : 0.1})
```
## Examples

See plasmidin_example.py for code to run and produce output:
```python3 plasmidin_examples.py```

This will use the example fasta files in the ```data/``` directory and produce output into ```plasmidin_output```