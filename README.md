# PlamidInsertChecker
Code for SWbioDTP Data Science and Machine Learning module

This is an installable python package that can be used to interrogate DNA sequences for restriction enzyme cutting sites. It can produce:  
    i. a table of cut sites with useful information on the number of cut sites, the cut site locations, the enzyme cut overhang, and the commercial suppliers  
    ii. a figure showing the cut locations on a plasmid or linear section of DNA

## Installation

Create a virtual environment using your favourite tool, clone the repsotiory to obtain the data and install the package using pip

```
virtualenv env_plasmidin
source env_plasmidin/bin/activate
git clone git@github.com:bmm514/PlamidInsertChecker.git
pip install plasmidin
```

## Basic useage
1. Import required plasmidin Classes
```
from plasmidin.plasmidin import RSInserter
from plasmidin.plasmid_diagrams import PlasmidDrawer
```

2. Input fasta file or bio python Seq() for plasmid and insert
```
plasmid_seq = 'path/to/plasmid_fasta.fa'
plasmid_linear = False
insert_seq = 'path/to/insert_fasta.fa'
insert_linear = True
remove_ambiguous_enzymes = True
```

3. Run analysis on these to find restriction enzyme cut sites and 
show the enzymes of interest i.e. single cut enzymes
```
rsinserter = RSInserter(plasmid_seq, plasmid_linear, insert_seq, insert_linear, remove_ambiguous_enzymes)

rsinserter_XbaI_BamHI.shared_single_enzymes
```

4. Select from the list appropriate restriction enzymes to cut
```
plasmid_cut_enzymes = ('EnzymeA', 'EnzymeB')
insert_cut_enzymes = ('EnzymeA', 'EnzymeB')
```

5. Integrate insert sequence into plasmid sequence 
```
rsinserter.integrate_seq(plasmid_cut_enzymes, insert_cut_enzymes)
```

6. Analyse the output restriction sites and save a table of the restriction sites present
```
integrated_table = '/path/to/restriction_enzymes.csv'

rsinserter.integrated_rsfinder.save_enzyme_table(integrated_table, delimiter = ',')
```

7. Create plasmid map for the integrated sequence
```
input_seq = rsinserter.integrated_rsfinder.input_seq
feature_info = rsinserter.integrated_rsfinder.feature_info
integrated_figure = '/path/to/integrated_restriction_map.pdf'

plasmid_drawer = PlasmidDrawer(input_seq, 'IntegratedSeq', feature_info)
plasmid_drawer.draw_gd_diagram(integrated_figure, 'circular', {'pagesize' : 'A4', 'circle_core' : 0.5, 'track_size' : 0.1})
```
## Examples

See plasmidin_example.py for code to run and produce output from PlasmidInsertChecker:
```
python3 scripts/plasmidin_examples.py
```

This will use the example fasta files in ```data/``` directory and produce output into ```plasmidin_output```