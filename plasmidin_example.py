from os import path
from os import makedirs
from plasmidin.plasmidin import RSInserter
from plasmidin.plasmid_diagrams import PlasmidDrawer

#Quick functions
def add_dirpath(outdir, filename):
    return path.join(outdir, filename)

#1) Set up the folder if necessary:
outdir = 'plasmidin_output'
if not path.isdir(outdir):
    makedirs(outdir)

#2) Set up the files to run
plasmid_seq = 'data/pUC19_plasmid.fa'
plasmid_linear = False
insert_seq = 'data/insert_XbaI_BamHI.fa'
insert_linear = True
remove_ambiguous_enzymes = True

#3) Create the RSInserter
rsinserter_XbaI_BamHI = RSInserter(plasmid_seq, insert_seq, plasmid_linear, insert_linear, remove_ambiguous = remove_ambiguous_enzymes)

#Look at enzymes shared
print(rsinserter_XbaI_BamHI.shared_single_enzymes)
#Chose enzymes to cut
plasmid_cut_enzymes = ('XbaI', 'BamHI')
insert_cut_enzymes = ('XbaI', 'BamHI')

#4) Insert sequence into plasmid
rsinserter_XbaI_BamHI.integrate_seq(plasmid_cut_enzymes, insert_cut_enzymes)

#5) Create output tables
#a) Make output files
plasmid_table = add_dirpath(outdir, 'pUC19_restriction_enzymes.csv')
insert_table = add_dirpath(outdir, 'myinsert_restriction_enzymes.csv')
integrated_table = add_dirpath(outdir, 'integrated_pUC19_myinsert_restriction_enzymes.csv')

#b) Create tables
rsinserter_XbaI_BamHI.backbone_rsfinder.save_enzyme_table(plasmid_table, delimiter = ',')
rsinserter_XbaI_BamHI.insert_rsfinder.save_enzyme_table(insert_table, delimiter = ',')
rsinserter_XbaI_BamHI.integrated_rsfinder.save_enzyme_table(integrated_table, delimiter = ',')

#6) Create plasmid maps
#a) Create the records to plot for the plasmid...
rsinserter_XbaI_BamHI.backbone_rsfinder.create_enzyme_records(max_n_cut_sites = 2)
#Plot them using the input_seq created previously and the feature_info generated in the backgroud...
pUC19_input_seq = rsinserter_XbaI_BamHI.backbone_rsfinder.input_seq
pUC19_feature_info = rsinserter_XbaI_BamHI.backbone_rsfinder.feature_info
pUC19_figure = add_dirpath(outdir, 'pUC19_restriction_map.pdf')

plasmid_drawer = PlasmidDrawer(pUC19_input_seq, 'pUC19', pUC19_feature_info)
plasmid_drawer.draw_gd_diagram(pUC19_figure, 'circular', {'pagesize' : 'A4', 'circle_core' : 0.5, 'track_size' : 0.1})

#b) Create the records to plot for the insert...
rsinserter_XbaI_BamHI.insert_rsfinder.create_enzyme_records(max_n_cut_sites = 2)
#Plot them using the input_seq created previously and the feature_info generated in the backgroud...
myinsert_input_seq = rsinserter_XbaI_BamHI.insert_rsfinder.input_seq
myinsert_feature_info = rsinserter_XbaI_BamHI.insert_rsfinder.feature_info
myinsert_figure = add_dirpath(outdir, 'myinsert_restriction_map.pdf')

plasmid_drawer = PlasmidDrawer(myinsert_input_seq, 'insert_XbaI_BamHI', myinsert_feature_info)
# plasmid_drawer.draw_gd_diagram(myinsert_figure, 'linear', {'pagesize' : 'A4', 'circle_core' : 0.5, 'track_size' : 1.0})
plasmid_drawer.draw_gd_diagram(myinsert_figure, 'linear', {'pagesize' : 'A4', 'track_size' : 1.0})

#c) Create the records to plot for the integrated sequence...
rsinserter_XbaI_BamHI.integrated_rsfinder.create_enzyme_records(max_n_cut_sites = 2)
#Plot them using the input_seq created previously and the feature_info generated in the backgroud...
integrated_input_seq = rsinserter_XbaI_BamHI.integrated_rsfinder.input_seq
integrated_feature_info = rsinserter_XbaI_BamHI.integrated_rsfinder.feature_info
integrated_figure = add_dirpath(outdir, 'integrated_pUC19_myinsert_restriction_map.pdf')

plasmid_drawer = PlasmidDrawer(integrated_input_seq, 'plasmid_insert_XbaI_BamHI', integrated_feature_info)
plasmid_drawer.draw_gd_diagram(integrated_figure, 'circular', {'pagesize' : 'A4', 'circle_core' : 0.5, 'track_size' : 0.1})

#Running for two XbaI in insert
# insert_seq = 'examples/insert_XbaI_XbaI.fa'
# insert_linear = True

# rsinserter_XbaI_XbaI = RSInserter(plasmid_seq, insert_seq, plasmid_linear, insert_linear, remove_ambiguous = remove_ambiguous_enzymes)

# #Look at enzymes shared
# print(rsinserter_XbaI_XbaI.shared_single_enzymes)
# #Chose enzymes to cut
# plasmid_cut_enzymes = ('XbaI', 'XbaI')
# insert_cut_enzymes = ('XbaI', 'XbaI')

# #Insert sequence into plasmid
# rsinserter_XbaI_XbaI.integrate_seq(plasmid_cut_enzymes, insert_cut_enzymes,1,2)
