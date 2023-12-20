from plasmidin.plasmidin import RSInserter
from plasmidin.plasmid_diagrams import PlasmidDrawer

#Set up the files to run
plasmid_seq = 'examples/plasmid.fa'
plasmid_linear = False
insert_seq = 'examples/insert_XbaI_BamHI.fa'
insert_linear = True
remove_ambiguous_enzymes = True

#Create the RSInserter
rsinserter_XbaI_BamHI = RSInserter(plasmid_seq, insert_seq, plasmid_linear, insert_linear, remove_ambiguous = remove_ambiguous_enzymes)

#Look at enzymes shared
print(rsinserter_XbaI_BamHI.shared_single_enzymes)
#Chose enzymes to cut
plasmid_cut_enzymes = ('XbaI', 'BamHI')
insert_cut_enzymes = ('XbaI', 'BamHI')

#Insert sequence into plasmid
rsinserter_XbaI_BamHI.integrate_seq(plasmid_cut_enzymes, insert_cut_enzymes)

rsinserter_XbaI_BamHI.backbone_rsfinder.save_enzyme_table('examples/plasmid_enzymes.csv', delimiter = ',')
rsinserter_XbaI_BamHI.insert_rsfinder.save_enzyme_table('examples/insert_enzymes.csv', delimiter = ',')
rsinserter_XbaI_BamHI.integrated_rsfinder.save_enzyme_table('examples/integrated_enzymes.csv', delimiter = ',')

#Create plasmid maps...

#Create the records to plot for the plasmid...
rsinserter_XbaI_BamHI.backbone_rsfinder.create_enzyme_records(max_n_cut_sites = 2)
#Plot them...
plasmid_drawer = PlasmidDrawer(rsinserter_XbaI_BamHI.backbone_rsfinder.input_seq, 'pBR322', rsinserter_XbaI_BamHI.backbone_rsfinder.feature_info)
plasmid_drawer.draw_gd_diagram('examples/pBR332_restriction_map.pdf', 'circular', {'pagesize' : 'A4', 'circle_core' : 0.5, 'track_size' : 0.1})

#Create the records to plot for the insert...
rsinserter_XbaI_BamHI.insert_rsfinder.create_enzyme_records(max_n_cut_sites = 2)
#Plot them...
plasmid_drawer = PlasmidDrawer(rsinserter_XbaI_BamHI.insert_rsfinder.input_seq, 'insert_XbaI_BamHI', rsinserter_XbaI_BamHI.insert_rsfinder.feature_info)
plasmid_drawer.draw_gd_diagram('examples/insert_XbaI_BamHI_restriction_map.pdf', 'linear', {'pagesize' : 'A4', 'circle_core' : 0.5, 'track_size' : 1.0})

#Create the records to plot for the integrated...
rsinserter_XbaI_BamHI.integrated_rsfinder.create_enzyme_records(max_n_cut_sites = 2)
#Plot them...
plasmid_drawer = PlasmidDrawer(rsinserter_XbaI_BamHI.integrated_rsfinder.input_seq, 'plasmid_insert_XbaI_BamHI', rsinserter_XbaI_BamHI.integrated_rsfinder.feature_info)
plasmid_drawer.draw_gd_diagram('examples/plasmid_insert_XbaI_BamHI_restriction_map.pdf', 'circular', {'pagesize' : 'A4', 'circle_core' : 0.5, 'track_size' : 0.1})

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
