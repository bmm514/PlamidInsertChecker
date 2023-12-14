from plasmidin import RSInserter

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

#Create plasmid maps


#Running for two XbaI in insert
insert_seq = 'examples/insert_XbaI_XbaI.fa'
insert_linear = True

rsinserter_XbaI_XbaI = RSInserter(plasmid_seq, insert_seq, plasmid_linear, insert_linear, remove_ambiguous = remove_ambiguous_enzymes)

#Look at enzymes shared
print(rsinserter_XbaI_XbaI.shared_single_enzymes)
#Chose enzymes to cut
plasmid_cut_enzymes = ('XbaI', 'XbaI')
insert_cut_enzymes = ('XbaI', 'XbaI')

#Insert sequence into plasmid
rsinserter_XbaI_XbaI.integrate_seq(plasmid_cut_enzymes, insert_cut_enzymes,1,2)
