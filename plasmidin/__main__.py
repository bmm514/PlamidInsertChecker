from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, AllEnzymes, Analysis, CommOnly

from plasmidin import cut_and_insert, RSFinder, RSInserter

#Sometimes need to select the correct interpreter in vscode using >python: Select Interpreter then chosing the env
def test_RSFinder():
    input_seq = Seq('ATTTTCTGAATTCGCTAACGTTA')
    dna_seq = Seq('AAAAGAATTCNNNNNNAACGTTTAT')

    linear = True

    rsfinder = RSFinder(input_seq, linear)
    dna_rsfinder = RSFinder(dna_seq, linear)

    # print(rsfinder.single_cut_enzymes)
    # print(rsfinder.all_cut_enzymes)
    # print(rsfinder.enzyme_cut_sites('EcoRI'))
    # print(rsfinder.enzyme_cut_sites('AbaSI'))

    # new_rb = RestrictionBatch(['EcoRI'])
    # rsfinder.change_rb(new_rb, update = True)
    # print(rsfinder.single_cut_enzymes)
    # print(rsfinder.enzyme_cut_sites('EcoRI'))

    # rsfinder.shared_restriction_enzymes(external_rsfinder)
    # rsfinder.shared_restriction_enzymes('', 3)

    # print(rsfinder.filter_enzymes(['EcoRI', 'AbaSI']))
    # print(rsfinder.enzyme_table.head())
    rsfinder.save_enzyme_table('/home/bmm41/PhD_VH/SWbioDTP_taught/DataSciMachLearn/plasmid_info.csv', delimiter = ',')
    dna_rsfinder.save_enzyme_table('/home/bmm41/PhD_VH/SWbioDTP_taught/DataSciMachLearn/dna_info.csv', delimiter = ',')
    # print(rsfinder.supplier_filtered)
    # print(rsfinder.supplier_names)
    rsfinder.filter_supplier({'B'})
    # print(rsfinder.supplier_filtered)
    rsfinder.save_supplier_table('/home/bmm41/PhD_VH/SWbioDTP_taught/DataSciMachLearn/plasmid_info_thermo.csv', delimiter = ',')

    # print(rsfinder.supplier_names)
def test_RSInserter():
    backbone_seq = Seq('ATGTTTAAACTTTCTGAATTCGCTAACGTGTACTA')
    backbone_linear = False
    insert_seq = Seq('AATAGTACTAGATGGTTTAAACAGTCTAAGAATTCNNNNNNNNNAACGTTTAGTACTAT')
    insert_linear = True

    rsinserter = RSInserter(backbone_seq, insert_seq, backbone_linear, insert_linear)
    # print(rsinserter.shared_single_enzymes)
    # print(rsinserter.backbone_single_cut_sites)
    # print(rsinserter.insert_single_cut_sites)
    # print('#'*20)
    # rsinserter.inegrate_seq(('EcoRI', 'MaeII'), ('EcoRI', 'MaeII'))
    # print(rsinserter.integrated_rsfinder.input_seq)
    # print('#'*20)
    # rsinserter.inegrate_seq(('EcoRI', 'MaeII'), ('MaeII', 'EcoRI'))
    # print(rsinserter.integrated_rsfinder.input_seq)
    # print('#'*20)   
    # rsinserter.inegrate_seq(('MaeII', 'EcoRI'), ('EcoRI', 'MaeII'))
    # print(rsinserter.integrated_rsfinder.input_seq)
    # print('#'*20)
    # rsinserter.inegrate_seq(('MaeII', 'EcoRI'), ('MaeII', 'EcoRI'))
    # print(rsinserter.integrated_rsfinder.input_seq)
    #include test for same insert site, need to find one that is not ambiguous!
    print('#'*20)  
    print(rsinserter.insert_rsfinder.n_cut_sites(2))
    print(rsinserter.insert_rsfinder.single_cut_enzymes)
    print(rsinserter.backbone_rsfinder.single_cut_enzymes)
    print(rsinserter.backbone_rsfinder.single_cut_enzymes['MssI'])
    print(rsinserter.shared_single_enzymes)
    #The problem is with MssI not being shared between the backbone and the insert sequences!!
    #Need a fix for the above but also need to find better restriction enzymes to test (AfaI and ScaI have are shared)
    rsinserter.inegrate_seq(('MssI', 'AfaI'), ('ScaI', 'ScaI'))
    print(rsinserter.integrated_rsfinder.input_seq)

if __name__ == '__main__':
    # test_main()
    # test_RSFinder()
    test_RSInserter()