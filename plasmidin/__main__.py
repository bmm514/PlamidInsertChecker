
from plasmidin import filter_seqs, cut_and_insert, RSFinder

def test_RSFinder():
    from Bio.Seq import Seq
    from Bio.Restriction import RestrictionBatch, AllEnzymes

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

def test_main():
    from Bio.Seq import Seq
    from Bio.Restriction import Analysis, RestrictionBatch, CommOnly
    plasmid_seq = Seq('ATTTTCTGAATTCGCTAACGTTA')
    dna_seq = Seq('AAAAGAATTCNNNNNNAACGTTTAT')

    (analysis_plasmid, analysis_dna), (filter_plasmid_rs, filter_dna_rs), shared_enzymes = filter_seqs(plasmid_seq, False, dna_seq, True)
    
    inserted_seq = cut_and_insert(
        plasmid_seq, filter_plasmid_rs, ('EcoRI', 'AclI'),
        dna_seq, filter_dna_rs, ('EcoRI', 'AclI')
        )

    csvout_plasmid = '/home/bmm41/PhD_VH/SWbioDTP_taught/DataSciMachLearn/plasmid_info.csv'
    csvout_dna = '/home/bmm41/PhD_VH/SWbioDTP_taught/DataSciMachLearn/dna_info.csv'
    csvout_inserted = '/home/bmm41/PhD_VH/SWbioDTP_taught/DataSciMachLearn/inserted_info.csv'

    # make_restriction_enzyme_table(analysis_plasmid, csvout_plasmid, shared_enzymes)
    # make_restriction_enzyme_table(analysis_dna, csvout_dna, shared_enzymes)

    analysis_inserted = Analysis(RestrictionBatch(CommOnly), inserted_seq, True)
    # make_restriction_enzyme_table(analysis_inserted, csvout_inserted, {})

    return analysis_plasmid, analysis_dna

if __name__ == '__main__':
    # test_main()
    test_RSFinder()