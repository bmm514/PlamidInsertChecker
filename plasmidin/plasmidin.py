from Bio.Restriction import Analysis, RestrictionBatch, CommOnly
import Bio.Restriction
from Bio.Seq import Seq
def find_restriction_sites(sequence: Seq, rb = RestrictionBatch(CommOnly), linear = True):
    """
    Find restriction sites from a DNA sequence
    
    sequence - Bio.Seq.Seq
        an object that contains the DNA sequence to be searched
    rb - Bio.Restriction.RestrictionBatch
        an object holding the restriction enzyme sites to search.
            Default is commerically available enzymes (Bio.Restriction.CommOnly)
    linear - bool
        whether the DNA sequence should be treated as linear or circular
    
    return a dictioanry of {restriction enzyme : cut_locations} that have > 1 cut site
    """
    analysis = Analysis(rb, sequence, linear)

    return analysis
def filter_seqs(backbone_seq, backbone_linear, insertion_seq, insertion_linear, rb = RestrictionBatch(CommOnly)):
    analysis_backbone = Analysis(rb, backbone_seq, backbone_linear)
    analysis_insertion = Analysis(rb, insertion_seq, insertion_linear)
    shared_restriction_enzymes = shared_restriction_sites(analysis_backbone, analysis_insertion)

    filter_backbone_rs, filter_insertion_rs = [filter_restriction_sites(analysis.with_N_sites(1), shared_restriction_enzymes) for analysis in [analysis_backbone, analysis_insertion]]

    return filter_backbone_rs, filter_insertion_rs

def shared_restriction_sites(analysis_1, analysis_2):
    restriction_sites_1 = analysis_1.with_sites()
    restriction_sites_2 = analysis_2.with_sites()
    shared_keys = set(restriction_sites_1.keys()) & set(restriction_sites_2.keys())

    return shared_keys

def filter_restriction_sites(restriction_sites, restriction_enzymes):
    filtered_sites = {key.__name__: value for key, value in restriction_sites.items() if key in restriction_enzymes}
    return filtered_sites

def cut_and_insert(backbone_seq, backbone_rs, backbone_enzymes, insertion_seq, insertion_rs, insertion_enzymes):
    backbone_lhs, _, backbone_rhs = cut_enzymes(backbone_seq, backbone_rs, backbone_enzymes)
    _, insertion_middle, _ = cut_enzymes(insertion_seq, insertion_rs, insertion_enzymes)

    return backbone_lhs + insertion_middle + backbone_rhs


def cut_enzymes(seq: Seq, restriction_sites: dict, enzymes: tuple):
    """
    Cut a Seq object by sites for two enzymes specified in (5' cut, 3' cut) order.
    There should only be one cut site for the restriction enzyme otherwise throw an error
    """
    if len(enzymes) == 2:
        lhs_enzyme = enzymes[0]
        rhs_enzyme = enzymes[1]
    elif len(enzymes) == 1:
        lhs_enzyme = rhs_enzyme = enzymes[0]
    else:
        raise TypeError('enzymes must have one or two inputs')

    lhs_sites = restriction_sites[lhs_enzyme]
    rhs_sites = restriction_sites[rhs_enzyme]

    if (len(lhs_sites) > 1) or (len(rhs_sites) > 1):
        raise TypeError('The enzyme selected has two cut sites to choose from. This will not create an unambiguous solution')
    
    lhs_cut = lhs_sites[0]
    rhs_cut = rhs_sites[0]
    
    lhs_seq, middle_seq, rhs_seq = seq[:lhs_cut], seq[lhs_cut:rhs_cut], seq[rhs_cut:]
    return lhs_seq, middle_seq, rhs_seq

def main():
    plasmid_seq = Seq('ATTTTCTGAATTCGCTAACGTTA')
    dna_seq = Seq('AAAAGAATTCNNNNNNAACGTTTAT')

    filter_plasmid_rs, filter_dna_rs = filter_seqs(plasmid_seq, False, dna_seq, True)
    # lhs_seq, middle_seq, rhs_seq = cut_enzymes(plasmid_seq, filtered_plasmid_rs, ('EcoRI', 'AclI'))
    seq = cut_and_insert(
        plasmid_seq, filter_plasmid_rs, ('EcoRI', 'AclI'),
        dna_seq, filter_dna_rs, ('EcoRI', 'AclI')
        )

    print(plasmid_seq, dna_seq)
    print(seq)

if __name__ == '__main__':
    main()