import pandas
from Bio.Restriction import Analysis, RestrictionBatch, CommOnly, AllEnzymes
from Bio.Seq import Seq

#Should I put this in RSFinder?
def enzyme_dict_to_string(n_cut_enzymes: dict):
    """Convert an analysis dictionary enzyme objects to the string name"""
    new_n_cut_enzymes = {}
    for enzyme, values in n_cut_enzymes.items():
        new_n_cut_enzymes[enzyme.__name__] = values
    
    return new_n_cut_enzymes

class RSFinder():
    """
    A class to find restriction enzyme sites within an input sequence
    """
    #Output can be used to compare common restriction sites
    def __init__(self, input_seq, linear: bool, rb = RestrictionBatch(CommOnly)):
        """
        input_seq - a Bio.Seq.Seq object
        linear_seq - boolean for whether the sequence is treated as linear or circular
        rb - the Bio.Restriction.RestrictionBatch to use. Defaults to commercially availably restriction enzymes
        """
        self._input_seq = input_seq
        self._linear = linear #do not want to be able to change in the class
        self._rb = rb

        self._analysis = self.restriction_site_analysis()
        self._single_cut_enzymes = self.single_cut_site()
        self._all_cut_enzymes = self.any_cut_sites()
        self._enzyme_table = self.create_restriction_enzyme_table()

        self._supplier_names = set()
        self._supplier_codes = set()
        self._supplier_filtered = None
        self._supplier_table = None

        self._supplier_codes_dict = {
            'B': 'Thermo Fisher Scientific',
            'C': 'Minotech Biotechnology',
            'E': 'Agilent Technologies',
            'F': 'Thermo Fisher Scientific',
            'I': 'SibEnzyme Ltd.',
            'J': 'Nippon Gene Co., Ltd.',
            'K': 'Takara Bio Inc.',
            'M': 'Roche Applied Science',
            'N': 'New England Biolabs',
            'O': 'Toyobo Biochemicals',
            'Q': 'Molecular Biology Resources - CHIMERx',
            'R': 'Promega Corporation',
            'S': 'Sigma Chemical Corporation',
            'V': 'Vivantis Technologies',
            'X': 'EURx Ltd.',
            'Y': 'SinaClon BioScience Co.'
        }

    @property
    def input_seq(self):
        return self._input_seq
    
    @property
    def linear(self):
        return self._linear
    
    @property
    def rb(self):
        return self._rb
    
    @rb.setter
    def rb(self, rb):
        if isinstance(rb, RestrictionBatch):
            self._rb = rb
        else:
            raise ValueError(f'rb is not a Bio.Restriction.RestrictionBatch object so not updating')
    
    @property
    def analysis(self):
        return self._analysis
    
    @property
    def single_cut_enzymes(self):
        return self._single_cut_enzymes
    
    @property
    def all_cut_enzymes(self):
        return self._all_cut_enzymes
    
    @property
    def enzyme_table(self):
        return self._enzyme_table
    
    @property
    def supplier_filtered(self):
        if self._supplier_filtered is None:
            print('There is no current RSFinder.supplier_filtered')
        return self._supplier_filtered
    
    @property
    def supplier_table(self):
        if self._supplier_table is None:
            print("There is no current RSFinder.supplier_table")
        return self._supplier_table
    
    @property
    def supplier_names(self):
        if self._supplier_names is set():
            print('There are no current suppliers that have been selected to be filtered')
        return self._supplier_names
    
    def change_rb(self, rb, update = True):
        """
        Change the RestrictionBatch held within RSFinder and update RSFinder._analysis,
        RSFinder._single_cut_enzymes and RSFinder.all_cut_enzymes
        """
        self.rb = rb
        if update:
            input_seq = self.input_seq
            linear = self.linear
            rb = self.rb
            self.__init__(input_seq, linear, rb)

    def restriction_site_analysis(self):
        """
        Run the Bio.Restriction.Analysis on self.input_seq"""
        rb = self.rb
        input_seq = self.input_seq
        linear = self.linear

        return Analysis(rb, input_seq, linear)
    
    def any_cut_sites(self):
        """Return the enzymes with any number of cuts in the input_seq"""
        any_cut_enzymes = self._analysis.with_sites()
        new_all_cut_enzymes = enzyme_dict_to_string(any_cut_enzymes)

        return new_all_cut_enzymes
    
    def n_cut_sites(self, n_sites):
        """Return the ezymes with n_sites number of cuts in the input_seq"""
        analysis = self._analysis
        n_cut_enzymes = analysis.with_N_sites(n_sites)
        new_n_cut_enzymes = enzyme_dict_to_string(n_cut_enzymes)
        
        return new_n_cut_enzymes

    def single_cut_site(self):
        """Return enzymes with a single cut site"""
        return self.n_cut_sites(1)
    
    def enzyme_cut_sites(self, restriction_enzyme):
        """Return the cut sites for the enzyme specified"""
        all_cut_enzymes = self._all_cut_enzymes
        try:
            return all_cut_enzymes[restriction_enzyme]
        except KeyError:
            print(f'Could not find {restriction_enzyme} in dictionary. Returning []')
            return []
    
    def _select_enzymes(self, n_cut_sites):
        if n_cut_sites is None:
            cut_enzymes = self.all_cut_enzymes
        elif n_cut_sites == 1: #saving computation
            cut_enzymes = self.single_cut_enzymes
        else:
            cut_enzymes = self.n_cut_sites(n_cut_sites)
        
        return cut_enzymes

    
    def filter_enzymes(self, restriction_enzymes, n_cut_sites = None):
        """
        Return a dictionary of {restriction_enzyme : cut_sites} for every restriction enzyme in restriction_enzymes.
        Use n_cut_sites to limit to specified number of cut sites, if None then do all cut sites
        Any missing keys will be skipped

        restriction_enzymes can be any iterable
        """
        cut_enzymes = self._select_enzymes(n_cut_sites)
        filtered_enzymes = {}
        for restriction_enzyme in restriction_enzymes:
            try:
                cut_sites = cut_enzymes[restriction_enzyme]
                filtered_enzymes[restriction_enzyme] = cut_sites
            except KeyError:
                print(f'Could not find {restriction_enzyme} in dictionary. Skipping {restriction_enzyme}')

        return filtered_enzymes

    def shared_restriction_enzymes(self, rsfinder, n_cut_sites = 1):
        """Extract the restriction enzymes in rsfinder that share the same sites with the current RSFinder"""
        if not isinstance(rsfinder, RSFinder):
            raise TypeError(f'rsfinder is not an RSFinder class')
        
        internal_enzymes = self._select_enzymes(n_cut_sites)
        external_enzymes = rsfinder._select_enzymes(n_cut_sites)

        shared_enzymes = set(internal_enzymes.keys()) & set(external_enzymes.keys())

        return shared_enzymes
    
    def _make_table(self, enzyme_dict):
        """Extract useful information from the restriction enzymes in enzyme_dict and turn into a dataframe"""
        data = []
        for enzyme_name, values in enzyme_dict.items():
            enzyme = self.rb.get(enzyme_name)
            n_sites = len(values)
            cut_locations = values
            #Can add to the dictionary if something would be useful
            data.append({
                'Name' : enzyme_name,
                'N_sites' : n_sites,
                'Cut_Locations' : '; '.join(map(str,cut_locations)),
                'Cut_type' : enzyme.overhang(),
                'CommerciallyAvailable' : enzyme.is_comm(),
                'Suppliers' : '; '.join(enzyme.supplier_list()),
            }
            )

        enzyme_df = pandas.DataFrame(columns = data[0].keys())
        for row in data:
            enzyme_df = enzyme_df._append(row, ignore_index = True)
        
        return enzyme_df


    def create_restriction_enzyme_table(self, n_cut_sites = None):
        """Take an Analysis object and create a table containing information on the Restriction Sites"""
        cut_enzymes = self._select_enzymes(n_cut_sites)
        enzyme_df = self._make_table(cut_enzymes)
        
        return enzyme_df
    
    def _save_table(self, df, table_out, delimiter = '\t'):
        df.to_csv(table_out, sep = delimiter, index = False)
    
    def save_enzyme_table(self, table_out, delimiter = '\t'):
        df = self.enzyme_table
        self._save_table(df, table_out, delimiter)
    
    def save_supplier_table(self, table_out, delimiter = '\t'):
        df = self.supplier_table
        if df is not None:
            self._save_table(df, table_out, delimiter)
        else:
            raise TypeError(f'There is no RSFinder.supplier_table present. Make one with RSFinder.filter_supplier')

    def filter_supplier(self, supplier_codes, n_cut_sites = None):
        """
        Select from a supplier code from below to filter(s) out enzyme that are present:

        'B': 'Thermo Fisher Scientific',
        'C': 'Minotech Biotechnology',
        'E': 'Agilent Technologies',
        'I': 'SibEnzyme Ltd.',
        'J': 'Nippon Gene Co., Ltd.',
        'K': 'Takara Bio Inc.',
        'M': 'Roche Applied Science',
        'N': 'New England Biolabs',
        'O': 'Toyobo Biochemicals',
        'Q': 'Molecular Biology Resources - CHIMERx',
        'R': 'Promega Corporation',
        'S': 'Sigma Chemical Corporation',
        'V': 'Vivantis Technologies',
        'X': 'EURx Ltd.',
        'Y': 'SinaClon BioScience Co.'

        n_cut_sites can be None or an integer to select enzymes with specific number of cut sites
        """
        cut_enzymes = self._select_enzymes(n_cut_sites)
        supplier_filtered = {}
        self._supplier_names = {self._supplier_codes_dict[code] for code in supplier_codes}
        self._supplier_codes = set(supplier_codes)
        
        for enzyme_name, values in cut_enzymes.items():
            enzyme_suppliers = set(self.rb.get(enzyme_name).suppl)
            retain_codes = set(supplier_codes)
            intersection = enzyme_suppliers.intersection(retain_codes)
            if len(intersection) > 0:
                supplier_filtered[enzyme_name] = values
        
        self._supplier_table = self._make_table(supplier_filtered)
        self._supplier_filtered = supplier_filtered

#This should be part of a seperate class that takes in 2 sequences i.e. plasmid and insert
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