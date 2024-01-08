from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.Seq import Seq

from reportlab.lib import colors
from reportlab.lib.units import cm

class PlasmidDrawer():
    def __init__(self, seq: Seq, seq_id: str, feature_info, coodinate_step = 500):
        self._seq = seq
        self._seq_id = seq_id
        self._seq_length = len(seq)
        self._coordinate_step = coodinate_step
        self._coordinate_track = 3
        self._feature_info = feature_info

        self._init_gd_diagram()
        self._init_features()

    @property
    def seq(self):
        return self._seq
    
    @property
    def seq_id(self):
        return self._seq_id
    
    @property
    def seq_length(self):
        return self._seq_length
    
    @property
    def feature_info(self):
        return self._feature_info
    
    @property
    def gd_diagram(self):
        return self._gd_diagram
    
    @property
    def gd_feature_set(self):
        return self._gd_feature_set
    
    def _init_gd_diagram(self):
        seq_id = self.seq_id
        seq_length = self.seq_length
        coord_step = self._coordinate_step

        self._gd_diagram = GenomeDiagram.Diagram(seq_id)

        self._gd_track_for_coordinates = self._gd_diagram.new_track(self._coordinate_track, name = 'Coodinates')
        self._gd_coordinate_set = self._gd_track_for_coordinates.new_set()
        for start_coord in range(0, seq_length, coord_step):
            feature = SeqFeature(SimpleLocation(start_coord, start_coord + 1))
            if start_coord == 0:
                coord_name = f'1 / {seq_length + 1}'
            else:
                coord_name = str(start_coord)
            self._gd_coordinate_set.add_feature(
                feature, 
                name = coord_name, 
                sigil = 'BOX',
                color = colors.black,
                label = True,
                label_size = 11,
                label_angle = 0
                )

        self._gd_track_for_features = self._gd_diagram.new_track(1, name = 'Restriction enzymes')
        self._gd_feature_set = self._gd_track_for_features.new_set()            

    def _init_features(self):
        feature_info = self._feature_info

        for feature, info in feature_info.values():
            # print(feature)
            # print(info)
            self.add_gd_feature(feature, info)

    def add_gd_feature(self, feature, info):
        feature_name = info.get('feature_name', '')
        sigil = info.get('sigil', 'BOX')
        color = info.get('color', colors.black)
        label = info.get('label', True)
        label_size = info.get('label_size', 11)
        label_angle = info.get('label_angle', 0)

        self.gd_feature_set.add_feature(
            feature,
            name = feature_name,
            sigil = sigil,
            color = color,
            label = label,
            label_size = label_size,
            label_angle = label_angle,
            )

    def remove_gd_feature(self, cut_site):
        raise NotImplementedError
    
    def draw_gd_diagram(self, diagram_file, diagram_format, draw_settings, filetype = 'PDF'):
        pagesize = draw_settings.get('pagesize', 'A4')
        start = 0
        end = self.seq_length

        if diagram_format == 'linear':
            self.gd_diagram.del_track(self._coordinate_track) #Remove the coordinates track
            self.gd_diagram.draw(
                format = diagram_format, 
                circular = False, 
                pagesize = pagesize,
                fragments = (self.seq_length // self._coordinate_step) + 1,
                start = start,
                end = end
                )
        elif diagram_format == 'circular':
            # pagesize = draw_settings.get('pagesize', 'A4')
            circle_core = draw_settings.get('circle_core', 0.5)
            track_size = draw_settings.get('track_size', 0.5)
            # start = 0
            # end = self.seq_length
            self.gd_diagram.draw(
                format = diagram_format, 
                circular = True, 
                pagesize = pagesize, 
                circle_core = circle_core,
                track_size = track_size,
                start = start,
                end = end
                )
        else:
            raise NotImplementedError('Please choose from linear or circular')
        
        self.gd_diagram.write(diagram_file, filetype)

def main():
    feature_info = (
        (
            SeqFeature(SimpleLocation(100,101)), {
            'feature_name' : 'feature1',
            'sigil' : 'BOX',
            'color' : colors.pink,
            'label' : True,
            'label_size' : 14,
            'label_angle' : 45
            }
        ),
        (
            SeqFeature(SimpleLocation(200,225)), {
            'feature_name' : 'feature2',
            'sigil' : 'BOX',
            'color' : colors.lightblue,
            'label' : True,
            'label_size' : 14,
            'label_angle' : 45
            }
        ),
        # (
        #     SeqFeature(SimpleLocation(20,27)), {
        #     'feature_name' : 'feature4',
        #     'sigil' : 'BOX',
        #     'color' : colors.purple,
        #     'label' : True,
        #     'label_size' : 14,
        #     'label_angle' : 0
        #     }
        # ),
        (
            SeqFeature(SimpleLocation(320, 325)), {
            'feature_name' : 'feature3',
            'sigil' : 'BOX',
            'color' : colors.darkblue,
            'label' : True,
            'label_size' : 14,
            'label_angle' : 45
            }
        ),
    )
    
    seq = 'TTCTCATGTTTGACAGCTTATCATCGATAAGCTTTAATGCGGTAGTTTATCACAGTTAAATTGCTAACGCAGTCAGGCACCGTGTATGAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTCACCCTGGATGCTGTAGGCATAGGCTTGGTTATGCCGGTACTGCCGGGCCTCTTGCGGGATATCGTCCATTCCGACAGCATCGCCAGTCACTATGGCGTGCTGCTAGCGCTATATGCGTTGATGCAATTTCTATGCGCACCCGTTCTCGGAGCACTGTCCGACCGCTTTGGCCGCCGCCCAGTCCTGCTCGCTTCGCTACTTGGAGCCACTATCGACTACGCGATCATGGCGACCACACCCGTCCTGTGGATCCTCTACGCCGGACGCATCGTGGCCGGCATCACCGGCGCCACAGGTGCGGTTGCTGGCGCCTATATCGCCGACATCACCGATGGGGAAGATCGGGCTCGCCACTTCGGGCTCATGAGCGCTTGTTTCGGCGTGGGTATGGTGGCAGGCCCCGTGGCCGGGGGACTGTTGGGCGCCATCTCCTTGCATGCACCATTCCTTGCGGCGGCGGTGCTCAACGGCCTCAACCTACTACTGGGCTGCTTCCTAATGCAGGAGTCGCATAAGGGAGAGCGTCGACCGATGCCCTTGAGAGCCTTCAACCCAGTCAGCTCCTTCCGGTGGGCGCGGGGCATGACTATCGTCGCCGCACTTATGACTGTCTTCTTTATCATGCAACTCGTAGGACAGGTGCCGGCAGCGCTCTGGGTCATTTTCGGCGAGGACCGCTTTCGCTGGAGCGCGACGATGATCGGCCTGTCGCTTGCGGTATTCGGAATCTTGCACGCCCTCGCTCAAGCCTTCGTCACTGGTCCCGCCACCAAACGTTTCGGCGAGAAGCAGGCCATTATCGCCGGCATGGCGGCCGACGCGCTGGGCTACGTCTTGCTGGCGTTCGCGACGCGAGGCTGGATGGCCTTCCCCATTATGATTCTTCTCGCTTCCGGCGGCATCGGGATGCCCGCGTTGCAGGCCATGCTGTCCAGGCAGGTAGATGACGACCATCAGGGACAGCTTCAAGGATCGCTCGCGGCTCTTACCAGCCTAACTTCGATCACTGGACCGCTGATCGTCACGGCGATTTATGCCGCCTCGGCGAGCACATGGAACGGGTTGGCATGGATTGTAGGCGCCGCCCTATACCTTGTCTGCCTCCCCGCGTTGCGTCGCGGTGCATGGAGCCGGGCCACCTCGACCTGAATGGAAGCCGGCGGCACCTCGCTAACGGATTCACCACTCCAAGAATTGGAGCCAATCAATTCTTGCGGAGAACTGTGAATGCGCAAACCAACCCTTGGCAGAACATATCCATCGCGTCCGCCATCTCCAGCAGCCGCACGCGGCGCATCTCGGGCAGCGTTGGGTCCTGGCCACGGGTGCGCATGATCGTGCTCCTGTCGTTGAGGACCCGGCTAGGCTGGCGGGGTTGCCTTACTGGTTAGCAGAATGAATCACCGATACGCGAGCGAACGTGAAGCGACTGCTGCTGCAAAACGTCTGCGACCTGAGCAACAACATGAATGGTCTTCGGTTTCCGTGTTTCGTAAAGTCTGGAAACGCGGAAGTCAGCGCCCTGCACCATTATGTTCCGGATCTGCATCGCAGGATGCTGCTGGCTACCCTGTGGAACACCTACATCTGTATTAACGAAGCGCTGGCATTGACCCTGAGTGATTTTTCTCTGGTCCCGCCGCATCCATACCGCCAGTTGTTTACCCTCACAACGTTCCAGTAACCGGGCATGTTCATCATCAGTAACCCGTATCGTGAGCATCCTCTCTCGTTTCATCGGTATCATTACCCCCATGAACAGAAATCCCCCTTACACGGAGGCATCAGTGACCAAACAGGAAAAAACCGCCCTTAACATGGCCCGCTTTATCAGAAGCCAGACATTAACGCTTCTGGAGAAACTCAACGAGCTGGACGCGGATGAACAGGCAGACATCTGTGAATCGCTTCACGACCACGCTGATGAGCTTTACCGCAGCTGCCTCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCGCAGCCATGACCCAGTCACGTAGCGATAGCGGAGTGTATACTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTGCAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAACACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTCTTCAAGAA'
    seq_id = 'TEST'

    plasmid_drawer = PlasmidDrawer(Seq(seq), seq_id, feature_info)
    plasmid_drawer.draw_gd_diagram('test.pdf', 'circular', {})

    return plasmid_drawer
        
if __name__ == '__main__':
    main()