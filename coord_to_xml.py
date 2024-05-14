import glob
import os
import numpy as np
from lxml import etree as ET
import re
import fire

from util.xml_parsing import create_xml_tree

class XmlFile:
    def __init__(self, file: str, output_path: str, full: bool = False):
        self.output_base_path = output_path
        self.file = file
        self.colors = {'tumorbuds': '#73d216', 'lymphocytes': '#ffaa00', 'hotspot': '#3465a4'}
        self.full_coord = full
        self.xmls = self.create_xml_trees()

    def read_txt_file(self, group, id=''):
        file_path = f'{self.file}{id}_coordinates_{group}.txt'
        if os.path.isfile(file_path):
            # load the file
            coordinates = np.loadtxt(file_path)
            return coordinates
        else:
            print(f'File {file_path} does not exist.')
            return None

    @property
    def data(self):
        data = []
        if self.full_coord:
            data = [{
                'lymphocytes': self.read_txt_file('lymphocytes'),
                'tumorbuds': self.read_txt_file('tumorbuds'),
                'output_file': os.path.join(self.output_base_path,
                                            '{}_full_asap.xml'.format(os.path.basename(self.file)))
            }]
        else:
            hotspots = self.read_txt_file('hotspot')
            if hotspots is None:
                return data
            elif len(hotspots.shape) == 1:  # we have one hotspot
                data = [{
                    'hotspot': hotspots,
                    'lymphocytes': self.read_txt_file('lymphocytes'),
                    'tumorbuds': self.read_txt_file('tumorbuds'),
                    'output_file': os.path.join(self.output_base_path,
                                                f'{os.path.basename(self.file)}_asap.xml')
                }]
            elif len(hotspots.shape) > 1:  # we have more than one hotspot present
                data = [{
                    'hotspot': h,
                    'lymphocytes': self.read_txt_file('lymphocytes', id=f'_hotspot{i}'),
                    'tumorbuds': self.read_txt_file('tumorbuds', id=f'_hotspot{i}'),
                    'output_file': os.path.join(self.output_base_path,
                                                f'{os.path.basename(self.file)}_hotspot{i}_asap.xml')
                } for i, h in enumerate(hotspots)]
        return data

    def create_xml_trees(self):
        xmls = {}
        for coord_dict in self.data:
            out_file = coord_dict.pop('output_file')
            xmls[out_file] = create_xml_tree(coord_dict, self.colors)
        return xmls

    def save_xml(self):
        if self.xmls is None:
            return
        for output_file, xml_string in self.xmls.items():
            with open(output_file, 'w') as f:
                f.write(xml_string)
            # e = ET.dump(xml_tree)

def create_asap_xml(coordinates_txt_files_folder: str, output_folder: str, separator_id: str = '_coord',
                    full: bool = False):
    # create output folder if it does not exist
    output_base_path = output_folder
    if not os.path.isdir(output_base_path):
        os.makedirs(output_base_path)

    # get a list of all the txt files to process
    input_path = os.path.join(coordinates_txt_files_folder, r'*_coordinates_*.txt')
    all_files = glob.glob(input_path)
    files_to_process = list(set([re.search(repr(f'(.*){separator_id}')[1:-1], f).group(1) for f in all_files]))

    for file in files_to_process:
        XmlFile(file=file, output_path=output_base_path, full=full).save_xml()


if __name__ == '__main__':
    fire.Fire(create_asap_xml)
