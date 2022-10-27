import glob
import os
import numpy as np
from lxml import etree as ET
import re
import fire


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
            xml_tree = ET.Element('ASAP_Annotations')
            # Make the annotations and coordinates
            xml_annotations = ET.SubElement(xml_tree, 'Annotations')
            # Make the groups
            xml_annotation_groups = ET.SubElement(xml_tree, 'AnnotationGroups')

            for group, coordinates in coord_dict.items():
                # check if file is not empty
                if coordinates is not None and len(coordinates) > 0:
                    # make the group
                    xml_group = ET.SubElement(xml_annotation_groups, 'Group',
                                              attrib={'Name': group, 'PartOfGroup': 'None',
                                                      'Color': self.colors[group]})
                    xml_group_attrib = ET.SubElement(xml_group, 'Attributes')

                    # Make the annotations and coordinates
                    if len(coordinates.shape) == 1:
                        coordinates = np.reshape(coordinates, (1, coordinates.shape[0]))
                    if coordinates.shape[1] == 2:
                        annotation_type = 'Dot'
                    else:
                        annotation_type = 'Rectangle'
                    # iterate over the list
                    for i, line in enumerate(coordinates):
                        attrib = {'Name': 'Annotation {}'.format(i), 'PartOfGroup': group,
                                  'Color': self.colors[group],
                                  'Type': annotation_type}
                        if annotation_type == 'Dot':
                            xml_annotation = ET.SubElement(xml_annotations, 'Annotation', attrib)
                            xml_coordinates = ET.SubElement(xml_annotation, 'Coordinates')
                            coor_attrib = {'Order': '0', 'X': str(line[0]), 'Y': str(line[1])}
                            xml_coordinate = ET.SubElement(xml_coordinates, 'Coordinate', attrib=coor_attrib)
                        else:
                            xml_annotation = ET.SubElement(xml_annotations, 'Annotation', attrib)
                            xml_coordinates = ET.SubElement(xml_annotation, 'Coordinates')
                            coor = [(line[i], line[i + 1]) for i in range(0, len(line), 2)]
                            for j, tup in enumerate(coor):
                                coor_attrib = {'Order': str(j), 'X': str(tup[0]), 'Y': str(tup[1])}
                                xml_coordinate = ET.SubElement(xml_coordinates, 'Coordinate', attrib=coor_attrib)
            xmls[out_file] = ET.ElementTree(xml_tree)
        return xmls

    def save_xml(self):
        if self.xmls is None:
            return
        for output_file, xml_tree in self.xmls.items():
            # e = ET.dump(xml_tree)
            e = xml_tree.write(output_file, pretty_print=True)


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
