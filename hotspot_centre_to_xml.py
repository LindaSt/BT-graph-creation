#TODO: add this to the readme

import fire
import pandas as pd
import os
import glob
import math
from lxml import etree as ET


class HotspotsCsvToXml:
    def __init__(self, hotspot_csvs_path: str, output_path: str, spacing: int = 0.242797397769517,
                 hotspot_size: int = 0.785, top: int = 10):
        self.spacing = spacing  # mikrometer per pixel
        self.hotspot_size = hotspot_size  # in mm2
        self.top = top  # top X hotspots extracted (sorted top-down by # buds). If set to 0, all hotspots are extracted.
        self.csv_files = hotspot_csvs_path
        self.output_path = output_path

    @property
    def csv_files(self):
        return self._csv_files

    @csv_files.setter
    def csv_files(self, path):
        if os.path.isdir(path):
            files = glob.glob(os.path.join(path, '*.csv'))
            if len(files) > 0:
                self._csv_files = files
            else:
                print(f"Folder {path} does not contain csv files.")
                exit(-1)
        else:
            print(f"Folder {path} does not exist.")
            exit(-1)

    @property
    def output_path(self):
        return self._output_path

    @output_path.setter
    def output_path(self, output_path):
        # make the output folder if it does not exist
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        self._output_path = output_path

    def process_files(self):
        # size of the hotspot square
        pixel_shift = math.sqrt(self.hotspot_size)*1000/self.spacing/2
        # process the files
        for csv_file in self.csv_files:
            file_id = os.path.basename(csv_file).split('.csv')[0]
            output_file = os.path.join(self.output_path, f'{file_id}_autom_hotspots.xml')
            hotspots = pd.read_csv(csv_file, index_col='hotspot_id')[['x', 'y']]
            hotspots = hotspots.head(self.top)
            xml_tree = self.coord_list_to_xml(hotspots.to_dict(orient='index'), pixel_shift)
            e = ET.ElementTree(xml_tree).write(output_file, pretty_print=True)

    @staticmethod
    def coord_list_to_xml(id_coord, pixel_shift):
        # initiate the tree
        xml_tree = ET.Element('ASAP_Annotations')
        # Make the annotations and coordinates
        xml_annotations = ET.SubElement(xml_tree, 'Annotations')
        # Make the groups
        xml_annotation_groups = ET.SubElement(xml_tree, 'AnnotationGroups')
        # make the group
        xml_group = ET.SubElement(xml_annotation_groups, 'Group',
                                  attrib={'Name': 'hotspot', 'PartOfGroup': 'None', 'Color': '#73d216'})
        xml_group_attrib = ET.SubElement(xml_group, 'Attributes')

        # Make the annotations and coordinates
        for annotation_id, centre_coord in id_coord.items():
            attrib = {'Name': 'Annotation {}'.format(annotation_id.split('_')[-1]), 'PartOfGroup': 'hotspot',
                      'Color': '#73d216', 'Type': 'Rectangle'}
            xml_annotation = ET.SubElement(xml_annotations, 'Annotation', attrib)
            xml_coordinates = ET.SubElement(xml_annotation, 'Coordinates')
            x, y = centre_coord['x']*4, centre_coord['y']*4  # because of level
            coor = [[x - pixel_shift, y + pixel_shift], [x + pixel_shift, y + pixel_shift],
                    [x + pixel_shift, y - pixel_shift], [x - pixel_shift, y - pixel_shift]]
            for j, tup in enumerate(coor):
                coor_attrib = {'Order': str(j), 'X': str(tup[0]), 'Y': str(tup[1])}
                xml_coordinate = ET.SubElement(xml_coordinates, 'Coordinate', attrib=coor_attrib)

        # e = ET.dump(xml_tree)
        return xml_tree


if __name__ == '__main__':
    fire.Fire(HotspotsCsvToXml).process_files()