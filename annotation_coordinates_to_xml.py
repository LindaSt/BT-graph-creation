import glob
import os
import numpy as np
from lxml import etree as ET
#import xml.etree.ElementTree as ET
import argparse
import re

#%%
def create_asap_xml(files_to_process, output_folder):
    """
    This function takes a coordinate txt file as input and creates an xml file that can be loaded in ASAP
    """
    colours = {'tumorbuds': '#73d216', 'lymphocytes': '#ffaa00', 'hotspot': '#3465a4'}

    for file in files_to_process:
        output_file = os.path.join(output_folder, '{}_asap.xml'.format(os.path.basename(file)))
        # initiate the tree
        xml_tree = ET.Element('ASAP_Annotations')
        # Make the annotations and coordinates
        xml_annotations = ET.SubElement(xml_tree, 'Annotations')
        # Make the groups
        xml_annotation_groups = ET.SubElement(xml_tree, 'AnnotationGroups')

        for group in ['lymphocytes', 'tumorbuds', 'hotspot']:
            # load the file
            file_path = '{}_coordinates_{}.txt'.format(file, group)

            if os.path.isfile(file_path):
                coordinates = np.loadtxt(file_path)
                # check if file is not empty
                if len(coordinates) > 0:
                    # make the group
                    xml_group = ET.SubElement(xml_annotation_groups, 'Group',
                                              attrib={'Name': group, 'PartOfGroup': 'None', 'Color': colours[group]})
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
                        attrib = {'Name': 'Annotation {}'.format(i), 'PartOfGroup': group, 'Color': colours[group], 'Type': annotation_type}
                        if annotation_type == 'Dot':
                            xml_annotation = ET.SubElement(xml_annotations, 'Annotation', attrib)
                            xml_coordinates = ET.SubElement(xml_annotation, 'Coordinates')
                            coor_attrib = {'Order': '0', 'X': str(line[0]), 'Y': str(line[1])}
                            xml_coordinate = ET.SubElement(xml_coordinates, 'Coordinate', attrib=coor_attrib)
                        else:
                            xml_annotation = ET.SubElement(xml_annotations, 'Annotation', attrib)
                            xml_coordinates = ET.SubElement(xml_annotation, 'Coordinates')
                            coor = [(line[i], line[i+1]) for i in range(0, len(line), 2)]
                            for j, tup in enumerate(coor):
                                coor_attrib = {'Order': str(j), 'X': str(tup[0]), 'Y': str(tup[1])}
                                xml_coordinate = ET.SubElement(xml_coordinates, 'Coordinate', attrib=coor_attrib)

            else:
                print('File {} does not exist. Continuing...'.format(file_path))
        e = ET.dump(xml_tree)
        e = ET.ElementTree(xml_tree).write(output_file, pretty_print=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates-txt-files-folder", type=str, required=True)
    parser.add_argument("--output-folder", type=str, required=True)
    args = parser.parse_args()

    # get a list of all the txt files to process
    input_path = os.path.join(args.coordinates_txt_files_folder, r'*_coordinates_*.txt')
    all_files = glob.glob(input_path)
    files_to_process = list(set([re.search(r'(.*)_coordinates', f).group(1) for f in all_files]))

    # create output folder if it does not exist
    output_base_path = args.output_folder
    if not os.path.isdir(output_base_path):
        os.mkdir(output_base_path)

    # create the xml files
    create_asap_xml(files_to_process, output_base_path)