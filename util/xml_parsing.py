import os
import xml.etree.ElementTree as ET
from xml.dom import minidom
import numpy as np


def parse_xml(file_path) -> dict:
    tree = ET.parse(file_path)
    root = tree.getroot()

    groups_colours = {i.attrib['Name']: i.attrib['Color'] for i in root.iter('Group')}
    groups = ['hotspot', 'lymphocytes', 'tumorbuds', 'lymphocytesR', 'tumorbudsR']
    annotations_elements = {g: [] for g in groups}

    for i in root.iter('Annotation'):
        annotations_elements[i.attrib['PartOfGroup']].append(i)

    annotations = {g: [] for g in groups}
    for group, element_list in annotations_elements.items():
        for element in element_list:
            if element.attrib['Type'] == 'Dot':
                annotations[group].append(
                    np.array([[float(i.attrib['X']), float(i.attrib['Y'])] for i in element.iter('Coordinate')][0]))
            else:
                if group in ['lymphocytes', 'tumorbuds']:
                    group = 'rectangles_' + group
                annotations[group].append(
                   np.array([[float(i.attrib['X']), float(i.attrib['Y'])] for i in element.iter('Coordinate')]))

    return annotations


def create_xml_tree(coord_dict, colors: dict = None):
    if colors is None:
        colors = {'hotspot': '#73d216', 'lymphocytes': '#3465a4', 'tumorbuds': '#cc0000',
                  'lymphocytesR': '#3465a4', 'tumorbudsR': '#cc0000'}

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
                                              'Color': colors[group]})
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
                          'Color': colors[group],
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

    xml_string = ET.tostring(xml_tree, encoding='unicode', method='xml')
    xml_pretty_string = minidom.parseString(xml_string).toprettyxml(indent="  ")
    return xml_pretty_string

