import os
import xml.etree.ElementTree as ET


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
                    [[float(i.attrib['X']), float(i.attrib['Y'])] for i in element.iter('Coordinate')][0])
            else:
                if group in ['lymphocytes', 'tumorbuds']:
                    group = 'rectangles_' + group
                annotations[group].append(
                    [[float(i.attrib['X']), float(i.attrib['Y'])] for i in element.iter('Coordinate')])

    return annotations

