import glob
import os
import numpy as np
from lxml import etree as ET
import argparse
import re
import shutil

from coord_to_xml import create_asap_xml
from xml_to_txt_file import process_xml_files


def setup_output_folders(output_path):
    # set / create the output
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    # make the output folders
    xml_output = os.path.join(output_path, 'asap_xml')
    if not os.path.isdir(xml_output):
        os.mkdir(xml_output)
    txt_output = os.path.join(output_path, 'coordinates_txt')
    if not os.path.isdir(txt_output):
        os.mkdir(txt_output)

    return xml_output, txt_output


def in_square(square_coordinates, point):
    assert len(square_coordinates) == 4 and len(point) == 2

    x_min, y_min = tuple(square_coordinates[0])
    x_max, y_max = tuple(square_coordinates[2])

    x_dim = x_min <= point[0] <= x_max
    y_dim = y_min <= point[1] <= y_max

    return x_dim and y_dim


def read_hotspot_xmls(hotspot_xmls):
    all_hotspots = {}
    for file_path in hotspot_xmls:
        tree = ET.parse(file_path)
        root = tree.getroot()
        filename = os.path.basename(os.path.splitext(file_path)[0])
        filename = re.sub(' ', '_', filename)

        group = 'hotspot'
        annotations_elements = [i for i in root.iter('Annotation') if i.attrib['PartOfGroup'] == group]

        annotations = [[[float(i.attrib['X']), float(i.attrib['Y'])] for i in element.iter('Coordinate')] for element in
                       annotations_elements]

        all_hotspots[filename] = annotations

    return all_hotspots


def parse_hotspot_xml(hotspot_path, txt_output):
    hotspot_xmls = glob.glob(os.path.join(hotspot_path, r'*.xml'))

    # make txt files
    process_xml_files(hotspot_xmls, txt_output)
    # rename the hotspot txt files to add '_output' into them so all the files have the same name
    hotspot_txts = glob.glob(os.path.join(txt_output, r'*hotspot*'))
    for ht in hotspot_txts:
        if '_output' not in ht:
            regex = re.search('(.*)(_coordinates.*)', ht)
            shutil.move(ht, '{}_output{}'.format(regex.group(1), regex.group(2)))

    # read in the hotspot xml
    hotspots = read_hotspot_xmls(hotspot_xmls)
    return hotspots


def create_asap_xmls(all_txt_files, xml_output):
    # create the asap xml files
    all_txt_file_names = [os.path.join(txt_output, os.path.basename(f)) for f in all_txt_files]
    files_to_process = list(set([re.search(r'(.*)_coordinates', f).group(1) for f in all_txt_file_names]))
    create_asap_xml(files_to_process, xml_output)


def check_output(txt_output, xml_output, all_hotspots):
    # make sure all the hotspots were processed and that there are three txt files for each hotspot
    txt_files = [os.path.basename(f).split('_CD8')[0] for f in glob.glob(txt_output+'/*.txt')]
    txt_files_ids = set(txt_files)
    xml_files_ids = set([os.path.basename(f).split('_CD8')[0] for f in glob.glob(xml_output+'/*.xml')])
    hotspot_ids = set([os.path.basename(f).split('_CD8')[0] for f in all_hotspots])

    # check that there are three files for each hotspot
    text_files_unique = np.unique(txt_files, return_counts=True)
    not_three =', '.join([text_files_unique[0][i] for i, count in enumerate(text_files_unique[1]) if count != 3])
    print(f'Not three output text files for files {not_three}')

    # check for missing files
    missing_txt_files = ', '.join(list(hotspot_ids - txt_files_ids))
    missing_xml_files = ', '.join(list(hotspot_ids - xml_files_ids))

    print(f'# hotspots: {len(hotspot_ids)} | # text files: {len(txt_files_ids)} | # xml files: {len(xml_files_ids)}')
    print(f'missing text files for hotspot(s) {missing_txt_files}')
    print(f'missing xml files for hotspot(s) {missing_xml_files}')


def create_hotspot_only_txt_files(coor_txt_files_path, xml_output, txt_output, all_hotspots):
    coor_txt_files_path = os.path.join(coor_txt_files_path, r'*_coordinates_*.txt')
    all_txt_files = glob.glob(coor_txt_files_path)
    txt_files_to_process = list(set([re.search(r'(.*)_output_coordinates', f).group(1) for f in all_txt_files]))

    output_text_files = []

    for file in txt_files_to_process:
        print('Processing file {}'.format(os.path.basename(file)))
        coord_in_hotspot = {}
        for group in ['lymphocytes', 'tumorbuds']:
            if os.path.basename(file) in all_hotspots.keys():
                hotspots = all_hotspots[os.path.basename(file)]

                # load the file
                file_path = '{}_output_coordinates_{}.txt'.format(file, group)
                output_txt_file = os.path.join(txt_output, os.path.basename(file_path))

                if os.path.isfile(file_path):
                    coordinates = np.loadtxt(file_path)

                    # iterate over hotspot files
                    for h in hotspots:
                        in_hotspot = [in_square(h, i) for i in coordinates]
                        coord_in_hotspot[group] = coordinates[in_hotspot]
                        to_save = coordinates[in_hotspot]
                        output_text_files.append(output_txt_file)
                        if not os.path.isfile(output_txt_file):
                            np.savetxt(output_txt_file, to_save, fmt='%.4f')
                        else:
                            print('The coordinates file {} already exists'.format(output_txt_file))
                else:
                    print('File {} does not exist. Continuing...'.format(file_path))
            else:
                print(f'No hotspots xml for file {os.path.basename(file)}')

    create_asap_xmls(output_text_files, xml_output)

    # make sure all the hotspots were processed and that there are three txt files for each hotspot
    check_output(txt_output, xml_output, all_hotspots.keys())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--xml-hotspots", type=str, required=True)
    parser.add_argument("--output-folder", type=str, required=True)
    parser.add_argument("--coordinate-txt-files", type=str, required=True)
    args = parser.parse_args()

    hotspot_path = args.xml_hotspots

    output_path = args.output_folder
    xml_output, txt_output = setup_output_folders(output_path)

    # get the hotspots (dict)
    hotspots = parse_hotspot_xml(hotspot_path, txt_output)

    coordinate_txt_files = args.coordinate_txt_files

    # create the hotspot asap and txt files
    create_hotspot_only_txt_files(coordinate_txt_files, xml_output, txt_output, hotspots)
