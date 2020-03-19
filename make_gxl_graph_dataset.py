
import glob
import os
import numpy as np
from lxml import etree as ET
#import xml.etree.ElementTree as ET
import argparse
import re
import json
# TODO: test and make more generic (radius_X, to_all, to_knn)
from util.edge_insertion_functions import to_closest_tb, connect_tumorbuds, radius_x

#%%

def make_gxl_tree(filename, node_dict, edge_dict):
    type_dict = {'str': 'string', 'int': 'int', 'float': 'float'}

    # initiate the tree
    xml_tree = ET.Element('gxl')

    # add the graph level info
    graph_attrib = {'id': filename, 'edgeids': 'false', 'edgemode': 'undirected'}
    graph_gxl = ET.SubElement(xml_tree, 'graph', graph_attrib)

    # add the nodes
    for node_id, node_attrib in node_dict.items():
        node_gxl = ET.SubElement(graph_gxl, 'node', {'id': '_{}'.format(node_id)})
        for attrib_name, attrib_value in node_attrib.items():
            attrib_gxl = ET.SubElement(node_gxl, 'attr', {'name': attrib_name})
            t = re.search(r'(\D*)', type(attrib_value).__name__).group(1)
            attrib_val_gxl = ET.SubElement(attrib_gxl, type_dict[t])
            attrib_val_gxl.text = str(attrib_value)

    # add the edges
    for edge_id, edge_attrib in edge_dict.items():
        from_, to_ = edge_attrib.pop('from_to')
        edge_gxl = ET.SubElement(graph_gxl, 'edge', {'from': '_{}'.format(from_), 'to': '_{}'.format(to_)})
        for attrib_name, attrib_value in edge_attrib.items():
            attrib_gxl = ET.SubElement(edge_gxl, 'attr', {'name': attrib_name})
            t = re.search(r'(\D*)', type(attrib_value).__name__).group(1)
            attrib_val_gxl = ET.SubElement(attrib_gxl, type_dict[t])
            attrib_val_gxl.text = str(attrib_value)

    # TODO: remove this savely
    e = ET.dump(xml_tree)
    return xml_tree

# node_dict = {0: {'type': 'tumorbud', 'X': 1.1, 'Y': 2.2}, 1: {'type': 'tumorbud', 'X': 6.1, 'Y': 4.2},
#              3: {'type': 'lymphocyte', 'X': 8.0, 'Y': 5.0}, 2: {'type': 'lymphocyte', 'X': 1.7, 'Y': 2.5}}
#
# edge_dict = {0: {'from_to': (0, 1), 'distance': 12.36}, 1: {'from_to': (1, 2), 'distance': 1.36}, 2: {'from_to': (0, 2), 'distance': 58.36}}

# tree = make_gxl_tree('bla', node_dict, edge_dict)
# e = ET.ElementTree(tree).write(r"test.gxl", pretty_print=True)

def create_gxl(files_to_process, output_folder, spacings, edge_definition, radius=None, connect_tbs=False):
    """
    This function takes a coordinate txt file as input and creates an xml file that can be loaded in ASAP
    """
    for file in files_to_process:
        filename = re.search(r'(.*)_output', os.path.basename(file)).group(1)
        print('Processing file {}...'.format(filename))
        output_file = os.path.join(output_folder, '{}.gxl'.format(os.path.basename(file)))
        assert spacings[filename][0] == spacings[filename][1]
        spacing = spacings[filename][0]

        # get all the nodes
        node_dict = {}
        for group in ['lymphocytes', 'tumorbuds']:
            node_name = group[:-1]
            # load the file
            file_path = '{}_coordinates_{}.txt'.format(file, group)

            if os.path.isfile(file_path):
                coordinates = np.loadtxt(file_path)
                # check if file actually contains coordinates
                if len(coordinates) > 0:
                    # Make the annotations and coordinates
                    if len(coordinates.shape) == 1:
                        coordinates = coordinates.reshape((1, 2))

                    assert coordinates.shape[1] == 2

                    # iterate over the list
                    # multiply by spacing to get the actual coordinates in mikro-meters
                    for i, line in enumerate(coordinates):
                        node_dict[i] = {'type': node_name, 'X': line[0]*spacing, 'Y': line[1]*spacing}

        # get all the edges
        if edge_definition is not None:
            if edge_definition == 'radius_x':
                edge_dict = eval(edge_definition)(node_dict, connect_tbs, radius)
            else:
                edge_dict = eval(edge_definition)(node_dict, connect_tbs)

        tree = make_gxl_tree(filename, node_dict, edge_dict)
        e = ET.ElementTree(tree).write(os.path.join(output_path, filename+'.gxl'), pretty_print=True)


def get_appendix(edge_definition, radius, connect_tbs):
    output_appendix = re.sub('x', str(radius), 'radius_x') if radius else edge_definition
    if connect_tbs:
        output_appendix = '{}_tb_connected'.format(output_appendix)
    return output_appendix


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates-txt-files-folder", type=str, required=True)
    parser.add_argument("--output-folder", type=str, required=True)
    parser.add_argument("--spacing-json", type=str, required=True)
    parser.add_argument("--edge-definition", type=str, required=False, default=None)
    parser.add_argument("--radius", type=str, required=False, default='50')
    parser.add_argument("--connect-tbs", required=False, action='store_true')

    args = parser.parse_args()

    edge_definition = args.edge_definition
    radius = float(args.radius) if 'radius_' in edge_definition else None

    connect_tbs = args.connect_tbs

    # get a list of all the txt files to process
    input_path = os.path.join(args.coordinates_txt_files_folder, r'*_coordinates_*.txt')
    all_files = glob.glob(input_path)
    files_to_process = list(set([re.search(r'(.*)_coordinates', f).group(1) for f in all_files]))

    # read the spacing json
    spacing_json = r'{}'.format(args.spacing_json)
    with open(spacing_json) as data_file:
        spacings = json.load(data_file)

    # create output folder if it does not exist
    output_path = '{}_{}'.format(args.output_folder, get_appendix(edge_definition, radius, connect_tbs))
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # create the xml files
    create_gxl(files_to_process, output_path, spacings, edge_definition, radius, connect_tbs)

