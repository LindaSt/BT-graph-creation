
import glob
import os
import numpy as np
from lxml import etree as ET
#import xml.etree.ElementTree as ET
import argparse
import re
import json
from scipy.spatial import distance

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

def split_nodes_dict(node_dict):
    # split nodes into the two classes
    all_buds = {}
    all_lymphs = {}
    for node_id, node_attrib in node_dict.items():
        type = node_attrib['type']
        x = node_attrib['X']
        y = node_attrib['Y']
        if type == 'tumorbud':
            all_buds[node_id] = (x, y)
        if type == 'lymphocyte':
            all_lymphs[node_id] = (x, y)

    return all_buds, all_lymphs


def get_edge_dict(coo_matrix, edge_features, edge_dict=None):
    # make the dictionary
    edge_dict = {} if edge_dict is None else edge_dict
    for i, (edge, dist) in enumerate(zip(coo_matrix, edge_features)):
        edge_dict[i] = {'from_to': tuple(edge), 'distance': dist}
    return edge_dict


def connect_tumorbuds(coo_matrix, edge_features, all_buds):
    all_bud_ids = list(all_buds.keys())
    if len(all_bud_ids) > 0:
        for i in range(len(all_buds)):
            for j in range(i+1, len(all_buds)):
                tb1 = all_buds[all_bud_ids[i]]
                tb2 = all_buds[all_bud_ids[j]]
                d = distance.euclidean(tb1, tb2)
                coo_matrix.append([all_bud_ids[i], all_bud_ids[j]])
                edge_features.append(d)
    return coo_matrix, edge_features


def radius_x(node_dict, connect_tbs, x):
    # split nodes into the two classes
    all_buds, all_lymphs = split_nodes_dict(node_dict)

    # calculate the distances
    coo_matrix = []
    edge_features = []
    # only insert edges if we have buds and lymphocytes
    if len(all_buds) > 0 and len(all_lymphs) > 0:
        for l_id, l_xy in all_lymphs.items():
            distances = np.empty((0, 2))
            for tb_id, tb_xy in all_buds.items():
                d = distance.euclidean(l_xy, tb_xy)
                # if d < x add edge
                if d <= x:
                    coo_matrix.append([l_id, tb_id])
                    edge_features.append(d)

        # connect the tbs if necessary
        if connect_tbs:
            coo_matrix, edge_features = connect_tumorbuds(coo_matrix, edge_features, all_buds)
    # get the the dictionary
    edge_dict = get_edge_dict(coo_matrix, edge_features)

    return edge_dict

def to_closest_tb(node_dict, connect_tbs):
    # split nodes into the two classes
    all_buds, all_lymphs = split_nodes_dict(node_dict)

    # calculate the distances
    coo_matrix = []
    edge_features = []
    # only insert edges if we have buds and lymphocytes
    if len(all_buds) > 0 and len(all_lymphs) > 0:
        for l_id, l_xy in all_lymphs.items():
            distances = np.empty((0, 2))
            for tb_id, tb_xy in all_buds.items():
                d = distance.euclidean(l_xy, tb_xy)
                distances = np.append(distances, np.array([tb_id, d]).reshape((1, 2)), axis=0)

            # find the closest bud
            ind = np.argmin(distances[:, 1])
            coo_matrix.append([l_id, int(distances[ind, 0])])
            edge_features.append(distances[ind, 1])

        # connect the tbs if necessary
        if connect_tbs:
            coo_matrix, edge_features = connect_tumorbuds(coo_matrix, edge_features, all_buds)

    # get the the dictionary
    edge_dict = get_edge_dict(coo_matrix, edge_features)

    return edge_dict


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

