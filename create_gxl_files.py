import glob
import os
import numpy as np
from lxml import etree as ET
import sys
import re
import json
from scipy.spatial import distance
from sklearn.neighbors import NearestNeighbors
import fire


class EdgeConfig:
    """
    This class decodes the edge definition arguments
    edge_def_tb_to_l: "radius-X" or "to-X-nn" where X is a integer --> how should tumor buds and lymphocytes be connected
    edge_def_tb_to_tb: "radius-X" or "to-X-nn" where X is a integer --> how should tumor buds be connected
    fully_connected: specify either 'all', 'tumorbuds' or 'lymphocytes' ('all' supersedes the other --edge-def* arguments)
    """

    def __init__(self, edge_def_tb_to_l=None, edge_def_tb_to_tb=None, fully_connected=None):
        self.edge_def_tb_to_l = self.decode(edge_def_tb_to_l)
        self.edge_def_tb_to_tb = self.decode(edge_def_tb_to_tb)
        self.fully_connected = fully_connected

    @property
    def fully_connected(self):
        return self._fully_connected

    @fully_connected.setter
    def fully_connected(self, fully_connected):
        if fully_connected is not None:
            assert fully_connected in ['all', 'tumorbuds', 'lymphocytes']
        self._fully_connected = fully_connected

    @property
    def edge_definitions(self) -> dict:
        """
        sets up a dictionary with the edge definitions
        {edge_type: [params]}
        """
        edge_def = {}
        if self.fully_connected:
            edge_def['fully_connected'] = self.fully_connected
        # fully connected 'all' supersedes the other edge definitions
        if self.fully_connected != 'all':
            if self.edge_def_tb_to_tb:
                edge_def['tb_to_tb'] = self.edge_def_tb_to_tb
            if self.edge_def_tb_to_l:
                edge_def['tb_to_l'] = self.edge_def_tb_to_l
        return edge_def

    @staticmethod
    def decode(edge_def):
        """
        Decodes the edge definition string from the command line input
        """
        if edge_def:
            if 'radius' in edge_def:
                return ['radius', int(edge_def.split('-')[-1])]
            elif '-nn' in edge_def:
                return ['kNN', int(edge_def.split('-')[1])]
            else:
                print(
                    f'Invalid input. Choose from "radius-X", "to-X-nn" and "fully-connected" (specify number instead of X)')
                sys.exit()
        else:
            return None

    def __str__(self):
        d = {k: ''.join([str(i) for i in v]) for k, v in self.edge_definitions.items()}
        return '-'.join([f'{k}_{v}' for k, v in d.items()])


class Graph:
    """
    Creates a graph object from a list of text files that all need to have the same ID

    file_id: slide_name
    file_path: path to the files (without the ending, e.g. folder/slide_name_output)
    spacing: spacing from ASAP (list, e.g [0.24, 0.24])
    edge_config: EdgeConfig object
    """

    def __init__(self, file_id, file_path, spacing, edge_config=None):
        # print(f'Creating graph for id {file_id}.')
        self.file_path = file_path
        self.file_id = file_id
        self.spacing = spacing

        # set the feature names
        self.node_feature_names = ['type', 'x', 'y']  # TODO: get this automatically?
        self.edge_feature_names = []  # will get added during self.add_edges()

        # get the node dict splits
        self.xy_tb_nodes, self.xy_lymph_nodes = self._get_split_node_dicts()

        # set up the edges
        self.edge_config = edge_config.edge_definitions
        self.edge_dict = {}
        self.add_edges()

    # *********** properties ***********
    @property
    def xy_all_nodes(self) -> dict:
        return {node_id: (node_attrib['x'], node_attrib['y']) for node_id, node_attrib in self.node_dict.items()}

    @property
    def spacing(self) -> float:
        return self._spacing

    @spacing.setter
    def spacing(self, spacing):
        assert spacing[0] == spacing[1]
        self._spacing = spacing[0]

    @property
    def hotspot_coordinates(self) -> list:
        # get the hotspot coordinates
        hotspot_path = f'{self.file_path}_coordinates_hotspot.txt'
        hotspot_coordinates = np.loadtxt(hotspot_path)
        assert len(hotspot_coordinates) == 8
        assert isinstance(hotspot_coordinates[0], float)
        hotspot_coordinates = [tuple([hotspot_coordinates[i], hotspot_coordinates[i + 1]]) for i in range(0, 8, 2)]
        return hotspot_coordinates

    @property
    def node_dict(self) -> dict:
        """
        returns a dict with all the nodes {node_id: [node_attributes]}
        """
        # get all the nodes
        node_dict = {}
        for group in ['lymphocytes', 'tumorbuds']:
            node_name = group[:-1]
            # load the file
            file_path = f'{self.file_path}_coordinates_{group}.txt'

            if os.path.isfile(file_path):
                coordinates = np.loadtxt(file_path)
                # check if file actually contains coordinates
                if len(coordinates) > 0:
                    # Make the annotations and coordinates
                    if len(coordinates.shape) == 1:
                        coordinates = coordinates.reshape((1, 2))

                    assert coordinates.shape[1] == 2

                    # multiply by spacing to get the actual coordinates in mikro-meters
                    for i, line in enumerate(coordinates):
                        # adjust x and y to make relative to the hot-spot coordinates.
                        line -= np.array(self.hotspot_coordinates[0])
                        assert min(line) >= 0
                        node_dict[i] = {'type': node_name, 'x': line[0] * self.spacing, 'y': line[1] * self.spacing}
            else:
                print(f'File {file_path} not found. Continuing...')

        return node_dict

    # *********** adding edges ***********
    def add_edges(self):
        """
        Adds all the edge to self.edge_dict in this format {edge_id: {'feature name 1': feature_value1, ...}, ...}
        The edge ID is the sorted string concatenation of the two node ids.
        The edge id of node 5 and node 0 is therefore '05'
        """
        # get all the edges
        if self.edge_config is not None and len(self.edge_config) > 0:
            for edge_type, param_list in self.edge_config.items():
                # edge_fct can be {'fully_connected', 'tb_to_tb', 'tb_to_l'}
                eval(f'self.{edge_type}')(param_list)

    def update_edge_dict(self, coo_matrix, edge_features, feature_name):
        """
        Updates self.edge_dict based on the coo_matrix and the edge_features list for a specific features (feature_name)

        The edge ID is the sorted string concatenation of the two node ids.
        The edge id of node 5 and node 0 is therefore '0_5'
        """
        # update the dictionary
        for edge, feature in zip(coo_matrix, edge_features):
            edge_id = [str(i) for i in sorted(edge)]
            edge_id_str = '_'.join(edge_id)
            # if the edge already exists, just add the edge features
            if edge_id_str in self.edge_dict.keys():
                assert feature_name not in self.edge_dict[edge_id_str]  # fix this for tb to tb
                self.edge_dict[edge_id_str][feature_name] = feature
            # if the edge does not exist, add it plus the feature
            else:
                self.edge_dict[edge_id_str] = {feature_name: feature}

    # *********** edge insertion functions ***********
    def tb_to_l(self, param_list):
        edge_fct = param_list[0]
        params = param_list[1:]
        # add the edges
        eval(f'self.{edge_fct}')(self.xy_tb_nodes, self.xy_lymph_nodes, params)

    def tb_to_tb(self, param_list):
        edge_fct = param_list[0]
        params = param_list[1:]
        # add the edges
        eval(f'self.{edge_fct}')(self.xy_tb_nodes, self.xy_tb_nodes, params)

    def fully_connected(self, params):
        # params should either be 'all', 'tumorbuds' or 'lymphocytes'
        assert params in ['all', 'tumorbuds', 'lymphocytes']
        node_dict = {'all': self.xy_all_nodes, 'tumorbuds': self.xy_tb_nodes, 'lymphocytes': self.xy_lymph_nodes}
        self.fully_connect(node_dict[params])

    def fully_connect(self, node_dict):
        # calculate the distances
        coo_matrix = []
        edge_features = []

        dict_ids = list(node_dict.keys())
        for i in range(len(node_dict)):
            for j in range(i + 1, len(node_dict)):
                n1 = node_dict[dict_ids[i]]
                n2 = node_dict[dict_ids[j]]
                d = distance.euclidean(n1, n2)
                coo_matrix.append([dict_ids[i], dict_ids[j]])
                edge_features.append(d)
        # update the dictionary
        self.update_edge_dict(coo_matrix, edge_features, feature_name='distance')

    def radius(self, center_dict, perimeter_dict, x):
        assert len(x) == 1
        x = x.pop()
        # calculate the distances
        coo_matrix = []
        edge_features = []

        if len(center_dict) > 0 and len(perimeter_dict) > 0:
            for id_c, xy_c in center_dict.items():
                for id_p, xy_p in perimeter_dict.items():
                    d = distance.euclidean(xy_c, xy_p)
                    # if d < x add edge
                    if d <= x:
                        coo_matrix.append([id_c, id_p])
                        edge_features.append(d)
        # update the dictionary
        self.update_edge_dict(coo_matrix, edge_features, feature_name='distance')

    def kNN(self, center_dict, perimeter_dict, k, distance_metric='euclidean'):
        assert len(k) == 1
        k = orig_k = k.pop()
        # calculate the distances
        coo_matrix = []
        edge_features = []
        # set-up in format for NearestNeighbors
        perimeter_keys = sorted(perimeter_dict.keys())
        center_keys = sorted(center_dict.keys())
        training_set = [perimeter_dict[i] for i in perimeter_keys]
        test_set = [center_dict[i] for i in center_keys]

        # only insert edges if we have elements in the lists
        if len(training_set) > 0 and len(test_set) > 0:
            neigh = NearestNeighbors(k, metric=distance_metric)
            neigh.fit(training_set)
            # if we are compare the same two sets, the first match will always be the point itself --> k += 1
            if training_set == test_set:
                k += 1
            # if #samples > k, set k to number of samples
            if k > len(training_set):
                k = len(training_set)

            distances_list, match_list = neigh.kneighbors(test_set, k, return_distance=True)
            for ind1, (indices, distances) in enumerate(zip(match_list, distances_list)):
                for ind2, d in zip(indices, distances):
                    # ignore self matches and check for duplicates
                    # get the actual node ids
                    node_ind1 = center_keys[ind1]
                    node_ind2 = perimeter_keys[ind2]
                    if node_ind1 != node_ind2 and sorted([node_ind1, node_ind2]) not in coo_matrix:
                        coo_matrix.append(sorted([node_ind1, node_ind2]))
                        edge_features.append(d)

        # update the dictionary
        self.update_edge_dict(coo_matrix, edge_features, feature_name='distance')

    # *********** gxl creation ***********
    def sanity_check(self):
        node_ids = self.node_dict.keys()
        edge_ids = self.edge_dict.keys()

        # make sure all edges point to an existing node
        nodes_in_edges = [i.split('_') for i in edge_ids]
        nodes_in_edges = set([int(val) for sublist in nodes_in_edges for val in sublist])
        assert len(nodes_in_edges - set(node_ids)) == 0

    def get_gxl(self):
        """
        returns the xml-tree for the gxl file
        """
        print(f'Creating gxl tree for {self.file_id}.')
        self.sanity_check()
        type_dict = {'str': 'string', 'int': 'int', 'float': 'float'}

        # initiate the tree
        xml_tree = ET.Element('gxl')

        # add the graph level info
        graph_attrib = {'id': self.file_id, 'edgeids': 'false', 'edgemode': 'undirected'}
        graph_gxl = ET.SubElement(xml_tree, 'graph', graph_attrib)

        # add hot-spot coordinates to gxl
        hotspot_gxl = ET.SubElement(xml_tree, 'hotspot-coordinates')
        for j, tup in enumerate(self.hotspot_coordinates):
            coor_attrib = {'Order': str(j), 'X': str(tup[0]), 'Y': str(tup[1])}
            xml_coordinate = ET.SubElement(hotspot_gxl, 'Coordinate', attrib=coor_attrib)

        # add the nodes
        for node_id, node_attrib in self.node_dict.items():
            node_gxl = ET.SubElement(graph_gxl, 'node', {'id': '_{}'.format(node_id)})
            for attrib_name, attrib_value in node_attrib.items():
                attrib_gxl = ET.SubElement(node_gxl, 'attr', {'name': attrib_name})
                t = re.search(r'(\D*)', type(attrib_value).__name__).group(1)
                attrib_val_gxl = ET.SubElement(attrib_gxl, type_dict[t])
                attrib_val_gxl.text = str(attrib_value)

        # add the edges
        for edge_id, (edge_name, edge_attrib) in enumerate(self.edge_dict.items()):
            from_, to_ = edge_name.split('_')
            edge_gxl = ET.SubElement(graph_gxl, 'edge', {'from': f'_{from_}', 'to': f'_{to_}'})
            for attrib_name, attrib_value in edge_attrib.items():
                attrib_gxl = ET.SubElement(edge_gxl, 'attr', {'name': attrib_name})
                t = re.search(r'(\D*)', type(attrib_value).__name__).group(1)
                attrib_val_gxl = ET.SubElement(attrib_gxl, type_dict[t])
                attrib_val_gxl.text = str(attrib_value)

        e = ET.dump(xml_tree)
        return xml_tree

    # *********** helper functions ***********
    def _get_split_node_dicts(self) -> tuple:
        # splits the nodes dict into the two classes
        all_buds = {}
        all_lymphs = {}
        for node_id, node_attrib in self.node_dict.items():
            type = node_attrib['type']
            x = node_attrib['x']
            y = node_attrib['y']
            if type == 'tumorbud':
                all_buds[node_id] = (x, y)
            if type == 'lymphocyte':
                all_lymphs[node_id] = (x, y)

        return all_buds, all_lymphs


class GxlFilesCreator:
    """
    Creates the xml trees from the text files with the coordinates
    """

    def __init__(self, files_to_process, spacings, edge_config, normalize=False):
        """
        files_to_process: list of paths to the files that should be processed
        spacings: dictionary that contains the spacing for each WSI (read from the spacing.json)
        edge_config: EdgeConfig object
        """
        self.files_to_process = files_to_process
        self.edge_config = edge_config
        self.spacings = spacings
        self.normalize = normalize

    @property
    def graphs(self) -> list:
        files_dict = {os.path.basename(f)[:-7]: f for f in self.files_to_process}  # get rid of '_output' at the end
        return [Graph(file_id, file_paths, self.spacings[file_id], self.edge_config) for file_id, file_paths in
                files_dict.items()]

    @property
    def gxl_trees(self) -> dict:
        """
        creates dictionary {file_id: xml-tree}
        """
        return {graph.file_id: graph.get_gxl() for graph in self.graphs}

    def save(self, output_folder):
        # create output folder if it does not exist
        output_path = os.path.join(output_folder, str(self.edge_config))
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        # save the xml trees
        print(f'Saving gxl files to {output_path}')
        for file_id, tree in self.gxl_trees.items():
            ET.ElementTree(tree).write(os.path.join(output_path, file_id + '.gxl'), pretty_print=True)


def make_gxl_dataset(coord_txt_files_folder, spacing_json, output_folder, edge_def_tb_to_l=None, edge_def_tb_to_tb=None,
                     fully_connected=None):
    """
    INPUT
    --coord-txt-files-folder: path to the folder with the coordinates text files
    --spacing-json: Path to json file that contains the spacing for each whole slide image. It is needed to compute the distance between elements.
    --edge-def-tb-to-l (optional):
      - radius-x: connect elements in radius X (in mikrometer)
      - to-X-nn: connect to k closest elements where X is the number of neighbours
    --edge-def-tb-to-tb (optional): same options as edge-def-tb-to-l
    --fully-connected: (optional) specify 'all', 'tumorbuds' or 'lymphocytes' ('all' supersedes the other --edge-def... arguments)
    --output-folder: path to where output folder should be created

    OUTPUT
    One gxl file per hotspot, which contains the graph (same structure as the gxl files from the IAM Graph Databse)
    The distances are centered (xy - xy(average)).
    """
    # get the edge definitions
    edge_def_config = EdgeConfig(edge_def_tb_to_l, edge_def_tb_to_tb, fully_connected)

    # read the spacing json
    spacing_json = r'{}'.format(spacing_json)
    with open(spacing_json) as data_file:
        spacings = json.load(data_file)

    # get a list of all the txt files to process
    input_path = os.path.join(coord_txt_files_folder, r'*_coordinates_*.txt')
    all_files = glob.glob(input_path)
    files_to_process = list(set([re.search(r'(.*)_coordinates', f).group(1) for f in all_files]))

    # Create the gxl files
    gxl_files = GxlFilesCreator(files_to_process, spacings, edge_def_config)

    # save the gxl files
    gxl_files.save(output_folder)


if __name__ == '__main__':
    fire.Fire(make_gxl_dataset)
