from scipy.spatial import distance
import numpy as np

# util
def _split_nodes_dict(node_dict):
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


def _get_edge_dict(coo_matrix, edge_features, edge_dict=None):
    # make the dictionary
    edge_dict = {} if edge_dict is None else edge_dict
    for i, (edge, dist) in enumerate(zip(coo_matrix, edge_features)):
        edge_dict[i] = {'from_to': tuple(edge), 'distance': dist}
    return edge_dict

