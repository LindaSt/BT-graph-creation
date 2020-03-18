# EDGE DEFINITION FUNCTIONS
from scipy.spatial import distance
import numpy as np


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
    all_buds, all_lymphs = _split_nodes_dict(node_dict)

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
    edge_dict = _get_edge_dict(coo_matrix, edge_features)

    return edge_dict


def to_closest_tb(node_dict, connect_tbs):
    # split nodes into the two classes
    all_buds, all_lymphs = _split_nodes_dict(node_dict)

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
    edge_dict = _get_edge_dict(coo_matrix, edge_features)

    return edge_dict


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

