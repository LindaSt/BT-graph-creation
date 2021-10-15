import fire
from lxml import etree as ET
import os
import numpy as np
import json


def get_indices(nb_elements, split=0.4):
    """
    returns a list of indices
    """
    rand = np.random.permutation(nb_elements)
    assert len(rand) > 3
    test_val_size = max(int(nb_elements*split/2), 1) # ensure that at least 1 element is present

    test_ind = rand[-test_val_size:]
    val_ind = rand[-2*test_val_size:-test_val_size]
    train_ind = rand[:-2*test_val_size]

    return train_ind, val_ind, test_ind


def get_file_class_pairs(patients, endpoint_var):
    """
    Creates the [{file=filename, class=class}] dict from the data
    """
    pairs = []
    for patient_id, patient_dict in patients.items():
        for filename, info in patient_dict.items():
            pairs.append({'file': f'{filename}.gxl', 'class': info[endpoint_var]})
    return pairs


def get_patient_id_class_pairs(patients, endpoint_var):
    """
    Creates the [{patient_id=patient id, class=class}] dict from the data
    """
    pairs = []
    for patient_id, patient_dict in patients.items():
        for filename, info in patient_dict.items():
            pairs.append({'patient_id': patient_id, 'class': info[endpoint_var]})
    return pairs


def get_splits_dict(patient_id_class_pair, split=0.4):
    """
    Makes a dataset split (train = 1-split, valid = 1-0.5*split, test = 1-0.5*split)
    Ensures an even distribution of all classes in all test sets and that each patient is only
    present in one subset

    """
    classes, patients_per_class = np.unique(np.array([i['class'] for i in patient_id_class_pair]), return_counts=True)
    class_patients = {c: [i['patient_id'] for i in patient_id_class_pair if i['class'] == c] for c in classes}

    # get the indices per class on how to split into train, test and val
    splits_per_class = {c: [[patients[i] for i in subset] for subset in get_indices(len(patients), split=split)] for c, patients in class_patients.items()}

    return splits_per_class


def get_file_class_pairs_per_split(data, patients_splits_dict, endpoint):
    """
    makes the filename and class list for each split
    """
    all_pairs = {'train': [], 'val': [], 'test': []}

    for c, patient_id_tup in patients_splits_dict.items():
        data_subsets = [{pat_id: data[pat_id] for pat_id in pat_id_list} for pat_id_list in patient_id_tup]
        for subset_name, subset in zip(['train', 'val', 'test'], data_subsets):
            all_pairs[subset_name] += get_file_class_pairs(subset, endpoint)

    return all_pairs['train'], all_pairs['val'], all_pairs['test']


def get_xml_tree(file_class_list, dataset_name, endpoint):
    """
    input: output of get_file_class_pairs_per_split() --> list of dict [{'file': filename, 'class': class}]
    """
    # initiate the tree
    xml_tree = ET.Element('GraphCollection')

    # add the dataset info
    cxl = ET.SubElement(xml_tree, dataset_name, {'counts': str(len(file_class_list)), 'endpoint': endpoint})

    # add the elements
    for dic in file_class_list:
        element_cxl = ET.SubElement(cxl, 'print', {'file': dic['file'], 'class': dic['class']})

    # ET.dump(xml_tree)
    return xml_tree


def json_to_cxl(json_path, output_folder, endpoint, dataset_name, split=0.4, seed=42):
    """
    Main function, input documentation see __main__
    """
    np.random.seed(seed)

    # load the json file
    with open(json_path) as json_file:
        data = json.load(json_file)

    # get the data set splits, equal distribution per class and separated by patient
    patient_id_class_pairs = get_patient_id_class_pairs(data, endpoint)
    patients_splits_dict = get_splits_dict(patient_id_class_pairs, split=split)

    # get the list of the filename and the class for each split
    train_split, val_split, test_split = get_file_class_pairs_per_split(data, patients_splits_dict, endpoint)

    # get the xml trees
    train_cxl = get_xml_tree(train_split, dataset_name, endpoint)
    valid_cxl = get_xml_tree(val_split, dataset_name, endpoint)
    test_cxl = get_xml_tree(test_split, dataset_name, endpoint)

    # save the xml trees
    for filename, tree in zip(['train.cxl', 'val.cxl', 'test.cxl'], [train_cxl, valid_cxl, test_cxl]):
        ET.ElementTree(tree).write(os.path.join(output_folder, filename), pretty_print=True)


if __name__ == '__main__':
    """
    INPUT:
    command line arguments:
    --output-path: where the cxl files should be saved to
    --json-path: path to the json file with the endpoint(s) data
    --split: (optional, default 0.4) how much should be split off for train and test (will be split in half for test valid) 
    --seed: (optional, default 42) set seed for split
    --endpoint: name of the variable (in the dict) that encodes the end-point
    --dataset-name: name that should be given to the dataset in the xml files
        
    json file has the structure
    patient-id:
        filename:
            folder: folder name
            endpoint: endpoint value
    """
    fire.Fire(json_to_cxl)
