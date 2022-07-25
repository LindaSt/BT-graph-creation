import fire
import pandas as pd
import os
import json
import numpy as np


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


class SplitJson:
    def __init__(self, excel_path: str, output_folder: str, endpoint_name: str, sheet_name: str = None,
                    json_path: str = None, seed: str = 42, split: float = 0.4, cross_val: int = 1):
        np.random.seed(seed)
        self.output_folder = output_folder
        self.endpoint_name = endpoint_name
        self.json_path = json_path
        self.input_filename = os.path.splitext(os.path.basename(excel_path))[0]
        self.split = split
        self.cross_val = cross_val

        self.endpoints_df = (excel_path, sheet_name)

        self.save_endpoint_jsons()


    @property
    def endpoints_df(self):
        return self._endpoints

    @endpoints_df.setter
    def endpoints_df(self, path_sheet):
        excel_path, sheet_name = path_sheet
        df = pd.read_excel(excel_path, sheet_name=sheet_name, engine='openpyxl')

        df.drop(df[df.Excluded == True].index, inplace=True)
        df.drop(df[df['Exclude for BTS'] == 'x'].index, inplace=True)
        df.drop(df[df['CD8 ID'] == 'na'].index, inplace=True)

        df.set_index('CD8 filename', inplace=True)
        self._endpoints = df

    @property
    def endpoints_dict(self):
        # set up the dictionary
        # load the json file if provided
        if self.json_path:
            with open(self.json_path) as json_file:
                endpoints_dict = json.load(json_file)
        else:
            endpoints_dict = {}
        for filename in self.endpoints_df.index:
            patient_id = self.endpoints_df.at[filename, 'CD8 ID']
            patient_nr = self.endpoints_df.at[filename, 'Patient-Nr']
            d = {self.endpoint_name: int(self.endpoints_df.at[filename, self.endpoint_name]),
                 'CD8 folder': self.endpoints_df.at[filename, 'CD8 folder'],
                 'patient-nr': patient_nr}
            if patient_id in endpoints_dict:
                endpoints_dict[str(patient_id)].update({filename: d})
            else:
                endpoints_dict[str(patient_id)] = {filename: d}

        return endpoints_dict

    @property
    def split_dict(self):
        # get the dataset split json
        # get the data set splits, equal distribution per class and separated by patient

        # Creates the [{patient_id=patient id, class=class}] dict from the data
        patient_id_class_pairs = [{'patient_id': patient_id, 'class': info[self.endpoint_name]} for patient_id, patient_dict in
                 self.endpoints_dict.items() for filename, info in patient_dict.items()]

        # Makes a dataset split (train = 1-split, valid = 1-0.5*split, test = 1-0.5*split)
        # Ensures an even distribution of all classes in all test sets and that each patient is only present in one subset
        classes, patients_per_class = np.unique(np.array([i['class'] for i in patient_id_class_pairs]),
                                                return_counts=True)
        class_patients = {c: [i['patient_id'] for i in patient_id_class_pairs if i['class'] == c] for c in classes}

        # get the indices per class on how to split into train, test and val
        splits = {c: {subset: indices for subset, indices in self.get_indices(len(patients), split=self.split).items()} for
                            c, patients in class_patients.items()}
        splits_per_class = {c: {subset: [patients[i] for i in indices] for subset, indices in self.get_indices(len(patients), split=self.split).items()} for
                            c, patients in class_patients.items()}

        # get the list of the filename and the class for each split
        split_dict = {'train': {}, 'val': {}, 'test': {}}
        for class_id, subset_list in splits_per_class.items():
            for subset, file_list in subset_list.items():
                split_dict[subset][str(class_id)] = file_list
                # d = {f: class_id for f in file_list}
                # split_dict[subset] = {**split_dict[subset], **d}

        return split_dict

    @staticmethod
    def get_indices(nb_elements, split=0.4):
        """
        returns a list of indices
        """
        rand = np.random.permutation(nb_elements)
        assert len(rand) > 3
        test_val_size = max(int(nb_elements * split / 2), 1)  # ensure that at least 1 element is present

        test_ind = rand[-test_val_size:]
        val_ind = rand[-2 * test_val_size:-test_val_size]
        train_ind = rand[:-2 * test_val_size]

        return {'train': train_ind, 'val': val_ind, 'test': test_ind}

    def get_file_class_pairs_per_split(self, data, patients_splits_dict, endpoint):
        """
        makes the filename and class list for each split
        """
        all_pairs = {'train': [], 'val': [], 'test': []}

        for c, patient_id_tup in patients_splits_dict.items():
            data_subsets = [{pat_id: data[pat_id] for pat_id in pat_id_list} for pat_id_list in patient_id_tup]
            for subset_name, subset in zip(['train', 'val', 'test'], data_subsets):
                all_pairs[subset_name] += self.get_file_class_pairs(subset, endpoint)

        return all_pairs

    @staticmethod
    def get_file_class_pairs(patients, endpoint_var):
        """
        Creates the [{file=filename, class=class}] dict from the data
        """
        # pairs = []
        # for patient_id, patient_dict in patients.items():
        #     for filename, info in patient_dict.items():
        #         pairs.append({'file': f'{filename}.gxl', 'class': info[endpoint_var]})
        return {f'{filename}.gxl': info[endpoint_var] for patient_id, patient_dict in patients.items() for filename, info in patient_dict.items()}

    def save_endpoint_jsons(self):
        # save the json
        with open(os.path.join(self.output_folder, f'{self.input_filename}-all.json'), 'w') as fp:
            json.dump(self.endpoints_dict, fp, indent=4, cls=NpEncoder)

        with open(os.path.join(self.output_folder, f'{self.input_filename}-split.json'), 'w') as fp:
            json.dump(self.split_dict, fp, indent=4, cls=NpEncoder)


if __name__ == '__main__':
    """
    INPUT:
    command line arguments:
    --output-path: where the json files should be saved to
    --excel-path: path to the excel with the data
    --sheet-name: (optional) name of the excel sheet that contains the data
    --endpoint-name: name of the column that should be used as an end-point --> comma separated if multiple
    --json-path: (optional) json file that should be extended
    --cross-val: (option) how many cross validation splits should be made (default is 1)
    
    Excel is expected to have the following columns:
    - 'CD8 filename': name of the slide file (e.g. patient1_I_AE1_AE3_CD8)
    - 'CD8 folder: folder where the slide is saved
    - 'CD8 ID': first part of the filename (e.g patient1)
    - 'Patient-Nr': anonymized patient number
    - A column named after the --endpoint argument, e.g.
        'Need resection?'
           0 : No (0 in TNM stage / follow up (>= 2 years), no recurrence)
           1 : True (1 in TNM stage / follow up (>= 2 years) with recurrence)
           
    OUTPUT:
    json file with all files with structure
        Filename-ID:
            CD8 filename:
                folder: CD8 Folder
                patient-nr: Patient-Nr
                endpoint: endpoint value
    
    json file dataset split with structure
        train:
            endpoint value: [Filename-ID, ...]
            ...
        val:
            endpoint value: [Filename-ID, ...]
            ...
        test:
            endpoint value: [Filename-ID, ...]
            ...
    """
    #TODO: add asserts to make sure no data is lost
    fire.Fire(SplitJson)

