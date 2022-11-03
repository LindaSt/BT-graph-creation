import fire
import pandas as pd
import os
import json
import numpy as np
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold


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
    def __init__(self, excel_path: str, output_folder: str, endpoint_name: str = 'Need resection?', sheet_name: str = None,
                    json_path: str = None, seed: str = 42, split: float = 0.4, cross_val: int = 1, multiple_hotspots=None):
        np.random.seed(seed)
        self.output_folder = output_folder
        self.endpoint_name = endpoint_name
        self.json_path = json_path
        self.input_filename = os.path.splitext(os.path.basename(excel_path))[0]
        self.split = split
        self.cross_val = cross_val
        self.seed = seed
        self.multiple_hotspots = multiple_hotspots

        self.endpoints_df = (excel_path, sheet_name)

        self.save_endpoint_jsons()


    @property
    def endpoints_df(self):
        return self._endpoints

    @endpoints_df.setter
    def endpoints_df(self, path_sheet):
        excel_path, sheet_name = path_sheet
        df = pd.read_excel(excel_path, sheet_name=sheet_name, engine='openpyxl')

        df.drop(df[df.Excluded==True].index, inplace=True)
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

        if self.multiple_hotspots is None:
            return endpoints_dict
        else:
            multiple_hs_dict = {}
            for file_id, fileinfo_d in endpoints_dict.items():
                multiple_hs_dict[file_id] = {f'{fname}_hotspot{i}': info for fname, info in fileinfo_d.items() for i in range(self.multiple_hotspots)}
            return multiple_hs_dict

    @property
    def split_dict(self):
        # TODO implement
        return NotImplementedError

    @property
    def split_dict_cv(self):
        # get the dataset split json
        # get the data set splits, equal distribution per class and separated by patient
        # split is either per patient or we only have one entry per patient
        pid = np.array(self.endpoints_df['Patient-Nr'])
        X = np.array(self.endpoints_df.index)
        y = np.array(self.endpoints_df[self.endpoint_name], dtype=int)
        split_dict = {int(i): {} for i in range(self.cross_val)}

        for i, (train_index, test_index) in enumerate(self._get_split_generator(pid, X, y)):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            filename_per_class = {int(c): [] for c in np.unique(y_test)}
            for c, filename in zip(y_test, X_test):
                if self.multiple_hotspots is None:
                    filename_per_class[c].append(filename)
                else:
                    filename_per_class[c] = filename_per_class[c] + [f'{filename}_hotspot{i}' for i in range(self.multiple_hotspots)]
            split_dict[i] = filename_per_class
        return split_dict

    def _get_split_generator(self, pid, X, y):
        if np.unique(pid).size == len(pid):
            skf = StratifiedKFold(n_splits=self.cross_val, random_state=self.seed, shuffle=True)
            split_generator = skf.split(X, y)
        else:
            gkf = StratifiedGroupKFold(n_splits=self.cross_val, random_state=self.seed, shuffle=True)
            split_generator = gkf.split(X, y, groups=pid)
        return split_generator

    def save_endpoint_jsons(self):
        # save the json
        if self.multiple_hotspots is None:
            json_filename = self.input_filename
        else:
            json_filename = f'{self.input_filename}-hotspot-top{self.multiple_hotspots}'
        with open(os.path.join(self.output_folder, f'{json_filename}-all.json'), 'w') as fp:
            json.dump(self.endpoints_dict, fp, indent=4, cls=NpEncoder)

        with open(os.path.join(self.output_folder, f'{json_filename}-cv{self.cross_val}.json'), 'w') as fp:
            if self.cross_val > 1:
                json.dump(self.split_dict_cv, fp, indent=4, cls=NpEncoder)
            else:
                json.dump(self.split_dict, fp, indent=4, cls=NpEncoder)


if __name__ == '__main__':
    """
    INPUT:
    command line arguments:
    --output-path: where the json files should be saved to
    --excel-path: path to the excel with the data
    --sheet-name: (optional) name of the excel sheet that contains the data
    --endpoint-name: name of the column that should be used as an end-point
    --json-path: (optional) json file that should be extended
    --cross-val: (option) how many cross validation splits should be made (default is 1)
    --seed: (option) set a seed (default is 42)
    --multiple_hotspots: (optional) set number, if multiple hotspots are present per slide
    
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
            
    or if with cross validation            
        0:
            endpoint value: [Filename-ID, ...]
            ...
        1:
            endpoint value: [Filename-ID, ...]
            ...
        ...
    """
    #TODO: add asserts to make sure no data is lost
    fire.Fire(SplitJson)

