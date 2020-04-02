import fire
import pandas as pd
import os
import json


def make_json(excel_path, output_folder, endpoint_encodings, sheet_name=None, json_path=None):
    endpoint_encodings = endpoint_encodings.split(',')
    input_filename = os.path.splitext(os.path.basename(excel_path))[0]
    # read the excel as pandas frame
    endpoints = pd.read_excel(excel_path, sheet_name=sheet_name)
    endpoints.set_index('filename', inplace=True)

    # set up the dictionary
    # load the json file if provided
    if json_path:
        with open(json_path) as json_file:
            endpoints_dict = json.load(json_file)
    else:
        endpoints_dict = {}
    for filename in endpoints.index:
        patient_id = endpoints.at[filename, 'Patient-ID']
        folder = endpoints.at[filename, 'folder']

        d = {encoding: str(endpoints.at[filename, encoding]) for encoding in endpoint_encodings}
        d['folder'] = folder
        if patient_id in endpoints_dict:
            endpoints_dict[patient_id].update({filename: d})
        else:
            endpoints_dict[patient_id] = {filename: d}

    # save the json
    with open(os.path.join(output_folder, f'{input_filename}.json'), 'w') as fp:
        json.dump(endpoints_dict, fp, sort_keys=True, indent=4)


if __name__ == '__main__':
    """
    INPUT:
    command line arguments:
    --output-path: where the json files should be saved to
    --excel-path: path to the excel with the data
    --sheet-name: (optional) name of the excel sheet that contains the data
    --endpoints: name(s) of the column that should be used as an end-point --> comma separated if multiple
    --json-path: (optional) json file that should be extended
    
    Excel is expected to have the following columns:
    - filename: name of the slide file (e.g. patient1_I_AE1_AE3_CD8)
    - folder: folder where the slide is saved
    - Patient-ID: first part of the filename (e.g patient1)
    - A column named after the --endpoint argument, e.g.
        N: N group
           0 : 0 in TNM stage
           1 : 1 in TNM stage
           2 : only follow up (>= 2 years), no recurrence
           3 : only follow up (>= 2 years) with recurrence
           
    OUTPUT:
    json file with structure
    patient-id:
        filename:
            folder: folder name
            endpoint: endpoint value
    """
    fire.Fire(make_json)

