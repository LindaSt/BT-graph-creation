import json
import os
import glob
import shutil
import fire


def split_dataset(split_json: str, gxl_folder_path: str, output_path: str = None):
    with open(split_json) as data_file:
        datasplit_dict = json.load(data_file)
    file_id_to_folder = {file_id: [split, cls] for split, d in datasplit_dict.items() for cls, file_ids in d.items()
                         for file_id in file_ids}

    # set up the folder structure
    folders_to_create = [os.path.join(output_path, split, cls) for split, d in datasplit_dict.items() for cls in
                         d.keys()]
    _ = [os.makedirs(p) for p in folders_to_create if not os.path.isdir(p)]

    # get file list
    gxl_files = glob.glob(os.path.join(gxl_folder_path, '*'))

    # move files
    for gxl in gxl_files:
        file_id = os.path.basename(gxl).split('.gxl')[0]
        try:
            if output_path is None:
                outfolder = os.path.join(gxl_folder_path, "/".join(file_id_to_folder[file_id]))
                shutil.move(src=gxl, dst=os.path.join(os.path.join(outfolder, os.path.basename(gxl))))
            else:
                outfolder = os.path.join(output_path, "/".join(file_id_to_folder[file_id]))
                shutil.copy2(src=gxl, dst=os.path.join(os.path.join(outfolder, os.path.basename(gxl))))
        except KeyError:
            print(f'Patient ID {file_id} not found in json file')
            continue


if __name__ == '__main__':
    fire.Fire(split_dataset)

