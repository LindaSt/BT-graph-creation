import glob
import os
import numpy as np
import xml.etree.ElementTree as ET
import argparse
from util.file_parsing import parse_xml


def process_xml_files(files_to_process, output_path):
    """
    Convert the xml files into numpy text files (one for each annotation group)
    """
    all_annotations = {}
    existing_files = set([os.path.basename(s).split('_output')[0] for s in glob.glob(os.path.join(output_path, '*.txt'))])

    for file_path in files_to_process:
        if os.path.basename(file_path).split('_output')[0] in existing_files:
            print('File {} already exists!'.format(file_path))
        else:
            try:
                print('Processing file {}'.format(file_path))
                filename = os.path.basename(os.path.splitext(file_path)[0])

                all_annotations[filename] = parse_xml(file_path)
            except:
                print(f'Something went wrong with file {filename}. Skipping...')

    # save the annotations as numpy readable files
    for filename, annotation_dict in all_annotations.items():
        for group, coordinates in annotation_dict.items():
            # check if we have coordinates for that group
            if len(coordinates) > 0:
                output_file = '{}_coordinates_{}.txt'.format(filename, group)
                coordinates = np.array(coordinates)
                output_file_path = os.path.join(output_path, output_file)

                if not os.path.isfile(output_file):
                    # dot annotations can be directly saved
                    if len(coordinates.shape) < 3:
                        to_save = coordinates
                    elif len(coordinates.shape) == 3:
                        # reformat np array to 2D by flattening the corner coordinates to one list [bl, br, tr, tl]
                        to_save = np.empty((0, 8), float)
                        for hotspot_coord in coordinates:
                            to_save = np.append(to_save,
                                      np.array([item for sublist in hotspot_coord for item in sublist]).reshape(1, 8),
                                      axis=0)
                    np.savetxt(output_file_path, to_save, fmt='%.4f')
                else:
                    print('The coordinates file {} already exists'.format(output_file_path))

# %%


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--xml-files-folder", type=str, required=True)
    parser.add_argument("--output-folder", type=str, required=True)
    args = parser.parse_args()

    # get the file list
    home_path = args.xml_files_folder
    files_to_process = glob.glob(os.path.join(home_path, r'*.xml'))

    # set / create the output
    output_path = args.output_folder
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # process the files
    process_xml_files(files_to_process, output_path)
