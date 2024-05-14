import glob
import os
import numpy as np
import xml.etree.ElementTree as ET
import argparse

from util.xml_parsing import parse_xml, create_xml_tree
from reduce_coord_to_core import in_circle


def reduce_spot_size(xml_data, diameter):
    """Reduce hotspot size in xml data and remove coordinates outside"""
    tl, tr, br, bl = tuple(xml_data['hotspot'][0])
    center_x = (tl[0] + br[0]) / 2
    center_y = (tl[1] + br[1]) / 2
    radius = diameter / 2
    hotspot = [center_x - radius, center_y + radius, center_x + radius, center_y + radius, center_x + radius,
               center_y - radius, center_x - radius, center_y - radius]

    lymph = in_circle(center_x, center_y, radius, xml_data['lymphocytes'])
    buds = in_circle(center_x, center_y, radius, xml_data['tumorbuds'])

    return {'hotspot': np.array([hotspot]), 'lymphocytes': lymph, 'tumorbuds': buds}


def process_xmls(files_to_process, output_path, diameter, overwrite=False):
    existing_files = set([os.path.basename(s) for s in glob.glob(os.path.join(output_path, '*.xml'))])

    for file_path in files_to_process:
        if os.path.basename(file_path).split('.xml')[0] in existing_files and not overwrite:
            print('File {} already exists!'.format(file_path))
        else:
            print('Processing file {}'.format(file_path))
            filename = os.path.basename(os.path.splitext(file_path)[0])
            xml_data = parse_xml(file_path)
            reduced_data = reduce_spot_size(xml_data, diameter)
            xml_out = create_xml_tree(reduced_data)
            if xml_out is not None:
                output_file_path = os.path.join(output_path, os.path.basename(file_path))
                with open(output_file_path, 'w') as f:
                    f.write(xml_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--xml-files-folder", type=str, required=True)
    parser.add_argument("--output-folder", type=str, required=True)
    parser.add_argument("--diameter", type=int, default=1024,
                        help="New size of the core spots in pixels (diameter)")
    parser.add_argument("--overwrite", action='store_true')
    args = parser.parse_args()

    # get the file list
    home_path = args.xml_files_folder
    files_to_process = glob.glob(os.path.join(home_path, r'*.xml'))

    # set / create the output
    output_path = args.output_folder
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # process the files
    process_xmls(files_to_process, output_path, args.diameter, args.overwrite)
