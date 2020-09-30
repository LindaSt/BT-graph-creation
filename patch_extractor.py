import os
import numpy as np
import glob
import fire
import xml.etree.ElementTree as ET
from PIL import Image
from openslide import open_slide


class BTPatchExtractor:
    def __init__(self, file_path: str, output_path: str, asap_xml_path: str, overwrite: bool = False, hotspot: bool = False, level: int = 0,
                 lymph_patch_size: int = 70, tb_patch_size: int = 250):
        """
        This Object extracts (patches of) an mrxs file to a png format.

        :param file_path: string
            path to the mrxs single file or folder of files.
        :param output_path: string
            path to the output folder. The output format is the same name as the mrxs file,
            with an appendix if multiple patches are extracted.
        :param asap_xml_path: string (optional)
            Path to the coordinate xml files (created with ASAP) single file or folder of files
            If not provided, the full image is converted into a png.
        :param overwrite: bool (optional)
            overides exisiting extracted patches (default is False)
        :param hotspot: bool (optional)
            set if hotspot should also be extracted (default False)
        :param level: int (optional)
            Level of the mrxs file that should be used for the conversion (default is 0).
        """
        # initiate the mandatory elements
        self.file_path = file_path
        self.output_path = output_path
        # instantiate optional parameters
        self.coord_path = asap_xml_path
        self.overwrite = overwrite
        self.staining = 'CD8'
        self.level = level
        self.extract_hotspot = hotspot
        self.lymph_patch_size = lymph_patch_size
        self.tb_patch_size = tb_patch_size

        self.groups = ['tumorbuds', 'lymphocytes', 'hotspot'] if self.extract_hotspot else ['tumorbuds', 'lymphocytes']

        self.process_files()

    @property
    def output_path(self):
        return self._output_path

    @output_path.setter
    def output_path(self, output_path):
        # make the output folder if it does not exist
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        self._output_path = output_path

    @property
    def wsi_files(self):
        if os.path.isdir(self.file_path):
            files = glob.glob(os.path.join(self.file_path, f'*{self.staining}.mrxs')) + glob.glob(os.path.join(self.file_path, f'*{self.staining}.ndpi'))
            return files
        else:
            return [self.file_path]

    @property
    def coord_files(self):
        return glob.glob(os.path.join(self.coord_path, f'*{self.staining}*asap.xml')) if os.path.isdir(self.coord_path) else [self.coord_path]

    @property
    def files_to_process(self):
        # we only have one file to process
        if len(self.wsi_files) == 1:
            filename = os.path.splitext(os.path.basename(self.file_path))[0]
            output_folder = os.path.join(self.output_path, filename)
            # skip if overwrite = False and folder exists
            if not self.overwrite and os.path.isdir(output_folder):
                print(f'Folder {output_folder} already exists. Output saving is skipped. To overwrite add --overwrite.')
            else:
                return [(output_folder, self.file_path, self.coord_path)]

        # we have multiple files to process
        else:
            # create a list of the paired mrxs and coordinate files
            # only take files that have a corresponding coordinates file
            files_to_process = []
            for wsi_path in self.wsi_files:
                filename = os.path.splitext(os.path.basename(wsi_path))[0]
                output_folder = os.path.join(self.output_path, filename)
                # skip if overwrite = False and folder exists
                if not self.overwrite and os.path.isdir(output_folder):
                    print(
                        f'Folder {output_folder} already exists. Output saving is skipped. To overwrite add --overwrite.')
                    continue

                checked = []
                for coord_file in self.coord_files:
                    if filename in coord_file:
                        checked.append(coord_file)
                if len(checked) != 1:
                    print(
                        f'File {filename}.mrxs does not have a / too many corresponding xml file/s. File will be skipped.')
                else:
                    files_to_process.append((output_folder, wsi_path, checked.pop()))

            return files_to_process

    def process_files(self):
        # process the files with coordinates
        if self.coord_path and ((os.path.isdir(self.file_path) and os.path.isdir(self.coord_path)) or (
                os.path.isfile(self.file_path) and os.path.isfile(self.coord_path))):

            for output_folder_path, wsi_path, coord_path in self.files_to_process:
                # make the output folder if it does not exist
                if not os.path.isdir(output_folder_path):
                    os.makedirs(output_folder_path)
                # open the wsi and get the coordinates
                wsi_img = open_slide(wsi_path)
                group_coordinates = self.parse_xml(coord_path)
                # iterate over the objects
                for group, coords in group_coordinates.items():
                    for id, coord in coords:
                        output_file_path = os.path.join(output_folder_path, f'{os.path.basename(output_folder_path)}-{group}-{id}.png')
                        # extract the patch
                        top_left_coord, size = self.get_rectangle_info(coord, group)
                        png = self.extract_crop(wsi_img, top_left_coord, size)
                        # save the image
                        print(f'Saving image {output_file_path}')
                        Image.fromarray(png[:, :, :3]).save(output_file_path)
        else:
            # Something went wrong
            print('Path(s) are invalid.')

    def get_rectangle_info(self, asap_coord, group):
        if group == 'hotspot':
            top_left_coord = [int(i) for i in asap_coord[0]]
            size = asap_coord[2][0] - asap_coord[0][0]

        elif group == 'lymphocytes':
            top_left_coord = [int(i-self.lymph_patch_size/2) for i in asap_coord]
            size = self.lymph_patch_size

        elif group == 'tumorbuds':
            top_left_coord = [int(i-self.tb_patch_size/2) for i in asap_coord]
            size = self.tb_patch_size

        else:
            print('Invalid group')
            return

        return top_left_coord, size

    def parse_xml(self, file_path):
        # reads the xml files and retrieves the coordinates of all elements with the coord_annotation_tag
        tree = ET.parse(file_path)
        root = tree.getroot()

        annotations_elements = {g: [] for g in self.groups}

        for i in root.iter('Annotation'):
            if i.attrib['PartOfGroup'] in annotations_elements:
                annotations_elements[i.attrib['PartOfGroup']].append(i)

        annotations = {g: [] for g in self.groups}
        for group, element_list in annotations_elements.items():
            for element in element_list:
                if element.attrib['Type'] == 'Dot':
                    annotation = [[float(i.attrib['X']), float(i.attrib['Y'])] for i in element.iter('Coordinate')][0]
                else:
                    annotation = [[float(i.attrib['X']), float(i.attrib['Y'])] for i in element.iter('Coordinate')]

                # get the id (used as node id later)
                annot_id = int(element.attrib['Name'].split(' ')[-1])
                annotations[group].append((annot_id, annotation))

        return annotations

    def extract_crop(self, wsi_img, top_left_coord, size):
        # crop the region of interest from the mrxs file on the specified level
        # get the level and the dimensions
        id_level = np.argmax(np.array(wsi_img.level_downsamples) == self.level)
        dims = wsi_img.level_dimensions[id_level]

        # TODO make sure the dimension we want to crop are within the image dimensions

        # extract the region of interest
        img = wsi_img.read_region(top_left_coord, id_level, (size, size))

        # Convert to img
        img = np.array(img)
        img[img[:, :, 3] != 255] = 255
        return img


if __name__ == '__main__':
    fire.Fire(BTPatchExtractor)