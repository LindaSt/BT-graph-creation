import os
import numpy as np
import glob
import fire
import xml.etree.ElementTree as ET
from PIL import Image
from openslide import open_slide
import pandas as pd
from multiprocessing import Process


class BTPatchExtractor:
    def __init__(self, file_path: str, output_path: str, asap_xml_path: str, overwrite: bool = False,
                 hotspot: bool = False, level: int = 0, lymph_patch_size: int = 200, tb_patch_size: int = 200,
                 matched_files_excel: str = None, n_threads: int = 6, no_multi_thread: bool = False,
                 staining: str = 'CD8'):
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
        :param lymph_patch_size: int (optional, default is 200 pixel)
            size of the patch around the lymphocyte coordinates
        :param tb_patch_size: int (optional, default is 200 pixel)
            size of the patch around the tumor bud coordinates
        :param level: int (optional, default is 0)
            Level of the mrxs file that should be used for the conversion.
        :param matched_files_excel: str
            Optional. If provided, then this file will be used to match the xmls to the mrxs file names
            (specify info in MATCHED_EXEL_INFO)
        :param staining: str (optional, default is 'CD8')
            Staining ID that is searched for when generating the list of WSI files
        """
        # initiate the mandatory elements
        self.file_path = file_path
        self.output_path = output_path
        # instantiate optional parameters
        self.coord_path = asap_xml_path
        self.overwrite = overwrite
        self.staining = staining
        self.level = level
        self.matched_files_excel = matched_files_excel
        self.extract_hotspot = hotspot
        self.lymph_patch_size = lymph_patch_size
        self.tb_patch_size = tb_patch_size
        self.n_threads = n_threads

        self.groups = ['tumorbuds', 'lymphocytes', 'hotspot'] if self.extract_hotspot else ['tumorbuds', 'lymphocytes']

        self.matched_excel_info = {'wsi_col': 'CD8 Filename', 'file_id_col': 'Algo coordinates text file ID', 'sheet_name': 'Masterfile',
                              'folder_col': 'Folder'}

        self.no_multi_thread = no_multi_thread

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
            # if we have a matched excel
            if self.matched_files_excel:
                files = self.file_path
            else:
                files = glob.glob(os.path.join(self.file_path, f'*{self.staining}*.mrxs')) + glob.glob(os.path.join(self.file_path, f'*{self.staining}.ndpi'))
            if len(files) == 0:
                print(f'No WSIs found in folder {self.file_path}!')
                exit(-1)
            return files
        # if we have just a single file
        elif os.path.isfile(self.file_path):
            return [self.file_path]
        else:
            print(f'Folder {self.file_path} not found.')
            exit(-1)

    @property
    def coord_files(self):
        if os.path.isdir(self.file_path):
            # if we have a matched excel
            if self.matched_files_excel:
                files = self.coord_path
            else:
                files = glob.glob(os.path.join(self.coord_path, f'*{self.staining}*asap.xml'))
            return files
        # if we have just a single file
        elif os.path.isfile(self.file_path):
            return [self.coord_path]

    @property
    def files_to_process(self):
        if self.matched_files_excel:
            return self._get_matched_files_excel()

        else:
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

    def _get_matched_files_excel(self):
        files_to_process = []
        df = self.parse_matched_files_excel()
        error = []
        for wsi_file, wsi_folder, xml_name in zip(df[self.matched_excel_info['wsi_col']],
                                                  df[self.matched_excel_info['folder_col']],
                                                  df[self.matched_excel_info['file_id_col']]):

            output_files_folder_path = os.path.join(self.output_path, f'{xml_name}-level{self.level}')
            wsi_path = os.path.join(self.wsi_files, os.path.join(wsi_folder, wsi_file))
            xml_coord_path = os.path.join(self.coord_files, f'{xml_name}_output_asap.xml')
            # check if files listed in excel actually exist
            if not os.path.isfile(wsi_path):
                print(f'WSI {wsi_path} not found (skipping file)')
                error.append(wsi_path)
                continue
            if not os.path.isfile(xml_coord_path):
                print(f'XML {xml_coord_path} not found (skipping file)')
                error.append(xml_coord_path)
                continue

            # skip if output foler exists if overwrite = False
            if not self.overwrite and os.path.ispath(output_files_folder_path):
                print(
                    f'File {output_files_folder_path} already exists. Output saving is skipped. To overwrite add --overwrite.')
                continue

            files_to_process.append((output_files_folder_path, wsi_path, xml_coord_path))

        return files_to_process

    def process_files(self):
        if self.no_multi_thread:
            for (output_folder_path, wsi_path, coord_path) in self.files_to_process:
                self.process_file(output_folder_path, wsi_path, coord_path)
        else:
            # process the files with coordinates
            chunks = np.array_split(self.files_to_process, self.n_threads)
            prcs = []
            for c in chunks:
                p = Process(target=self.process_chunk, args=(c,))
                p.start()
                prcs.append(p)
            [pr.join() for pr in prcs]

    def process_chunk(self, chunk):
        for c in chunk:
            output_folder_path, wsi_path, coord_path = tuple(c)
            self.process_file(output_folder_path, wsi_path, coord_path)

    def process_file(self, output_folder_path, wsi_path, coord_path):
        # make the output folder if it does not exist
        if not os.path.isdir(output_folder_path):
            os.makedirs(output_folder_path)
        # open the wsi and get the coordinates
        wsi_img = open_slide(wsi_path)
        group_coordinates = self.parse_xml(coord_path)
        # iterate over the objects
        offset_tb = len(group_coordinates['lymphocytes'])
        for group, coords in group_coordinates.items():
            for id, coord in coords:
                # to ensure that enumeration is continous (in ASAP xml it starts at 0 for each group)
                if group == 'tumorbuds':
                    id += offset_tb
                output_file_path = os.path.join(output_folder_path,
                                                f'{os.path.basename(output_folder_path)}_{group}_{id}_{"-".join([str(i) for i in coord])}.png')
                # extract the patch
                top_left_coord, size = self.get_rectangle_info(coord, group)
                png = self.extract_crop(wsi_img, top_left_coord, size)
                # save the image
                print(f'Saving image {output_file_path}')
                Image.fromarray(png[:, :, :3]).save(output_file_path)

    def get_rectangle_info(self, asap_coord, group):
        if group == 'hotspot':
            top_left_coord = [int(i) for i in asap_coord[0]]
            size = int(asap_coord[2][0] - asap_coord[0][0])
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

    def parse_matched_files_excel(self) -> pd.DataFrame:
        # TODO: this is probably not working anymore
        df = pd.read_excel(self.matched_files_excel, sheet_name=self.matched_excel_info['sheet_name'], engine='openpyxl')
        # remove two empty top lines and set third line to header
        df.columns = df.iloc[2]
        df = df.iloc[3:]
        # drop all rows that do not contain 0 or 1 in column "Need resection?" (excluded because no data available)
        df = df.drop(df[~df["Need resection?"].isin([0, 1])].index)
        # drop all rows that do not contain a file name
        # TODO: make this neater
        df = df[df[self.matched_excel_info['wsi_col']].notna()]
        df = df[df[self.matched_excel_info['file_id_col']].notna()]
        df = df.drop(df[df[self.matched_excel_info['wsi_col']].isin(["tbd", "na"])].index)
        df = df.drop(df[df[self.matched_excel_info['file_id_col']].isin(["tbd", "na"])].index)
        return df


if __name__ == '__main__':
    fire.Fire(BTPatchExtractor).process_files()
