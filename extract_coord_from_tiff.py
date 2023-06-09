import numpy as np
import os
import multiresolutionimageinterface as mir
import glob
from scipy import ndimage
from multiprocessing import Process
import json
import fire


# export PYTHONPATH="${PYTHONPATH}:/opt/ASAP/bin"
# import sys
# sys.path.append(r"C:\Program Files\ASAP 2.0\bin")

class CoordinatesFromTiffExtractor:
    def __init__(self, tif_files_folder: str, output_folder: str, window_size: int = 1024,
                 spacing_json_filepath: str = None, lymph_only: bool = False, buds_only: bool = False,
                 overwrite: bool = False,
                 multi_thread: bool = True, n_threads: int = 8):

        self.lymph_only = lymph_only
        self.buds_only = buds_only
        self.overwrite = overwrite
        self.all_spacing = spacing_json_filepath
        self.multi_thread = multi_thread
        self.n_threads = n_threads
        # get a list of all the tif files to process
        self.files_to_process = glob.glob(os.path.join(tif_files_folder, r'*.tif'))
        #

        # set the indices for the bud and lymphocyte annotations
        # bud output: 0 for background, 1 for foreground, >1 for buds (mostly 2, but overlaps can create 3 or more)
        self.bud_indx = 2
        self.lymp_indx = 1

        # create output folder if it does not exist
        self.output_base_path = output_folder

        # sets the window size of the parsed patch
        self.step_size = window_size

    @property
    def all_spacing(self):
        return self._spacing

    @all_spacing.setter
    def all_spacing(self, json_path):
        # if a spacing json file is provided load it, else create a new dict
        if json_path:
            all_spacing = json.load(json_path)
        else:
            all_spacing = {}
        self._spacing = all_spacing

    @property
    def output_base_path(self):
        return self._output_base_path

    @output_base_path.setter
    def output_base_path(self, output_path):
        # make the output folder if it does not exist
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        self._output_base_path = output_path

    @staticmethod
    def get_tile(img, level0_coords, level):
        # snip region in world coords
        ul_corner = level0_coords[:2]
        dr_corner = level0_coords[2:]

        # convert coords to new level coords
        scale_ratio = img.getLevelDownsample(level)
        dr_corner_x = int(np.round((dr_corner[0] - ul_corner[0]) / scale_ratio))
        dr_corner_y = int(np.round((dr_corner[1] - ul_corner[1]) / scale_ratio))

        # acquire new image
        tile = img.getUCharPatch(int(ul_corner[0]), int(ul_corner[1]), dr_corner_x, dr_corner_y, level)

        return tile

    @staticmethod
    def convert_xy(img, level, coords, newlevel):
        """
        """
        ratio = img.getLevelDownsample(level) / img.getLevelDownsample(newlevel)
        new_coords = np.array(np.array(coords) * ratio, dtype=np.int)

        return new_coords, ratio

    def get_bud_coords(self, img_patch, x_indx, y_indx, ratio):
        """
        """
        patch = np.zeros(img_patch.shape, np.uint8)
        patch[img_patch >= self.bud_indx] = 1
        # do some erosions and dillations to close gaps
        # kernel = np.ones((5, 5), np.uint8)
        # bud_patch = cv2.erode(bud_patch, kernel, iterations=1)
        # bud_patch = cv2.dilate(bud_patch, kernel, iterations=1)

        # Get all center coords of every get_obj
        labeled_img, number_of_objects = ndimage.label(patch)
        if number_of_objects > 0:
            labels = [i for i in np.unique(labeled_img) if i > 0]
            # patch_size = np.sqrt(np.unique(labeled_img == labels[0], return_counts=True)[1][1])
            center_of_masses = ndimage.measurements.center_of_mass(patch, labeled_img, labels)
            # cm_level_0 = [tuple(np.array(cm) * (patch_size * ratio -1) / 2) for cm in center_of_masses]

            # make sure all the buds are actually buds, because some of the non-buds have an artefact that
            # gives them a "frame" of the same encoding as the buds
            # obj_labels = [[img_patch[tuple(t)] for t in [[int(i) for i in com] for com in center_of_masses]]]
            center_of_masses = [cm for cm in center_of_masses if
                                img_patch[tuple([int(i) for i in cm])] == self.bud_indx]
            coords = np.array([[(cm[1] + x_indx) * ratio, (cm[0] + y_indx) * ratio] for cm in center_of_masses])
            if len(coords) == 0:
                coords = None
        else:
            coords = None
        return coords

    @staticmethod
    def get_lymph_coords(patch, x_indx, y_indx, ratio):
        """
        """
        # Get all center coords of every get_obj
        labeled_img, number_of_objects = ndimage.label(patch)
        if number_of_objects > 0:
            labels = [i for i in np.unique(labeled_img) if i > 0]
            # patch_size = np.sqrt(np.unique(labeled_img == labels[0], return_counts=True)[1][1])
            center_of_masses = ndimage.measurements.center_of_mass(patch, labeled_img, labels)
            # cm_level_0 = [tuple(np.array(cm) * (patch_size * ratio -1) / 2) for cm in center_of_masses]

            coords = np.array([[(cm[1] + x_indx) * ratio, (cm[0] + y_indx) * ratio] for cm in center_of_masses])
        else:
            coords = None
        return coords

    def process_files(self):
        """
        Processes a list of tif files with lymphocyte and tumor bud detections (from JMB)
        """
        if self.overwrite:
            print('Existing files will be overwritten!')

        if self.multi_thread:
            chunks = np.array_split(self.files_to_process, self.n_threads)
            prcs = []
            for c in chunks:
                p = Process(target=self.process_chunk, args=(c,))
                p.start()
                prcs.append(p)
            [pr.join() for pr in prcs]
        else:
            for file in self.files_to_process:
                self.process_file(file)

        # Save spacing as json to later calculate the distance
        with open(os.path.join(self.output_base_path, 'spacing_cd8_files.json'), 'w') as fp:
            json.dump(self.all_spacing, fp)

    def process_chunk(self, chunk):
        for file in chunk:
            self.process_file(file)

    def process_file(self, file):
        if "combined" in file:
            file_name = os.path.splitext(os.path.basename(file))[0].split("_combined")[0]
        else:
            file_name = f'{os.path.splitext(os.path.basename(file))[0].split(".tif")[0]}'

        # save the coordinates
        all_bud_coords = np.empty((0, 2), int)
        output_file_bud = os.path.join(self.output_base_path, file_name + "_coordinates_tumorbuds.txt")

        output_file_lymp = os.path.join(self.output_base_path, file_name + "_coordinates_lymphocytes.txt")
        all_lymph_coords = np.empty((0, 2), int)

        # the annotations correspond to level 1
        level = 1

        print("Processing: {}".format(file_name))

        # check if the files are already present
        if (not os.path.isfile(output_file_lymp) or not os.path.isfile(output_file_bud)) or self.overwrite:
            img_obj = mir.MultiResolutionImageReader().open(file)
            assert file_name not in self.all_spacing
            self.all_spacing[file_name] = img_obj.getSpacing()[0]
            dim_x, dim_y = img_obj.getLevelDimensions(level)

            # sliding window over the whole image
            for y_indx in range(0, dim_y, self.step_size):
                for x_indx in range(0, dim_x, self.step_size):
                    coords_level_0, ratio = self.convert_xy(img_obj, level,
                                                            [x_indx, y_indx, x_indx + self.step_size,
                                                             y_indx + self.step_size],
                                                            0)

                    # get the patch
                    img_patch = self.get_tile(img_obj, coords_level_0, level)
                    img_patch = img_patch.squeeze()
                    # print(np.unique(img_patch, return_counts=True))
                    # get the patch
                    # img_patch_0 = get_tile(img_obj, coords_level_0, 0)
                    # img_patch_0 = img_patch.squeeze()

                    # get the coordinates of the center of each annotation
                    # Coordinates are additionally multiplied by 2 because TIF annotations level 0 is actually level 1 on the WSI
                    # buds
                    if not self.lymph_only:
                        bud_coords = self.get_bud_coords(img_patch, x_indx, y_indx, ratio * 2)

                        if bud_coords is not None:
                            all_bud_coords = np.append(all_bud_coords, bud_coords, axis=0)

                    # lymphocytes
                    if not self.buds_only:
                        lymph_patch = np.zeros(img_patch.shape, np.uint8)
                        lymph_patch[img_patch == self.lymp_indx] = 1
                        lymp_coords = self.get_lymph_coords(lymph_patch, x_indx, y_indx, ratio * 2)
                        if lymp_coords is not None:
                            all_lymph_coords = np.append(all_lymph_coords, lymp_coords, axis=0)

            if not self.buds_only and len(all_lymph_coords) > 0:
                np.savetxt(output_file_lymp, all_lymph_coords, fmt='%.3f')
            if not self.lymph_only and len(all_bud_coords) > 0:
                np.savetxt(output_file_bud, all_bud_coords, fmt='%.3f')
        else:
            print('The coordinates files {} or {} already exist. To overwrite use --overwrite.'.format(output_file_lymp,
                                                                                                       output_file_bud))


if __name__ == '__main__':
    fire.Fire(CoordinatesFromTiffExtractor).process_files()

