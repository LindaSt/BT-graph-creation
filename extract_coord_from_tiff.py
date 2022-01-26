import numpy as np
import os
import multiresolutionimageinterface as mir
import glob
from scipy import ndimage
import argparse
import json
import cv2

# export PYTHONPATH="${PYTHONPATH}:/opt/ASAP/bin"

def get_tile(img, level0_coords, level):
    """
    """
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


def convert_xy(img, level, coords, newlevel):
    """
    """
    ratio = img.getLevelDownsample(level) / img.getLevelDownsample(newlevel)
    new_coords = np.array(np.array(coords) * ratio, dtype=np.int)

    return new_coords, ratio


def get_obj_coords(patch, x_indx, y_indx, ratio):
    """
    """
    # Get all center coords of every get_obj
    labeled_img, number_of_objects = ndimage.label(patch)
    if number_of_objects > 0:
        labels = [i for i in np.unique(labeled_img) if i > 0]
        # patch_size = np.sqrt(np.unique(labeled_img == labels[0], return_counts=True)[1][1])
        center_of_masses = ndimage.measurements.center_of_mass(patch, labeled_img, labels)
        #cm_level_0 = [tuple(np.array(cm) * (patch_size * ratio -1) / 2) for cm in center_of_masses]

        coords = np.array([[(cm[1] + x_indx) * ratio, (cm[0] + y_indx) * ratio] for cm in center_of_masses])
    else:
        coords = None
    return coords


def process_files(files_to_process, output_base_path, step_size, bud_indx, lymp_indx, spacing_json_filepath=None,
                  lymph_only=False, buds_only=False, overwrite=False):
    """
    Processes a list of tif files with lymphocyte and tumor bud detections (from JMB)
    """

    # if a spacing json file is provided load it, else create a new dict
    if spacing_json_filepath:
        all_spacing = json.load(spacing_json_filepath)
    else:
        all_spacing = {}

    print('Existing files will be overwritten!')
    for file in files_to_process:
        if "combined" in file:
            file_name = os.path.splitext(os.path.basename(file))[0].split("_combined")[0]
        else:
            file_name = f'{os.path.splitext(os.path.basename(file))[0].split(".tif")[0]}'

        # save the coordinates
        all_bud_coords = np.empty((0, 2), int)
        output_file_bud = os.path.join(output_base_path, file_name + "_coordinates_tumorbuds.txt")

        output_file_lymp = os.path.join(output_base_path, file_name + "_coordinates_lymphocytes.txt")
        all_lymph_coords = np.empty((0, 2), int)

        # the annotations correspond to level 1
        level = 1

        print("Processing: {}".format(file_name))

        # check if the files are already present
        if (not os.path.isfile(output_file_lymp) or not os.path.isfile(output_file_bud)) or overwrite:
            img_obj = mir.MultiResolutionImageReader().open(file)
            assert file_name not in all_spacing
            all_spacing[file_name] = img_obj.getSpacing()[0]
            dim_x, dim_y = img_obj.getLevelDimensions(level)

            # sliding window over the whole image
            for y_indx in range(0, dim_y, step_size):
                for x_indx in range(0, dim_x, step_size):
                    coords_level_0, ratio = convert_xy(img_obj, level,
                                                       [x_indx, y_indx, x_indx + step_size, y_indx + step_size],
                                                       0)

                    # get the patch
                    img_patch = get_tile(img_obj, coords_level_0, level)
                    img_patch = img_patch.squeeze()
                    # print(np.unique(img_patch, return_counts=True))
                    # get the patch
                    # img_patch_0 = get_tile(img_obj, coords_level_0, 0)
                    # img_patch_0 = img_patch.squeeze()

                    # get the coordinates of the center of each annotation
                    # Coordinates are additionally multiplied by 2 because TIF annotations level 0 is actually level 1 on the WSI

                    # buds
                    if not lymph_only:
                        bud_patch = np.zeros(img_patch.shape, np.uint8)
                        bud_patch[img_patch == bud_indx] = 1
                        # do some erosions and dillations to close gaps
                        kernel = np.ones((5, 5), np.uint8)
                        bud_patch = cv2.erode(bud_patch, kernel, iterations=1)
                        bud_patch = cv2.dilate(bud_patch, kernel, iterations=1)

                        bud_coords = get_obj_coords(bud_patch, x_indx, y_indx, ratio * 2)
                        if bud_coords is not None:
                            bud_coords = bud_coords
                            all_bud_coords = np.append(all_bud_coords, bud_coords, axis=0)

                    # lymphocytes
                    if not buds_only:
                        lymph_patch = np.zeros(img_patch.shape, np.uint8)
                        lymph_patch[img_patch == lymp_indx] = 1
                        lymp_coords = get_obj_coords(lymph_patch, x_indx, y_indx, ratio * 2)
                        if lymp_coords is not None:
                            lymp_coords = lymp_coords
                            all_lymph_coords = np.append(all_lymph_coords, lymp_coords, axis=0)

            if not buds_only and len(all_lymph_coords) > 0:
                np.savetxt(output_file_lymp, all_lymph_coords, fmt='%.3f')
            if not lymph_only and len(all_bud_coords) > 0:
                np.savetxt(output_file_bud, all_bud_coords, fmt='%.3f')
        else:
            print('The coordinates files {} or {} already exist. To overwrite use --overwrite.'.format(output_file_lymp, output_file_bud))

    # Save spacing as json to later calculate the distance
    with open(os.path.join(output_base_path, 'spacing_cd8_files.json'), 'w') as fp:
        json.dump(all_spacing, fp)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tif-files-folder", type=str, required=True)
    parser.add_argument("--output-folder", type=str, required=True)
    parser.add_argument("--window-size", type=int, default=1024)
    parser.add_argument("--spacing-json", type=str, default=None)
    parser.add_argument("--lymph-only", action='store_true', default=False)
    parser.add_argument("--buds-only", action='store_true', default=False)
    parser.add_argument("--overwrite", action='store_true', default=False)
    args = parser.parse_args()

    spacing_json = args.spacing_json

    # set the indices for the bud and lymphocyte annotations
    # bud output: 0 for background, 1 for foreground, >1 for buds (mostly 2, but overlaps can create 3 or more)
    bud_indx = 2
    lymp_indx = 1

    # get a list of all the tif files to process
    tif_files = os.path.join(args.tif_files_folder, r'*.tif')
    files_to_process = glob.glob(tif_files)

    # create output folder if it does not exist
    output_base_path = args.output_folder
    if not os.path.isdir(output_base_path):
        os.mkdir(output_base_path)

    # sets the window size of the parsed patch
    step_size = args.window_size

    # read the annotations and create the coordinate files
    process_files(files_to_process=files_to_process, output_base_path=output_base_path, step_size=step_size,
                  bud_indx=bud_indx, lymp_indx=lymp_indx, spacing_json_filepath=spacing_json,
                  lymph_only=args.lymph_only, buds_only=args.buds_only, overwrite=args.overwrite)
