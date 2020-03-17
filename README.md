# BT-graph-creation
This repo contains scripts to parse the [ASAP](https://computationalpathologygroup.github.io/ASAP/) 
tif-overlays to maunal annotation xml-files, numpy files with the coordinates and gxl files (graphs).

##Why?
whole slide image (WSI)


##Files and Functionalities
1. `tif_annotations_to_numpy.py`
   - **Input**:
     - `--tif-files-folder`: path to folder with the tif files that contain the overlay. It expects the detected buds to
     have an index of 3 and the lymphocytes to have an index of 9 (value in the image). The detections are marked as rectangles.
     -  `--output-folder`: folder where the output is saved to.
     - `--window-size`: optional. Size (in pixels) of the window that is processed at once.
     - `--spacing-json`: optional. Path to json file that contains the spacing for each whole slide image that should
     be updated. If not provided a new file is created. This is later used to compute the distance between elements.
   - **Output**:
     - 2 text files per WSI: one file each for the tumor buds and the lymphocytes with the coordinates.
     - Json file: contains the spacing for each WSI.


1. `annotation_coordinates_to_xml.py` 
   - **Input**:
   - **Output**:

##Installation requirements
Install the conda environment with `conda env create -f environment.yml`. 

Additionally you need to install [ASAP](https://github.com/computationalpathologygroup/ASAP).
