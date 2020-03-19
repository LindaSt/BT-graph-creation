# BT-graph-creation
This repo contains scripts to parse the [ASAP](https://computationalpathologygroup.github.io/ASAP/) 
tif-overlays to maunal annotation xml-files, numpy files with the coordinates and gxl files (graphs).

## Why?
The annotations of the tumor bud and lymphocyte detections are in tif format, where each of them are represented as 
small squares (value 3 or 9). The tif image can be loaded as an overlay in ASAP.
 
The coordinates on the whole slide image (WSI) of all these elements have to be extracted. We also want to be able
to reduce the elements to a certain area (e.g. hotspot). The coordinates are saved as text files as well as xml files,
which are compatible to be loaded in ASAP.

Finally, the graph-representations are saved in gxl format.

## Files and Functionalities
1. `tif_annotations_to_numpy.py`: Extracts the detection output
   - **Input**:
     - `--tif-files-folder`: path to folder with the tif files that contain the overlay. It expects the detected buds to
     have an index of 3 and the lymphocytes to have an index of 9 (value in the image). The detections are marked as rectangles.
     - `--output-folder`: folder where the output is saved to.
     - `--window-size`: optional. Size (in pixels) of the window that is processed at once.
     - `--spacing-json`: optional. Path to json file that contains the spacing for each whole slide image that should
     be updated. If not provided a new file is created. This is later used to compute the distance between elements.
   - **Output**:
     - Two text files per WSI: one file each for the tumor buds and the lymphocytes with the coordinates.
     - Json file: contains the spacing for each WSI.

1. `annotation_coordinates_to_xml.py`: Creates an ASAP-xml file from the coordinate text files
   - **Input**:
     - `--coordinates-txt-files-folder`: folder with the text files of the coordinates for 
      'lymphocytes', 'tumorbuds' or 'hotspot'(one file per annotation, 
        e.g. output of `tif_annotations_to_numpy.py`).
     - `--output-folder`: folder where the output is saved to.
   - **Output**:
     - One xml file per WSI that can be loaded in ASAP.

1. `xml_to_txt_file.py`: convert ASAP-xml files to coordinate text files
   - **Input**:
     - `--xml-files-folder`: folder with the xml files to be converted
     - `--output-folder`: folder where the output is saved to.
   - **Output**:
     - Up to three text files per WSI (depending on how many annotations there are in the xml).

1. `create_hotspot_graphs.py`: retrieves only the elements within a certain area.
   - **Input**:
     - `---xml-hotspots`: folder with the xml files of the hotspots
     - `--coordinate-txt-files`: folder with the text files of the coordinates of the lymphocytes and tumor buds
     - `--output-folder`: folder where the output is saved to
   - **Output**: Contains only the detected elements within the hotspot.
     - Folder with coordinate text files: three files per WSI (one for the hotspot, one for the lymphocytes, one
       for the tumor buds)
     - Folder with xml files: one pre WSI, can be loaded in ASAP. Useful for quality control and manual correction of the
       annotations.
       
1. TODO: `make_gxl_graph_dataset.py`: creates the gxl files (graphs)
   - **Input**:
     - `--coordinates-txt-files-folder`: path to the folder with the coordinates text files
     - `--spacing-json`: Path to json file that contains the spacing for each whole slide image. 
     It is needed to compute the distance between elements.
     - `--edge-definition-tb-to-l` and `--edge-definition-tb-to-tb` have the following options: 
       - `radius-x`: connect elements in radius X
       - `to-knn`: connect to k closest elements
       - `to-all`: connect to all elements
     - `--fully-connected`: make a fully connected graph (supersedes the `--edge-definition...`)
     - `--output-folder`: folder where the output is saved to
   - **Output**: 
     - one gxl file per hotspot, which contains the graph (same structure as the gxl files from the IAM Graph Databse)

## Installation requirements
Install the conda environment with `conda env create -f environment.yml`. 

Additionally you need to install [ASAP](https://github.com/computationalpathologygroup/ASAP).
