# BT-graph-creation
This repo contains scripts to parse the [ASAP](https://computationalpathologygroup.github.io/ASAP/) 
tif-overlays to manual annotation xml-files, numpy files with the coordinates and gxl files (graphs).

## Why?
The annotations of the tumor bud and lymphocyte detections are in tif format, where each of them are represented as 
small squares (value 3 or 9). The tif image can be loaded as an overlay in ASAP.
 
The coordinates on the whole slide image (WSI) of all these elements have to be extracted. We also want to be able
to reduce the elements to a certain area (e.g. hotspot). The coordinates are saved as text files as well as xml files,
which are compatible to be loaded in ASAP.

Finally, the graph-representations are saved in gxl format.

## Files and Functionalities
### Convert detections to graphs (gxl files)
1. `extract_coord_from_tiff`: Extracts the detection output
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

1. `coord_to_xml.py`: Creates an ASAP-xml file from the coordinate text files
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

1. `reduce_coord_to_hotspot.py`: retrieves only the elements within a certain area. It expects the hotspot files to have 
   the same name as the matching coordinate files (run `check_hotspot_xml.py` first)
   - **Input**:
     - `--xml-hotspot-folder`: folder with the xml files of the hotspots
     - `--coordinate-txt-files-folder`: folder with the text files of the coordinates of the lymphocytes and tumor buds
     - `--output-folder`: folder where the output is saved to
   - **Output**: Contains only the detected elements within the selection.
     - Folder with coordinate text files: three files per WSI (one for the hotspots, one for the lymphocytes, one
       for the tumor buds)
     - Folder with xml files: one per WSI, can be loaded in ASAP. Useful for quality control and manual correction of the
       annotations.
       
1. `create_gxl_files.py`: creates the graphs as gxl files
   - **Input**:
     - `--asap_xml_files_folder`: path to the folder with the coordinates text files
     - `--edge-definition-tb-to-l` and `--edge-definition-tb-to-tb` have the following options: 
       - `radius-x`: connect elements in radius X (in mikrometer)
       - `to-X-nn`: connect to k closest elements where X is the number of neighbours
     - `--fully-connected`: supersedes the `--edge-definition*`. Following options:
       - `all`: fully connected graph
       - `lymphocytes`: only the lymphocytes are fully connected
       - `tumorbuds`: only the tumorbuds are fully connected
       - `--output-folder`: path to where output folder should be created
     - `--node-feature-csvs`: optional. Path to folder with csv files that contain additional node features.
       The first column needs to have the node index number. The headers will be used as the feature name.
       If there is a column named "filename", it will be dropped.
     - `--spacing-json`: optional. Path to json file that contains the spacing for each whole slide image. 
     It is needed to compute the distance between elements. (default is 0.242797397769517, which corresponds to level 0 
       for the slide scanner used in this project)
   - **Output**: 
     - Folder with one gxl file per hotspot, which contains the graph (same structure as the gxl files 
       from the IAM Graph Databse). The x,y coordinates and edge distance labels are in mikro-meter.
       The graphs are not pre-processed (features are not normalized, x,y coordinates are not centered)!

### Extract image-based features
1. `patch_extractor.py`: let's you extract patches from a single mrxs file or a folder based on annotations in
an ASAP xml file (expects the following annotations groups: `lymphocytes`, `tumorbuds` and `hotspot`).
    - **Input**:
        - `file_path`: path to the mrxs single file or folder of files. 
        - `output_path`: path to the output folder. The output format is the same name as the mrxs file, with an appendix if multiple patches are extracted.
        - `asap_xml_path`: Path to the coordinate xml files (created with ASAP) single file or folder of files
        - `overwrite`: optional. Overrides existing extracted patches (default is False)
        - `hotspot`: optional. Set if hotspot should also be extracted (default False)
        - `lymph_patch_size`: optional. Size of the patch around the lymphocyte coordinates in pixels (default is 300))
        - `tb_patch_size`: optional. Size of the patch around the tumorbud coordinates in pixels (default is 300))
        - `level`: optional. Level of the mrxs file that should be used for the conversion (default is 0).
        - `matched_files_excel`: optional. If provided, then this file will be used to match the xmls to the mrxs file names
            (specify info in MATCHED_EXEL_INFO)
    - **Output**:
        - Folder with one sub-folder per mrxs file containing the corresponding cropped patches.


(### Convert endpoints to dataset split (cxl files))

TODO: update this
1. `make_endpoint_json.py`: sets up json dictionary based on an excel file
   - **Input**:
     - `--output-path`: where the json files should be saved to
     - `--excel-path`: path to the excel with the data
     - `--sheet-name`: (optional) name of the excel sheet that contains the data
     - `--endpoints`: name(s) of the column that should be used as an end-point --> comma separated if multiple
     - `--json-path`: (optional) json file that should be extended
    
     The excel file is expected to have the following columns:
     - filename: name of the slide file (e.g. patient1.mrxs)
     - folder: folder where the slide is saved
     - Patient-ID: first part of the filename (e.g patient1)
     - A column(s) named after the `--endpoints` argument, e.g. N group:<br />
       >  0 : N0 in TNM stage<br />
          1 : N1 in TNM stage<br />
          2 : only follow up (>= 2 years), no recurrence<br />
          3 : only follow up (>= 2 years) with recurrence
           
   - **Input**:
     - json file with structure: <br />
       >patient-id: patient id<br />
       >>filename:  file name<br />
         folder: folder name<br />
         endpoint1: endpoint value<br />
         endpoint2: endpoint value<br />
         ...

1. `endpoints_json_to_cxl.py`: splits the data from the json file into train, vaild and test cxl files. 
    It is split based on the patient-id, so if your slides are very unequally distributed amongh the patients,
    you will have make sure you still have roughly the percentage split you want.
   - **Input**:
     - `--output-path`: where the cxl files should be saved to
     - `--json-path`: path to the json file with the endpoint(s) data
     - `--endpoint`: name of the variable (in the dict) that encodes the end-point (only 1!)
     - `--dataset-name`: name that should be given to the dataset in the xml files
     - `--split`: (optional, default 0.4) how much should be split off for train and test (will be split in half for test valid) 
     - `--seed`: (optional, default 42) set seed for split
        
     json file needs to be in the following structure (output of `make-endpoint-json.py`): <br />
       >patient-id: patient-id<br />
       >>filename: file name<br />
         folder: folder name<br />
         endpoint1: endpoint value<br />
         endpoint2: endpoint value<br />
         ...

   - **Output**:
      - 3 xml files (train.cxl, vaild.cxl and test.cxl) as in the IAMDB structure:
          > \<GraphCollection><br />
          \<BT-Hotspots counts="49" endpoint="N"><br />
                \<print file="file1.gxl" class="0"/><br />
                \<print file="file2.gxl" class="0"/><br />
                ...<br />
                \<print file="file49.gxl" class="2"/><br />
          </BT-Hotspots><br />
          </GraphCollection>
          
          Where `BT-Hotspots` is the `--dataset-name`, `counts` is the number of elements in the subset,
          and `enpoint` is the `--endpoint`.

### Other utility scripts
- `check_hotspot_xmls`: Checks if hotspot xml(s) have the right format and corrects them if possible 
  (for an example of the expected structure see the script file). It also changes the filename according to the provided 
  excel master-file to match the output text file names.
  - Input:
    - `--input-path`: path to the xml hotspot file/folder
    - `--output-path`: path to the folder, where the corrected xml files should be saved to
    - `--overwrite`: overwrites existing xml files (default is False)
    - `--matched-files-excel`: optional (default False). If specified, the files are matched based on this excel
      (check script to change column headers)
  - Output:
    - File / folder of corrected hotspot xml files

## Installation requirements
Install the conda environment with `conda env create -f environment.yml`. 

Additionally you need to install [ASAP](https://github.com/computationalpathologygroup/ASAP), if you want
to use the `extract_coord_from_tiff.py` script.
