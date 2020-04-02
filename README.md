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
### Convert detections to graphs (gxl files)
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
       
1. In development: creating the gxl files (graphs)

### Convert endpoints to dataset split (cxl files)
1. `make-endpoint-json.py`: sets up json dictionary based on excel file
   - **Input**:
     - `--output-path`: where the json files should be saved to
     - `--excel-path`: path to the excel with the data
     - `--sheet-name`: (optional) name of the excel sheet that contains the data
     - `--endpoints`: name(s) of the column that should be used as an end-point --> comma separated if multiple
     - `--json-path`: (optional) json file that should be extended
    
     The excel file is expected to have the following columns:
     - filename: name of the slide file (e.g. patient1_I_AE1_AE3_CD8)
     - folder: folder where the slide is saved
     - Patient-ID: first part of the filename (e.g patient1)
     - A column(s) named after the `--endpoints` argument, e.g. N group:<br />
       >  0 : 0 in TNM stage<br />
          1 : 1 in TNM stage<br />
          2 : only follow up (>= 2 years), no recurrence<br />
          3 : only follow up (>= 2 years) with recurrence
           
   - **Input**:
     - json file with structure: <br />
       >patient-id:<br />
       >>filename:<br />
         folder: folder name<br />
         endpoint1: endpoint value<br />
         endpoint2: endpoint value<br />
         ...

1. `endpoints-json-to-cxl.py`: sets up json dictionary based on excel file
   - **Input**:
     - `--output-path`: where the cxl files should be saved to
     - `--json-path`: path to the json file with the endpoint(s) data
     - `--endpoint`: name of the variable (in the dict) that encodes the end-point (only 1!)
     - `--dataset-name`: name that should be given to the dataset in the xml files
     - `--split`: (optional, default 0.4) how much should be split off for train and test (will be split in half for test valid) 
     - `--seed`: (optional, default 42) set seed for split
        
     json file needs to be in the following structure (output of `make-endpoint-json.py`): <br />
       >patient-id:<br />
       >>filename:<br />
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


## Installation requirements
Install the conda environment with `conda env create -f environment.yml`. 

Additionally you need to install [ASAP](https://github.com/computationalpathologygroup/ASAP).
