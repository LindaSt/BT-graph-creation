# example:
# <?xml version="1.0"?>
# <ASAP_Annotations>
# 	<Annotations>
# 		<Annotation Name="Annotation 0" Type="Rectangle" PartOfGroup="hotspot" Color="#F4FA58">
# 			<Coordinates>
# 				<Coordinate Order="0" X="21676.3555" Y="72592.5781" />
# 				<Coordinate Order="1" X="25325.4883" Y="72592.5781" />
# 				<Coordinate Order="2" X="25325.4883" Y="76241.7109" />
# 				<Coordinate Order="3" X="21676.3555" Y="76241.7109" />
# 			</Coordinates>
# 		</Annotation>
# 	</Annotations>
# 	<AnnotationGroups>
# 		<Group Name="hotspot" PartOfGroup="None" Color="#64FE2E">
# 			<Attributes />
# 		</Group>
# 	</AnnotationGroups>
# </ASAP_Annotations>

import glob
import fire
import os
import xmltodict
import pandas as pd

# file_path = "Z:\\GRP Dawson\\Hotspot Annotations\\B02.2702_ID_AE1_AE3_CD8.xml"
# output_path = "Z:\\GRP Dawson\\Hotspot Annotations Checked\\B02.2702_ID_AE1_AE3_CD8.xml"


def check_xml(file_path):
    with open(file_path) as fd:
        doc = xmltodict.parse(fd.read())

    if len(doc) == 0:
        print(f'File {file_path} is empty (skipped)!')
        return

    # check and if necessary correct the structure
    if doc["ASAP_Annotations"]["Annotations"]["Annotation"]["@Name"] != "Annotation 0":
        doc["ASAP_Annotations"]["Annotations"]["Annotation"]["@Name"] = "Annotation 0"

    if doc["ASAP_Annotations"]["Annotations"]["Annotation"]["@PartOfGroup"] != "hotspot":
        doc["ASAP_Annotations"]["Annotations"]["Annotation"]["@PartOfGroup"] = "hotspot"

    if not doc["ASAP_Annotations"]["AnnotationGroups"]:
        doc["ASAP_Annotations"]["AnnotationGroups"] = {"Group": {"@Name": "hotspot", "@PartOfGroup": "None",
                                                                 "@Color": "#64FE2E", "@Attributes": None}}
    xml = xmltodict.unparse(doc, pretty=True)
    return xml


def parse_matched_files_excel(matched_files_excel):
    df = pd.read_excel(matched_files_excel, sheet_name='Masterfile', engine='openpyxl')
    # drop all rows that do not contain 0 or 1 in column "Need resection?" (excluded because no data available)
    # df = df.drop(df[~df["Need resection?"].isin([0, 1])].index)
    # drop all rows that do not contain a file name
    # TODO: make this neater
    df = df[df['Hotspot filename'].notna()]
    df = df[df['Algo coordinates text file ID'].notna()]
    df = df.drop(df[df['Hotspot filename'].isin(["tbd", "na"])].index)
    df = df.drop(df[df['Algo coordinates text file ID'].isin(["tbd", "na"])].index)
    # drop all the other columns
    files_to_process = df[['Hotspot filename', 'Algo coordinates text file ID']]
    return list(files_to_process.itertuples(index=False, name=None))


def check_hotspots(input_path: str, output_path: str, overwrite: bool = False, matched_files_excel: str = None):
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    if os.path.isdir(input_path):
        file_list = [(i, os.path.basename(i).split('.')[0]) for i in glob.glob(os.path.join(input_path, "*.xml"))]
    else:
        file_list = [(input_path, os.path.basename(input_path).split('.')[0])]

    if matched_files_excel:
        file_list = [(os.path.join(input_path, hotspot_name), outfile_name) for hotspot_name, outfile_name in parse_matched_files_excel(matched_files_excel)]

    for hotspot_file_path, outfile_name in file_list:
        try:
            filename = os.path.basename(hotspot_file_path)
            if overwrite and os.path.isfile(os.path.join(output_path, filename)):
                print(f'File {filename} already exists (skipping file).')
            else:
                xml = check_xml(hotspot_file_path)
                with open(os.path.join(output_path, f"{outfile_name}.xml"), "w") as text_file:
                    text_file.write(xml)

        except FileNotFoundError:
            print(f'File {hotspot_file_path} not found (skipping file).')


if __name__ == '__main__':
    """
    INPUT:
    command line arguments:
    --input-path: path to the xml hotspot file/folder
    --output-path: path to the folder, where the corrected xml files should be saved to
    --overwrite: overwrites existing xml files (default is False)
    --matched-files-excel: excel file which matches up the hotspot names with the text file id
    """
    fire.Fire(check_hotspots)

