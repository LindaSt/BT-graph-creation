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


def check_hotspots(input_path: str, output_path: str, overwrite: bool = False):
    if os.path.isdir(input_path):
        file_list = glob.glob(os.path.join(input_path, "*.xml"))
    else:
        file_list = [input_path]

    for file_path in file_list:
        filename = os.path.basename(file_path)
        if overwrite and os.path.isfile(os.path.join(output_path, filename)):
            print(f'File {filename} already exists (skipping file).')
        else:
            # print(file_path)
            xml = check_xml(file_path)
            with open(os.path.join(output_path, filename), "w") as text_file:
                text_file.write(xml)


if __name__ == '__main__':
    """
    INPUT:
    command line arguments:
    --input-path: path to the xml hotspot file/folder
    --output-path: path to the folder, where the corrected xml files should be saved to
    --overwrite: overwrites existing xml files (default is False)
    """
    fire.Fire(check_hotspots)

