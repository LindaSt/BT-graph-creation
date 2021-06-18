# correct Heather's wrong hotspot xml structure
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

import os
import numpy as np
from lxml import etree as ET
from collections import defaultdict
import xmltodict

file_path = "Z:\\GRP Dawson\\Hotspot Annotations\\B02.2702_ID_AE1_AE3_CD8.xml"
output_path = "Z:\\GRP Dawson\\Hotspot Annotations Checked\\B02.2702_ID_AE1_AE3_CD8.xml"

tree = ET.parse(file_path)
root = tree.getroot()
filename = os.path.basename(os.path.splitext(file_path)[0])

with open(file_path) as fd:
    doc = xmltodict.parse(fd.read())


if doc["ASAP_Annotations"]["Annotations"]["Annotation"]["@Name"] != "Annotation 0":
    doc["ASAP_Annotations"]["Annotations"]["Annotation"]["@Name"] = "Annotation 0"

if doc["ASAP_Annotations"]["Annotations"]["Annotation"]["@PartOfGroup"] != "hotspot":
    doc["ASAP_Annotations"]["Annotations"]["Annotation"]["@PartOfGroup"] = "hotspot"

if not doc["ASAP_Annotations"]["AnnotationGroups"]:
    doc["ASAP_Annotations"]["AnnotationGroups"] = {"Group": {"@Name": "hotspot", "@PartOfGroup": "None",
                                                             "@Color": "#64FE2E", "@Attributes": None}}

#%%

xml = xmltodict.unparse(doc, pretty=True)
print(xml)
#%%
with open(output_path, "w") as text_file:
    text_file.write(xml)