# RAS Commander Arc Hydro Tools - Automated Tool Metadata Updater
#
# This script updates ArcGIS Pro tool metadata XML files from the Python source code.
#
# -----------------------------------------------------------------------------------
# How to Use:
#
# 1. Run this script from within ArcGIS Pro's Python environment:
#       python update_tool_metadata.py
#    or double-click update_metadata.bat (Windows).
#
# 2. The script will:
#    - Read tool info from Scripts/archydro/
#    - Update .pyt.xml files in toolboxes/
#
# 3. View updated metadata in ArcGIS Pro by right-clicking a tool and selecting "Item Description".
#
# For details, see METADATA_README.md.
#
# -----------------------------------------------------------------------------------
# What Gets Updated:
#   - Tool titles, descriptions, parameter docs, usage, examples, tags, credits
#
# When to Update:
#   - After adding/modifying tools or parameter help
#
# Customization:
#   - Edit tool class docstrings/parameter descriptions, then rerun this script
#
# Troubleshooting:
#   - If changes don't appear, remove/re-add toolbox in ArcGIS Pro and clear cache
#
# For more info, see METADATA_README.md in the repo.
# -----------------------------------------------------------------------------------

# -*- coding: utf-8 -*-
"""
update_tool_metadata.py

Automatically populates ArcGIS Pro tool metadata XML files from Python tool definitions.
This script reads the tool descriptions and parameter information from the Python
tool classes and updates the corresponding .pyt.xml metadata files.
"""

import os
import sys
import xml.etree.ElementTree as ET
import re
import datetime
import importlib.util
import base64

# Add the Scripts directory to path
scripts_dir = os.path.join(os.path.dirname(__file__), 'Scripts', 'archydro')
sys.path.insert(0, scripts_dir)

# Tool modules to process
TOOL_MODULES = [
    ('rc_load_hecras_1d_geometry', 'LoadHECRAS1DGeometry'),
    ('rc_load_hecras_2d_geometry', 'LoadHECRAS2DGeometry'),
    ('rc_load_hecras_2d_results', 'LoadHECRAS2DResults'),
    ('rc_load_ras_terrain', 'LoadRASTerrain'),
    ('rc_organize_ras_project', 'OrganizeRASProject')  # Fixed module name
]

# Common tags for all tools
COMMON_TAGS = [
    "HEC-RAS",
    "hydraulic modeling", 
    "RAS Commander",
    "Arc Hydro",
    "flood modeling",
    "2D modeling",
    "HDF5"
]

# Credits text
CREDITS_TEXT = """Development of the RAS Commander Arc Hydro Tools was sponsored by CLB Engineering (https://clbengineering.com/) in cooperation with ESRI.

Based on the ras-commander library: https://github.com/gpt-cmdr/ras-commander"""

def format_for_arcgis_html(text):
    """Convert text to ArcGIS-compatible HTML format."""
    # Bold text between **
    text = re.sub(r'\*\*(.*?)\*\*', r'<B>\1</B>', text)
    
    # Italic text between *
    text = re.sub(r'\*(.*?)\*', r'<I>\1</I>', text)
    
    # Code formatting between backticks
    text = re.sub(r'`(.*?)`', r'<SPAN STYLE="font-family: monospace">\1</SPAN>', text)
    
    return text

def clean_description_text(text):
    """Clean and format description text for XML with HTML-like formatting."""
    # Check if text is already formatted with HTML tags
    if '<DIV>' in text or '<P>' in text or '<UL>' in text:
        # Text is already formatted, return as-is
        return text
    
    # Remove extra whitespace and newlines
    text = re.sub(r'\n\s+', '\n', text)
    text = re.sub(r'\n{3,}', '\n\n', text)
    
    # Convert bullet points to HTML list items
    lines = text.split('\n')
    formatted_lines = []
    in_list = False
    
    for line in lines:
        line = line.strip()
        if line.startswith('•') or line.startswith('*'):
            if not in_list:
                formatted_lines.append('<DIV><UL>')
                in_list = True
            # Remove bullet and add as list item
            item_text = line[1:].strip()
            formatted_lines.append(f'<LI>{item_text}</LI>')
        else:
            if in_list and line:
                formatted_lines.append('</UL></DIV>')
                in_list = False
            if line:
                # Wrap paragraphs in DIV tags
                formatted_lines.append(f'<DIV><P>{line}</P></DIV>')
    
    # Close any open list
    if in_list:
        formatted_lines.append('</UL></DIV>')
    
    return '\n'.join(formatted_lines)

def encode_image_for_metadata(image_path):
    """Encode an image file to base64 for inclusion in metadata XML."""
    try:
        with open(image_path, 'rb') as img_file:
            img_data = img_file.read()
            img_base64 = base64.b64encode(img_data).decode('utf-8')
            return img_base64
    except Exception as e:
        print(f"Warning: Could not encode image {image_path}: {e}")
        return None

def extract_tool_info(module_name, class_name):
    """Extract tool information from Python class."""
    # Import the module
    module_path = os.path.join(scripts_dir, f"{module_name}.py")
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    # Get the tool class
    tool_class = getattr(module, class_name)
    tool_instance = tool_class()
    
    # Extract basic info
    info = {
        'label': tool_instance.label,
        'description': tool_instance.description,  # Don't clean here, will be cleaned when used
        'parameters': []
    }
    
    # Get parameter info
    params = tool_instance.getParameterInfo()
    for param in params:
        param_info = {
            'name': param.name,
            'displayName': param.displayName,
            'datatype': str(param.datatype),
            'parameterType': param.parameterType,
            'direction': param.direction,
            'description': param.description if hasattr(param, 'description') else ""
        }
        
        # Handle special datatypes
        if hasattr(param.datatype, '__iter__') and not isinstance(param.datatype, str):
            param_info['datatype'] = ' or '.join(str(dt) for dt in param.datatype)
        
        info['parameters'].append(param_info)
    
    return info

def create_usage_text(tool_info):
    """Create usage text from tool information with HTML formatting."""
    usage_parts = []
    
    # Add main description with HTML formatting
    desc_html = clean_description_text(tool_info['description'])
    usage_parts.append(desc_html)
    
    # Add parameter usage
    usage_parts.append('\n<DIV><P><B>Parameters:</B></P></DIV>')
    
    for param in tool_info['parameters']:
        param_type = "Required" if param['parameterType'] == 'Required' else "Optional"
        
        # Create parameter entry with HTML formatting
        param_html = f'<DIV><P><B>{param["displayName"]} ({param_type}):</B></P></DIV>'
        usage_parts.append(param_html)
        
        if param['description']:
            # Indent parameter description
            desc_html = f'<DIV STYLE="margin-left:20px;"><P>{param["description"]}</P></DIV>'
            usage_parts.append(desc_html)
    
    return '\n'.join(usage_parts)

def create_code_sample(class_name, tool_info):
    """Create a code sample for the tool."""
    # Extract required parameters
    required_params = [p for p in tool_info['parameters'] if p['parameterType'] == 'Required']
    
    # The code content remains the same, but we need to ensure it's properly escaped
    sample = f"""# Example: Using {tool_info['label']} in a Python script

import arcpy

# Set input parameters
"""
    
    # Add parameter examples based on the tool
    if class_name == 'LoadHECRAS1DGeometry':
        sample += """hdf_file = r"C:\\RASModels\\MyProject.g01.hdf"
geometry_elements = ["Cross Sections", "River Centerlines"]

# Run the tool
arcpy.RASCommander.LoadHECRAS1DGeometry(
    input_hdf=hdf_file,
    geometry_elements=geometry_elements,
    output_cross_sections=r"memory\\CrossSections",
    output_centerlines=r"memory\\RiverCenterlines"
)

print("1D geometry loaded successfully!")"""
    
    elif class_name == 'LoadHECRAS2DGeometry':
        sample += """hdf_file = r"C:\\RASModels\\MyProject.p01.hdf"
geometry_elements = ["Mesh Area Perimeters", "Mesh Cell Centers"]

# Run the tool
arcpy.RASCommander.LoadHECRAS2DGeometry(
    input_hdf=hdf_file,
    geometry_elements=geometry_elements,
    output_perimeters=r"C:\\Output\\Output.gdb\\MeshPerimeters",
    output_cell_points=r"C:\\Output\\Output.gdb\\CellCenters"
)

print("2D geometry loaded successfully!")"""
    
    elif class_name == 'LoadHECRAS2DResults':
        sample += """hdf_file = r"C:\\RASModels\\MyProject.p01.hdf"  # Must contain results
results_elements = ["Max WSE at Cell Centers", "Max Vel at Cell Faces"]

# Run the tool
arcpy.RASCommander.LoadHECRAS2DResults(
    input_hdf=hdf_file,
    results_elements=results_elements,
    output_max_wse=r"C:\\Output\\Output.gdb\\MaxWSE",
    output_max_face_vel=r"C:\\Output\\Output.gdb\\MaxVelocity"
)

print("Results loaded successfully!")"""
    
    elif class_name == 'LoadRASTerrain':
        sample += """prj_file = r"C:\\RASModels\\MyProject.prj"
terrains_to_load = ["Existing Terrain", "Proposed Terrain"]

# Run the tool
arcpy.RASCommander.LoadRASTerrain(
    in_ras_project=prj_file,
    import_all=False,
    terrains_to_load=terrains_to_load
)

print("Terrain layers added to map!")"""
    
    elif class_name == 'OrganizeRASProject':
        sample += """project_folder = r"C:\\RASModels\\MyProject"
output_gdb = r"C:\\Output\\MyProject_Organized.gdb"

# Run the tool - processes all plan files in the folder
arcpy.RASCommander.OrganizeRASProject(
    input_path=project_folder,
    output_gdb=output_gdb,
    include_cell_polygons=False,  # Set True for smaller models
    extract_all_results=True
)

print("Project organized successfully!")"""
    
    # Make sure to escape any special XML characters
    sample = sample.replace('&', '&amp;')
    sample = sample.replace('<', '&lt;')
    sample = sample.replace('>', '&gt;')
    
    return sample

def update_metadata_xml(xml_path, tool_info, class_name):
    """Update the metadata XML file with tool information."""
    # Parse existing XML
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    # Update or create dataIdInfo section
    data_id_info = root.find('.//dataIdInfo')
    if data_id_info is None:
        data_id_info = ET.SubElement(root, 'dataIdInfo')
    
    # Update title
    id_citation = data_id_info.find('idCitation')
    if id_citation is None:
        id_citation = ET.SubElement(data_id_info, 'idCitation')
    
    res_title = id_citation.find('resTitle')
    if res_title is None:
        res_title = ET.SubElement(id_citation, 'resTitle')
    res_title.text = tool_info['label']
    
    # Update abstract/description
    id_abs = data_id_info.find('idAbs')
    if id_abs is None:
        id_abs = ET.SubElement(data_id_info, 'idAbs')
    id_abs.text = clean_description_text(tool_info['description'])
    
    # Update purpose
    id_purp = data_id_info.find('idPurp')
    if id_purp is None:
        id_purp = ET.SubElement(data_id_info, 'idPurp')
    id_purp.text = create_usage_text(tool_info)
    
    # Update credits
    id_credit = data_id_info.find('idCredit')
    if id_credit is None:
        id_credit = ET.SubElement(data_id_info, 'idCredit')
    id_credit.text = CREDITS_TEXT
    
    # Add thumbnail image as Binary/Enclosure (for individual tool XMLs only)
    binary = root.find('Binary')
    if binary is None:
        # Binary should be added after dataIdInfo but before mdHrLv
        elements = list(root)
        insert_pos = len(elements)  # Default to end
        for i, elem in enumerate(elements):
            if elem.tag in ['mdHrLv', 'mdDateSt', 'distInfo']:
                insert_pos = i
                break
        binary = ET.Element('Binary')
        root.insert(insert_pos, binary)
    
    # Clear any existing enclosures
    for enclosure in binary.findall('Enclosure'):
        binary.remove(enclosure)
    
    # Get the image path relative to the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    image_path = os.path.join(script_dir, 'Images', 'ras-commander-archydro-revised.png')
    
    # Encode the image
    img_base64 = encode_image_for_metadata(image_path)
    
    if img_base64:
        # Create Enclosure element
        enclosure = ET.SubElement(binary, 'Enclosure')
        enclosure.set('rel', 'side-panel-help')
        
        # Create Data element with image
        data = ET.SubElement(enclosure, 'Data')
        data.set('EsriPropertyType', 'Image')
        data.set('OriginalFileName', 'ras-commander-archydro-revised.png')
        data.text = img_base64
    
    # Update search keywords
    search_keys = data_id_info.find('searchKeys')
    if search_keys is None:
        search_keys = ET.SubElement(data_id_info, 'searchKeys')
    
    # Clear existing keywords and add new ones
    for keyword in search_keys.findall('keyword'):
        search_keys.remove(keyword)
    
    # Add tool-specific keywords
    for tag in COMMON_TAGS:
        keyword = ET.SubElement(search_keys, 'keyword')
        keyword.text = tag
    
    # Add tool-specific keywords based on class
    if '1D' in class_name:
        for tag in ['1D modeling', 'cross sections', 'river geometry']:
            keyword = ET.SubElement(search_keys, 'keyword')
            keyword.text = tag
    elif '2D' in class_name:
        for tag in ['2D modeling', 'mesh', 'computational mesh']:
            keyword = ET.SubElement(search_keys, 'keyword')
            keyword.text = tag
    
    # Update tool section
    tool_elem = root.find('.//tool')
    if tool_elem is None:
        tool_elem = ET.SubElement(root, 'tool')
    
    # Update tool description
    tool_elem.set('displayname', tool_info['label'])
    tool_elem.set('name', class_name)
    tool_elem.set('toolboxalias', 'RASCommander')
    
    # Update parameters
    parameters = tool_elem.find('parameters')
    if parameters is None:
        parameters = ET.SubElement(tool_elem, 'parameters')
    
    # Clear existing parameters
    for param in parameters.findall('param'):
        parameters.remove(param)
    
    # Add updated parameters
    for param_info in tool_info['parameters']:
        param = ET.SubElement(parameters, 'param')
        param.set('name', param_info['name'])
        param.set('displayname', param_info['displayName'])
        param.set('type', param_info['parameterType'])
        param.set('direction', param_info['direction'])
        param.set('datatype', param_info['datatype'])
        
        # Create expression
        if param_info['parameterType'] == 'Required':
            param.set('expression', param_info['name'])
        else:
            param.set('expression', '{' + param_info['name'] + '}')
        
        # Add dialog reference with HTML formatting
        dialog_ref = ET.SubElement(param, 'dialogReference')
        if param_info['description']:
            dialog_ref.text = f'<DIV><P>{param_info["description"]}</P></DIV>'
        else:
            dialog_ref.text = f'<DIV><P>Parameter: {param_info["displayName"]}</P></DIV>'

        # Add python reference with HTML formatting
        python_ref = ET.SubElement(param, 'pythonReference')
        if param_info['description']:
            python_ref.text = f'<DIV><P>{param_info["description"]}</P></DIV>'
        else:
            python_ref.text = f'<DIV><P>Python parameter: {param_info["name"]}</P></DIV>'
    
    # Add summary
    summary = tool_elem.find('summary')
    if summary is None:
        summary = ET.SubElement(tool_elem, 'summary')
    summary.text = clean_description_text(tool_info['description'])
    
    # Add usage
    usage = tool_elem.find('usage')
    if usage is None:
        usage = ET.SubElement(tool_elem, 'usage')
    usage.text = create_usage_text(tool_info)
    
    # Add code example
    script_examples = tool_elem.find('scriptExamples')
    if script_examples is None:
        script_examples = ET.SubElement(tool_elem, 'scriptExamples')
    
    # Clear existing examples
    for example in script_examples.findall('scriptExample'):
        script_examples.remove(example)
    
    # Add new example
    script_example = ET.SubElement(script_examples, 'scriptExample')
    example_title = ET.SubElement(script_example, 'title')
    example_title.text = f"Using {tool_info['label']}"
    
    example_para = ET.SubElement(script_example, 'para')
    example_para.text = f"Example of using {tool_info['label']} in a Python script."
    
    example_code = ET.SubElement(script_example, 'code')
    example_code.text = create_code_sample(class_name, tool_info)
    
    # Update timestamps
    esri = root.find('.//Esri')
    if esri is not None:
        mod_date = esri.find('ModDate')
        if mod_date is None:
            mod_date = ET.SubElement(esri, 'ModDate')
        mod_date.text = datetime.datetime.now().strftime('%Y%m%d')
        
        mod_time = esri.find('ModTime')
        if mod_time is None:
            mod_time = ET.SubElement(esri, 'ModTime')
        mod_time.text = datetime.datetime.now().strftime('%H%M%S00')
    
    # Write updated XML
    tree.write(xml_path, encoding='UTF-8', xml_declaration=True)
    print(f"Updated metadata for {class_name}")

def main():
    """Main function to update all tool metadata."""
    toolboxes_dir = os.path.join(os.path.dirname(__file__), 'toolboxes')
    
    print("RAS Commander Metadata Updater")
    print("=" * 50)
    
    # Process each tool
    for module_name, class_name in TOOL_MODULES:
        try:
            print(f"\nProcessing {class_name}...")
            
            # Extract tool information
            tool_info = extract_tool_info(module_name, class_name)
            
            # Find corresponding XML file
            xml_filename = f"RAS Commander.{class_name}.pyt.xml"
            xml_path = os.path.join(toolboxes_dir, xml_filename)
            
            if not os.path.exists(xml_path):
                print(f"  Warning: XML file not found: {xml_filename}")
                continue
            
            # Update the XML
            update_metadata_xml(xml_path, tool_info, class_name)
            
        except Exception as e:
            print(f"  Error processing {class_name}: {str(e)}")
            import traceback
            traceback.print_exc()
    
    # Also update the main toolbox metadata
    print("\nUpdating main toolbox metadata...")
    main_xml_path = os.path.join(toolboxes_dir, "RAS Commander.pyt.xml")
    if os.path.exists(main_xml_path):
        try:
            update_main_toolbox_metadata(main_xml_path)
            print("Main toolbox metadata updated successfully!")
        except Exception as e:
            print(f"Error updating main toolbox metadata: {str(e)}")
    
    print("\nMetadata update complete!")

def update_main_toolbox_metadata(xml_path):
    """Update the main toolbox metadata."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    # Update toolbox info
    toolbox = root.find('.//toolbox')
    if toolbox is not None:
        toolbox.set('name', 'RAS Commander')
        toolbox.set('alias', 'RASCommander')
    
    # Update dataIdInfo
    data_id_info = root.find('.//dataIdInfo')
    if data_id_info is None:
        data_id_info = ET.SubElement(root, 'dataIdInfo')
    
    # Update title
    id_citation = data_id_info.find('idCitation')
    if id_citation is None:
        id_citation = ET.SubElement(data_id_info, 'idCitation')
    
    res_title = id_citation.find('resTitle')
    if res_title is None:
        res_title = ET.SubElement(id_citation, 'resTitle')
    res_title.text = "Arc Hydro RAS Commander Tools"
    
    # Update description
    id_abs = data_id_info.find('idAbs')
    if id_abs is None:
        id_abs = ET.SubElement(data_id_info, 'idAbs')
    id_abs.text = """<DIV><P>The RAS Commander Toolbox provides tools for loading and visualizing HEC-RAS 1D and 2D geometry, terrain, and results data from HDF5 files directly within ArcGIS Pro.</P></DIV>

<DIV><P>This toolbox is a direct port of the HDF5 data extraction logic from the ras-commander library, sponsored by CLB Engineering in cooperation with ESRI.</P></DIV>

<DIV><P>Tools included:</P></DIV>
<DIV><UL>
<LI>Load HEC-RAS 1D Geometry Layers - Extract cross sections, river centerlines, and structures</LI>
<LI>Load HEC-RAS 2D Geometry Layers - Extract mesh elements, breaklines, and boundary conditions</LI>
<LI>Load HEC-RAS 2D Results Summary Layers - Extract maximum WSE and velocity results</LI>
<LI>Load HEC-RAS Terrain - Load terrain VRT files from RAS Mapper</LI>
<LI>Organize HEC-RAS Project - Process entire projects into organized geodatabases</LI>
</UL></DIV>"""
    
    # Update credits
    id_credit = data_id_info.find('idCredit')
    if id_credit is None:
        id_credit = ET.SubElement(data_id_info, 'idCredit')
    id_credit.text = CREDITS_TEXT
    
    # Do NOT add thumbnail image for the main toolbox .pyt.xml
    # (The main toolbox XML cannot accept the image data in the same format as the individual tools.)
    # So, skip any Binary/Enclosure/image logic here.
    
    # Write updated XML
    tree.write(xml_path, encoding='UTF-8', xml_declaration=True)

if __name__ == "__main__":
    main()