To assist with packaging of this toolbox with Arc Hydro Tools, the scripts will be located in /Scripts/archydro

All .py scripts will be prefixed by "rc_" (for ras commander)

Due to the fact that __init__.py already exists in the Arc Hydro folder, we will need to manage our toolbox without using __init__.py, and be aware that in production a different __init__.py exists that we cannot control.  

The Arc Hydro Team prefers that the RAS Commander Arc Hydro Tools conform to the existing flat file structure under /Scripts/archydro.  No other files exist in this folder with an rc_ prefix.  use lower case file names 




Folder Mapping for ArcHydro Production: 

/Images/ras-commander-archydro-revised.png
This is used in the Help file and documentation

/Scripts/archydro 
All supporting .py files, organized to fit into the existing flat archydro file structure.  All .py files should start with rc_ prefix, use lower case file names and should not use __init__.py, as one already exists in this folder.  All imports should be relative, as the repository mirrors the final production folder structure.  

/Templates/Layers/archydro
Any layer lyrx files in this folder should be named lower case, with an "rc_" prefix.  I need a full guide on how to incorporate the use the .lyrx files into my toolbox's operations.  I understand I may need to calculate the max/min of a dataset to set the color ramp limits before inserting.  I need to figure out how to do this by looking up the ArcGIS Toolbox documentation and confirming how it is done elsewhere.  

/toolboxes
RAS Commander.pyt and .pyt.xml files for each tool box are located in this folder. 