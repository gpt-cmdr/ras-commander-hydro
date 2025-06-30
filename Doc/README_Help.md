# RAS Commander Help Documentation

This directory contains help documentation for the RAS Commander Toolbox.

## Files

- **RASCommander_Help.html** - Main help documentation file that displays when users click the help button in ArcGIS Pro

## Help Integration

The help system is integrated into each tool through the `getHelp()` method in each tool class. When users click the help button (?) in the tool dialog:

1. ArcGIS Pro calls the tool's `getHelp()` method
2. The method returns a file:/// URL pointing to the HTML help file
3. The URL includes an anchor (#) to jump to the specific tool's section

## Updating Help

To update the help documentation:

1. Edit `RASCommander_Help.html` 
2. Make sure to maintain the anchor IDs for each section:
   - `#load-hec-ras-2d-geometry-layers`
   - `#load-hec-ras-2d-results-summary-layers`
   - `#load-hec-ras-terrain`

## VRT Limitations Note

The Load HEC-RAS Terrain tool includes important warnings about VRT limitations:
- The tool only loads base terrain VRT files
- Vector terrain modifications from RAS Mapper are NOT included
- This limitation is documented in:
  - The tool's class docstring
  - The tool's description
  - Warning messages during execution
  - The help documentation

## Styling

The HTML help file uses inline CSS for portability. The styling includes:
- Responsive design that works in various browser windows
- Color-coded sections for warnings and notes
- Clear parameter descriptions
- Professional appearance matching ArcGIS Pro's style