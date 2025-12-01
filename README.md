# FI-radii
**Short script to automate the determination of fluid inclusion radii from images stored in a WITec datafile using WITio in MATLAB.**

This script loops through an Excel table of analyzed fluid inclusion (FI) point analysis and pulls the corresponding image files labeled with "_100x, "_50x" or "_20x". For the loop to work, the object names within the WIN file must be consistent and unique. 

This MATLAB script requires the plugin [WITio: A MATLAB data evaluation toolbox plugin for MATLAB](https://github.com/ElsevierSoftwareX/SOFTX-D-20-00088). Users should also credit: 
>Holmi, J.T., and  Lipsanen, H., 2022, WITio: A MATLAB data evaluation toolbox to script broader insights into big data from WITec microscopes: SoftwareX, v. 18, p. 101009,  doi:10.1016/j.softx.2022.101009.
