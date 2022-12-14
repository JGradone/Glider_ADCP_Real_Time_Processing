# Glider_ADCP_Real_Time_Processing

Glider_ADCP_Real_Time_Processing is a package designed for processing Teledyne RDI Pathfinder data collected in real-time. The initial development of this code was for Teledyne Webb Research Slocum gliders, however it is intended to be platform and (eventually) sensor agnostic. The Python script developed here runs on a RaspberryPi and does the following:

1) Read a Teledyne RDI Pathfinder .PD0 file
2) QAQC the ADCP data
3) Perform a least squares linear inversion to extract a horizontal velocity profile
4) Save output as a NetCDF

This script can be run on a user's local machine to perform the above steps as well.

Acknowledgments
----------------------
Please cite this package using this [DOI]. The intention for publishing this code on GitHub is twofold: 1) to make glider based acoustic current profiler data easier to work with and 2) to improve upon the processing workflow and make it more transparent. For these reasons, I encourage comments, questions, concerns, and for users of this code to report any potential bugs. Please do not hesitate to reach out to me at jgradone@marine.rutgers.edu.

The PI of this project is Dr. Travis Miles. Joe Gradone, Eli Hunter, and Julia Engdahl contributed to the code development.

Citations for the methods used in this package are included in the code where relevant.
