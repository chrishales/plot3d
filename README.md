plot3d
======

CASA task to quickly inspect a MS for RFI by plotting time vs frequency vs amplitude.

Current version: 1.0

Written by Christopher A. Hales. Correspondence regarding plot3d is always welcome.

plot3d is released under a BSD 3-Clause Licence; refer to the licence in this repository or the header of ```task_plot3d.py``` for details.

An example screenshot of the 3D-rotatable view produced by plot3d is shown below:
![screenshot](./screenshot.png)

## Installation

Download the source files into a directory containing your measurement set. Without changing directories, open CASA and type
```
os.system('buildmytasks')
```
then exit CASA. A number of files should have been produced, including ```mytasks.py```. Reopen CASA and type
```
execfile('mytasks.py')
```
and then
```
inp plot3d
```
Set some parameters and press go!

For a more permanent installation, go to the hidden directory ```.casa``` which resides in your home directory and create a file called ```init.py```. In this file, put the line
```
execfile('/<path_to_task_directory>/mytasks.py')
```

## Acknowledging the use of plot3d

plot3d is provided in the hope that it (or elements of its code) will be useful for your work. If you find that it is, I would appreciate an acknowledgement in your work, but this is not required. Happy plot3d-ing!
