# procPharm

### General Install
Calcium Imaging Analysis
This Package allows you tot view calcium imaging data along with visualize the cells. This package has not been released on CRAN yet, nor have the python components been released on PYPI.

To install the package on your home computer follow the code below,
````
# First install devtools to gain access to the function install_github
install.packages("devtools")

# load up all the devtools
require(devtools)

# Now install the package
install_github("https://github.com/leeleavitt/procPharm.git")

````
### Mac OS setup
If you are working on a Mac, please follow these [instructions](./extras/procPharm_MacOS_setup.md) to set your computer up.

### Lab Computers
If installing R on computers in the lab (these need to be connected to the network drive see [here](./extras/Z_drive_Mounting_Information_1.docx)). From there you need to follow the instructions [here](./extras/README.R)

### Python Installation
This R package depends on some python processing. To install
1. Install the [anaconda distribution](https://www.anaconda.com/distribution/). From there you will need to use `pip` to install our software. 
2. Now opn the *Anaconda Prompt*, from your programs
3. Install git with `conda install git`
4. Install the python package with, 

````
pip install -e "git+https://github.com/leeleavitt/procPharm/#egg=pkg&subdirectory=python_packages/python_pharmer"
````
5. Computer will need a restart after this

