
# Mac OS setup
If you are working on a Mac. First follow these [instructions](./extras/procPharm_MacOS_setup.md) to set up your computer. Come back to this page when you have completed it.


# General install
This Package allows you to view calcium imaging data along with visualize the cells. This package has not been released on CRAN yet, nor have the python components been released on PYPI.

To install the package on your home computer. First open `R`. If R is not installed make sure to install now [here](https://cran.r-project.org/bin/windows/base/old/3.5.3/). While install select **custom installation** and select the **multiple windows** option.
````
# First install devtools to gain access to the function 
install.packages(c("devtools", "backports", "fs"))
````
On MAC, make sure to answer no, if prompted.

````
# Now install the package
devtools::install_github("https://github.com/leeleavitt/procPharm.git")

````

Make Sure to read the console if any ERROR's arise. If so Close R and paste the commands again.

# Lab Computers
If installing R on computers in the lab (these need to be connected to the network drive see [here](./extras/Z_drive_Mounting_Information_1.docx)). From there you need to follow the instructions [here]("../labDriveInstall/README.R")

# Python Installation
**Skip this if;**
1. You don't need to process the raw files 
2. If you are installing this software to look at traces.

This R package depends on some python processing. To install
1. Install the [anaconda distribution](https://www.anaconda.com/distribution/). Make sure to select for a single user during this installation. From there you will need to use `pip` to install our software. 
2. Now open the *Anaconda Prompt*, from your programs. This means click the windows logo and type `anaconda prompt`. Double click this and a terminal window will open.
3. Install git with 
````
conda install git
````
4. Install the python package with, 

````
pip install -e "git+https://github.com/leeleavitt/procPharm/#egg=pkg&subdirectory=python_packages/python_pharmer"
````

5. 
````
conda install scikit-image
````
6. Computer will need a restart after this

