# procPharm
Calcium Imaging Analysis
This Package allows you tot view calcium imaging data along with visualize the cells. This package has not been released on CRAN yet, nor have the python components been released on PYPI.

To install the package on your home computer follow the code below,
````
# First install devtools to gain access to the fucntion install_github
install.packages("devtools")

# load up all the devtools
require(devtools)

# Now install the package
install_github("https://github.com/leeleavitt/procPharm.git")
````

This R package depends on some python processing. After installing the package you need install the python packages in order for the data processing to work. We recommend you install [anaconda distribution](https://www.anaconda.com/distribution/). From there you will need to use `pip` to install our software. 

````
pip install PATH_TO_WHERE_PROCPHARM_IS\procPharm\python_packages\python_pharmer
````
`PATH_TO_WHERE_PROCPHARM_IS` is most likely located in `C:/Users/username/Documents/R/win-library/3.6` To find out exactly where r is install your packages use,
````
.libpaths()
````
For example this is how I installed the package
````
pip install C:\Users\leele\Documents\R\win-library\3.6\python_packages\python_pharmer\
````