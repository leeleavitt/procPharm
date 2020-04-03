# The location of this file is in,
"Z:/Computer setup/R"

# This is how we install our software!
# Where this should be located is in
# "Z:/Computer setup/R"
# First double click "R-3.5.3-win.exe"
# MAKE SURE TO SELECT CUSTOM INSTALL
# SELECT SEPERATE WINDOWS

# First copy the Rprofile.site file to the desktop
# Now run the following code

# Gives me the location of the current R
homeDir <- Sys.getenv("R_HOME")

# Create a file.path for the etc to replace the rprofile.site
etcPath <- file.path(homeDir, 'etc')

file.copy("~/Desktop/Rprofile.site", etcPath, TRUE, TRUE)

