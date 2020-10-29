# The location of this file is in,
".extras/labDriveInstall/README.R"


# First copy the Rprofile.site file to the desktop
# Now run the following code

# Gives me the location of the current R
homeDir <- Sys.getenv("R_HOME")

# Create a file.path for the etc to replace the rprofile.site
etcPath <- file.path(homeDir, 'etc')

setwd("..")
file.copy("./Desktop/Rprofile.site", etcPath, TRUE, TRUE)

