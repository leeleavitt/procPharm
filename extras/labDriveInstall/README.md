This is how we install our software for somputers that want access to the software on the network drives.

Open up R as an administrator. To do this right click the R icon and select run as administrator.

Now copy this [file](./Rprofile.site) to the desktop of your computer.

Copy this series of code to the R Console.
````
# Gives me the location of the current R
homeDir <- Sys.getenv("R_HOME")

# Create a file.path for the etc to replace the rprofile.site
etcPath <- file.path(homeDir, 'etc')

setwd("..")
file.copy("./Desktop/Rprofile.site", etcPath, TRUE, TRUE)
````

If you see `TRUE` written out, you have succeeded. if not review the instructions and contact Lee.
