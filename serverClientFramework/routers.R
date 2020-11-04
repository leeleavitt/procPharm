# router.R
# Identifies incomming requests and routes them to R functions.

#* @param msg The message to echo back.
#* @get /echo
function(msg="") {
    list(msg = paste0("The test message is: '", msg, "'"))
}

#* @post /import
function(req) {
    print(req)
}