library(plumber)

pr("routers.R") %>%
    pr_run(port=8000)