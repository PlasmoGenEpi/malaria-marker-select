if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop(
    "Package 'rsconnect' is required to deploy.\n",
    "Install it with: install.packages('rsconnect')"
  )
}

app_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

args <- list(
  appDir = app_dir,
  appName = "malaria-marker-select",
  forceUpdate = TRUE
)

account <- Sys.getenv("SHINYAPPS_ACCOUNT", unset = "")
if (nzchar(account)) {
  args$account <- account
}

# If the user hasn't run setAccountInfo() yet, rsconnect will error with a clear message.
do.call(rsconnect::deployApp, args)

