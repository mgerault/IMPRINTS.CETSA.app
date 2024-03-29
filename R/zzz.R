.onAttach <- function(libname, pkgname) {
  if (!("IMPRINTS.CETSA" %in% installed.packages())) {
    packageStartupMessage(
      "Please install the last version of `IMPRINTS.CETSA` by",
      " `devtools::install_github('nkdailingyun/IMPRINTS.CETSA')`",
      " or by procuring you the corresponding tar.gz file."
    )
  }
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    packageStartupMessage(
        "Please install `STRINGdb` package by",
        " `BiocManager::install('STRINGdb')`"
    )
  }

  packageStartupMessage(
    "\n",
    "Welcome to IMPRINTS.CETSA.app package! To launch the app, run runCETSAapp() function.\n",
    "To access the documentation type browseVignettes(package = 'IMPRINTS.CETSA.app') or \n",
    "with vignette('IMPRINTS.CETSA.app_doc', package = 'IMPRINTS.CETSA.app')"
    )

  WD <<- getwd()
  packageStartupMessage(
    "\n",
    "A variable named 'WD' has been created to make it easier to browse files\n",
    "when using the app. You can change its value but it must be a valid file path."
  )
}
