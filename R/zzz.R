#check if package mineCETSA has been downloaded
.onAttach <- function(libname, pkgname) {

  if (!("mineCETSA" %in% rownames(installed.packages()))) {
    packageStartupMessage(
      paste0(
        "Please install the last version of `mineCETSA` by",
        " `devtools::install_github('nkdailingyun/mineCETSA')`",
        " or by procuring you the tar.gz file."
      )
    )
  }
  if (!("STRINGdb" %in% rownames(installed.packages()))) {
    packageStartupMessage(
      paste0(
        "Please install `STRINGdb` package by",
        " `BiocManager::install('STRINGdb')`"
      )
    )
  }
  if (!("EBImage" %in% rownames(installed.packages()))) {
    packageStartupMessage(
      paste0(
        "Please install `EBImage` package by",
        " `BiocManager::install('EBImage')`"
      )
    )
  }

}

#save the working directory of the user in order to save file when using the app
.onLoad <- function(...){
  WD <<- getwd()
  message(
    paste0("Your working directory has been save to the variable WD as ",
           WD, " Please don't change this variable unless you want to change
           your saving directory when using the app.")
  )
}
