.onAttach=function(libname,pkgname){
   packageStartupMessage("Loaded glmnet ", as.character(packageDescription("glmnet")[["Version"]]),"\n")
}
