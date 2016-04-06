
.onAttach <- function(libname, pkgname) {
   proj <- getCloudProject();
   if (length(proj)==0) {
       cat ("Option cloudproject is not set\n");
   } else {
       cat ("Option cloudproject is ", proj, "\n");
   }
}

