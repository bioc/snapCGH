.First.lib <-
    function(libname, pkgname)
    library.dynam("snapCGH", pkgname, libname)

.Last.lib <-
    function(libpath)
    dyn.unload(file.path(libpath,
                         "libs",
                         paste("snapCGH",
                               .Platform$"dynlib.ext",
                               sep = "")))
