
getCloudProject <- function() {
   getOption('cloudproject')
}

#' @export
cloudServer <- function (name, intAddr, extAddr, ...) {
    chmCreateServer (name,
                     sprintf ("http://%s:8080/chm", intAddr),
                     list (serverProtocol="manager",
                           deployServer=sprintf("http://%s:18080/chm/manager/rest", intAddr),
                           serviceName="default",
                           viewServer=sprintf("http://%s:8080/chm", extAddr),
                           ...))
}
