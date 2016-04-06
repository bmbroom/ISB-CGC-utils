
getCloudProject <- function() {
   getOption('cloudproject')
}

#' @export
cloudServer <- function (name, intAddr, extAddr, chmPort=80, managerPort=18080, ...) {
    chmCreateServer (name,
                     sprintf ("http://%s:%d/chm", intAddr, chmPort),
                     list (serverProtocol="manager",
                           deployServer=sprintf("http://%s:%d/chm/manager/rest", intAddr, managerPort),
                           serviceName="default",
                           viewServer=sprintf("http://%s:%d/chm", extAddr, chmPort),
                           ...))
}
