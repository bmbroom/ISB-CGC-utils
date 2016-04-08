# ISB-CGC-utils
This repository provides a small R package of example functions for simplifying the creation of
[Next-Generation Clustered Heat Maps (NG-CHMs)](http://bioinformatics.mdanderson.org/main/NG-CHM:Overview)
using the ISB CGC cloud pilot.

The [rstudio-isb-ngchm-demo repository](https://github.com/bmbroom/rstudio-isb-ngchm-demo) provides
an example RStudio system extended with this package (and the
[NG-CHM R package](https://github.com/bmbroom/NGCHMR)).

## Usage
Load the NGCHM package and create a connection to the desired NG-CHM server(s):

```R
> library(NGCHM)
> chmCreateManagedServer('cloud','name-of-ngchm-instance','external-ip-of-ngchm-instance')
```
The example shows how to connect to the managed NG-CHM server created by the
[ngchm-system-one system](https://github.com/bmbroom/ngchm-system-one) in a cloud environment.

Specify the name of the google cloud project you will use and load the ISBCHM package:
```R
> options(cloudproject='name-of-your-google-cloud-project')
> library(ISBCHM)
```

Obtain data for approximately 410 cancer-related genes and 545 prostate cancer samples
from the ISB CGC cloud and create an NG-CHM:
```R
> chm <- demoCHM()
```

Display the NG-CHM in the viewer pane of your RStudio console:
```R
> plot(chm)
```

You can expand the viewer pane to a whole browser tab/window.  You can also display the
NG-CHM directly in a browser tab/window:
```R
> browseURL(chmGetURL(chm))
```
