Individual scripts for computing pseudotime trajectories for each dataset. FTP downloads are included. You need `Seurat`, `SeuratDisk`, `SingleCellExperiment`, `tidyverse`, and `monocle3` to run the files.

Monocle3 is pretty tricky to install. We were unable to install it on a RHEL/CentOS 7 server, but can get it installed on MacOS and RHEL/CentOS8. 

Each dataset has an individual script for runnig because there are unique paramters for computing the trajecory for each dataset. 


After running all the individual dataset scripts you run the extract_MSTs.r script to get the relevant summary stastitics needed for running downstream analysis in python

