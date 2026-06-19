CytoSOM: a package for easy use of FlowSOM
==========================================
Introduction
==============
CytoSOM helps to apply [FlowSOM](https://github.com/SofieVG/FlowSOM) to flow cytometry data, and perform statistical analysis. You can use the [userguide](https://github.com/gautierstoll/CytoSOM/blob/master/CytoSOM%20User%20Guide.pdf) or read the "HowTo" below.

## Contact
cytosom.package@gmail.com

## Installation

Some dependent packages may not be available directly. They have to be installed first, within bioconductor: [flowCore](https://www.bioconductor.org/packages/release/bioc/html/flowCore.html), [CytoML](https://www.bioconductor.org/packages/release/bioc/html/CytoML.html), [flowWorkspace](https://www.bioconductor.org/packages/release/bioc/html/flowWorkspace.html), [FlowSOM](https://bioconductor.org/packages/release/bioc/html/FlowSOM.html).

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("flowCore","flowWorkspace","CytoML","FlowSOM"))
```

Then, you can install CytoSOM:

```R
devtools::install_github(repo ="gautierstoll/CytoSOM")
```

## HowTo

### Using gated data from FlowJo

`CytoSOM` is applied when data has been gated within FlowJo software. A gating can be constructed manually, see below. The data can also be extracted from previously contructed meta-clusters, see below. 

`CytoSOM` needs to be installed. In a folder (eg `MyData`), put your FlowJo sesssion file (eg `file.wsp`) and a sub-folder (eg `FCSdata`) that contains the `.fcs` files (eg `Tube`...`.fcs`)

1. Launch `RStudio`

2. Set your working directory as the folder where your FlowJo session and your data are stored, using the full folder name:
```R
setwd("full_name_before_myData/MyData")
```
Alternatively, it can be done directly in `RStudio`: Session -> Set Working Directory -> Choose directory

3. Create an R object that collects your data (eg `CytoData`), indicating the gate (i.e. cell population) of interest drawn in FlowJo (eg `CD45`):
```R
CytoData <- CytoSOM::DownLoadCytoData(dirFCS="FCSdata","CD45",fcsPattern = "Tube",compensate=FALSE)
```
Note in this case, data have already been compensated during acquisition on the flow cytometer.

4. Build a clustering tree (eg `CytoTree`), indicating the list of makers used for clustering (eg `c("CD3","CD4","CD8","CD11b","FOXP3","CD19")`), the size of the cluster grid (eg `7`), the number of meta-clusters (eg `25`), and the seed of the random generator (eg `0`):
```R
CytoTree <- CytoSOM::buildFSOMTree(CytoData,c("CD3","CD4","CD8","CD11b","FOXP3","CD19"),7,25,0)
```
Note that `buildFSOMTree` need to have the compensated names of the markers, that can be extracted from the object `CytoData` by the commande `as.vector(gsub(" <.*>","",CytoData$fSOMData$prettyColnames))`

5. Plot the tree
```R
CytoSOM::PlotStarsMSTRm(CytoTree$fSOMTree,CytoTree$metaCl,"Title Name",0)
```
If the smallest clusters make the image diffcult to interpret, you can plot the tree by removing the smallest clusters (eg remove the 3 smallest clusters):
```R
CytoSOM::PlotStarsMSTRm(CytoTree$fSOMTree,CytoTree$metaCl,"Title Name",3)
```
It is possible to represent the clustering in term of a heatmap, meta-cluster vs median of marker, eg:
```R
CytoSOM::HeatMapTree(CytoTree,c("CD3","CD4","CD8","CD11b","FOXP3","CD19"))
```
It is also possible to examine the histogram of a given meta-cluster for a given marker (meta-cluster in red compared the whole data in black), eg (marker CD3 for meta-cluster 8):
```R
CytoSOM::HistMarkerMetaCl(CytoTree,"CD3","8")
```

If the tree looks satisfactory, you can move on to the next step. Otherwise, try to rebuild the tree with different parameters (ie size of the cluster grid and/or number of meta-clusters)

6. Rename the meta-clusters, switching from a meta-cluster number to an explicit phenotype. This can be done manually by providing 
a data frame with a column "oldName" and a column "newName" (names in column "newName" must be all different); this renaming is partial if not all meta-cluster names are in "oldName" column. For instance, for renaming meta-cluster 2 and 3:
```R
renameDF <- data.frame(oldName = c(2,3),newName = c("Type_two","Type_three"))
CytoTreeMRn <- CytoSOM::TreeMetaManualRenaming(CytoTree,renameDF)
```
This can be done automatically, by providing the set of markers and the number of levels used for each markers (from 2 to 4), eg:
```R
CytoTreeRn <- CytoSOM::TreeMetaRenaming(CytoTree,c("CD4","CD8","CD11b","FOXP3", "CD19"),4)
```

7. Download an annotation table in `.csv` format (eg `AnnotationTable.csv`), that contains a column 'files', a column 'Treatment' and an optional column 'NormalizationFactor', indicating the separator (eg `;`):
```R
tableTreatmentFCS <- read.csv("AnnotationTable.csv",sep=";",dec=".")
```

8. Plot trees in `.pdf` files, for different population sizes and different markers (eg `c("MHCII","PD1","PDL1","PDL2")`):
```R
CytoSOM::plotTreeSet(CytoTreeRn ,c("MHCII","PD1","PDL1","PDL2"),"TitleTrees",rmClNb=0,tableTreatmentFCS,globalScale=T)
```
This command generates two `.pdf` files for visual comparison of the clusters (i.e. cell subsets) based on their size and on the expression profile of markers to define (eg `c("MHCII","PD1","PDL1","PDL2")`): `TitleTrees_ClusterTree.pdf` and `TitleTrees_MarkerTree.pdf`.

9. Perform statistical analysis of meta-cluster sizes:
```R
StatAnalysisSizes <- CytoSOM::BoxPlotMetaClust(CytoTreeRn,Title="MyTitle",tableTreatmentFCS,ControlTreatmen="PBS",
BottomMargin=3,yLab="CD45",Norm=FALSE,Robust = TRUE,ClustHeat=TRUE)
```
If `Norm` is set to `TRUE`, the column 'NormalizationFactor' of the `.csv` table is used to normalize the meta-cluster sizes (be sure that this column provide real numbers). Otherwise, the analysis is perfomed on relative size (percentage). The control treatment is used for statistical annotation of population size heatmap. A file `MyTitle_BoxPlotPercentMetaClust.pdf` is produced. Two files containing p-values are also generated: `MyTitle_LmPvalPercentMetacl.csv` and `MyTitle_PairwisePvalPercentMetacl.csv`

10. Perform statistical analysis of a given marker (eg PD1) MFI, across the different meta-clusters:
```R
StatAnalysisPD1 <- CytoSOM::BoxPlotMarkerMetaClust(CytoTreeRn,Title="MyTitle",tableTreatmentFCS,ControlTreatmen="PBS",
BottomMargin=3,"PD1",Robust = TRUE,ClustHeat=TRUE)
```
A file `MyTitle_BoxPlotPD1MetaClust.pdf` is produced. Two files containing p-values are also generated: `MyTitle_LmPvalPD1Metacl.csv` and `MyTitle_PairwisePvalPD1Metacl.csv`

### Gating without FlowJo

`CytoSOM` needs to be installed. In a folder (eg `MyData`), put a sub-folder (eg `FCSdata`) that contains the `.fcs` files (eg `Tube`...`.fcs`)

1. Launch `RStudio`

2. Set your working directory as the folder where your data are stored, using the full folder name:
```R
setwd("full_name_before_myData/MyData")
```
Alternatively, it can be done directly in `RStudio`: Session -> Set Working Directory -> Choose directory

3. Download the data:

```R
RawData <- FlowSOM::ReadInput(input = "FCSdata",pattern = "Tube",compensate = F)
```
Note: in this case, data have already been compensated during acquisition on the flow cytometer.

4. Create two polygon gates, eg the first one within the 2D plot "FSC-A" x "SSC-A" using `.fcs` files 1 and 3, the second one within the 2D plot "FSC-A" x "Livedead" (from 0 to 10000) using `.fcs` files 2 and 4 (order of files are the order that can be seen in `RawData$metaData`):
```R
Poly1 <- CytoSOM::InteractivePolyGate(RawData,marker1 = "FSC-A",marker2 = "SSC-A",fcsFiles = c(1,3))
Poly2 <- CytoSOM::InteractivePolyGate(RawData,marker1 = "FSC-A",marker2 = "Livedead",fcsFiles = c(2,4),ylim=c(0,10000))
```
Any list of files can be used to contruct a polygon gate. The order of files are the order that can be seen in `RawData$metaData`. Instead, file names can be provided.


5. Create an R object (dataset) with the intersection of the two polygon gates above, with a chosen name (eg "myGate"):
```R
CytoData <- CytoSOM::PolygonGatingRawData(RawData,Polygons = list(Poly1,Poly2),gatingName = "myGate”)
```

6. A new polygon can be applied to this gated dataset:

```R
Poly3 <- CytoSOM::InteractivePolyGate(CytoData$fSOMData,marker1 = "SSC-H",marker2 = "Livedead",fcsFiles = c(8,10))
CytoData <- CytoSOM::PolygonGatingGatedData(CytoData,Polygons = list(Poly3),gatingName = "CD45”)
```

Then the analysis can be continued at point 4 above ("Using gated data from FlowJo").

### Extracting data from meta-clusters

Suppose that data has been downloaded (eg `CytoData`), and a cluster tree has been constructed (eg `CytoTree`). A sub-dataset can be extracted from a list of metaclusters (eg `c(1,3)`):
```R
SubCytoData <- CytoSOM::DataFromMetaClust(CytoData$fSOMData,CytoTree,c(1,3))
```
Then a new tree can be constructed, continuing at point 4 of HowTo above ("Using gated data from FlowJo").

If meta-clusters have been renamed, the list of exact names should be provided. For that, the function `FindMetaClustNames` can be useful. For instance, 
```R
clusterNames <- FindMetaClustNames(subNames = c("CD4high","FOXP3-"),TreeMetaClust = CytoTree)
``` 
finds all meta-clusters having "CD4+" and "FOXP3-" in their name. In a similar way, if meta-clusters have been renamed by `TreeMetaRenaming`, 
```R
clusterNames <- FindMetaClustNames(subNames = "4_",TreeMetaClust = CytoTree,start=T)
```
finds the name of the meta-cluster number 4.
