## III. Analysis of simulated data

In this R software tutorial we review key concepts of weighted gene co-expression network analysis (WGCNA). The tutorial also serves as a small introduction to clustering procedures in R. We use simulated gene expression data to evaluate different module detection methods and gene screening approaches.

### Data description and download

Although the tutorial uses simulated data, in Section 2 it also demonstrates loading of summary, expression and trait data. Data files used for that section are generated in Section 7; however, we also provide them here for download:

- [Simulated data](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/SimulatedData.zip)

Please download and unzip them in a folder of your choice. We recommend a folder separate from the mouse analyses above, but same folder will work as well. Note the name of the folder; when you start an R session, the first command should be to change the R working directory into this folder.

### R Tutorial

Individual sections of the tutorial can be viewed in PDF format by clicking on the links below. We recommend starting at the top working through the sections in the order they are presented here. Each section saves its results on disk and the results needed as input for the subsequent parts can be loaded from disk, so repeated execution of any of the sections does not require re-working previous sections again.

1. Simulation of expression and trait data: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-01-dataSimulation.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-01-dataSimulation.R)
2. Loading of expression data, an alternative to data simulation, provided to illustrate data loading of real data: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-02-dataLoading.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-02-dataLoading.R)
3. Basic data preprocessing illustrates rudimentary techniques for handling missing data and removing outliers: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-03-Preprocessing.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-03-Preprocessing.R)
4. Standard gene screening illustrates gene selection based on Pearson correlation and shows that the results are not satisfactory: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-04-StandardScreening.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-04-StandardScreening.R)
5. Construction of a weighted gene co-expression network and network modules illustrated step-by-step; includes a discussion of alternate clustering techniques: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-05-NetworkConstruction.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-05-NetworkConstruction.R)
6. Relating modules and module eigengenes to external data illustrates methods for relating modules to external microarray sample traits: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-06-RelatingToExt.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-06-RelatingToExt.R)
7. Module membership, intramodular connectivity, and screening for intramodular hub genes illustrates using the intramodular connectivity to define measures of module membership and to screen for genes based on network information: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-07-Membership.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-07-Membership.R)
8. Visualization of gene networks: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-08-Visualization.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-08-Visualization.R)  

