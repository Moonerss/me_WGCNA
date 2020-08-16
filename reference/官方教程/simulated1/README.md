## I. Network analysis of liver expression data from female mice: finding modules related to body weight

### Data description and download

This tutorial guides the reader through the analysis of an empirical data set. The data are gene expression measurements from livers of female mouse of a specific F2 intercross. For a detailed description of the data and the biological implications we refer the reader to Ghazalpour et al (2006), *Integrating Genetics and Network Analysis to Characterize Genes Related to Mouse Weight* ([link to paper](http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.0020130); [link to additional information](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/MouseWeight/)). We note that the data set contains 3600 measured expression profiles. These were filtered from the original over 20,000 profiles by keeping only the most variant and most connected probes. In addition to the expression data, several physiological quantitative traits were measured for the mice. Please download the following

- [zipped data sets,](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-Data.zip)

and unzip them in a folder of your choice, preferably a new folder created specifically for this tutorial. Note the name of the folder; when you start an R session, the first command should be to change the R working directory into this folder.

### R Tutorial

The flowchart of the tutorial is shown below.

![img](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/overview-femaleLiver.png)


Individual sections can be viewed in PDF format by clicking on the links below. Plain text R code from each section is also available by clicking on the corresponding R script link. For first-time users we recommend starting at the top of the list and working down. Each section of the tutorial saves results on disk and the results needed as input for the subsequent sections can be loaded from disk, so repeated execution of any of the sections does not require re-working previous sections again.

1. Data input and cleaning: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.R)
2. Network construction and module detection
   1. Automatic, one-step network construction and module detection: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.R)
   2. Step-by-step network construction and module detection: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.R)
   3. Dealing with large datasets: block-wise network construction and module detection: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-blockwise.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-blockwise.R)
3. Relating modules to external clinical traits and identifying important genes: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.R)
4. Interfacing network analysis with other data such as functional annotation and gene ontology [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-04-Interfacing.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-04-Interfacing.R)
5. Network visualization using WGCNA functions: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-05-Visualization.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-05-Visualization.R)
6. Export of networks to external software: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-06-ExportNetwork.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-06-ExportNetwork.R)