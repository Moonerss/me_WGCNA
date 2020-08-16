## II. Consensus analysis of female and male liver expression data

### Data description and download

In this tutorial we illustrate a consensus network analysis on the example of two expression data sets, the female liver analyzed in Tutorial I, and a corresponding expression data set from livers of male mice. The two sets are biologically very similar, but significant differences exists as well. The consensus analysis parallels the female data analysis very closely and some sections are copied almost verbatim. We concentrate on the parts of the analysis that illustrate the idea behind a consensus analysis, and we leave out parts such as functional enrichment analysis for which the analysis and code would be exactly or nearly the same.

To run the tutorials, the following two zip bundles of data sets are necessary:

- [Female data](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-Data.zip) (same as in Part I - no need to download again)
- [Male data](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/MaleLiver-Data.zip)

Please download and unzip them in a folder of your choice, for example in the same folder as the female data (file names in the female and consensus analyses do not conflict). Note the name of the folder; when you start an R session, the first command should be to change the R working directory into this folder.

### R Tutorial

The flowchart of the tutorial is shown below.

![img](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/overview-Consensus.png)


Individual sections can be viewed in PDF format by clicking on the links below. We highly recommend that the user first works through the female expression data analysis, because it explains many of the same basic analysis techniques on a simpler example, without the additional complications of analyzing two sets at the same time. We recommend starting at the top working through the sections in the order they are presented here. Each section saves its results on disk and the results needed as input for the subsequent parts can be loaded from disk, so repeated execution of any of the sections does not require re-working previous sections again.

1. Data input and cleaning, including re-formatting the data for consensus analysis: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-DataInput.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-DataInput.R)
2. Network construction and consensus module detection
   1. Automatic, one-step network construction and consensus module detection: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-auto.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-auto.R)
   2. Step-by-step network construction and module detection, including scaling of Topological Overlap Matrices: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-man.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-man.R)
   3. Dealing with large datasets: block-wise network construction and consensus module detection, including comparing the block-wise approach to the standard single-block method: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-blockwise.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-blockwise.R)
3. Relating the consensus modules to female set-specific modules (this section requires the results of Section 2.a of the female turorial): [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-RelateToFemMods.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-RelateToFemMods.R)
4. Relating consensus module to external microarray sample traits and exporting the results of network analysis: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-RelateModsToTraits.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-RelateModsToTraits.R)
5. Studying and comparing the relationships among modules and traits between the two data sets, including the visualization of consensus eigengene networks and the results of the differential analysis: [PDF document](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-EigengeneNetworks.pdf), [R script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-EigengeneNetworks.R)  

