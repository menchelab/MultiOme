# Network analysis reveals rare disease signatures across multiple levels of biological organization

Supplementary codes, reproducible walk-through reports, and Shiny app complementing Buphamalai. et.al.


## This repository
The structue of the project Github repo is as described below:
```
menchelab/Multiome/
├── Explorer ----> Shiny application for disease modularity inspection
├── data --------> Data used in the analyses
├── functions ---> Common functions called in main analyses
├── report ------> Reproducible, walkthrough, R-markdown powered report
├── source ------> Main analysis files
├── cache -------> Pre-computed results for heavier tasks. Must be downloaded (see link below) 
└── .gitignore
```

## The Explorer
The Explorer is a Shiny app made for exploring results included in the manuscript in details, and can be used as additional resources for investigating gene connectivity of a disease of interest.
![Overview](https://github.com/menchelab/MultiOme/blob/main/Explorer/Figs/github/github_readme_overview1.png?raw=true)

### Differential Modularity
![Differential Modularity](https://github.com/menchelab/MultiOme/blob/main/Explorer/Figs/github/github_readme_DiffMod.png?raw=true)


### Disease-Network Landscape
![Disease-Network Landscape](https://github.com/menchelab/MultiOme/blob/main/Explorer/Figs/github/github_readme_Landscape.png?raw=true)


### Detailed Inspection
![Detailed Inspection](https://github.com/menchelab/MultiOme/blob/main/Explorer/Figs/github/github_readme_Inspection.png?raw=true)

The Explorer can be launched via [this link](http://lem.westeurope.cloudapp.azure.com:8880/app/MultiOmeExplorer)

## The Repoducible Walkthrough Guides

[This supplementary report](/report/Supplementary_report.html) is aimed to be a reproducible walk-through guide for figures and analyses complementing the manuscript: Buphamalai et.al., *Network analysis reveals rare disease signatures across multiple levels of biological organization*, submitted to Nature Communications. (Almost) all of the figures and analyses can be reproduced via walk-through reports organized in the same order as the manuscript in RMarkdown format, and can be found in 

PLease find the following guidelines:

1. The appearance of sections in this document is at the same order as in the manuscript, and can be navigated using the Table of Content (ToC) appeared on the top left corner of this document, with subsections corresponding to exact figures/statistics 
2. The corresponding code chunks to each figures/analyses are provided above the output, and are hidden by default. To expand each code chunk, click `Show code`.
3. This report mainly contains visualization and post-processing of major analyses. Heavier computations were pre-computed, and corresponding `R`, `sh`, or `py` scripts required for each analysis are mentioned for each section. These files are available in [/source](/source) folder.
4. The pre-computed results are saved in `./cache` folder and can be downloaded from: [link](https://drive.google.com/file/d/1T7tJMojIbELeT-aLOD_Pv639eUgdGqUh/view). Unzip the folder under the main directory (`./cache`). 
5. The corresponding `Rmd` files used to produce this report can be found in [/report](/report) folder.

