# Network analysis reveals rare disease signatures across multiple levels of biological organization

Supplementary codes, repoducible walk-through reports, and Shiny app complementing Buphamalai. et.al.


## This repository
The structue of the project Github repo is as described below:
```
menchelab/Multiome/
├── Explorer ----> Shiny application for disease modularity inspection
├── data --------> Data used in the analyses
├── functions ---> Common functions called in main analyses
├── report ------> Reproducible, walkthrough, R-markdown powered report
├── source ------> Main analysis files
├── cache -------> Pre-computed results for heavier tasks. Must be download (see link below) 
└── .gitignore
```

## The Explorer
The Explorer is a Shiny app made for exploring results included in the manuscript in details, and can be used as additional resources for investigating gene connectivity of a disease of interest.
![Overview](https://github.com/menchelab/MultiOme/blob/main/Explorer/Figs/github/github_readme_overview.png?raw=true)

### Differential Modularity
![Differential Modularity](https://github.com/menchelab/MultiOme/blob/main/Explorer/Figs/github/github_readme_DiffMod.png?raw=true)


### Disease-Network Landscape
![Disease-Network Landscape](https://github.com/menchelab/MultiOme/blob/main/Explorer/Figs/github/github_readme_Landscape.png?raw=true)


### Detailed Inspection
![Detailed Inspection](https://github.com/menchelab/MultiOme/blob/main/Explorer/Figs/github/github_readme_Inspection.png?raw=true)

## The Repoducible Walkthrough Guides
(Almost) all of the figures and analyses can be reproduced via walk-through reports organized in the same order as the manuscript in RMarkdown format, and can be found in [/report folder](report).
The html report can be found [here](/report/Supplementary_report.html)
