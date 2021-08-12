# Network analysis reveals rare disease signatures across multiple levels of biological organization

Supplementary codes, reproducible walk-through reports, and Shiny app complementing Buphamalai. et.al.

- [Repository structure](#this-repository)
- [The Explorer introduction](#the-explorer)
- [Reproducible report](#the-reproducible-walkthrough-guides)
- [Session information](#session-information)
- [System requirements](#system-requirements)

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

The Explorer can be launched via [http://menchelab.com/MultiOmeExplorer](http://menchelab.com/MultiOmeExplorer)

## The Reproducible Walkthrough Guides

[This supplementary report](/report/Supplementary_report.html) is aimed to be a reproducible walk-through guide for figures and analyses complementing the manuscript: Buphamalai et.al., *Network analysis reveals rare disease signatures across multiple levels of biological organization*, submitted to Nature Communications. (Almost) all of the figures and analyses can be reproduced via walk-through reports organized in the same order as the manuscript in RMarkdown format, and can be found in [report/](report) folder. 

Please find the following guidelines:

1. The appearance of sections in this document is at the same order as in the manuscript, and can be navigated using the Table of Content (ToC) appeared on the top left corner of this document, with subsections corresponding to exact figures/statistics 
2. The corresponding code chunks to each figures/analyses are provided above the output, and are hidden by default. To expand each code chunk, click `Show code`.
3. This report mainly contains visualization and post-processing of major analyses. Heavier computations were pre-computed, and corresponding `R`, `sh`, or `py` scripts required for each analysis are mentioned for each section. These files are available in [/source](/source) folder.
4. The pre-computed results are saved in `./cache` folder and can be downloaded from: [link](https://drive.google.com/file/d/1IpH8v2RLDRRXdqFvLZq-rUTwyqhGAezP/view?usp=sharing). Unzip the folder under the main directory (`./cache`). 
5. The corresponding `Rmd` files used to produce this report can be found in [/report](/report) folder.

## Session information

Analyses were conducted using the R Statistical language (version 3.6.3; R Core Team, 2020) on macOS 10.16, using the packages voronoiTreemap (version 0.2.1; Alexander Kowarik et al., 2021), cowplot (version 1.1.1; Claus Wilke, 2020), igraph (version 1.2.6; Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. https://igraph.org), RColorBrewer (version 1.1.2; Erich Neuwirth, 2014), ggplot2 (version 3.3.3; Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.), stringr (version 1.4.0; Hadley Wickham, 2019), tidyr (version 1.1.2; Hadley Wickham, 2020), forcats (version 0.5.1; Hadley Wickham, 2021), scales (version 1.1.1; Hadley Wickham and Dana Seidel, 2020), readr (version 1.4.0; Hadley Wickham and Jim Hester, 2020), dplyr (version 1.0.4; Hadley Wickham et al., 2021), ggforestplot (version 0.1.0; Ilari Scheinin et al., 2021), rmarkdown (version 2.7; JJ Allaire and Yihui Xie and Jonathan McPherson and Javier Luraschi and Kevin Ushey and Aron Atkins and Hadley Wickham and Joe Cheng and Winston Chang and Richard Iannone, 2021), ggrepel (version 0.9.1; Kamil Slowikowski, 2021), tibble (version 3.1.0; Kirill Müller and Hadley Wickham, 2021), purrr (version 0.3.4; Lionel Henry and Hadley Wickham, 2020), report (version 0.3.0; Makowski et al., 2020), treemap (version 2.4.2; Martijn Tennekes, 2017), ggstatsplot (version 0.7.0; Patil, 2018), pacman (version 0.5.1; Rinker et al., 2017), ggraph (version 2.0.4; Thomas Lin Pedersen, 2020), patchwork (version 1.1.1; Thomas Lin Pedersen, 2020), tidygraph (version 1.2.0; Thomas Lin Pedersen, 2020), MASS (version 7.3.53; Venables et al., 2002), tidyverse (version 1.3.0; Wickham et al., 2019), pROC (version 1.17.0.1; Xavier Robin et al., 2011) and knitr (version 1.31; Yihui Xie, 2021).

### References
  - Alexander Kowarik, Bernhard Meindl, Malaver Vojvodic, Mike Bostock and Franck Lebeau (2021). voronoiTreemap: Voronoi Treemaps with Added Interactivity by Shiny. R package version 0.2.1. https://github.com/uRosConf/voronoiTreemap
  - Claus O. Wilke (2020). cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'. R package version 1.1.1. https://CRAN.R-project.org/package=cowplot
  - Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. https://igraph.org
  - Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer
  - H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
  - Hadley Wickham (2019). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.4.0. https://CRAN.R-project.org/package=stringr
  - Hadley Wickham (2020). tidyr: Tidy Messy Data. R package version 1.1.2. https://CRAN.R-project.org/package=tidyr
  - Hadley Wickham (2021). forcats: Tools for Working with Categorical Variables (Factors). R package version 0.5.1. https://CRAN.R-project.org/package=forcats
  - Hadley Wickham and Dana Seidel (2020). scales: Scale Functions for Visualization. R package version 1.1.1. https://CRAN.R-project.org/package=scales
  - Hadley Wickham and Jim Hester (2020). readr: Read Rectangular Text Data. R package version 1.4.0. https://CRAN.R-project.org/package=readr
  - Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.4. https://CRAN.R-project.org/package=dplyr
  - Ilari Scheinin, Maria Kalimeri, Vilma Jagerroos, Juuso Parkkinen, Emmi Tikkanen, Peter Würtz and Antti Kangas (2021). ggforestplot: Forestplots of Measures of Effects and Their Confidence Intervals. https://nightingalehealth.github.io/ggforestplot/index.html, https://github.com/nightingalehealth/ggforestplot.
  - JJ Allaire and Yihui Xie and Jonathan McPherson and Javier Luraschi and Kevin Ushey and Aron Atkins and Hadley Wickham and Joe Cheng and Winston Chang and Richard Iannone (2021). rmarkdown: Dynamic Documents for R. R package version 2.7. URL https://rmarkdown.rstudio.com.
  - Kamil Slowikowski (2021). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. R package version 0.9.1. https://CRAN.R-project.org/package=ggrepel
  - Kirill Müller and Hadley Wickham (2021). tibble: Simple Data Frames. R package version 3.1.0. https://CRAN.R-project.org/package=tibble
  - Lionel Henry and Hadley Wickham (2020). purrr: Functional Programming Tools. R package version 0.3.4. https://CRAN.R-project.org/package=purrr
  - Makowski, D., Ben-Shachar, M.S., Patil, I. & Lüdecke, D. (2020). Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption. CRAN. Available from https://github.com/easystats/report. doi: .
  - Martijn Tennekes (2017). treemap: Treemap Visualization. R package version 2.4-2. https://CRAN.R-project.org/package=treemap
  - Patil, I. (2018). ggstatsplot: 'ggplot2' Based Plots with Statistical Details. CRAN. Retrieved from https://cran.r-project.org/web/packages/ggstatsplot/index.html
  - R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
  - Rinker, T. W. & Kurkiewicz, D. (2017). pacman: Package Management for R. version 0.5.0. Buffalo, New York. http://github.com/trinker/pacman
  - Thomas Lin Pedersen (2020). ggraph: An Implementation of Grammar of Graphics for Graphs and Networks. R package version 2.0.4. https://CRAN.R-project.org/package=ggraph
  - Thomas Lin Pedersen (2020). patchwork: The Composer of Plots. R package version 1.1.1. https://CRAN.R-project.org/package=patchwork
  - Thomas Lin Pedersen (2020). tidygraph: A Tidy API for Graph Manipulation. R package version 1.2.0. https://CRAN.R-project.org/package=tidygraph
  - Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0
  - Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
  - Xavier Robin, Natacha Turck, Alexandre Hainard, Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez and Markus Müller (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics, 12, p. 77. DOI: 10.1186/1471-2105-12-77 <http://www.biomedcentral.com/1471-2105/12/77/>
  - Yihui Xie (2021). knitr: A General-Purpose Package for Dynamic Report Generation in R. R package version 1.31.


The reference list of softwares listed above can also be found [here](/report/report_session.md).

### System requirements

Heavier computing tasks were pre-computing on a cluster and stored as cache, available for download in the link listed in the previous section. With pre-computed cached data, the [Report](/report/Supplementary_report.Rmd) should be executed in a local machine (quad-core CPU, 8GB RAM) within 5 minutes. A minimum storage space of ~2.5GB is required.




