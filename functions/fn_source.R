############################
# load tables from the database
############################
# define function: fetch from database

fetchfromDB = function(dbname, tablename, column = "*", additional = ""){
    require(DBI)
    db <- dbConnect(RMySQL::MySQL(),
                    dbname = dbname,
                    user  = 'root',
                    password = 'cqsr4h',
                    host = 'menchelabdb.int.cemm.at')

    command = paste("SELECT",column, "FROM", tablename, additional, sep = " ")

    rs <- dbSendQuery(db, command)
    rd = dbFetch(rs, n = -1)
    dbClearResult(rs)
    dbDisconnect(db)
    return(rd)
}



###########################
# producing graphical summary of a data
############################

f_summary <- function(data_to_plot)
{
  ## univariate data summary
  require(nortest)
  #data <- as.numeric(scan ("data.txt")) #commenting out by mike
  data <- na.omit(as.numeric(as.character(data_to_plot))) #added by mike
  dataFull <- as.numeric(as.character(data_to_plot))

  # first job is to save the graphics parameters currently used
  def.par <- par(no.readonly = TRUE)
  par("plt" = c(.2,.95,.2,.8))
  layout( matrix(c(1,1,2,2,1,1,2,2,4,5,8,8,6,7,9,10,3,3,9,10), 5, 4, byrow = TRUE))

  #histogram on the top left
  h <- hist(data, breaks = "Sturges", plot = FALSE)
  xfit<-seq(min(data),max(data),length=100)
  yfit<-yfit<-dnorm(xfit,mean=mean(data),sd=sd(data))
  yfit <- yfit*diff(h$mids[1:2])*length(data)
  plot (h, axes = TRUE, main = paste(deparse(substitute(data_to_plot))), cex.main=2, xlab=NA)
  lines(xfit, yfit, col="blue", lwd=2)
  leg1 <- paste("mean = ", round(mean(data), digits = 4))
  leg2 <- paste("sd = ", round(sd(data),digits = 4))
  count <- paste("count = ", sum(!is.na(dataFull)))
  missing <- paste("missing = ", sum(is.na(dataFull)))
  legend(x = "topright", c(leg1,leg2,count,missing), bty = "n")

  ## normal qq plot
  qqnorm(data, bty = "n", pch = 20)
  qqline(data)
  p <- ad.test(data)
  leg <- paste("Anderson-Darling p = ", round(as.numeric(p[2]), digits = 4))
  legend(x = "topleft", leg, bty = "n")

  ## boxplot (bottom left)
  boxplot(data, horizontal = TRUE)
  leg1 <- paste("median = ", round(median(data), digits = 4))
  lq <- quantile(data, 0.25)
  leg2 <- paste("25th percentile =  ", round(lq,digits = 4))
  uq <- quantile(data, 0.75)
  leg3 <- paste("75th percentile = ", round(uq,digits = 4))
  legend(x = "top", leg1, bty = "n")
  legend(x = "bottom", paste(leg2, leg3, sep = "; "), bty = "n")

  ## the various histograms with different bins
  h2 <- hist(data,  breaks = (0:20 * (max(data) - min (data))/20)+min(data), plot = FALSE)
  plot (h2, axes = TRUE, main = "20 bins")

  h3 <- hist(data,  breaks = (0:10 * (max(data) - min (data))/10)+min(data), plot = FALSE)
  plot (h3, axes = TRUE, main = "10 bins")

  h4 <- hist(data,  breaks = (0:8 * (max(data) - min (data))/8)+min(data), plot = FALSE)
  plot (h4, axes = TRUE, main = "8 bins")

  h5 <- hist(data,  breaks = (0:6 * (max(data) - min (data))/6)+min(data), plot = FALSE)
  plot (h5, axes = TRUE,main = "6 bins")

  ## the time series, ACF and PACF
  plot (data, main = "Time series", pch = 20, ylab = paste(deparse(substitute(data_to_plot))))
  acf(data, lag.max = 20)
  pacf(data, lag.max = 20)

  ## reset the graphics display to default
  par(def.par)

  #original code for f_summary by respiratoryclub

}


############################################
## create a function to extract nodes of interest out of the network and get the component size
############################################
lcc_subtract = function(g, gene_set){
  g_exclude = g - setdiff(V(g)$name, gene_set)
  lcc = components(g_exclude)$csize
  return(lcc)
}


##########################################
## Translational table from entrez to gene name
##########################################

# This function can take any of the columns(org.Hs.eg.db) as type and keys as long as the row names are in the format of the keys argument
IDconvert <- function(genelist, from="ENTREZID", to="SYMBOL"){
  #note from and to is taken from Affymatrix Human Genome U95 Set Annotation data; library('hgu95av2.db')
  # columns(hgu95av2.db) shows all possibility for to and from
  #  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
  # [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
  # [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"        
  # [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
  # [21] "PROBEID"      "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
  # [26] "UNIGENE"      "UNIPROT"   
  
  require("AnnotationDbi")
  require("org.Hs.eg.db")
  
  geneSymbols <- mapIds(org.Hs.eg.db, keys=genelist, column=to, keytype=from, multiVals="first")
  
  # get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
  #inds <- which(!is.na(geneSymbols))
  #found_genes <- geneSymbols[inds]
  
  return(geneSymbols)
}

# same as previous but lead with matrix row names
# This function can take any of the columns(org.Hs.eg.db) as type and keys as long as the row names are in the format of the keys argument
IDconvert_matrix <- function(df,  from="ENTREZID", to="SYMBOL"){
  require("AnnotationDbi")
  require("org.Hs.eg.db")
  
  geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(df), column=to, keytype=from, multiVals="first")
  
  # get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
  inds <- which(!is.na(geneSymbols))
  found_genes <- geneSymbols[inds]
  
  # subset your data frame based on the found_genes
  df2 <- df[names(found_genes), ]
  rownames(df2) <- found_genes
  return(df2)
}


############################
# dlapply: same as lapply but recursive so the function applies to all sublists (level 2)
############################
## define a function for applying on two levels of list
dlapply = function(list, func){
  nested_output = lapply(list, function(l1) 
    lapply(l1, function(l2) func(l2)))
  return(nested_output)
}

#############################
# GTEx: tissue mapped to gganatogram tissues
#############################



##############################
# create ggplot2 custom theme
###############################

theme_new <- function(base_size = 12, base_family = "Helvetica", xtilt = 0){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      #line = element_line(colour="black"),
      #text = element_text(colour="black"),
      axis.title = element_text(size = 14),
      #axis.text = element_text(colour="black", size=8),
      #strip.text = element_text(size=12),
      legend.key=element_rect(colour=NA, fill =NA),
      panel.grid = element_blank(),   
      panel.border = element_rect(fill = NA, colour = "black", size=1),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = NA)
    )
}
