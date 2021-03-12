## 1. from a list of network, provide a list of file names with weights
## 2. create the supraadjacency list of those networks - with common gene index
## 3. 


library(AnnotationDbi)
library(Matrix)
library(readMM)
gene_allnet = dlapply(g, function(x) V(x)$name)
gene_allnet = unique(unlist(gene_allnet))
gene_allnet = gene_allnet[!is.na(gene_allnet) & gene_allnet != ""]

gene_allnet_id = 1:length(gene_allnet); names(gene_allnet_id) = gene_allnet


for(network in names(g)){
  for(subnetwork in names(g[[network]])){
    print(subnetwork)
    
    g_current = g[[network]][[subnetwork]]
    el = as_data_frame(as_edgelist(g_current))
    colnames(el) = c("to", "from")
    el$val = 1
    el = el[el$to != el$from,] #get rid of loop
    el = el[!is.na(el$from) & !is.na(el$to),] #remove NA
    
    missing_node = setdiff(gene_allnet, V(g_current)$name)
    el_zero = data.frame(from = missing_node, to = missing_node, val = 0)
    
    el_full = rbind(el, el_zero)
    
    el_full$pos1 = gene_allnet_id[el_full$from]
    el_full$pos2 = gene_allnet_id[el_full$to]
    el_full = el_full[!is.na(el_full$pos1) & !is.na(el_full$pos2),] #remove NA
    
    el_full$swap = el_full$pos1 > el_full$pos2
    el_full$pos1_final = el_full$pos1
    el_full$pos2_final = el_full$pos2
    el_full$pos1_final[el_full$swap] = el_full$pos2[el_full$swap]
    el_full$pos2_final[el_full$swap] = el_full$pos1[el_full$swap]
    
    #g_current = graph_from_data_frame(d = el_full, directed = F, )
    
    A = sparseMatrix(i=el_full$pos1_final,j=el_full$pos2_final,x=el_full$val, dimnames = list(gene_allnet, gene_allnet), symmetric = T)
    diag(A) = 1 #prevent matrix being singular
    degrees = rowSums(A)
    M = A/degrees
    # Minv = Matrix::solve(M)
    # colnames(Ainv) = rownames(Ainv) = V(g[[network]][[subnetwork]])$name
    #colnames(A) = rownames(A) = V(g[[network]][[subnetwork]])$name
    saveRDS(M, paste0("~/Documents/projects/Interactome/Adjacency_matrices_transition//",network,"_",subnetwork,".RDS"))
    # writeMM(M,file=paste0('~/Documents/projects/Interactome/Adjacency_matrices_MM/',network,"_",subnetwork,".txt"))
    #saveRDS(Minv, paste0("~/Documents/projects/Interactome/inverse_adj_",network,"_",subnetwork,".RDS"))
  }
}

