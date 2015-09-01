
rmTable <- function(idlist,
                      idtype=c("SYMBOL","ENTREZ","REFSEQ",
                     "ENSEMBL","ACCNUM","UNIPROT","PMID"),
                      symbol = TRUE,
                      lib, library.loc=NULL,
                      n.upp=20, n.low=5,
                      sets=c("GO","KEGG"),
                      nupstart=10, by=1, alpha=0.01, gamma=0.02, p, nburn=1000000, 
                      ngen=10000000, sub=1000, penalty=5, initial="random",
                      annotate=TRUE){

    sets <- match.arg(sets)
    idtype <- match.arg(idtype)

    out <- gs2edge(idlist=idlist, idtype=idtype, lib=lib,
                   library.loc=library.loc,
                   n.upp=n.upp, n.low=n.low,
                   sets=sets)
    
    if(missing(p)) p <- estP(out$I, out$y, alpha=alpha, gamma=gamma)
    
    message("Running MAP analysis ...")
    result.rm <- sequentialRM(out$I, out$y, nupstart=nupstart, by=by,
                              alpha=alpha, gamma=gamma, p=p)
    message("Running MCMC analysis ...")
    
    result.bp <- bp(whole=colnames(out$I), part=out$y, edge=out$edge,
                    alpha=alpha, gamma=gamma, p=p, nburn=nburn, ngen=ngen, 
                    sub=sub, penalty=penalty, initial=initial)

    pr <- result.bp[result.bp$Name %in% result.rm$onwholes, "ActiveProbability"]
    pord <- order(pr,decreasing=TRUE)

    require(paste(lib, "db", sep = "."), character.only = TRUE, 
        lib.loc = library.loc)
    if(symbol){
    	  hh <- get(paste(lib, "SYMBOL", sep = ""))
      mapped_genes <- mappedLkeys(hh)
      xx <- as.list(hh[mapped_genes])
      symbols <- unlist(xx[rownames(out$I)])
    }
    else 
      symbols <- rownames(out$I)
    set.size <- apply(out$I[,as.character(result.rm$onwholes)], 2, sum)
    number.associated.genes <- as.vector( out$y%*%out$I[,as.character(result.rm$onwholes)] )
    genes <- apply(out$I[,as.character(result.rm$onwholes)], 2, function(x) paste(symbols[x*out$y==1], collapse=';'))
    result.out <- data.frame(MAP=result.rm$onwholes, P.MFA=pr, set.size=set.size, 
    number.associated.genes=number.associated.genes, genes=genes)[pord,]

    if(annotate){
        message("Labelling table ...")
        require(paste(sets,"db",sep="."), character.only=TRUE)
        main <- switch(sets,
               GO={
                   gterms <- toTable(GOTERM)
                   data.frame(gterms[match(result.out$MAP,gterms$go_id),
                                     c("Term","Ontology")], result.out)
               },
               KEGG={
                   kterms <- toTable(KEGGPATHID2NAME)
                   data.frame(kterms[match(result.out$MAP, kterms$path_id),
                                     "path_name", drop=FALSE], result.out)
               })
    } else main <- result.out
    
    return(list(rm.table=main, rm.sol=result.rm, bp.sol=result.bp, incidence.mat=out))
}
