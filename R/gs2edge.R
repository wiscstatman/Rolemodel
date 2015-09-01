## Given an org library, create I and y for input into ILP ##
## idlist = character vector of gene symbols
## lib = db library

gs2edge <- function(idlist,
                 idtype=c("SYMBOL","ENTREZID","REFSEQ",
                     "ENSEMBL","ACCNUM","UNIPROT","PMID"),
                 lib, library.loc=NULL,
                 n.upp=50,n.low=5,
                 sets=c("GO","KEGG")){
    require(paste(lib,"db",sep="."),character.only=TRUE,
            lib.loc=library.loc)

    idtype <- match.arg(idtype)
    sets <- match.arg(sets)
    gs <- c(GO="GO2ALLEGS",KEGG="PATH2EG")

    ## idlist to egs ##
    egs <- if(idtype=="ENTREZID") idlist else
    toTable(revmap(get(paste(lib,idtype,sep="")))[idlist])[,1]
    
    gset <- get(paste(lib,gs[sets],sep=""))
    tgset <- unique(toTable(gset)[,1:2])
    ngset <- table(tgset[,2])

    ngset.sub <- names(ngset)[ngset>=n.low & ngset<=n.upp]
    egs.sub <- mappedLkeys(gset[ngset.sub])
    tgset.sub <- tgset[tgset[,2] %in% ngset.sub,]

    imat <- do.call(cbind,tapply(tgset.sub[,1],tgset.sub[,2],
                   function(i) as.numeric(egs.sub %in% i)))
    yvec <- as.numeric(egs.sub %in% egs)
    rownames(imat) <- names(yvec) <- egs.sub

    list(I=imat,y=yvec,edge=tgset.sub[,2:1])
}

I2edge <- function(I)
     do.call(rbind,lapply(colnames(I), function(i)
         cbind(rep(i,sum(I[,i])),rownames(I)[as.logical(I[,i])])))
