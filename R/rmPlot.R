rmPlot <- function(I, y, onwholes, sets=c("GO", "KEGG")){
	
   sets <- match.arg(sets)
   require(paste(sets,"db",sep="."), character.only=TRUE)
   terms <- switch(sets,
      GO={
            gterms <- toTable(GOTERM)
            gterms[match(onwholes,gterms$go_id), "Term"]

          },
      KEGG={
            kterms <- toTable(KEGGPATHID2NAME)
            kterms[match(onwholes, kterms$path_id),
                                     "path_name", drop=FALSE]
          })
    
   aMat <- I[,onwholes]
   ord <- order(y%*%aMat, decreasing=T)
   aOrd <- aMat[,ord]
   tmp <- aOrd*y
   aOrd <- tmp[apply(tmp,1,sum)>0,]

   gOrd <- rownames(aOrd)[aOrd[,1]==1]
   for(i in 2:ncol(aOrd)){
	  aa <- rownames(aOrd)[aOrd[,i]==1]
	  gOrd <- c(gOrd,setdiff(aa,gOrd))
   }
   aOrd <- aOrd[gOrd,]

   term <- terms[ord]
   id <- onwholes[ord]
   p <- par(no.readonly=TRUE)
   glab="none"
   slab="none"
   xlabs<-rownames(aOrd)
   zlabs<-1:ncol(aOrd)

   ## Space for text in xlim ##
   xpos <- apply(aOrd,2,function(x)
       (1:length(x))[x>0][which.max((1:length(x))[x>0])])
   mwidth1 <- strwidth(xlabs,units="inches")
   mwidth2 <- strwidth(id,units="inches")
   mwidth4 <- strwidth(zlabs,units="inches")

   ## Check twidth <= par("pin")[1] ##  
   twidth <- strwidth(term, units="inches") + 0.5*par("cin")[1]
   pwidth <- (xpos+1)*par("cin")[2]+twidth ## plot width, inches ##
   mspace <- p$mgp[2]*par("cin")[2]+0.5*p$mai[4]

   par(pin=c(if(glab=="none") p$pin[1] else max(pwidth),p$pin[2]),
      mai=c(if(glab=="none") p$mai[1] else max(mwidth1)+mspace,
        max(mwidth2)+mspace,p$mai[3],
        if(slab=="none") p$mai[4] else max(mwidth4)+mspace))

   image(1:nrow(aOrd),1:ncol(aOrd), aOrd[,ncol(aOrd):1,drop=FALSE],
        xlab="",ylab="", yaxt="n",
        xaxt=ifelse(glab=="none","s","n"),
        col=gray(seq(1,0.5, length=64)),
        xlim=c(0.5,max((xpos+0.5)/(1-twidth/par("pin")[1]))))
   abline(v=1.5:(nrow(aOrd)+0.5), col='white')
  
   #if(glab!="none") 
   #axis(side=1,at=1:nrow(aOrd),labels=xlabs,las=2)
   axis(side=2,at=1:ncol(aOrd), labels=colnames(aOrd)[ncol(aOrd):1],las=1)
   text(xpos+0.5,ncol(aOrd):1, term,pos=4, cex=0.7)
   par(p)

}