### R code from vignette source 'Rolemodel.Rnw'

###################################################
### code chunk number 1: Rolemodel.Rnw:53-55
###################################################
 options(width=77, continue="  ")
 pdf.options(pointsize=8)


###################################################
### code chunk number 2: genesamp
###################################################
library(org.Hs.eg.db)

genes <- sample(mappedRkeys(org.Hs.egSYMBOL),100)
head(genes)


###################################################
### code chunk number 3: <loadlib
###################################################
library(Rolemodel)


###################################################
### code chunk number 4: gs2edge
###################################################
gs <- gs2edge(genes, n.upp=30, idtype="SYMBOL", lib="org.Hs.eg")

I <- gs$I
I[1:6,1:6]

y <- gs$y
head(y)

edge <- gs$edge
head(edge)


###################################################
### code chunk number 5: t2d
###################################################
data(t2d)
dim(t2d$I)
sum(t2d$y)
newI <- subRM(t2d$I, 5, 20)
newy <- t2d$y[rownames(newI)]
dim(newI)
sum(newy)


###################################################
### code chunk number 6: params
###################################################
alpha <- 0.00019
gamma <- 0.02279
pest <- estP(I=newI, y=newy, alpha, gamma)
pest


###################################################
### code chunk number 7: ILP
###################################################
 alpha <- 0.00019
 gamma <- 0.02279
 p <- 0.00331
# do not run
# res <- ILP(newI, newy, alpha, gamma, p)


###################################################
### code chunk number 8: sequentialRM
###################################################
res <- sequentialRM(newI, newy, nupstart=10, by=1, alpha, gamma, p)
str(res)


###################################################
### code chunk number 9: sequentialRM
###################################################
names(res$sol$solution)[res$sol$solution == 1][1:8]


###################################################
### code chunk number 10: bp
###################################################
eout <- I2edge(newI)
 bp.out <- bp(whole=colnames(newI), part=newy, edge=eout, alpha, gamma, p, 
nburn=1000000, ngen=10000000, sub=1000, penalty=5, initial="random")
 str(bp.out)


###################################################
### code chunk number 11: bp
###################################################
bp.out$Name[bp.out$ActiveProbability>0.5]


###################################################
### code chunk number 12: summarize the result
###################################################
tmp <- bp.out$ActiveProbability[bp.out$Name%in%res$onwholes]
ord <- order(tmp, decreasing = T)
output<- data.frame(MAP = res$onwholes[ord], P.MFA = tmp[ord])
print(output)


###################################################
### code chunk number 13: do ILP and MCMC analysis simultaneously
###################################################
data(T2D_genelist)
idlist <- idlist[-57] ## no Entrez id for "KLHDC5"
res <- rmTable (idlist, lib="org.Hs.eg", n.upp=20, n.low=5, nupstart=10, by=1, alpha=0.00019, gamma=0.02279, p=0.00331)
str(res)


###################################################
### code chunk number 14: fig1
###################################################
rmPlot(res[[4]]$I, res[[4]]$y, res[[2]]$onwholes)


