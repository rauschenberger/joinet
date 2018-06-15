
#' @export
#' @title
#' Prepare data matrices
#' 
#' @description
#' This function removes duplicate samples from each matrix,
#' only retains samples appearing in all matrices,
#' and brings them into the same order.
#' 
#' @param ...
#' matrices with samples in the rows and variables in the columns
#' (separated by commas)
#' 
#' @param message
#' display messages\strong{:} logical
#' 
#' @examples
#' NA
#' 
match.samples <- function(...,message=TRUE){
    
    list <- list(...)
    if(length(list)==1 & is.list(list[[1]])){list <- list[[1]]}
    if(is.null(names(list))){
        names(list) <- sapply(substitute(list(...))[-1],deparse)
    }
    names <- names(list)
    

    # check input
    cond <- sapply(list,function(x) !is.matrix(x))
    if(any(cond)){
        stop("Provide matrices!",call.=FALSE)
    }
    cond <- sapply(list,function(x) is.null(rownames(x)))
    if(any(cond)){
        stop("Provide row names!",call.=FALSE)
    }
    
    # remove duplicated samples
    duplic <- lapply(list,function(x) duplicated(x))
    for(i in seq_along(list)){
        percent <- round(100*mean(duplic[[i]]))
        if(message){message(percent,"% duplicates in \"",names[i],"\"")}
        list[[i]] <- list[[i]][!duplic[[i]],]
    }
    
    # retain overlapping samples
    all <- Reduce(f=intersect,x=lapply(list,rownames))
    for(i in seq_along(list)){
        percent <- round(100*mean(rownames(list[[i]]) %in% all))
        if(message){message(percent,"% overlap in \"",names[i],"\"")}
        list[[i]] <- list[[i]][all,]
    }
    
    # check output
    cond <- sapply(list,function(x) any(duplicated(rownames(x))))
    if(any(cond)){
        stop("Duplicate samples!",call.=FALSE)
    }
    cond <- sapply(list,function(x) nrow(x)!=nrow(list[[1]]))
    if(any(cond)){
        stop("Different sample sizes!",call.=FALSE)
    }
    cond <- sapply(list,function(x) any(rownames(x)!=rownames(list[[1]])))
    if(any(cond)){
        stop("Different sample names!",call.=FALSE)
    }
    
    return(list)
}

#' @export
#' @title
#' Adjust library sizes
#' 
#' @description
#' This function adjusts RNA-seq expression data for different library sizes.
#' 
#' @param x
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (variables)
#' 
#' @examples
#' NA
#' 
adjust.samples <- function(x){
    if(!is.matrix(x)){
        stop("no matrix argument",call.=FALSE)
    }
    if(!is.numeric(x)){
        stop("no numeric argument",call.=FALSE)
    }
    if(!is.integer(x)){
        warning("non-integer values",call.=FALSE)
    }
    if(any(x)<0){
        warning("negative values",call.=FALSE)
    }
    n <- nrow(x); p <- ncol(x)
    lib.size <- rowSums(x)
    norm.factors <- edgeR::calcNormFactors(object=t(x),lib.size=lib.size)
    gamma <- norm.factors*lib.size/mean(lib.size)
    gamma <- matrix(gamma,nrow=n,ncol=p,byrow=FALSE)
    x <- x/gamma
    return(x)
}

#' @export
#' @title
#' Adjust exon length
#' 
#' @description
#' This function adjusts exon expression data for different exon lengths.
#' 
#' @param x
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (exons)
#' 
#' @param group
#' gene (not exon) names\strong{:} vector of length \eqn{p}
#' 
#' @param offset
#' exon start positions\strong{:} vector of length \eqn{p}
#' exon end positions\strong{:} vector of length \eqn{p}
#' 
#' @details
#' No information on chromosomes required.
#' 
#' @examples
#' NA
#' 
adjust.covariates <- function(x,offset,group){
    if(!is.numeric(x)|!is.matrix(x)){
        stop("Argument \"x\" is no numeric matrix.",call.=FALSE)
    }
    if(!is.numeric(offset)|!is.vector(offset)){
        stop("Argument \"offset\" is no numeric vector.",call.=FALSE)
    }
    if(any(offset<0)){
        stop("Argument \"offset\" takes negative values",call.=FALSE)   
    }
    if(!is.character(group)|!is.vector(group)){
        stop("Argument \"group\" is no character vector.",call.=FALSE)
    }
    if(ncol(x)!=length(group)|ncol(x)!=length(offset)){
        stop("Contradictory dimensions.",call.=FALSE)
    }
    n <- nrow(x); p <- ncol(x); names <- dimnames(x)
    x <- as.numeric(x)
    offset <- rep(offset,each=n)
    group <- strsplit(group,split=",")
    group <- sapply(group,function(x) x[[1]][1])
    group <- rep(group,each=n)
    lmer <- lme4::lmer(x ~ offset + (1|group))
    x <- matrix(stats::residuals(lmer),nrow=n,ncol=p,dimnames=names)
    x <- x-min(x)
    return(x)
}

#' @export
#' @title
#' Search for genes
#' 
#' @description
#' This function retrieves all genes on a chromosome.
#' 
#' @param chr
#' chromosome\strong{:} integer 1-22
#' 
#' @param path
#' path to gene transfer format files (.gtf)
#' 
#' @param release
#' character "NCBI36", "GRCh37", or "GRCh38"
#' 
#' @param build
#' integer 49-91
#' 
#' @details
#' This functions ...
#' 
#' @examples
#' NA
#' 
map.genes <- function(chr,path=getwd(),release="GRCh37",build="71"){
    file <- paste0("Homo_sapiens.",release,".",build,".gtf")
    if(!file.exists(file.path(path,file))){
        url <- paste0("ftp://ftp.ensembl.org/pub/release-",build,
                      "/gtf/homo_sapiens/",file,".gz")
        destfile <- file.path(path,paste0(file,".gz"))
        utils::download.file(url=url,destfile=destfile,method="auto")
        R.utils::gunzip(filename=destfile,remove=FALSE,overwrite=TRUE)
    }
    object <- refGenome::ensemblGenome()
    refGenome::basedir(object) <- path
    refGenome::read.gtf(object,filename=file)
    genes <- refGenome::getGenePositions(object=object,by="gene_id")
    genes <- genes[genes$seqid==chr & genes$gene_biotype=="protein_coding",]
    genes <- genes[,c("gene_id","seqid","start","end")]
    rownames(genes) <- NULL
    colnames(genes)[colnames(genes)=="seqid"] <- "chr"
    return(genes)
}

#' @export
#' @title
#' Search for exons
#' 
#' @description
#' This function
#' 
#' @param gene_id
#' gene names\strong{:} vector with one entry per gene (gene names!)
#' 
#' @param exon_id
#' exon names\strong{:} vector with one entry per exon
#' (also gene names! separated by comma if multiple genes)
#' 
#' @details
#' The exon names should contain the gene names. For each gene, this function
#' returns the indices of the exons.
#' 
#' @examples
#' NA
#'
map.exons <- function(gene_id,exon_id){
    p <- length(gene_id)
    exons <- list()
    pb <- utils::txtProgressBar(min=0,max=p,style=3)
    for(i in seq_len(p)){
        utils::setTxtProgressBar(pb=pb,value=i)
        which <- as.integer(grep(pattern=gene_id[i],x=exon_id))
        exons[[i]] <- which
    }
    return(exons)
}

#' @export
#' @title
#' Search for SNPs
#' 
#' @description
#' This function
#' 
#' @param gene.chr
#' chromosome\strong{:}
#' numeric vector with one entry per gene
#' 
#' @param gene.start
#' start position\strong{:}
#' numeric vector with one entry per gene
#' 
#' @param gene.end
#' end position\strong{:}
#' numeric vector with one entry per gene
#' 
#' @param snp.chr
#' integer 1-22
#' 
#' @param snp.pos
#' chromosomal position of SNPs\strong{:}
#' numeric vector with one entry per SNP
#' 
#' @details
#' This functions ...
#' 
#' @examples
#' NA
#' 
map.snps <- function(gene.chr,gene.start,gene.end,snp.chr,snp.pos){
    if(length(gene.chr)!=length(gene.start)|length(gene.chr)!=length(gene.end)){
        stop("Invalid.",call.=FALSE)
    }
    p <- length(gene.start)
    snps <- data.frame(from=integer(length=p),to=integer(length=p))
    pb <- utils::txtProgressBar(min=0,max=p,style=3)
    for(i in seq_len(p)){ # 
        utils::setTxtProgressBar(pb=pb,value=i)
        chr <- snp.chr == gene.chr[i]
        if(!any(chr)){next}
        start <- snp.pos >= gene.start[i] - 1*10^3
        end <- snp.pos <= gene.end[i] + 0
        which <- as.integer(which(chr & start & end))
        if(length(which)==0){next}
        snps$from[i] <- min(which)
        snps$to[i] <- max(which)
        if(length(which)==1){next}
        if(!all(diff(which)==1)){stop("SNPs are in wrong order!")}
    }
    return(snps)
}

#' @export
#' @title
#' Drop "trivial" genes
#' 
#' @description
#' This function
#' 
#' @param map
#' list with names "genes", "exons", and "snps"
#' (output from \code{map.genes}, \code{map.exons}, and \code{map.snps})
#' 
#' @details
#' This functions drops genes without SNPs or with a single exon.
#' 
#' @examples
#' NA
#' 
drop.trivial.genes <- function(map){
    
    # check input
    if(length(map)!=3){
        stop("Unexpected argument length.",call.=FALSE)
    }
    if(any(names(map)!=c("genes","exons","snps"))){
        stop("Unexpected argument names.",call.=FALSE)
    }
    
    # search
    p <- nrow(map$genes)
    pass <- rep(NA,times=p)
    pb <- utils::txtProgressBar(min=0,max=p,style=3)
    for(i in seq_len(p)){
        utils::setTxtProgressBar(pb=pb,value=i)
        ys <- map$exons[[i]]
        check <- logical()
        # Exclude genes without SNPs:
        check[1] <- map$snps$from[i] > 0
        check[2] <- map$snps$to[i] > 0
        # Exclude genes with single exon:
        check[3] <- length(ys) > 1
        pass[i] <- all(check)
    }
    
    # check output
    if(any(pass[map$snps$to==0 & map$snps$from==0])){
        stop("Genes without any SNPs.",call.=FALSE)
    }
    if(any(pass[sapply(map$exons,length)<2])){
        stop("Genes without multiple exons.",call.=FALSE)
    }
    
    #map$genes <- map$genes[pass,]
    #map$exons <- map$exons[pass]
    #map$snps <- map$snps[pass,]
    return(pass) # temporary
}


#' @export
#' @title
#' Conduct tests
#' 
#' @description
#' This function
#' 
#' @param Y
#' exon expression\strong{:}
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (exons)
#' 
#' @param X
#' SNP genotype\strong{:}
#' matrix with \eqn{n} rows (samples) and \eqn{q} columns (SNPs)
#' 
#' @param map
#' list with names "genes", "exons", and "snps"
#' (output from \code{map.genes}, \code{map.exons}, and \code{map.snps})
#' 
#' @param i
#' gene index\strong{:}
#' integer between \eqn{1} and \code{nrow(map$genes)}
#' 
#' @param limit
#' cutoff for rounding \code{p}-values
#' 
#' @param steps
#' size of permutation chunks\strong{:}
#' integer vector
#' 
#' @param 
#' 
#' @details
#' The maximum number of permutations equals \code{sum(steps)}. Permutations is
#' interrupted if at least \code{limit} test statistics for the permuted data
#' are larger than the test statistic for the observed data.
#' 
#' @examples
#' NA
#' 
test.single <- function(Y,X,map,i,limit=NULL,steps=NULL,rho=c(0,0.5)){
    
    if(is.null(limit)){limit <- 5}
    if(is.null(steps)){steps <- c(10,20,20,50)}
    
    # check input
    if(!is.numeric(limit)){
        stop("Argument \"limit\" is not numeric.",call.=FALSE)
    }
    if(limit<1){
        stop("Argument \"limit\" is below one.",call.=FALSE)
    }
    if(!is.numeric(step)|!is.vector(steps)){
        stop("Argument \"steps\" is no numeric vector.",call.=FALSE)
    }
    if(sum(steps)<2){
        stop("Too few permutations \"sum(steps)\".",call.=FALSE)
    }
    
    # extract data
    ys <- map$exons[[i]]
    y <- Y[,ys,drop=FALSE]
    xs <- seq(from=map$snps$from[i],to=map$snps$to[i],by=1)
    x <- X[,xs,drop=FALSE]
    
    # test effects
    pvalue <- rep(x=NA,times=length(rho))
    for(j in seq_along(rho)){
        tstat <- spliceQTL:::G2.multin(
            dep.data=y,indep.data=x,nperm=steps[1]-1,rho=rho[j])$Sg
        for(nperm in steps[-1]){
            tstat <- c(tstat,spliceQTL:::G2.multin(
                dep.data=y,indep.data=x,nperm=nperm,rho=rho[j])$Sg[-1])
            if(sum(tstat >= tstat[1]) >= limit){break}
        }
        pvalue[j] <- mean(tstat >= tstat[1],na.rm=TRUE)
    }

    return(pvalue)
}


test.multiple <- function(Y,X,map){
    p <- nrow(map$genes)
    
    # permutations
    min <- 100
    max <- p/0.05+1
    limit <- ceiling(0.05*max/p)
    base <- 1.5 # adjust sequence
    from <- log(min,base=base)
    to <- log(max,base=base)
    steps <- c(min,diff(unique(round(base^(seq(from=from,to=to,length.out=20))))))
    if(max != sum(steps)){stop("Invalid combination?",call.=FALSE)}
    
    # parallel computation
    type <- ifelse(test=.Platform$OS.type=="windows",yes="PSOCK",no="FORK")
    cluster <- parallel::makeCluster(spec=8,type=type)
    parallel::clusterSetRNGStream(cl=cluster,iseed=1)
    test.single <- spliceQTL::test.single
    parallel::clusterExport(cl=cluster,varlist=c("test.single","Y","X","map","limit","steps"),envir=environment())
    parallel::clusterEvalQ(cl=cluster,expr=library(spliceQTL))
    start <- Sys.time()
    pvalue <- parallel::parLapply(cl=cluster,X=1:p,fun=function(i) test.single(Y=Y,X=X,map=map,i=i,limit=limit,steps=steps))
    end <- Sys.time()
    parallel::stopCluster(cluster)
    rm(cluster)
    
    # tyding up
    #names <- names(pvalue[[1]])
    #pvalue <- lapply(names,function(x) lapply(pvalue,function(y) y[[x]]))
    #names(pvalue) <- names
    #attributes(pvalue)$time <- end - start
    return(pvalue)
}




#--- spliceQTL test functions --------------------------------------------------

# Function: G2.multin
# This is to compute the G2 test statistic under the assumption that the response follows a multinomial distribution
### Input 
### dep data and indep data with samples on the rows and genes on the columns
### grouping: Either a logical value = F or a matrix with a single column and same number of rows as samples. 
###         Column name should be defined.
###         Contains clinical information of the samples. 
###         Should have two groups only. 
### nperm : number of permutations 
### rho: the null correlation between SNPs
### mu: the null correlation between observations corresponding to different exons and different individuals

### Output
### A list containing G2 p.values and G2 test statistics

### Example : G2T = G2(dep.data = cgh, indep.data = expr, grouping=F, stand=TRUE, nperm=1000)
### G2 p.values : G2T$G2p
### G2 TS : G2T$$Sg

G2.multin <- function(dep.data,indep.data,stand=TRUE,nperm=100,grouping=F,rho=0,mu=0){
    
    nperm = nperm
    ## check for the number of samples in dep and indep data
    
    
    if (nrow(dep.data)!=nrow(indep.data)){
        cat("number of samples not same in dep and indep data","\n")
    }
    
    if(any(abs(rho)>1)){
        cat("correlations rho larger than abs(1) are not allowed")
    }
    
    nresponses <- ncol(dep.data)
    ncovariates <- ncol(indep.data)
    ### centering and standardizing the data are not done in this case
    
    #  dep.data = scale(dep.data,center=T,scale=stand)
    #  indep.data = scale(indep.data,center=T,scale=stand)
    
    #### No  grouping of the samples.
    
    ## Calculate U=(I-H)Y and UU', where Y has observations on rows; also tau.mat=X*W.rho*X', 
    ##   where X has observations on rows and variables on columns
    ##  and W.rho = I + rho*(J-I), a square matrix with as many rows as columns in X
    ## NOTE: this formulation uses X with n obs on the rows and m covariates no the columns, so it is the transpose of the first calculations
    nsamples <- nrow(dep.data)
    n.persample <- rowSums(dep.data)
    n.all <- sum(dep.data)
    H <- (1/n.all)*matrix( rep(n.persample,each=nsamples),nrow=nsamples,byrow=T)
    U <- (diag(rep(1,nsamples)) - H) %*% dep.data
    ## Now we may have a vector of values for rho - so we define tau.mat as an array, with the 3rd index corresponding to the value of rho
    tau.mat <- array(0,dim=c(nsamples,nsamples,length(rho)))
    for(xk in 1:length(rho))  
    {  
        if (rho[xk]==0) { tau.mat[,,xk] <- tcrossprod(indep.data) } 
        else { w.rho <- diag(rep(1,ncovariates)) + rho[xk]*(tcrossprod(rep(1,ncovariates)) -diag(rep(1,ncovariates))  )
        tau.mat[,,xk] <- indep.data %*% w.rho %*% t(indep.data)}
        
    }
    ######################################
    ### NOTES ARMIN START ################
    # all(X %*% t(X) == tau.mat[,,1]) # rho = 0 -> TRUE
    # all(X %*% (t(X) %*% X) %*% t(X) == tau.mat[,,1]) # rho = 1
    # plot(as.numeric(X %*% (t(X) %*% X) %*% t(X)),as.numeric(tau.mat[,,1]))
    ### NOTES ARMIN END ##################
    ######################################
    samp_names = 1:nsamples ## this was rownames(indep.data), but I now do this so that rownames do not have to be added to the array tau.mat
    Sg = get.g2stat.multin(U,mu=mu,rho=rho,tau.mat)
    ### now we will have a vector as result, with one value per combination of values of rho and mu
    #
    ### G2 
    ### Permutations
    # When using permutations: only the rows of tau.mat are permuted
    # To check how the permutations can be efficiently applied, see tests_permutation_g2_multin.R
    
    
    perm_samp = matrix(0, nrow=nrow(indep.data), ncol=nperm)   ## generate the permutation matrix
    for(i in 1:ncol(perm_samp)){
        perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
    }
    
    ## permutation starts - recompute tau.mat  (or recompute U each time)
    for (perm in 1:nperm){
        tau.mat.perm = tau.mat[perm_samp[,perm],,,drop=FALSE]          # permute rows
        tau.mat.perm = tau.mat.perm[,perm_samp[,perm],,drop=FALSE]     # permute columns
        
        Sg = c(Sg,spliceQTL:::get.g2stat.multin(U, mu=mu,rho=rho,tau.mat.perm) )
    }
    
    
    ########################################################################
    
    #### G2 test statistic
    # *** recompute for a vector of values for each case - just reformat the result with as many rows as permutations + 1,
    # and as many columns as combinations of values of rho and mu
    Sg = matrix(Sg,nrow=nperm+1,ncol=length(mu)*length(rho))
    colnames(Sg) <- paste(rep("rho",ncol(Sg)),rep(1:length(rho),each=length(mu)),rep("mu",ncol(Sg)),rep(1:length(mu),length(rho)) )
    
    ### Calculte G2 pval
    G2p =  apply(Sg,2,spliceQTL:::get.pval.percol) 
    
    return (list(perm = perm_samp,G2p = G2p,Sg = Sg))
}

# Function: get.g2stat.multin
# Computes the G2 test statistic given two data matrices, under a multinomial distribution
# This is used internally by the G2 function
# Inputs: 
#  U = (I-H)Y, a n*K matrix where n=number obs and K=number multinomial responses possible
#  tau.mat = X' W.rho X, a n*n matrix : both square, symmetric matrices with an equal number of rows
# Output: test statistic (single value)
# 
get.g2stat.multin <- function(U, mu, rho, tau.mat)
{
    g2tstat <- NULL
    for(xk in 1:length(rho))
    {
        for(xj in 1:length(mu))
        {
            if(mu[xj]==0) { g2tstat <- c(g2tstat, sum( diag( tcrossprod(U) %*% tau.mat[,,xk] ) ) )
            } else {
                g2tstat <- c(g2tstat, (1-mu[xj])*sum(diag( tcrossprod(U) %*% tau.mat[,,xk] ) ) + mu[xj]*sum( t(U) %*% tau.mat[,,xk] %*% U )  )
            }
            
        }
    }
    g2tstat
}

# Function: get.pval.percol
# This function takes a vector containing the observed test stat as the first entry, followed by values generated by permutation,
# and computed the estimated p-value
# Input
# x: a vector with length nperm+1
# Output
# the pvalue computed
get.pval.percol <- function(x){
    pval = mean(x[1]<= c(Inf , x[2:length(x)]))
    pval
}

