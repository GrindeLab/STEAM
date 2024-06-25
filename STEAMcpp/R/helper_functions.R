# extract object from .RData file
get_obj <- function(Rdata){
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata,
                  "\nReturning only the first object"))
  }
  return(get(objname))
}

# convert genetic distance into recombination fraction
L_to_theta <- function(cM){
  L <- cM/100
  theta <- 0.5*(1-exp(-2*L))
  return(theta)
}

# function to get distance between pairs of SNPs
get_dist <- function(s1, s2, gen.map){
  pos1 <- gen.map$cM[match(s1, gen.map$snp.id)]
  pos2 <- gen.map$cM[match(s2, gen.map$snp.id)]
  return(abs(pos2-pos1))
}

# get correlation for one pair of SNPs
get_corr <- function(s1, s2, afr, eur, nam){
  # load ancestries at SNPs s1 and s2
  afr.anc <- snpgdsGetGeno(afr, snp.id = c(s1,s2), with.id = TRUE, verbose=F)
  eur.anc <- snpgdsGetGeno(eur, snp.id = c(s1,s2), with.id = TRUE, verbose=F)
  nam.anc <- snpgdsGetGeno(nam, snp.id = c(s1,s2), with.id = TRUE, verbose=F)
  # get indices for each SNP (make sure they're in order we expected)
  e1 <- which(eur.anc$snp.id==s1); e2 <- which(eur.anc$snp.id == s2)
  a1 <- which(afr.anc$snp.id==s1); a2 <- which(afr.anc$snp.id == s2)
  n1 <- which(nam.anc$snp.id==s1); n2 <- which(nam.anc$snp.id == s2)
  # get correlation between all possible pairs of ancestries
  c11 <- cor(afr.anc$g[,a1], afr.anc$g[,a2])
  c12 <- cor(afr.anc$g[,a1], eur.anc$g[,e2])
  c13 <- cor(afr.anc$g[,a1], nam.anc$g[,n2])
  c21 <- cor(eur.anc$g[,e1], afr.anc$g[,a2])
  c22 <- cor(eur.anc$g[,e1], eur.anc$g[,e2])
  c23 <- cor(eur.anc$g[,e1], nam.anc$g[,n2])
  c31 <- cor(nam.anc$g[,n1], afr.anc$g[,a2])
  c32 <- cor(nam.anc$g[,n1], eur.anc$g[,e2])
  c33 <- cor(nam.anc$g[,n1], nam.anc$g[,n2])
  return(c(c11,c12,c13,c21,c22,c23,c31,c32,c33))
}
