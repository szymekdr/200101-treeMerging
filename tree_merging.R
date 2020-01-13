## Szymek Drobniak
## 16/12/2019
## Individuality analyses for Linhart et al. individual vocal signals paper

BiocManager::install("ggtree")

required.pck <-
  c(
    "asremlPlus",
    "MCMCglmm",
    "gdata",
    "asreml",
    "ape",
    "geiger",
    "ggtree",
    "phangorn",
    "ggstance"
  )
for (p in required.pck) {
  status <- require(p, character.only = T)
  if (!status) {
    install.packages(p)
    require(p, character.only = T)
  }
}


#### trees loading and prep ----

trees_bird <- read.nexus("output_bird.nex")
str(trees_bird[[1]])
ggtree(trees_bird[[1]]) + geom_treescale(0,-1) +
  geom_tiplab() +
  xlim_tree(150)

trees_mammal <- read.nexus("output_mammal.nex")
ggtree(trees_mammal[[1]]) + geom_treescale(0,-1) +
  geom_tiplab() +
  xlim_tree(250)


#### function ----

bind.trees.mrca <-
  function (trees1,
            trees2,
            N = min(length(trees1), length(trees2)),
            MRCA,
            MRCA.SE = NULL) {
    require(ape)
    require(phytools)
    
    ## function used to join two trees provided as random sets
    ## from a posterior of trees
    ## function uses fixed MCRA timepoint or a CI of MRCA estimate
    ## both trees have to be rootes and have zero root edge
    
    ## trees1, trees2 - nexus lists of posterior tree samples
    ## N - total number of trees in the smaller set (if sizes differ)
    ## MRCA - estimate in tree time units (usually Mya) of the MCRA of two trees
    ## (usually fossil based estimate ect)
    ## MRCA.SE if available (as SE or 95% CI) uncertainty around the MRCA
    
    
    tree1 <- trees1[[sample(1:N, 1)]]
    tree2 <- trees2[[sample(1:N, 1)]]
    
    if (!is.rooted(tree1) |
        !is.rooted(tree2))
      stop("Trees have to be rooted")
    if (any(!is.null(tree1$root.edge), !is.null(tree2$root.edge)))
      stop("Root edges must be NULL")
    
    tree1_length <- nodeheight(tree1, 1)
    tree2_length <- nodeheight(tree2, 1)
    
    if (is.null(MRCA.SE))
      MRCA.SE <- 0
    mrca_current <- rnorm(1, mean = MRCA, sd = MRCA.SE)
    
    tree1$root.edge <- mrca_current - tree1_length
    tree2$root.edge <- mrca_current - tree2_length
    
    out_tree <-
      bind.tree(tree1,
                tree2,
                where = Ntip(tree1) + 1,
                position = mrca_current - tree1_length)
    
    ## make tree ultrameric using negative least squares method
    ## on a cophenetic tree scaffold
    
    out_tree <-
      nnls.tree(cophenetic(out_tree),
                out_tree,
                rooted = T,
                trace = 0)
    return(out_tree)
  }

newtree <-
  bind.trees.mrca(
    trees1 = trees_bird,
    trees2 = trees_mammal,
    N = 1000,
    MRCA = 320,
    MRCA.SE = 7.4
  )