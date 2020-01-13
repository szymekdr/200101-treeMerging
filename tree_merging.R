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

#### load data ----

data <- read.table("191216_individuality.csv", sep = ";", head = T)
summary(data)

plot_data <- data[, 4:20]

data$bodymass_S <- scale(data$bodymass, center = T, scale = F)
data$ds_S <- scale(data$ds, center = T, scale = F)
data$HSest_S <- scale(data$HSest, center = T, scale = F)
data$longevity_S <- scale(data$longevity, center = T, scale = F)

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


#### visualize trait ----

ggtree(newtree) +
  geom_treescale(0,-1) +
  geom_tiplab(cex = 2) +
  xlim_tree(400)

newtree <- groupClade(newtree, c(160, 161), group_name = "Group")

master_plot <- ggtree(newtree) +
  #groupClade(newtree, list(birds = newtree$tip.label[1:59], mammals = newtree$tip.label[60:159]), "Group") +
  # aes(color = Group) +
  #geom_treescale(0, -1) +
  geom_tiplab(cex = 2) +
  xlim_tree(400) # +
# scale_color_manual(values = c("darkred", "darkblue"))

plot(master_plot)

large_plot <-
  master_plot %<+% plot_data[, c("species", "context")] + geom_tippoint(aes(color = context)) +
  scale_color_manual(
    values = c(
      "lightgreen",
      "coral",
      "dodgerblue",
      "lightblue",
      "yellow3",
      "orange",
      "black",
      "turquoise3",
      "mediumorchid1",
      "slateblue1",
      "darkgreen"
    ),
    na.value = "gray93"
  )

plot(large_plot)

final_plot <- large_plot +
  geom_facet(
    panel = "Data",
    data = plot_data[, c("species", "ds")],
    geom = ggstance::geom_barh,
    aes(x = ds, color = context, fill = context),
    stat = "identity",
    width = 0.6
  ) +
  theme_tree2()
plot(final_plot)

final_plot <- large_plot +
  geom_facet(
    panel = "Individuality (DS)",
    data = plot_data[, c("species", "HSest")],
    geom = ggstance::geom_barh,
    aes(x = HSest, fill = context),
    stat = "identity",
    width = 0.6
  ) +
  theme_tree2() +
  scale_fill_manual(
    values = c(
      "lightgreen",
      "coral",
      "dodgerblue",
      "lightblue",
      "yellow3",
      "orange",
      "black",
      "turquoise3",
      "mediumorchid1",
      "slateblue1",
      "darkgreen"
    ),
    na.value = "gray93"
  )

plot(final_plot)


#### MCMCglmm analysis ----

data$species2 <- data$species
levels(data$groupliving2) <- c("z_notknown", "gr", "sol")
data$groupliving2 <- relevel(data$groupliving2, ref = "gr")

phylo_A <- inverseA(newtree, nodes = "ALL")

model1 <- MCMCglmm(ds_S ~ bodymass_S + longevity_S + context + groupliving2,
                   data = na.omit(data[,c("ds_S", "bodymass_S", "longevity_S", "context", "groupliving2", "species", "species2")]),
                   random = ~ species + species2,
                   ginverse = list(species2 = phylo_A$Ainv),
                   prior = list(R = list(V = 1, nu = 0.002),
                                G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4),
                                         G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4))),
                   nitt = 200000, thin = 150, burnin = 50000)

autocorr.plot(model1$Sol)
autocorr.plot(model1$VCV)
plot(model1$VCV)
summary(model1)


model2 <- MCMCglmm(ds_S ~ bodymass_S + longevity_S + context + groupliving2,
                   data = na.omit(data[,c("ds_S", "bodymass_S", "longevity_S", "context", "groupliving2", "species", "species2")]),
                   random = ~ species,
                   prior = list(R = list(V = 1, nu = 0.002),
                                G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4))),
                   nitt = 200000, thin = 150, burnin = 50000)

autocorr.plot(model2$Sol)
autocorr.plot(model2$VCV)

summary(model2)


model3 <- MCMCglmm(ds_S ~ bodymass_S + context + groupliving2,
                   data = na.omit(data[,c("ds_S", "bodymass_S", "context", "groupliving2", "species", "species2")]),
                   random = ~ species + species2,
                   ginverse = list(species2 = phylo_A$Ainv),
                   prior = list(R = list(V = 1, nu = 0.002),
                                G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4),
                                         G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4))),
                   nitt = 200000, thin = 150, burnin = 50000)

autocorr.plot(model3$Sol)
autocorr.plot(model3$VCV)
plot(model3$VCV)
summary(model3)


model4 <- MCMCglmm(ds_S ~ bodymass_S + context + groupliving2,
                   data = na.omit(data[,c("ds_S", "bodymass_S", "context", "groupliving2", "species", "species2")]),
                   random = ~ species,
                   prior = list(R = list(V = 1, nu = 0.002),
                                G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4))),
                   nitt = 200000, thin = 150, burnin = 50000)

autocorr.plot(model4$Sol)
autocorr.plot(model4$VCV)

summary(model4)


model3_HS <- MCMCglmm(HSest_S ~ bodymass_S + context + groupliving2,
                   data = na.omit(data[,c("HSest_S", "bodymass_S", "context", "groupliving2", "species", "species2")]),
                   random = ~ species + species2,
                   ginverse = list(species2 = phylo_A$Ainv),
                   prior = list(R = list(V = 1, nu = 0.002),
                                G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4),
                                         G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4))),
                   nitt = 200000, thin = 150, burnin = 50000)

autocorr.plot(model3_HS$Sol)
autocorr.plot(model3_HS$VCV)
plot(model3_HS$VCV)
summary(model3_HS)


model4_HS <- MCMCglmm(HSest_S ~ bodymass_S + context + groupliving2,
                   data = na.omit(data[,c("HSest_S", "bodymass_S", "context", "groupliving2", "species", "species2")]),
                   random = ~ species,
                   prior = list(R = list(V = 1, nu = 0.002),
                                G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1e4))),
                   nitt = 200000, thin = 150, burnin = 50000)

autocorr.plot(model4_HS$Sol)
autocorr.plot(model4_HS$VCV)

summary(model4_HS)
