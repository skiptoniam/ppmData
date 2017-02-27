# Species specific offsets - could deliniate the difference between collection method?
# rather than PA - PO data sets - you could generate method... is that meaningful? abundance - pa - po...
# need to generate the same number of observations as n_species * n_backgroundpoints.
n.species <- 20
y <- rep(0,nrow(x)*n.species)
offset <- rep(0,nrow(x)*n.species)
species <- 1:n.species
species.PO <- 1:n.species
species.PA <- sort(sample(1:n.species,5))
n.PA <- length(species.PA)
n.PO <- length(species.PO)
p.sdm <- 2#ncol(sdm.BG.model.matrix) - 1
p.bias <- 1#ncol(bias.BG.model.matrix)
PA.good.rows <- rownames(PA)
quadrat.size <- 1
region.size <- n.bg*10

for(k in 1:n.species) {
  yk <- rep(0,nrow(x))
  yk[1:n.species] <- 1*(1:n.species == k)     # sdm margin only active for species k
  yk[1 + n.species] <- 1*(1 == k)             # bias margin only active for species no. 1

  if(species[k] %in% species.PA) {
    yk[1 + n.species + (1:n.PA)] <- PA[PA.good.rows,species[k]]  # PA sites
  } else {
    yk[1 + n.species + (1:n.PA)] <- NA
  }
  if(species[k] %in% species.PO) {
    yk[1 + n.species + (1:n.bg)] <- 0             # BG sites
  } else {
    yk[1 + n.species + n.PA + (1:n.bg)] <- NA
  }
  yk[1 + n.species + n.sites + (1:(p.sdm + 2 + p.bias))] <- 0  # penalty

  y[(k-1)*nrow(x) + 1:nrow(x)] <- yk

  offk <- rep(0,nrow(x))
  offk[1 + n.species + (1:n.PA)] <- log(quadrat.size)
  offk[1 + n.species + n.PA + (1:n.bg)] <- log(region.size) - log(n.bg)  # BG sites
  offset[(k-1)*nrow(x) + 1:nrow(x)] <- offk
}
