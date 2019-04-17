#'  Simulate data from the prior/likelihood model assumed in SPASIBA
#'  
#' @param ploidy 1 or 2
#' @param n.loc integer nuber of loci
#' @param    flat domai
#' @param size.pop.ref  Matrix with one row per population and one column 
#' per locus giving haploid population sample size
#' @param coloc.ref  TRUE / FALSE indicating whether colocated samples are allowed
#' @param n.unknown  Number of individuls with unknown spatial coordinates
#' @param nx.grid  Number of pixels horizontally
#' @param ny.grid  Number of pixels vertically
#' @param win  c(x.min,x.max,y.min,y.max)
#' @param kappa.geno  spatial scale parameter, hidden Gaussian RF, genotypes
#' @param nu.geno   spatial smoothing parameter, hidden Gaussian RF, genotypes
#' @param sigma2.geno  marginal variance, hidden Gaussian RF, genotypes
#' @param n.quant  number of quantitative covariates
#' @param kappa.quant   spatial scale parameter, quant  variables
#' @param nu.quant  spatial smoothing parameter, quant  variables
#' @param sigma2.gf.quant  marginal variance, quant  variables
#' @param marginal.quant  distribution of quantitative variables ('norm' or 'lnorm')
#' @param tau.noise.quant   precision of distribution of quantitative variables
#' @param store.gf TRUE / FALSE indicating whether Gaussian random field 
#' values used to simulated allele frequencies should be stored
#' @param obsolete  TRUE / FALSE indicating whether random field 
#' simulations should be done with the obsolete function GaussRF or the more recent function RFsimulate
#' 
#' 
#' @return A list with components 
#' \itemize{
#' \item sim.method: The simulation method. The function currently
#'   implements the "geostat" model
#'  \item sphere: A logical indicating whether coordinates are handled as
#'   (Lat,Lon) on the sphere or Euclidean coordinates in a 2D flat domain
#'  \item geno.ref: A matrix containing the genetic data (individual or 
#'  \item population allele counts)  of the reference 
#' samples. This matrix has one line per population (or individual) and one column per locus
#'  \item coord.ref: A matrix containing the geographic coordinates of 
#'  \item reference populations. This matrix has one line per population (or individual) and two columns
#'  \item size.pop.ref: A matrix containing the haploid population size of reference populations. 
#' This matrix has one line per population (or individual) and one column per locus
#'  \item geno.unknown: Indivudual genotypes of individuals of unknown geographic origin
#'  \item true.coord.unknown: True geographic coordinates of individuals assumed here to be of 
#' unknown geographic origin
#'  \item ploidy: 1 or 2 
#'  \item kappa.geno: Scale parameter used in the simulation of the underlying Matern field 
#'  \item nu.geno: Smoothing parameter used in the simulation of the underlying Matern field
#'  \item sigma2.geno: Variance parameter used in the simulation of the underlying Matern field
#'  \item coord.grid=coord: Coordinates of the regular grid at which simulation of 
#' Gaussian fields and allele frequencies is performed
#' }
#' 
#' @export
#' @examples
#' dat <- SPASIBA.sim(ploidy=2,
#'                    n.loc=3, # nb of loci
#'                    sphere=FALSE,
#'                    size.pop.ref=rep(2,30), # haploid pop size at each site
#'                    n.unknown=20, # nb indiv with unknown origin
#'                    nx.grid=100, # nb hor. pixels on grid 
#'                    ny.grid=100, # nb vert.  pixels on grid 
#'                    win=c(0,10,40,50),
#'                    kappa.geno=1/10,# spatial scale parameter
#'                    nu.geno=1,# spatial smoothing parameter
#'                    sigma2.geno=1# marginal variance of hidden Gaussian RF
#'                    )
#'                    summary(dat)
SPASIBA.sim <-
function (ploidy = 2, n.loc = 10, n.quant = 0, sphere = FALSE, 
    size.pop.ref, coloc.ref = FALSE, n.unknown, nx.grid, ny.grid, 
    win, kappa.geno = NULL, nu.geno = 1, sigma2.geno = NULL, 
    kappa.quant = NULL, nu.quant = 1, sigma2.gf.quant = NULL, 
    marginal.quant = "lnorm", tau.noise.quant = NULL, store.gf = FALSE, 
    obsolete = FALSE) {
    cp = colorRampPalette(c("blue", "cyan", "yellow", "red"))
    n.site.ref <- length(size.pop.ref)
    if (!(marginal.quant %in% c("norm", "lnorm"))) 
        stop("marginal.quant should be norm or lnorm")
    if (!sphere) {
        print("Simulation in R^2")
        u <- rep(seq(win[1], win[2], l = nx.grid))
        v <- rep(seq(win[3], win[4], l = ny.grid))
        n.grid <- nx.grid * ny.grid
        coord <- matrix(nrow = n.grid, ncol = 2)
        coord[, 1] <- rep(seq(win[1], win[2], l = nx.grid), ny.grid)
        coord[, 2] <- as.vector(matrix(nrow = nx.grid, ncol = ny.grid, 
            byrow = TRUE, rep(seq(win[3], win[4], l = ny.grid), 
                nx.grid)))
        sub.ref.site <- sample(1:n.grid, size = n.site.ref, replace = coloc.ref)
        sub.unknown <- sample(setdiff(1:n.grid, sub.ref.site), 
            size = n.unknown, replace = FALSE)
        themodel = RMwhittle(nu = nu.geno, var = sigma2.geno, 
            scale = 1/kappa.geno)
        if (n.loc > 0) {
            if (obsolete) {
                gf.grid.geno <- GaussRF(x = u, y = v, grid = TRUE, 
                  model = "whittle", param = c(0, sigma2.geno, 
                    0, 1/kappa.geno, nu.geno), n = n.loc)
            }
            else {
                gf.grid.geno <- RFsimulate(x = u, y = v, grid = TRUE, 
                  model = themodel, n = n.loc)
            }
        }
        if (n.quant > 0) {
            if (obsolete) {
                gf.grid.quant <- GaussRF(x = u, y = v, grid = TRUE, 
                  model = "whittle", param = c(0, sigma2.gf.quant, 
                    0, 1/kappa.quant, nu.quant), n = n.quant)
            }
            else {
                gf.grid.quant <- RFsimulate(x = u, y = v, grid = TRUE, 
                  model = themodel, n = n.quant)
            }
        }
        if (n.loc > 0) {
            if (obsolete) {
                gf.geno <- matrix(nrow = n.grid, ncol = n.loc)
                for (l in 1:n.loc) {
                  gf.geno[, l] <- as.vector(gf.grid.geno[, , 
                    l])
                }
            }
            else {
                gf.geno <- gf.grid.geno@data
            }
        }
        if (n.quant > 0) {
            if (obsolete) {
                gf.quant <- matrix(nrow = n.grid, ncol = n.quant)
                for (l in 1:n.quant) {
                  gf.quant[, l] <- as.vector(gf.grid.quant[, 
                    , l])
                }
            }
            else {
                gf.quant <- gf.grid.quant@data
            }
        }
        y <- matrix(nrow = n.grid, ncol = n.loc)
        size.pop <- rep(ploidy, n.grid)
        size.pop[sub.ref.site] <- size.pop.ref
        size.pop[sub.unknown] <- ploidy
        size.pop <- rep(size.pop, n.loc)
        size.pop <- matrix(ncol = n.loc, size.pop)
        for (l in 1:n.loc) {
            y[, l] <- rbinom(n.grid, size = size.pop[, l], prob = 1/(1 + 
                exp(-gf.geno[, l])))
        }
        geno.ref <- y[sub.ref.site, ]
        coord.ref <- coord[sub.ref.site, ]
        size.pop.ref <- size.pop[sub.ref.site, ]
        geno.unknown <- y[sub.unknown, ]
        true.coord.unknown <- coord[sub.unknown, ]
        if (n.quant > 0) {
            quant.ref <- matrix(nrow = n.site.ref, ncol = n.quant, 
                rnorm(n = n.site.ref * n.quant, mean = as.vector(gf.quant[sub.ref.site, 
                  ]), sd = 1/sqrt(tau.noise.quant)))
            quant.unknown <- matrix(nrow = n.unknown, ncol = n.quant, 
                rnorm(n = n.unknown * n.quant, mean = as.vector(gf.quant[sub.unknown, 
                  ]), sd = 1/sqrt(tau.noise.quant)))
            if (marginal.quant == "lnorm") {
                quant.ref <- exp(quant.ref)
                quant.unknown <- exp(quant.unknown)
            }
        }
    }
    else {
        print("Simulation in S^2")
        n.tot <- n.site.ref + n.unknown
        coord.all <- cbind(runif(n.tot, win[1], win[2]), runif(n.tot, 
            win[3], win[4]))
        lonlat3D = function(lon, lat) {
            cbind(cos((lon/180) * pi) * cos((lat/180) * pi), 
                sin((lon/180) * pi) * cos((lat/180) * pi), sin((lat/180) * 
                  pi))
        }
        true.radius.of.earth = 6371
        coord.all.3 <- lonlat3D(coord.all[, 1], coord.all[, 2])
        mesh <- inla.mesh.create(loc = coord.all.3)
        omesh = old.mesh.class(mesh)
        proj <- inla.mesh.projector(mesh, dims = c(nx.grid, ny.grid), 
            xlim = c(win[1], win[2]), ylim = c(win[3], win[4]), 
            projection = "longlat")
        spde = inla.spde.create(mesh, model = "matern")
        tau.geno = 1/(4 * pi * kappa.geno^2 * sigma2.geno)^0.5
        gf.geno <- matrix(nrow = n.tot, ncol = n.loc)
        size.pop <- c(size.pop.ref, rep(ploidy, n.unknown))
        size.pop <- rep(size.pop, n.loc)
        size.pop <- matrix(ncol = n.loc, size.pop)
        y <- matrix(nrow = n.tot, ncol = n.loc)
        if (n.quant > 0) {
            gf.quant <- matrix(nrow = n.tot, ncol = n.quant)
        }
        for (iloc in 1:n.loc) {
            print(paste("locus ", iloc))
            x <- inla.spde.query(spde, sample = list(tau = tau.geno, 
                kappa2 = kappa.geno^2))$sample
            plotdata.x = inla.mesh.project(proj, x)
            xdat <- numeric(nrow(coord.all.3))
            for (ipop in 1:nrow(coord.all.3)) {
                xc <- coord.all[ipop, 1]
                yc <- coord.all[ipop, 2]
                ii <- which.min(abs(xc - seq(win[1], win[2], 
                  length = nx.grid)))
                jj <- which.min(abs(yc - seq(win[3], win[4], 
                  length = ny.grid)))
                gf.geno[ipop, iloc] <- plotdata.x[ii, jj]
            }
            y[, iloc] <- rbinom(n.tot, size = size.pop[, iloc], 
                prob = 1/(1 + exp(-gf.geno[, iloc])))
        }
        if (n.quant > 0) {
            tau.quant = 1/(4 * pi * kappa.quant^2 * sigma2.gf.quant)^0.5
            for (iquant in 1:n.quant) {
                print(paste("locus ", iquant))
                x <- inla.spde.query(spde, sample = list(tau = tau.quant, 
                  kappa2 = kappa.quant^2))$sample
                plotdata.x = inla.mesh.project(proj, x)
                xdat <- numeric(nrow(coord.all.3))
                for (ipop in 1:nrow(coord.all.3)) {
                  xc <- coord.all[ipop, 1]
                  yc <- coord.all[ipop, 2]
                  ii <- which.min(abs(xc - seq(win[1], win[2], 
                    length = nx.grid)))
                  jj <- which.min(abs(yc - seq(win[3], win[4], 
                    length = ny.grid)))
                  gf.quant[ipop, iquant] <- plotdata.x[ii, jj]
                }
            }
        }
        coord.ref <- coord.all[1:n.site.ref, ]
        true.coord.unknown <- coord.all[-(1:n.site.ref), ]
        geno.ref <- y[1:n.site.ref, ]
        geno.unknown <- y[-(1:n.site.ref), ]
        size.pop.ref <- size.pop[1:n.site.ref, ]
        if (n.quant > 0) {
            if (marginal.quant == "norm") {
                quant.ref <- matrix(nrow = n.site.ref, ncol = n.quant, 
                  rnorm(n = n.site.ref * n.quant, mean = as.vector(gf.quant[1:n.site.ref, 
                    ]), sd = 1/sqrt(tau.noise.quant)))
                quant.unknown <- matrix(nrow = n.unknown, ncol = n.quant, 
                  rnorm(n = n.unknown * n.quant, mean = as.vector(gf.quant[-(1:n.site.ref), 
                    ]), sd = 1/sqrt(tau.noise.quant)))
            }
            if (marginal.quant == "lnorm") {
                quant.ref <- exp(quant.ref)
                quant.unknown <- exp(quant.unknown)
            }
        }
    }
    res <- list(sim.method = "geostat", geno.ref = geno.ref, 
        coord.ref = coord.ref, sphere = sphere, size.pop.ref = size.pop.ref, 
        geno.unknown = geno.unknown, true.coord.unknown = true.coord.unknown, 
        ploidy = ploidy, kappa.geno = kappa.geno, nu.geno = nu.geno, 
        sigma2.geno = sigma2.geno, coord.grid = coord)
    if (n.quant > 0) {
        res <- append(res, list(quant.ref = quant.ref, quant.unknown = quant.unknown, 
            kappa.quant = kappa.quant, nu.quant = nu.quant, sigma2.gf.quant = sigma2.gf.quant, 
            marginal.quant = marginal.quant, tau.noise.quant = tau.noise.quant))
    }
    if (!sphere & store.gf) {
        res <- append(res, list(u = u, v = v, gf.geno = gf.geno, 
            gf.grid.geno = gf.grid.geno, sub.ref.site = sub.ref.site))
        if (n.quant > 0) {
            res <- append(res, list(gf.quant = gf.quant))
        }
    }
    res
}
