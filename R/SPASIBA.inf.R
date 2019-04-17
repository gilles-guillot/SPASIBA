#' @title Inference, prediction for Continuous Spatial Assignment with INLA
#' 
#' 
#' @param geno.ref: Matrix of allele counts of reference data. One row per population, one column per locus.
#' @param ploidy: 1/2
#' @param coord.ref: coordinates of reference sampling sites. One row per sampling site, two columns (cartesian or lon-lat)
#' @param sphere: Logical (\code{TRUE} or \code{FALSE}) indicating whether samples are living on a sphere or a flat domain
#' @param size.pop.ref: Matrix with one row per population and one column per locus giving 
#' haploid population sample sizes of the various reference populations
#' @param geno.unknown: Genotypes  of individuals of unknwon geographic origin. 
#' One row per indivdual, one column per locus. This should contain allele counts of an 
#' arbitrary reference allele at each locus (hence 0,1, 2 or NA)
#' @param quant.ref: Matrix of quantitative covariates for reference individuals. One row per individual. Currently not implemented.
#' @param quant.unknown: Matrix of quantitative covariates for individuals of unknwon geographic origin. Currently not implemented.
#' @param make.inf: Logical (\code{TRUE} or \cite{FALSE}) indicating whether inference of covariance function (under the GMRF-SPDE model) should be carried out
#' @param loc.infcov: Subset of loci to be used to carry out inference of covariance function (GMRF-SPDE model)
#' @param kappa.geno.plug: Spatial scale parameter kappa of the hidden Gaussian random field 
#' (if known from inference carried out with the present function earlier)
#' @param tau.geno.plug: Precision (inverse of variance) of the hidden Gaussian random field 
#' (if known from inference carried out earlier with the present function)
#' @param window: A 2x2 matrix specifying the xmin, xmax, ymin, ymax of the rectangular window over which predicted allele frequencies 
#' maps will be computed. As default, the smallest rectangular window with edges parrallel to the axes will be used.
#' @param max.edge.mesh: Max length of an edge in triangulation for GMRF-SPDE model. Change default value only if you are familiar with 
#' the triangulation function for the SPDE model in INLA. 
#' @param offset.mesh: Offset in mesh for  GMRF-SPDE model. Change default value only if you are familiar with 
#' the triangulation function for the SPDE model in INLA.

#' @param marginal.quant: Distribution of quantitative variables (norm or lnorm). Currently not implemented. 
#' @param make.pred: TRUE/FALSE indicating whether prediction of allele  frequencies should be carried out
#' @param nx.pred: Number of pixels horizontally in allele frequencies prediction
#' @param ny.pred: Number of pixels vertically in allele frequencies prediction
#' @param make.assign: TRUE/FALSE indicating whether assignment of individuals of unknown origin should be made
#' @param proj.predref: Projection of grid used in prediction (typically as returned by a previous call to the present function when make.pred=TRUE)
#' @param freq.predref: Allele frequencies (returned by this function when make.pred=TRUE)  
#' @param verbose.inf: TRUE if you want to get the \code{inla()}  blurb
#' @param verbose.pred: TRUE if you want to get the \code{inla()} blurb
#'
#' @return A list with components below (some of them being only optionally returned). 
#' The most important in the list returned are \code{ll} and \code{coord.unknown.est}. 
#' \itemize{
#' \item mesh.ref: An object of class \code{inla.mesh} containing info
#'   about the mesh used under the inference of the GMRF-SPDE
#'   parameters. See \code{?inla.mesh.2d} for details.
#' \item cpu.used.infcov: CPU time used during the inference of the GMRF-SPDE
#'   parameters.
#' \item tau.geno.est: Precision (inverse of the variance) in the Matern
#'   field involved in the GMRF-SPDE model
#' \item kappa.geno.est: Scale parameter in the Matern
#'   field involved in the GMRF-SPDE model
#' \item make.inf: Logical (\code{TRUE} or \code{FALSE}) indicating
#'   whether the last call to the present function requested to estimate
#'   parameters of the GMRF-SPDE model.
#' \item coord.pred: Coordinates of the regular grid at which allele
#'   frequencies are predicted.
#' \item mesh.predref: Mesh (an object of class \code{inla.mesh}) used
#'   involved in the prediction of allele   frequencies. See \code{?inla.mesh.2d} for details.
#' \item proj.predref: An object typically returned by
#'   \code{inla.mesh.project} used to convert info between an irregular
#'   trangulation and a regular grid. See   \code{?inla.mesh.projector} for details.
#' \item u: The \code{x} component of a \code{proj.predref} object 
#' \item v: The \code{y} component of a \code{proj.predref} object
#' \item cpu.used.predfreq.perloc: CPU time spent computing allele
#'   frequency maps
#' \item freq.predref: An array with dimensions \code{(nx.pred,ny.pred,nloc)} containing
#'   predicted frequency maps, with \code{(nx.pred,ny.pred)} dimension of the
#'   grid, \code{nloc}: number of loci.
#' \item ll: An array with dimensions \code{(nx.pred,ny.pred,n.unknown)}
#'   containing the log-likelihood of each individual to have origin at the
#'   various nodes of a regular grid. \code{(nx.pred,ny.pred)} dimension of the
#'   grid, \code{n.unknown}: number of individuals of unknown origin.
#' \item coord.unknown.est: The estimated geographic coordinates of the individuals.
#' }
#'
#' @author Gilles Guillot
#' @export
#'
#' 
SPASIBA.inf <-
function (geno.ref, ploidy, coord.ref, sphere = FALSE, size.pop.ref, 
    geno.unknown, quant.ref = NULL, quant.unknown = NULL, make.inf = FALSE, 
    loc.infcov = NULL, kappa.geno.plug = NULL, tau.geno.plug = NULL, 
    window = NULL, max.edge.mesh = NULL, offset.mesh = -0.1, 
    marginal.quant = "none", make.pred = FALSE, nx.pred = 20, 
    ny.pred = 20, make.assign = FALSE, proj.predref = NULL, freq.predref = NULL, 
    verbose.inf = TRUE, verbose.pred = TRUE) 
{
    lonlat3D = function(lon, lat) {
        cbind(cos((lon/180) * pi) * cos((lat/180) * pi), sin((lon/180) * 
            pi) * cos((lat/180) * pi), sin((lat/180) * pi))
    }
    if (!(marginal.quant %in% c("none", "norm", "lnorm"))) 
        stop("marginal.quant should be none, norm or lnorm")
    mesh.predref <- omesh.predref <- NULL
    use.quant <- !is.null(quant.ref)
    use.geno <- !is.null(geno.ref)
    final.res <- NULL
    true.radius.of.earth = 6371
    if (is.null(window)) {
        xmin <- min(coord.ref[, 1])
        xmax <- max(coord.ref[, 1])
        ymin <- min(coord.ref[, 2])
        ymax <- max(coord.ref[, 2])
    }
    else {
        xmin <- min(window[, 1])
        xmax <- max(window[, 1])
        ymin <- min(window[, 2])
        ymax <- max(window[, 2])
    }
    if (is.null(max.edge.mesh)) {
        max.edge.mesh = 0.1 * max(xmax - xmin, ymax - ymin)
    }
    if (use.geno) {
        n.loc <- ncol(geno.ref)
        if (is.null(loc.infcov)) {
            loc.infcov <- 1:n.loc
        }
        n.loc.infcov <- length(loc.infcov)
        n.unknown <- ifelse(is.null(geno.unknown), 0, nrow(geno.unknown))
    }
    if (use.quant) {
        n.quant <- ncol(quant.ref)
        n.unknown <- ifelse(is.null(quant.unknown), 0, nrow(quant.unknown))
    }
    if (use.geno) {
        if (!make.inf) {
            print("PLUGGING COVARIANCE PARAMETERS")
            kappa.geno.est <- kappa.geno.plug
            tau.geno.est <- tau.geno.plug
        }
        else {
            if (sphere) {
                coord.tmp <- lonlat3D(coord.ref[, 1], coord.ref[, 
                  2])
                window <- lonlat3D(window[, 1], window[, 2])
            }
            else {
                coord.tmp <- coord.ref
            }
            print("")
            print("Estimating covariance parameters")
            mesh.ref <- inla.mesh.2d(loc = coord.tmp, loc.domain = window, 
                offset = offset.mesh, max.edge = max.edge.mesh, 
                cutoff = 0)
            spde.infcov = inla.spde.create(mesh.ref, model = "matern", 
                param = list(alpha = 2))
            if (!sphere) {
                data.df <- data.frame(xcoord = rep(coord.tmp[, 
                  1], n.loc.infcov), ycoord = rep(coord.tmp[, 
                  2], n.loc.infcov), allele.counts = as.vector(geno.ref[, 
                  loc.infcov]), size = as.vector(size.pop.ref[, 
                  loc.infcov]), spatial = rep(mesh.ref$idx$loc, 
                  n.loc.infcov))
            }
            else {
                data.df <- data.frame(xcoord = rep(coord.tmp[, 
                  1], n.loc.infcov), ycoord = rep(coord.tmp[, 
                  2], n.loc.infcov), zcoord = rep(coord.tmp[, 
                  3], n.loc.infcov), allele.counts = as.vector(geno.ref[, 
                  loc.infcov]), size = as.vector(size.pop.ref[, 
                  loc.infcov]), spatial = rep(mesh.ref$idx$loc, 
                  n.loc.infcov))
            }
            r = rep(1:n.loc.infcov, each = nrow(coord.ref))
            formula = (allele.counts ~ -1 + as.factor(r) + f(spatial, 
                model = spde.infcov, replicate = r))
            result.infcov <- (inla(formula, family = "binomial", 
                Ntrials = as.vector(size.pop.ref[, loc.infcov]), 
                data = data.df, control.family = list(link = "logit"), 
                control.compute = list(return.marginals = FALSE), 
                control.inla = list(tolerance = 1e-05), verbose = verbose.inf))
            tau.geno.est = exp(result.infcov$summary.hyperpar[1, 
                "mean"])
            kappa.geno.est = exp(result.infcov$summary.hyperpar[2, 
                "mean"]/2)
            final.res <- append(final.res, list(mesh.ref = mesh.ref, 
                cpu.used.infcov = result.infcov$cpu.used, tau.geno.est = tau.geno.est, 
                kappa.geno.est = kappa.geno.est, make.inf = make.inf))
        }
    }
    if (make.pred) {
        cpu.used.predfreq.perloc <- matrix(nrow = ncol(geno.ref), 
            ncol = 4)
        print("")
        print("Mapping allele frequencies")
        n.pred <- nx.pred * ny.pred
        coord.pred <- matrix(nrow = n.pred, ncol = 2)
        coord.pred[, 1] <- rep(seq(xmin, xmax, l = nx.pred), 
            ny.pred)
        coord.pred[, 2] <- as.vector(matrix(nrow = nx.pred, ncol = ny.pred, 
            byrow = TRUE, rep(seq(ymin, ymax, l = ny.pred), nx.pred)))
        final.res <- append(final.res, list(coord.pred = coord.pred))
        coord.predref <- rbind(coord.pred, coord.ref)
        if (!sphere) {
            coord.predref.tmp <- coord.predref
        }
        else {
            coord.predref.tmp <- lonlat3D(coord.predref[, 1], 
                coord.predref[, 2])
            window <- lonlat3D(window[, 1], window[, 2])
        }
        print("Starting to creat mesh")
        mesh.predref <- inla.mesh.2d(loc = coord.predref.tmp, 
            loc.domain = window, offset = offset.mesh, max.edge = max.edge.mesh, 
            cutoff = 0)
        print("Mesh creating completed")
        if (use.geno) {
            size.pop.pred <- matrix(nrow = n.pred, ncol = n.loc, 
                2)
            y.pred <- matrix(nrow = n.pred, ncol = n.loc)
            y.predref <- rbind(y.pred, geno.ref)
            size.pop.predref <- rbind(size.pop.pred, size.pop.ref)
            spde.predref <- inla.spde2.matern(mesh.predref, B.tau = matrix(log(tau.geno.est), 
                1, 1), B.kappa = matrix(log(kappa.geno.est), 
                1, 1))
            proj.predref <- inla.mesh.projector(mesh.predref, 
                dims = c(nx.pred, ny.pred), xlim = c(xmin, xmax), 
                ylim = c(ymin, ymax))
            u <- proj.predref$x
            v <- proj.predref$y
            formula = (allele.counts ~ 1 + f(spatial, model = spde.predref))
            for (iloc in 1:n.loc) {
                print(paste("locus ", iloc))
                if (!sphere) {
                  data.df <- data.frame(xcoord = coord.predref.tmp[, 
                    1], ycoord = coord.predref.tmp[, 2], allele.counts = y.predref[, 
                    iloc], size = size.pop.predref[, iloc], spatial = mesh.predref$idx$loc)
                }
                else {
                  data.df <- data.frame(xcoord = coord.predref.tmp[, 
                    1], ycoord = coord.predref.tmp[, 2], zcoord = coord.predref.tmp[, 
                    3], allele.counts = y.predref[, iloc], size = size.pop.predref[, 
                    iloc], spatial = mesh.predref$idx$loc)
                }
                result.predref <- (inla(formula, family = "binomial", 
                  Ntrials = size.pop.predref[, iloc], data = data.df, 
                  control.family = list(link = "logit"), control.compute = list(return.marginals = FALSE), 
                  control.inla = list(tolerance = 1e-05), verbose = verbose.pred))
                cpu.used.predfreq.perloc[iloc, ] <- result.predref$cpu.used
                n.mesh.spde.predref <- length(result.predref$summary.random$spatial[, 
                  "mean"])
                if (iloc == 1) {
                  tmp <- inla.mesh.project(proj.predref, result.predref$summary.random$spatial[, 
                    "mean"])
                  gf.geno.predref <- freq.predref <- array(dim = c(nrow(tmp), 
                    ncol(tmp), n.loc))
                  gf.geno.predref[, , 1] <- tmp
                }
                else {
                  gf.geno.predref[, , iloc] <- inla.mesh.project(proj.predref, 
                    result.predref$summary.random$spatial[, "mean"])
                }
                freq.predref[, , iloc] <- 1/(1 + exp(-(as.numeric(result.predref$summary.fixed[1]) + 
                  gf.geno.predref[, , iloc])))
            }
        }
        if (use.quant) {
            spde.predref <- inla.spde.create(mesh.predref, model = "matern", 
                param = list(alpha = 2))
            proj.predref <- inla.mesh.projector(mesh.predref, 
                dims = c(nx.pred, ny.pred), xlim = c(xmin, xmax), 
                ylim = c(ymin, ymax))
            u <- proj.predref$x
            v <- proj.predref$y
            event <- rep(1, nrow(coord.predref.tmp))
            if (marginal.quant == "norm") {
                formula = (quant ~ 1 + f(spatial, model = spde.predref))
                family.quant <- "gaussian"
            }
            if (marginal.quant == "lnorm") {
                formula <- (inla.surv(quant, event) ~ 1 + f(spatial, 
                  model = spde.predref))
                family.quant <- "lognormal"
            }
            mean.quant <- sd.noise.quant <- rep(NA, n.quant)
            for (iquant in 1:n.quant) {
                print(paste("quantitative variable ", iquant))
                if (!sphere) {
                  data.df <- data.frame(xcoord = coord.predref.tmp[, 
                    1], ycoord = coord.predref.tmp[, 2], quant = c(rep(NA, 
                    n.pred), quant.ref[, iquant]), event = event, 
                    spatial = mesh.predref$idx$loc)
                }
                else {
                  data.df <- data.frame(xcoord = coord.predref.tmp[, 
                    1], ycoord = coord.predref.tmp[, 2], zcoord = coord.predref.tmp[, 
                    3], quant = c(rep(NA, n.pred), quant.ref[, 
                    iquant]), event = event, spatial = mesh.predref$idx$loc)
                }
                result.predref <- (inla(formula, family = family.quant, 
                  data = data.df, verbose = verbose.pred))
                mean.quant[iquant] <- result.predref$summary.fixed[, 
                  "mean"]
                sd.noise.quant[iquant] <- inla.emarginal(function(x) sqrt(1/x), 
                  result.predref$marginals.hyperpar[[1]])
                if (iquant == 1) {
                  tmp <- inla.mesh.project(proj.predref, result.predref$summary.random$spatial[, 
                    "mean"])
                  gf.quant.predref <- array(dim = c(nrow(tmp), 
                    ncol(tmp), n.quant))
                  gf.quant.predref[, , 1] <- mean.quant[iquant] + 
                    tmp
                }
                else {
                  gf.quant.predref[, , iquant] <- (mean.quant[iquant] + 
                    inla.mesh.project(proj.predref, result.predref$summary.random$spatial[, 
                      "mean"]))
                }
            }
        }
        final.res <- append(final.res, list(mesh.predref = mesh.predref, 
            proj.predref = proj.predref, u = u, v = v, cpu.used.predfreq.perloc = cpu.used.predfreq.perloc))
        if (use.geno) {
            final.res <- append(final.res, list(freq.predref = freq.predref))
        }
        if (use.quant) {
            final.res <- append(final.res, list(gf.quant.predref = gf.quant.predref, 
                mean.quant = mean.quant, sd.noise.quant = sd.noise.quant))
        }
    }
    if (make.assign) {
        print("")
        print(" Assigning individuals of unknown origin")
        node.opt <- matrix(nrow = n.unknown, ncol = 2)
        ll <- array(dim = c(nrow(freq.predref[, , 1]), ncol(freq.predref[, 
            , 1]), n.unknown), data = 0)
        for (iindiv in 1:n.unknown) {
            print(paste("indiv ", iindiv))
            if (use.geno) {
                for (iloc in 1:n.loc) {
                  gg <- geno.unknown[iindiv, iloc]
                  if (!is.na(gg)) {
                    ll[, , iindiv] <- ll[, , iindiv] + dbinom(x = gg, 
                      size = ploidy, prob = as.vector(freq.predref[, 
                        , iloc]), log = TRUE)
                  }
                }
            }
            if (use.quant) {
                for (iquant in 1:n.quant) {
                  qq <- quant.unknown[iindiv, iquant]
                  if (marginal.quant == "norm") {
                    if (!is.na(q)) {
                      ll[, , iindiv] <- ll[, , iindiv] + dnorm(x = qq, 
                        mean = gf.quant.predref[, , iquant], 
                        sd = sd.noise.quant[iquant], log = TRUE)
                    }
                  }
                  if (marginal.quant == "lnorm") {
                    if (!is.na(q)) {
                      ll[, , iindiv] <- ll[, , iindiv] + dlnorm(x = qq, 
                        meanlog = gf.quant.predref[, , iquant], 
                        sdlog = sd.noise.quant[iquant], log = TRUE)
                    }
                  }
                }
            }
            mm <- max(ll[, , iindiv])
            node.opt[iindiv, ] <- which(ll[, , iindiv] == mm, 
                arr.ind = TRUE)[1, ]
        }
        coord.unknown.est <- cbind(proj.predref$x[node.opt[, 
            1]], proj.predref$y[node.opt[, 2]])
        final.res <- append(final.res, list(ll = ll, coord.unknown.est = coord.unknown.est))
    }
    return(final.res)
}
