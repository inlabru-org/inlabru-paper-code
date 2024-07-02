# Script to reproduce the results in the
# inlabru paper.
#########################################

# Load required packages, and install missing packages
required_packages <- list(
  "inlabru" = "2.10.1.9012",
  "sf" = "1.0.16",
  "ggplot2" = "3.5.1",
  "patchwork" = "1.2.0",
  "dplyr" = "1.1.4",
  "spdep" = "1.3.5",
  "fmesher" = "0.1.6.9003",
  "INLA" = "24.06.27",
  "progress" = "1.2.3"
)
for (req_pkg in names(required_packages)) {
  if (!requireNamespace(req_pkg, quietly = TRUE) ||
      (utils::packageVersion(req_pkg) < required_packages[[req_pkg]])) {
    if (req_pkg == "INLA") {
      if (Sys.info()["sysname"] %in% c("Windows", "Darwin")) {
        install.packages("INLA",
          repos = c(
            INLA = "https://inla.r-inla-download.org/R/testing",
            getOption("repos")
          ),
          dependencies = TRUE,
          type = "binary"
        )
      } else {
        install.packages("INLA",
          repos = c(
            INLA = "https://inla.r-inla-download.org/R/testing",
            getOption("repos")
          ),
          dependencies = TRUE
        )
      }
    } else if ((req_pkg %in% c("inlabru", "fmesher")) &&
               !is.na(base::package_version(required_packages[[req_pkg]])[1, 4])) {
      install.packages(req_pkg,
                       repos = c(
                         inlabru_universe = "https://inlabru-org.r-universe.dev",
                         getOption("repos")
                       ),
                       dependencies = TRUE
      )
    } else {
      install.packages(req_pkg)
    }
  }
  library(package = req_pkg, character.only = TRUE)

  if (utils::packageVersion(req_pkg) < required_packages[[req_pkg]]) {
    warning(
      paste0(
        "Detected '",
        req_pkg,
        "' version ", utils::packageVersion(req_pkg),
        " which is older than the version used for code testing, version '",
        required_packages[[req_pkg]],
        "'."
      )
    )
  }
}


# create cache to save fitted model objects
# and method checking results
if (!dir.exists(here::here("R_cache"))) {
  dir.create(here::here("R_cache"), recursive = TRUE)
}

# create data if it does not already exist
if (!file.exists(here::here("Data", "aggregation_data.rds"))) {
  source(here::here("R", "spatial_example_data.R"))
}


# Section 1.1:  A Simple Example
################################

# Load data
data(toypoints, package = "inlabru")
point_data <- toypoints$points # sf point data
mesh <- toypoints$mesh # SPDE mesh
prediction_locs <- toypoints$pred_locs # sf point prediction locations

# Define model components

# Spatially structured random effect
matern <- inla.spde2.pcmatern(mesh,
  prior.range = c(2, 0.05),
  prior.sigma = c(2, 0.05)
)

cmp <- ~ Intercept(1) +
  grf(
    main = geometry,
    model = matern
  )

# Construct likelihood
lik <- like(
  formula = z ~ .,
  data = point_data,
  family = "gaussian"
)

# Fit model
fit <- bru(lik,
  components = cmp
)

# Make predictions
predictions <- predict(fit,
  prediction_locs,
  formula = ~ Intercept + grf
)

# Plot data and predictions
g1 <- ggplot() +
  geom_sf(
    data = point_data,
    aes(colour = z)
  ) +
  scale_colour_viridis_c() +
  theme_classic()

g2 <- ggplot() +
  gg(
    data = predictions,
    aes(fill = mean),
    geom = "tile"
  ) +
  scale_fill_viridis_c() +
  theme_classic()

# Figure 1:
############

# \textwidth as returned by the following:
# \usepackage{layouts}
# ...
# \printinunitsof{cm}\prntlen{\textwidth}
# 15.55528 cm
tw <- 15.55528
tw_in <- tw / 2.54 # In inches, for pdf()

png(
  filename = here::here("Figures", "intro_example.png"),
  width = tw,
  height = tw / 2,
  units = "cm",
  res = 300
)

g1 + g2 +
  plot_annotation(tag_levels = "A")

dev.off()

pdf(
  file = here::here("Figures", "intro_example.pdf"),
  width = tw_in,
  height = tw_in / 2 # ,
  #  units = "cm",
  #  res = 300
)

g1 + g2 +
  plot_annotation(tag_levels = "A")

dev.off()

g1 + g2 +
  plot_annotation(tag_levels = "A")

# Section 4: Approximate Bayesian Method Checking
#################################################

# Helper function for the KS plot checking for uniformity:
plot_ks <- function(values, mapping = aes(), diff = TRUE) {
  df <- data.frame(values = values) %>%
    dplyr::mutate(values_cdf = (rank(values) - 0.5) / length(values))
  circle_df <- data.frame(x = seq(0, 1, length = 2000))
  if (diff) {
    circle_df <-
      circle_df %>%
      dplyr::mutate(
        lower = -2 * sqrt(x * (1 - x)),
        upper = +2 * sqrt(x * (1 - x))
      )
  } else {
    circle_df <-
      circle_df %>%
      dplyr::mutate(
        lower = x - 2 * sqrt(x * (1 - x) / nrow(df)),
        upper = x + 2 * sqrt(x * (1 - x) / nrow(df))
      )
  }
  test <- ks.test(df$values, punif)
  if (diff) {
    geom <-
      list(
        ggplot2::geom_abline(slope = 0, intercept = 0, col = "grey", data = NULL),
        ggplot2::geom_line(
          data = circle_df,
          mapping = utils::modifyList(aes(x, lower), mapping),
          col = "gray", linewidth = 1
        ),
        ggplot2::geom_line(
          data = circle_df,
          mapping = utils::modifyList(aes(x, upper), mapping),
          col = "gray", linewidth = 1
        ),
        ggplot2::geom_line(
          mapping = utils::modifyList(
            aes(
              x = values,
              y = (values_cdf - values) * sqrt(nrow(df)),
              lty = "ECDF"
            ),
            mapping
          ),
          data = df
        ),
        ggplot2::geom_line(
          mapping = utils::modifyList(
            aes(
              x = values,
              y = -max(abs(values_cdf - values)) * sqrt(nrow(df)),
              lty = "KS statistic"
            ),
            mapping
          ),
          data = df
        ),
        ggplot2::geom_line(
          mapping = utils::modifyList(
            aes(
              x = values,
              y = max(abs(values_cdf - values)) * sqrt(nrow(df)),
              lty = "KS statistic"
            ),
            mapping
          ),
          data = df
        ),
        ggplot2::xlab("Values"),
        ggplot2::ylab("Scaled ECDF values")
        # ggplot2::ggtitle(paste0("K-S p-value = ", signif(test$p.value, 4)))
      )
  } else {
    geom <-
      list(
        ggplot2::geom_abline(slope = 1, intercept = 0, col = "grey", data = NULL),
        ggplot2::geom_line(
          data = circle_df,
          mapping = utils::modifyList(aes(x, lower), mapping),
          col = "gray", linewidth = 1
        ),
        ggplot2::geom_line(
          data = circle_df,
          mapping = utils::modifyList(aes(x, upper), mapping),
          col = "gray", linewidth = 1
        ),
        ggplot2::geom_line(
          mapping = utils::modifyList(
            aes(
              x = values,
              y = values_cdf,
              lty = "ECDF"
            ),
            mapping
          ),
          data = df
        ),
        ggplot2::geom_line(
          mapping = utils::modifyList(
            aes(
              x = values,
              y = values - max(abs(values_cdf - values)),
              lty = "KS statistic"
            ),
            mapping
          ),
          data = df
        ),
        ggplot2::geom_line(
          mapping = utils::modifyList(
            aes(
              x = values,
              y = values + max(abs(values_cdf - values)),
              lty = "KS statistic"
            ),
            mapping
          ),
          data = df
        ),
        ggplot2::xlab("Values"),
        ggplot2::ylab("ECDF")
        # ggplot2::ggtitle(paste0("K-S p-value = ", signif(test$p.value, 4)))
      )
  }
  geom
}

# Simulation based calibration:
###############################

n <- 100 # number of data points in each simulation
gam <- 0.5 # gamma value fixed

# Define the non-linear predictor function
lam_param <- function(u, gam) {
  bru_forward_transformation(qexp, u, rate = gam)
}

# number of simulations
n_sims <- 500

# matrix to store results from the simulation
res <- matrix(NA,
  nrow = n_sims,
  ncol = 4
)
colnames(res) <- c(
  "lam",
  "max_y",
  "n_iter",
  "e_cdf"
)

# seeds for reproducibility
set.seed(1042839)
gen_seeds <- sample((1:n_sims * 2), n_sims) # seeds for generate()

if (!file.exists(here::here("R_cache", "method_checking_results.rds")) ||
  !file.exists(here::here("R_cache", "method_checking_fit.rds"))) {
  prog_bar <- progress::progress_bar$new(
    format = "Computing simulations: [:bar] :percent eta: :eta",
    total = n_sims
  )

  for (i in seq_len(n_sims)) {
    prog_bar$tick()

    # simulate data
    u <- rnorm(1)
    lam <- lam_param(u, gam)
    y <- rpois(n, lam)
    df <- data.frame(
      y = y,
      id = seq_len(n)
    )

    # model components and formula
    u_prior <- list(prec = list(initial = 0, fixed = TRUE))
    # Manual transformation:
    #    cmp <- ~ -1 + u(rep(1, n), model = "iid", hyper = u_prior)
    #    fml <- y ~ log(lam_param(u, gam))
    # Automatic transformation:
    cmp <- ~ 0 + lambda(rep(1, n),
      model = "iid", hyper = u_prior,
      marginal = bru_mapper_marginal(qexp, rate = gam)
    )
    fml <- y ~ log(lambda)

    # fit model
    fit <- bru(
      components = cmp,
      formula = fml,
      data = df,
      family = "poisson",
      options = list(
        # bru_verbose = 1,
        bru_max_iter = 20
      )
    )

    # sample from posterior
    n_mc <- 5000
    lam_post <- generate(
      fit,
      # Manual:
      #   formula = ~ lam_param(u, gam)
      # Automatic:
      formula = ~ lambda[1],
      n.samples = n_mc,
      seed = gen_seeds[i]
    )
    lam_post <- as.vector(lam_post)


    # calculate empirical CDF value
    e_cdf <- 1 / n_mc * sum(lam_post < lam) - 1 / (2 * n_mc)

    # things to save
    n_iter <- max(fit$bru_iinla$track$iteration)
    to_save <- c(
      lam,
      max(y),
      n_iter,
      e_cdf
    )

    res[i, ] <- to_save

    if (i == 1) {
      saveRDS(list(fit = fit,
                   df = list(y = y, n = n),
                   true_lambda = lam),
              file = here::here("R_cache", "method_checking_fit.rds"))
    }
  }

  res <- data.frame(res)

  # save result to R_cache
  saveRDS(res, file = here::here("R_cache", "method_checking_results.rds"))
}
res <- readRDS(file = here::here("R_cache", "method_checking_results.rds"))
fit_info <- readRDS(file = here::here("R_cache", "method_checking_fit.rds"))
fit <- fit_info[["fit"]]
fit_df <- fit_info[["df"]]
fit_true_lambda <- fit_info[["true_lambda"]]

# plot K-S test figure
g1 <- ggplot() +
  plot_ks(res$e_cdf,
    mapping = aes(col = "red"),
    diff = FALSE
  ) +
  scale_linetype_manual("",
    values = c(1, 2)
  ) +
  scale_colour_discrete(guide = "none") +
  theme_bw() +
  theme(aspect.ratio = 1)

g2 <- ggplot() +
  plot_ks(res$e_cdf,
    mapping = aes(col = "red"),
    diff = TRUE
  ) +
  scale_linetype_discrete(name = "") +
  scale_colour_discrete(guide = "none") +
  theme_bw() +
  theme(aspect.ratio = 1)

# histogram plot
g3 <- ggplot(res) +
  geom_histogram(
    aes(
      x = e_cdf,
      y = after_stat(density)
    ),
    bins = 15,
    binwidth = 1 / 15,
    boundary = 0,
    colour = "black",
    alpha = 0.8
  ) +
  geom_hline(aes(yintercept = 1),
    linetype = "dashed"
  ) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  xlab("ECDF values")

# Figure 3:
############

pdf(
  file = here::here("Figures", "toy_sim.pdf"),
  width = tw_in,
  height = tw_in / 2.2 # ,
  #  units = "cm",
  #  res = 300
)

g2 + g3 +
  plot_annotation(tag_levels = "A")
dev.off()

g2 + g3 +
  plot_annotation(tag_levels = "A")

# Plot posterior density for the first simulation
# compared to exact density
n.mc.samples <- 20000
lam_post <- generate(fit,
  formula = ~ lambda[1],
  n.samples = n.mc.samples,
  seed = 1232
)

lam_df <- data.frame(lambda = as.vector(lam_post))

# Compute bootstrapped uncertainty for the density estimate
density_boot <- function(x, probs = c(0.05, 0.95), n_boot = 1000,
                         method = c("t", "quantile")) {
  method <- match.arg(method)
  y <- stats::density(x,
                      from = min(x),
                      to = max(x))
  boot_res <- list()
  for (i in seq_len(n_boot)) {
    boot_res[[i]] <- stats::density(sample(x, replace = TRUE),
                                    from = min(x),
                                    to = max(x))
  }
  boot_y <- do.call(cbind, lapply(boot_res, function(z) z$y))
  res <- data.frame(
    x = y$x,
    ymin = numeric(nrow(boot_y)),
    ymid = numeric(nrow(boot_y)),
    ymax = numeric(nrow(boot_y))
  )
  for (i in seq_len(nrow(boot_y))) {
    # Derivation of bias-correcting bootstrap intervals:
    # Under an additive Bootstrap error principle, we have
    #   P( a < true - y < b) = P( a < y - boot < b) = probs gives a and b
    # Then
    #   P( y + a < true < y + b) = probs
    # fulfils the defining property of a confidence interval procedure.
    if (method == "t") {
      m <- mean(y$y[i] - boot_y[i, ])
      s <- sd(y$y[i] - boot_y[i, ])
      dof <- length(boot_y[i, ]) - 1
      res[i, "ymid"] <- y$y[i] + m
      res[i, c("ymin", "ymax")] <- y$y[i] + m + s * qt(probs, df = dof)
    } else {
      q <- quantile(y$y[i] - boot_y[i, ], c(0.5, probs))
      res[i, "ymid"] <- y$y[i] + q[1]
      res[i, c("ymin", "ymax")] <- y$y[i] + q[1] + q[c(2, 3)]
    }
  }
  res
}

set.seed(2435)
lam_boot <- density_boot(lam_df$lambda, probs = c(0.005, 0.995), n_boot = 1000,
                         method = "t")

g_single <- ggplot() +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax),
              data = lam_boot,
              alpha = 0.5, fill = "grey75") +
  geom_density(
    data = lam_df,
    aes(x = lambda, col = "Approximate", linetype = "Approximate"),
    linewidth = 1
  ) +
  geom_vline(
    aes(xintercept = fit_true_lambda, col = "True value", linetype = "True value"),
    linewidth = 1
  ) +
  geom_function(
    fun = dgamma,
    args = list(shape = 1 + sum(fit_df$y), rate = gam + length(fit_df$y)),
    mapping = aes(col = "Exact", linetype = "Exact"),
    n = 1000,
    linewidth = 1
  ) +
  scale_colour_manual(
    values = c("Approximate" = "blue", "Exact" = "red", "True value" = "black")
  ) +
  scale_linetype_manual(
    values = c("Approximate" = "dashed", "Exact" = "solid", "True value" = "solid")
  ) +
  guides(
    col = guide_legend(title = "Density"),
    fill = guide_legend(title = "Density"),
    linetype = guide_legend(title = "Density")
  ) +
  theme_classic() +
  ylab("Posterior density") +
  xlab(expression(lambda)) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

# Figure 2:
############

pdf(
  file = here::here("Figures", "toy_single.pdf"),
  width = tw_in,
  height = tw_in/2 # ,
  #  units = "cm",
  #  res = 300
)

print(g_single)

dev.off()

g_single

# Spatial modelling examples
#############################

# Load Glasgow area data + simulate counts #
############################################

# aggregation_data.rds consists of:
# glasgow  - glasgow intermediate zone polygons
# glasgow_bnd - glasgow boundary polygon
# mesh - SPDE mesh
# ips - integration points and weights for integrating the random field in
#       each intermediate zone

glasgow_data <- readRDS(here::here("Data", "aggregation_data.rds"))
glasgow <- glasgow_data$glasgow_iz
bnd <- glasgow_data$glasgow_bnd
mesh <- glasgow_data$mesh

# Simulate from Matern GRF on a regular grid
pts <- fm_pixels(mesh, mask = bnd)

# matern covariance parameters
rho <- 5
sigma <- 0.5

seed <- 3733
set.seed(seed)
samp <- fm_matern_sample(mesh, rho = rho, sigma = sigma)[, 1]
pts$z <- fm_evaluate(mesh, loc = pts, field = samp)

# linear predictor parameters
intercept <- 1.5
b_z <- 1

# intensity field
pts$lambda <- exp(intercept + b_z * pts$z)

# aggregate intensity to areal unit
pts$ID <- st_join(pts, glasgow)$ID
ips <- glasgow_data$ips
agg <- bru_mapper_logsumexp(
  rescale = FALSE,
  n_block = nrow(glasgow)
)

z_ips <- fm_evaluate(mesh, loc = ips, field = samp)
glasgow$lambda <- ibm_eval(
  agg,
  input = list(block = ips$.block, weights = ips$weight),
  state = intercept + b_z * z_ips,
  log = FALSE
)

# simulate counts in area units
set.seed(1234)
glasgow$count <- rpois(nrow(glasgow), glasgow$lambda)

# Example 5.1: BYM area model
#############################

# adjacency matrix
glasgow_sp <- as(glasgow, "Spatial")
glasgow_nb <- poly2nb(glasgow_sp)
B <- nb2mat(glasgow_nb, style = "B") # binary yes/no neighbours

# manually add neighbours with links over the Clyde river
clyde_nbs <- list(
  c(120, 15),
  c(95, 19),
  c(94, 49),
  c(93, 49),
  c(91, 49),
  c(50, 48),
  c(50, 47),
  c(52, 47)
)

for (j in 1:length(clyde_nbs)) {
  id <- clyde_nbs[[j]]
  B[id[1], id[2]] <- 1
  B[id[2], id[1]] <- 1
}

# inla.graph object of neighbourhood structure
x <- mat2listw(B, style = "W")
nb2INLA(here::here("Data", "glasgow.adj"), x$neighbours)
g <- inla.read.graph(here::here("Data", "glasgow.adj"))

# specify model components
cmp <- ~ 0 + beta(1) + w(ID,
  model = "bym",
  graph = g
)

# model formula
fml <- count ~ .

# construct likelihood
lik <- like(
  family = "poisson",
  formula = fml,
  data = glasgow
)

# fit model
fit <- bru(
  components = cmp,
  lik
)

# generate predictions at area level
pred <- predict(
  fit,
  formula = ~ exp(beta_latent + w_latent),
  n.samples = 500
)

glasgow$lambda_post_mean <- pred$mean[seq_len(nrow(glasgow))]

# plot prediction alongside observed count
vals <- c(glasgow$count, glasgow$lambda_post_mean)
lower <- min(vals)
upper <- max(vals)
scl <- scale_fill_viridis_c("rate",
  limits = c(lower, upper)
)

# posterior mean rate per area unit
p_pred <- ggplot() +
  geom_sf(
    data = glasgow,
    aes(fill = lambda_post_mean), alpha = 1
  ) +
  scl +
  theme_void() +
  coord_sf(datum = st_crs(glasgow)) +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("Posterior expected count") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# observed count per area unit
p_truth <- ggplot() +
  geom_sf(
    data = glasgow,
    aes(fill = count), alpha = 1
  ) +
  scl +
  theme_void() +
  coord_sf(datum = st_crs(glasgow)) +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("Observed count") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )


# Figure 4
############

p_spatial_example_1 <-
  p_pred + p_truth +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "right") &
    guides(fill = guide_colourbar(
      title = "",
      title.theme = element_text(size = 18),
      label.theme = element_text(size = 12)
    ))

pdf(
  file = here::here("Figures", "spatial_example_1.pdf"),
  width = tw_in,
  height = tw_in / 2 # ,
  #  units = "cm",
  #  res = 300,
  #  type = "cairo-png"
)

p_spatial_example_1

dev.off()

p_spatial_example_1


# Example 5.2: Aggregation of SPDE to areal unit ###
####################################################

# define model components
matern <- inla.spde2.pcmatern(mesh,
  prior.range = c(2, 0.1),
  prior.sigma = c(3, 0.1)
)

cmp <- ~ 0 + beta(1) + xi(main = geometry, model = matern)

# aggregation mapper
agg <- bru_mapper_logsumexp(
  rescale = FALSE,
  n_block = nrow(glasgow)
)

# model formula
fml <- count ~ ibm_eval(agg,
  input = list(block = .block, weights = weight),
  state = beta + xi,
  log = TRUE
)

# likelihood
lik <- like(
  formula = fml,
  family = "poisson",
  data = ips,
  response_data = glasgow
)

if (!file.exists(here::here("R_cache", "spatial_example_2_fit.rds"))) {
  # Fit model
  fit_2 <- bru(
    components = cmp,
    lik,
    options = list(
      bru_verbose = 2,
      bru_max_iter = 10
    )
  )

  saveRDS(fit_2, here::here("R_cache", "spatial_example_2_fit.rds"))
}
fit_2 <- readRDS(here::here("R_cache", "spatial_example_2_fit.rds"))


# Assess convergence

# Figure 4
###########

p_conv <- bru_convergence_plot(fit_2)

pdf(
  file = here::here("Figures", "spatial_example_fit2_convergence.pdf"),
  width = 1.5 * tw_in,
  height = tw_in # ,
  #  units = "cm",
  #  res = 300,
  #  type = "cairo-png"
)

p_conv

dev.off()

p_conv

# predictions at area level
pred_2 <- predict(
  fit_2,
  newdata = ips,
  formula = ~ ibm_eval(agg,
    input = list(block = .block, weights = weight),
    beta + xi,
    log = FALSE
  ),
  n.samples = 500
)
glasgow$lambda_post_mean_2 <- pred_2$mean

# colour scale
vals <- c(glasgow$lambda, glasgow$lambda_post_mean_2)
lower <- min(vals)
upper <- max(vals)
scl <- scale_fill_viridis_c("rate",
  limits = c(lower, upper)
)

# predictions at pixel level
pred_2_pt <- predict(fit_2,
  newdata = pts,
  formula = ~ exp(beta + xi),
  n.samples = 500
)

pts$lambda_post_mean_2 <- pred_2_pt$mean

# colour scale
vals <- c(pts$lambda, pts$lambda_post_mean_2)
lower <- min(vals)
upper <- max(vals)
scl_px <- scale_fill_viridis_c("rate",
  limits = c(lower, upper)
)

# plots

# posterior mean rate per area unit
p_pred <- ggplot() +
  geom_sf(
    data = glasgow,
    aes(fill = lambda_post_mean_2)
  ) +
  scl +
  theme_void() +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("Posterior expected count") +
  theme(plot.title = element_text(hjust = 0.5))

# true rate per area unit
p_truth <- ggplot() +
  geom_sf(
    data = glasgow,
    aes(fill = lambda)
  ) +
  scl +
  theme_void() +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("True expected count") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# posterior mean rate field
p_pred_pts <- ggplot() +
  gg(
    aes(
      fill = lambda_post_mean_2
    ),
    geom = "tile",
    data = pts,
  ) +
  coord_sf(crs = fm_crs(pts)) +
  scl_px +
  theme_void() +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("Posterior mean field") +
  theme(plot.title = element_text(hjust = 0.5))

# true rate field
p_truth_pts <- ggplot() +
  gg(
    aes(
      fill = lambda
    ),
    geom = "tile",
    data = pts,
  ) +
  coord_sf(crs = fm_crs(pts)) +
  scl_px +
  theme_void() +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("True field") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Figure 6:
############

p_spatial_example_2 <-
  (p_pred + p_truth +
    plot_layout(guides = "collect") &
    guides(
      fill = guide_colourbar(
        title = expression(lambda[i]),
        title.hjust = 0.12,
        title.theme = element_text(size = 14),
        label.theme = element_text(size = 12)
      )
    )) /
    (p_pred_pts + p_truth_pts +
      plot_layout(guides = "collect") &
      guides(
        fill = guide_colourbar(
          title = expression(lambda(s)),
          # title.hjust = 0.1,
          title.theme = element_text(size = 14),
          label.theme = element_text(size = 12)
        )
      )) +
    plot_annotation(tag_levels = "A") &
    theme(plot.title = element_text(hjust = 0.5))

pdf(
  file = here::here("Figures", "spatial_example_2.pdf"),
  width = tw_in,
  height = tw_in,
  #  units = "cm",
  #  res = 300,
  #  type = "cairo-png"
)

p_spatial_example_2

dev.off()

p_spatial_example_2


# Example 3) Aggregation to areal unit level with ####
############ latent field observed at points. ########
######################################################

# observe latent covariate (field realisation) at a set of points
set.seed(10301)
z_locs <- st_sf(
  geometry = st_sample(bnd,
    size = 100,
    type = "random"
  )
)

# Shortcut version:
z_locs$z <- fm_evaluate(mesh, loc = z_locs, field = samp)

# specify model components
cmp <- ~ 0 + xi(main = geometry, model = matern) +
  alpha_0(1) +
  beta_0(1) +
  beta_1(1)

# z formula
z_fml <- z ~ alpha_0 + xi

# count formula
count_fml <- count ~ ibm_eval(
  agg,
  input = list(weights = weight, block = .block),
  state = beta_0 + beta_1 * xi
)

# z likelihood
z_lik <- like("Gaussian",
  formula = z_fml,
  data = z_locs
)

# count likelihood
count_lik <- like("Poisson",
  formula = count_fml,
  data = ips,
  response_data = glasgow
)

if (!file.exists(here::here("R_cache", "spatial_example_3_fit.rds"))) {
  # Fit model
  fit_3 <- bru(
    components = cmp,
    z_lik,
    count_lik,
    options = list(
      bru_verbose = 1,
      bru_max_iter = 10#,
#      bru_initial = list(beta_0 = 1, beta_1 = 1, xi = 1)
    )
  )

  saveRDS(fit_3, here::here("R_cache", "spatial_example_3_fit.rds"))
}
fit_3 <- readRDS(here::here("R_cache", "spatial_example_3_fit.rds"))

# Assess convergence

# Figure 7
##########

p_conv <- bru_convergence_plot(fit_3)

pdf(
  file = here::here("Figures", "spatial_example_fit3_convergence.pdf"),
  width = 1.5 * tw_in,
  height = tw_in # ,
  #  units = "cm",
  #  res = 300,
  #  type = "cairo-png"
)

p_conv

dev.off()

p_conv

# model predictions

# rate field
pred_3_pts <- predict(fit_3,
  newdata = pts,
  formula = ~ exp(beta_0 + beta_1 * xi)
)

pts$lambda_post_mean_3 <- pred_3_pts$mean

# latent field predictions
pred_3_field <- predict(fit_3,
  newdata = pts,
  formula = ~xi
)

pts$z_post_mean <- pred_3_field$mean

# colour scales for plots
vals <- c(pts$lambda, pts$lambda_post_mean_3)
lower <- min(vals)
upper <- max(vals)
sfl <- scale_fill_viridis_c(limits = c(lower, upper))
scl <- scale_colour_viridis_c(limits = c(lower, upper))

# posterior mean rate field
p_pred_pts_3 <- ggplot() +
  geom_tile(
    data = pts,
    mapping = aes(geometry = geometry, fill = lambda_post_mean_3),
    stat = "sf_coordinates"
  ) +
  geom_sf(data = bnd, alpha = 0) +
  sfl +
  theme_void() +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("Posterior mean field") +
  theme(plot.title = element_text(hjust = 0.5))

# True rate field
p_truth_pts_3 <- ggplot() +
  geom_tile(
    data = pts,
    mapping = aes(geometry = geometry, fill = lambda),
    stat = "sf_coordinates"
  ) +
  geom_sf(data = bnd, alpha = 0) +
  sfl +
  theme_void() +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("True field") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )


# colour scales for plots of z
vals <- c(pts$z, pts$z_post_mean)
lower <- min(vals)
upper <- max(vals)
z_sfl <- scale_fill_viridis_c(limits = c(lower, upper))
z_scl <- scale_colour_viridis_c(limits = c(lower, upper))

# posterior mean for z field
p_z_pred <- ggplot() +
  geom_tile(
    data = pts,
    aes(geometry = geometry, fill = z_post_mean),
    stat = "sf_coordinates"
  ) +
  geom_sf(data = bnd, alpha = 0) +
  z_sfl +
  theme_void() +
  xlab("Easting") +
  ylab("Northing") +
  ggtitle("Covariate posterior mean") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# covariate observations
p_z <- ggplot() +
  geom_sf(
    data = z_locs,
    mapping = aes(colour = z),
    cex = 0.8
  ) +
  geom_sf(
    data = bnd,
    alpha = 0
  ) +
  scale_colour_viridis_c("z") +
  theme_void() +
  ggtitle("Covariate observations") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.15)
  ) +
  guides(colour = guide_colourbar(label.theme = element_text(size = 12)))

# add new scale to z observations plot
p_z2 <- p_z +
  z_scl +
  theme(legend.position = "none")

# Figure 8
##########

p_spatial_example_3b <-
  (p_pred_pts_3 + p_truth_pts_3 +
    plot_layout(guides = "collect") &
    guides(
      fill = guide_colourbar(
        title = expression(lambda(s)),
        # title.hjust = 0.1,
        title.theme = element_text(size = 14),
        label.theme = element_text(size = 12)
      )
    )) /
    (p_z_pred + p_z2 +
      plot_layout(guides = "collect") &
      guides(
        fill = guide_colourbar(
          title = "z",
          title.hjust = 0.1,
          title.theme = element_text(size = 14),
          label.theme = element_text(size = 12)
        )
      )) +
    plot_annotation(tag_levels = "A") &
    theme(plot.title = element_text(hjust = 0.5))

pdf(
  file = here::here("Figures", "spatial_example_3b.pdf"),
  width = tw_in,
  height = tw_in # ,
  #  units = "cm",
  #  res = 300,
  #  type = "cairo-png"
)

p_spatial_example_3b

dev.off()

p_spatial_example_3b

sessionInfo()
