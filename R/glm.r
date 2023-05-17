#' glmPairsBoot
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle Incremental claims triangle.
#' @param nboot DESCRIPTION.
#' @param nsim DESCRIPTION.
#' @param keep DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
glmPairsBoot <- function(triangle, nboot, nsim, keep) {
  ndev <- ncol(triangle)
  colnames(triangle) <- rownames(triangle) <- 1:ndev

  long <- as.data.frame(triangle)
  long <- as.data.frame(lapply(long, as.numeric))
  latest <- rev(long[long$dev + long$origin == ndev + 1, "value"])[-1]

  long.lower <- long[is.na(long$value), ]
  long.upper <- long[!is.na(long$value), ]
  nobs <- nrow(long.upper)
  npred <- ndev**2 - nobs

  if (keep == "corners") {
    long.keep <- long.upper[-c(ndev, nobs), ]
  } else if (keep == "rowcol") {
    long.keep <- rbind(long.upper[long.upper$origin == 1, ], long.upper[long.upper$dev == 1, ])
  }

  model <- glm(value ~ factor(origin) + factor(dev), quasipoisson(), long.upper)
  long.lower$value <- predict(model, long.lower, type = "response")
  disp <- summary(model)$dispersion

  combine <- function(list1, list2) {
    res <- list()
    res$a <- rbind(list1$a, list2$a)
    res$b <- rbind(list1$b, list2$b)
    res$intercept <- c(list1$intercept, list2$intercept)
    res$disp <- c(list1$disp, list2$disp)
    res$reserve <- rbind(list1$reserve, list2$reserve)
    res$pred.error <- rbind(list1$pred.error, list2$pred.error)
    res$pred.dist <- c(list1$pred.dist, list2$pred.dist)
    return(res)
  }

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)

  res.list <- foreach(
    iboot = 1:nboot,
    .combine = combine,
    .packages = "ChainLadder",
    .options.snow = opts) %dopar% {
    deficient <- TRUE
    while (deficient) {
      if (keep == "randomised") {
        long.keep <- data.frame(origin = numeric(), dev = numeric(), value = numeric())
        for (j in 1:ndev) {
          long.keep <- rbind(long.keep,
            long.upper[long.upper$dev == j & long.upper$origin == sample(1:(ndev + 1 - j), 1), ]
          )
        }
        for (i in 1:ndev) {
          long.keep <- rbind(long.keep,
            long.upper[long.upper$origin == i & long.upper$dev == sample(1:(ndev + 1 - i), 1), ]
          )
        }
      }
      nidxs <- nobs - nrow(long.keep)
      idxs.boot <- sample(seq_len(nrow(long.upper)), nidxs, replace = TRUE)
      long.boot <- rbind(long.upper[idxs.boot, ], long.keep)

      model.boot <- glm(value ~ factor(origin) + factor(dev), quasipoisson(), long.boot)
      deficient <- (length(names(model.boot$coefficients)) != 2 * ndev - 1) || (model.boot$qr$rank < 2 * ndev - 1)
    }
    a.boot <- model.boot$coefficients[2:ndev]
    b.boot <- model.boot$coefficients[(ndev + 1):(2 * ndev - 1)]
    intercept.boot <- model.boot$coefficients[1]
    disp.boot <- summary(model.boot)$dispersion

    long.proj <- long.lower
    long.proj$value <- predict(model.boot, long.lower, type = "response")
    reserve.proj <- tapply(long.proj$value, long.proj$origin, sum, na.rm = TRUE)

    pred.error <- matrix(rep(0, nsim * (ndev - 1)), nrow = nsim)
    pred.dist <- rep(0, nsim)
    for (isim in 1:nsim) {
      long.pred <- long.proj
      long.fit <- long.lower
      long.pred$value <- rnorm(npred, long.pred$value, sqrt(disp.boot * long.pred$value))
      long.fit$value <- rnorm(npred, long.fit$value, sqrt(disp * long.fit$value))

      reserve.pred <- tapply(long.pred$value, long.pred$origin, sum, na.rm = TRUE)
      reserve.fit <- tapply(long.fit$value, long.fit$origin, sum, na.rm = TRUE)
      pred.error[isim, ] <- reserve.proj - reserve.fit
      pred.dist[isim] <- sum(reserve.pred)
    }

    list(
      a = a.boot,
      b = b.boot,
      intercept = intercept.boot,
      disp = disp.boot,
      reserve = reserve.proj,
      pred.error = pred.error,
      pred.dist = pred.dist
    )
  }
  close(pgb)

  res <- list(
    table = data.frame(
      idx = 2:ndev,
      a = colMeans(res.list$a),
      b = colMeans(res.list$b),
      reserve = colMeans(res.list$reserve),
      prederror = sqrt(colMeans(res.list$pred.error**2))
    ),
    point = data.frame(
      intercept = mean(res.list$intercept),
      dispersion = mean(res.list$disp)
    ),
    pred.dist = res.list$pred.dist
  )
  return(res)
}

#' glmSemiParamBoot
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#' @param nboot DESCRIPTION.
#' @param nsim DESCRIPTION.
#' @param resids.type DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
glmSemiParamBoot <- function(triangle, nboot, nsim, resids.type) {
  ndev <- ncol(triangle)
  colnames(triangle) <- rownames(triangle) <- 1:ndev

  long <- as.data.frame(triangle)
  long <- as.data.frame(lapply(long, as.numeric))
  latest <- rev(long[long$dev + long$origin == ndev + 1, "value"])[-1]
  long.lower <- long[is.na(long$value), ]
  long.upper <- long[!is.na(long$value), ]
  nobs <- nrow(long.upper)
  npred <- ndev**2 - nobs

  model <- glm(value ~ factor(origin) + factor(dev), quasipoisson(), long.upper)
  long.lower$value <- predict(model, long.lower, type = "response")
  fitted <- model$fitted.values

  n <- (ndev**2 + ndev) / 2
  p <- 2 * ndev - 1

  if (resids.type == "pearson") {
    resids <- resid(model, type = "pearson") # not divided by dispersion
    leverages <- hatvalues(model)
    leverages <- leverages[-c(ndev, nobs)] # exclude corner points
    resids <- resids[-c(ndev, nobs)]
    resids <- resids / sqrt(1 - leverages)
    resids <- resids - mean(resids)
  } else if (resids.type == "deviance") {
    y <- model$y
    mu <- model$fitted.values
    resids <- sign(y - mu) * sqrt(2 * (y * log(y / mu) + mu - y))
    rm.idxs <- unique(c(which(!(resids > -sqrt(2 * min(fitted)))), ndev, nobs))
    resids <- resids[-rm.idxs]
    leverages <- hatvalues(model)
    leverages <- leverages[-rm.idxs]
    resids <- resids / sqrt(1 - leverages)
    resids <- resids - mean(resids)

    # define helper function to invert residuals
    devResInv <- function(mu, resid) {
      diffFunc <- function(x) {
        w <- resid**2 / (2 * mu) - 1
        res <- x * log(x) - x - w
        return(res)
      }

      dFun <- function(x) { log(x) }

      if (resid < 0) {
        lb <- 0
        ub <- 1
      } else {
        w <- resid**2 / (2 * mu) - 1
        if (w == -1) {
          lb <- ub <- 1
        } else if (-1 < w && w < exp(1)**2) {
          lb <- 1 + (w + 1) * (exp(1)**2 - 1) / (exp(1)**2 + 1)
          ub <- 1 + sqrt(w + 1) * (exp(1)**2 - 1) / sqrt(exp(1)**2 + 1)
        } else if (w == exp(1)**2) {
          lb <- ub <- exp(1)**2
        } else {
          lb <- w / log(w) + exp(1)**2 / 2
          ub <- 2 * w / log(w)
        }
      }
      res <- pracma::newtonRaphson(diffFunc, 1.5, dFun)
      return(res$root * mu)
    }
  }

  combine <- function(list1, list2) {
    res <- list()
    res$a <- rbind(list1$a, list2$a)
    res$b <- rbind(list1$b, list2$b)
    res$intercept <- c(list1$intercept, list2$intercept)
    res$disp <- c(list1$disp, list2$disp)
    res$reserve <- rbind(list1$reserve, list2$reserve)
    res$pred.error <- rbind(list1$pred.error, list2$pred.error)
    res$pred.dist <- c(list1$pred.dist, list2$pred.dist)
    return(res)
  }

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)

  res.list <- foreach(
    iboot = 1:nboot,
    .combine = combine,
    .options.snow = opts) %dopar% {
    long.boot <- long.upper
    if (resids.type == "pearson") {
      negative <- TRUE
      while (negative) {
        resids.boot <- sample(resids, nobs, replace = TRUE)
        long.boot$value <- fitted + sqrt(fitted) * resids.boot
        negative <- any(long.boot$value < 0)
      }
    } else if (resids.type == "deviance") {
      resids.boot <- sample(resids, nobs, replace = TRUE)
      for (i in 1:nobs) {
        long.boot$value[i] <- devResInv(fitted[i], resids.boot[i])
      }
    }

    model.boot <- glm(value ~ factor(origin) + factor(dev), quasipoisson(), long.boot)
    a.boot <- model.boot$coefficients[2:ndev]
    b.boot <- model.boot$coefficients[(ndev + 1):(2 * ndev - 1)]
    intercept.boot <- model.boot$coefficients[1]
    disp.boot <- summary(model.boot)$dispersion

    long.proj <- long.lower
    long.proj$value <- predict(model.boot, long.lower, type = "response")
    reserve.proj <- tapply(long.proj$value, long.proj$origin, sum, na.rm = TRUE)

    pred.error <- matrix(rep(0, nsim * (ndev - 1)), nrow = nsim)
    pred.dist <- rep(0, nsim)
    for (isim in 1:nsim) {
      long.pred <- long.proj
      long.fit <- long.lower
      if (resids.type == "pearson") {
        resids.sim <- sample(resids, npred, replace = TRUE)
        long.pred$value <- long.pred$value + sqrt(long.pred$value) * resids.sim
        long.fit$value <- long.fit$value + sqrt(long.fit$value) * resids.sim
      } else if (resids.type == "deviance") {
        resids.sim <- sample(resids, npred, replace = TRUE)
        for (i in 1:npred) {
          long.pred$value[i] <- devResInv(long.pred$value[i], resids.sim[i])
          long.fit$value[i] <- devResInv(long.fit$value[i], resids.sim[i])
        }
      }

      reserve.pred <- tapply(long.pred$value, long.pred$origin, sum, na.rm = TRUE)
      reserve.fit <- tapply(long.fit$value, long.fit$origin, sum, na.rm = TRUE)
      pred.error[isim, ] <- reserve.fit - reserve.proj
      pred.dist[isim] <- sum(reserve.pred)
    }

    list(
      a = a.boot,
      b = b.boot,
      intercept = intercept.boot,
      disp = disp.boot,
      reserve = reserve.proj,
      pred.error = pred.error,
      pred.dist = pred.dist
    )
  }
  close(pgb)

  res <- list(
    table = data.frame(
      idx = 2:ndev,
      a = colMeans(res.list$a),
      b = colMeans(res.list$b),
      reserve = colMeans(res.list$reserve),
      prederror = sqrt(colMeans(res.list$pred.error**2))
    ),
    point = data.frame(
      intercept = mean(res.list$intercept),
      disp = mean(res.list$dis)
    ),
    pred.dist = res.list$pred.dist
  )

  return(res)
}

#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#' @param nboot DESCRIPTION.
#' @param nsim DESCRIPTION.
#' @param dist DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
glmParamBoot <- function(triangle, nboot, nsim, dist) {
  ndev <- ncol(triangle)
  colnames(triangle) <- rownames(triangle) <- 1:ndev

  long <- as.data.frame(triangle)
  long <- as.data.frame(lapply(long, as.numeric))
  latest <- rev(long[long$dev + long$origin == ndev + 1, "value"])[-1]
  long.lower <- long[is.na(long$value), ]
  long.upper <- long[!is.na(long$value), ]
  nobs <- nrow(long.upper)
  npred <- ndev**2 - nobs

  model <- glm(value ~ factor(origin) + factor(dev), quasipoisson(), long.upper)
  long.lower$value <- predict(model, long.lower, type = "response")
  disp <- summary(model)$dispersion
  fitted <- model$fitted.values

  combine <- function(list1, list2) {
    res <- list()
    res$a <- rbind(list1$a, list2$a)
    res$b <- rbind(list1$b, list2$b)
    res$intercept <- c(list1$intercept, list2$intercept)
    res$disp <- c(list1$disp, list2$disp)
    res$reserve <- rbind(list1$reserve, list2$reserve)
    res$pred.error <- rbind(list1$pred.error, list2$pred.error)
    res$pred.dist <- c(list1$pred.dist, list2$pred.dist)
    return(res)
  }

  pgb <- txtProgressBar(max = nboot)
  progress <- function(n) setTxtProgressBar(pgb, n)
  opts <- list(progress = progress)

  res.list <- foreach(
    iboot = 1:nboot,
    .combine = combine,
    .packages = "ChainLadder",
    .options.snow = opts) %dopar% {
    long.boot <- long.upper

    negative <- TRUE
    while (negative) {
      if (dist == "normal") {
        long.boot$value <- rnorm(nobs, fitted, sqrt(disp * fitted))
      } else if (dist == "gamma") {
        alpha <- fitted**2 / (disp * fitted)
        beta <- fitted / (disp * fitted)
        long.boot$value <- rgamma(nobs, alpha, beta)
      } else {
        long.boot$value <- disp * rpois(nobs, fitted / disp)
      }
      negative <- any(long.boot$value < 0)
    }

    model.boot <- glm(value ~ factor(origin) + factor(dev), quasipoisson(), long.boot)
    a.boot <- model.boot$coefficients[2:ndev]
    b.boot <- model.boot$coefficients[(ndev + 1):(2 * ndev - 1)]
    intercept.boot <- model.boot$coefficients[1]
    disp.boot <- summary(model.boot)$dispersion

    long.proj <- long.lower
    long.proj$value <- predict(model.boot, long.lower, type = "response")
    reserve.proj <- tapply(long.proj$value, long.proj$origin, sum, na.rm = TRUE)

    pred.error <- matrix(rep(0, nsim * (ndev - 1)), nrow = nsim)
    pred.dist <- rep(0, nsim)
    for (isim in 1:nsim) {
      long.pred <- long.proj
      long.fit <- long.lower
      long.pred$value <- rnorm(npred, long.pred$value, sqrt(disp.boot * long.pred$value))
      long.fit$value <- rnorm(npred, long.fit$value, sqrt(disp * long.fit$value))

      reserve.pred <- tapply(long.pred$value, long.pred$origin, sum, na.rm = TRUE)
      reserve.fit <- tapply(long.fit$value, long.fit$origin, sum, na.rm = TRUE)
      pred.error[isim, ] <- reserve.proj - reserve.fit
      pred.dist[isim] <- sum(reserve.pred)
    }

    list(
      a = a.boot,
      b = b.boot,
      intercept = intercept.boot,
      disp = disp.boot,
      reserve = reserve.proj,
      pred.error = pred.error,
      pred.dist = pred.dist
    )
  }
  close(pgb)

  res <- list(
    table = data.frame(
      idx = 2:ndev,
      a = colMeans(res.list$a),
      b = colMeans(res.list$b),
      reserve = colMeans(res.list$reserve),
      prederror = sqrt(colMeans(res.list$pred.error**2))
    ),
    point = data.frame(
      intercept = mean(res.list$intercept),
      disp = mean(res.list$disp)
    ),
    pred.dist = res.list$pred.dist
  )
  return(res)
}

#' glmFit
#'
#' FUNCTION_DESCRIPTION
#'
#' @param triangle DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @export
glmFit <- function(triangle) {
  ndev <- ncol(triangle)

  long <- as.data.frame(triangle)
  long <- transform(long, origin = factor(origin, levels = dimnames(triangle)[[1]]))
  long <- long[order(long$origin, long$dev), ]

  long.lower <- long[is.na(long$value), ]
  long.upper <- long[!is.na(long$value), ]

  model <- glm(value ~ factor(origin) + factor(dev), quasipoisson(), long.upper, maxit = 1e3)
  a <- model$coefficients[2:ndev]
  b <- model$coefficients[(ndev + 1):(2 * ndev - 1)]
  intercept <- model$coefficients[1]
  disp <- summary(model)$dispersion

  long.lower$value <- predict(model, long.lower, type = "response")
  triangle.proj <- incr2cum(as.triangle(rbind(long.upper, long.lower)))

  latest <- triangle.proj[
    row(triangle.proj) + col(triangle.proj) == ndev + 1][1:(ndev - 1)
  ]
  reserve <- triangle.proj[2:ndev, ndev] - rev(latest)

  resids.pears <- resid(model, type = "pearson")
  resids.dev <- sign(model$y - model$fitted.values) *
    sqrt(zapsmall(2 * (model$y * log(model$y / model$fitted.values) + model$fitted.values - model$y)))

  param.var.quasi <- summary(model)$cov.scaled # (X^T X)^-1 * \phi
  param.var.pois <- summary(model)$cov.unscaled # (X^T X)^-1
  model.terms <- delete.response(terms(model))
  X <- model.matrix(model.terms, long.lower, xlev = model$xlevels)

  eta.vars.pois <- list()
  eta.vars.quasi <- list()
  idx <- 1
  for (i in 2:ndev) {
    ncols <- i - 1
    eta.vars.pois[[i - 1]] <- X[idx:(idx + ncols - 1), , drop = FALSE] %*% param.var.pois %*% t(X[idx:(idx + ncols - 1), , drop = FALSE]) # nolint
    eta.vars.quasi[[i - 1]] <- X[idx:(idx + ncols - 1), , drop = FALSE] %*% param.var.quasi %*% t(X[idx:(idx + ncols - 1), , drop = FALSE]) # nolint
    idx <- idx + ncols
  }

  triangle.lower <- as.triangle(long.lower)
  pred.error.pois <- rep(0, ndev - 1)
  pred.error.quasi <- rep(0, ndev - 1)
  for (i in 1:(ndev - 1)) {
    mu <- triangle.lower[i, !is.na(triangle.lower[i, ]), drop = FALSE]
    pred.error.pois[i] <- sqrt(sum(triangle.lower[i, ], na.rm = TRUE) + mu %*% eta.vars.pois[[i]] %*% t(mu))
    pred.error.quasi[i] <- sqrt(disp * sum(triangle.lower[i, ], na.rm = TRUE) + mu %*% eta.vars.quasi[[i]] %*% t(mu))
  }

  res <- list(
    table = data.frame(
      idx = 2:ndev,
      a = a,
      b = b,
      reserve = reserve,
      prederrorpois = pred.error.pois,
      prederrorquasi = pred.error.quasi
    ),
    point = data.frame(
      intercept = intercept,
      disp = disp
    ),
    resids = data.frame(
      resids.pears = resids.pears,
      resids.dev = resids.dev
    ),
    fitted = model$fitted.values,
    model = model
  )

  return(res)
}
