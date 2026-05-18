library(gamlss)

# ── Helpers de transformación ─────────────────────────────────────────────────
apply_transform <- function(y, transform) {
  switch(transform,
    "none"  = y,
    "log"   = { if (any(y <= 0, na.rm=TRUE))
                  stop(sprintf("transform='log' requiere y>0 (%d valores<=0)",
                               sum(y <= 0, na.rm=TRUE)))
                log(y) },
    "log1p" = log(y + 1),
    "sqrt"  = { if (any(y < 0, na.rm=TRUE))
                  stop(sprintf("transform='sqrt' requiere y>=0 (%d valores<0)",
                               sum(y < 0, na.rm=TRUE)))
                sqrt(y) },
    stop(sprintf("transform desconocido: '%s'", transform))
  )
}

apply_inv_transform <- function(y_t, transform) {
  switch(transform,
    "none"  = y_t,
    "log"   = exp(y_t),
    "log1p" = exp(y_t) - 1,
    "sqrt"  = y_t^2,
    stop(sprintf("inv_transform desconocido: '%s'", transform))
  )
}

# ── Función principal ─────────────────────────────────────────────────────────
fit_one_prefix <- function(prefix, DATA_DIR=".", N_AGE_GRID=500, FAMILY="NO") {

  input_file <- file.path(DATA_DIR, paste0(prefix, "_gamlss_data.csv"))
  if (!file.exists(input_file)) {
    cat(sprintf("[%s] archivo no encontrado, SKIP\n", prefix))
    return(invisible(NULL))
  }

  cat(sprintf("\n%s\nPREFIJO: %s\n%s\n", strrep("=",60), prefix, strrep("=",60)))

  # ── Leer transformación desde _gamlss_meta.csv (generado por export_for_gamlss.m) ──
  transform <- "none"
  meta_file <- file.path(DATA_DIR, paste0(prefix, "_gamlss_meta.csv"))
  if (file.exists(meta_file)) {
    meta     <- read.csv(meta_file, stringsAsFactors=FALSE)
    meta_map <- setNames(meta$value, meta$field)
    if ("transform" %in% names(meta_map))
      transform <- meta_map[["transform"]]
    cat(sprintf("  transform = '%s'  (leído de %s)\n", transform, basename(meta_file)))
  } else {
    cat(sprintf("  transform = 'none'  (no se encontró %s)\n", basename(meta_file)))
  }

  dat       <- read.csv(input_file, stringsAsFactors=FALSE, na.strings="NA")
  age       <- dat$age
  var_names <- setdiff(names(dat), c("age","equipment"))
  N <- nrow(dat); P <- length(var_names)
  cat(sprintf("  %d sujetos, %d variables\n", N, P))

  age_valid <- age[is.finite(age) & age > 0]
  age_grid  <- exp(seq(log(max(min(age_valid), 1)),
                       log(max(age_valid)),
                       length.out = N_AGE_GRID))

  results_centiles    <- data.frame(age = age_grid)
  results_diagnostics <- data.frame(
    variable       = var_names,
    pct_outside    = NA_real_,
    n_valid        = NA_integer_,
    family         = NA_character_,
    transform_used = NA_character_,
    stringsAsFactors = FALSE
  )
  coef_export    <- data.frame(age = age_grid)
  centile_levels <- c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975)
  centile_names  <- c("p025","p10","p25","p50","p75","p90","p975")
  n_ok <- 0

  for (p in seq_along(var_names)) {
    vname <- var_names[p]
    y     <- dat[[vname]]
    valid <- !is.na(y) & is.finite(y) & !is.na(age) & is.finite(age) & age > 0
    yy_orig <- y[valid]
    aa      <- age[valid]
    n_v     <- sum(valid)
    results_diagnostics$n_valid[p] <- n_v

    cat(sprintf("  Ajustando %s (%d/%d, n=%d)... ", vname, p, P, n_v))
    if (n_v < 30) { cat("SKIP (n<30)\n"); next }

    # ── Columna constante = 0: no hay nada que modelar ────────────────────
    # Ocurre con índices de asimetría cuando todos los sujetos tienen el mismo
    # valor (ej. electrodo faltante) o con variables que no aplican al prefijo.
    if (all(yy_orig == 0)) {
      cat("SKIP (todos cero)\n")
      results_diagnostics$family[p]         <- "ZERO"
      results_diagnostics$transform_used[p] <- "none"
      results_diagnostics$pct_outside[p]    <- 0
      # Rellenar curvas y centiles con ceros para que MATLAB pueda leer el CSV
      for (cn in centile_names) {
        results_centiles[[paste0(vname, "_", cn)]] <- rep(0, N_AGE_GRID)
      }
      results_centiles[[paste0(vname, "_mu")]]    <- rep(0, N_AGE_GRID)
      results_centiles[[paste0(vname, "_sigma")]] <- rep(0, N_AGE_GRID)
      coef_export[[paste0(vname, "_mu")]]    <- rep(0, N_AGE_GRID)
      coef_export[[paste0(vname, "_sigma")]] <- rep(0, N_AGE_GRID)
      next
    }

    # ── Selección de familia: evaluar SIEMPRE sobre datos CRUDOS (yy_orig) ──
    # Los criterios distribucionales se evalúan en escala original porque:
    # (1) BCT y BCPE requieren y > 0, que se cumple en datos crudos positivos
    # (2) Si los datos crudos son positivos y necesitan BCT, la transformación
    #     log es redundante — BCT maneja asimetría y colas con su parámetro nu.
    # (3) Si los datos crudos no son positivos (ej. Valence), BCT puede usarse
    #     directamente sobre ellos si n_crit >= 2.
    # La transformación log (u otra) solo se aplica si la familia elegida es NO.

    kurt_val    <- mean((yy_orig - mean(yy_orig))^4) / var(yy_orig)^2
    kurt_excess <- kurt_val - 3
    skew_val    <- mean((yy_orig - mean(yy_orig))^3) / var(yy_orig)^1.5
    sw_p        <- tryCatch({
      n_sw <- min(length(yy_orig), 2000)
      set.seed(42)
      samp <- if (length(yy_orig) > n_sw) sample(yy_orig, n_sw) else yy_orig
      shapiro.test(samp)$p.value
    }, error = function(e) 1.0)

    crit_kurt <- kurt_excess > 1
    crit_skew <- abs(skew_val) > 0.5
    crit_sw   <- sw_p < 0.05
    n_crit    <- sum(c(crit_kurt, crit_skew, crit_sw))

    all_positive_orig <- all(yy_orig > 0)

    if (n_crit >= 2 && all_positive_orig) {
      # Datos positivos + no-Normal → BCT/BCPE sobre datos CRUDOS
      # (transform log es redundante; BCT lo maneja con nu)
      fam_candidates <- c("BCT", "BCPE", "NO")
      yy             <- yy_orig   # usar datos crudos, ignorar transform
      transform_used <- "none"    # BCT absorbe la asimetría
      cat(sprintf("[kurt_ex=%.2f skew=%.2f SW_p=%.3f -> BCT(raw,%d/3)] ",
                  kurt_excess, skew_val, sw_p, n_crit))
    } else if (n_crit >= 2 && !all_positive_orig) {
      # Datos con negativos + no-Normal:
      # Empíricamente verificado que para distribuciones spike-and-slab
      # (ej. Valence/FAA: masa concentrada en 0 con colas largas),
      # las familias JSU/SHASHo sobreestiman sigma y empeoran la
      # calibración de centiles centrales respecto a NO.
      # NO es la opción más robusta: el pct_outside global (~5%) es correcto
      # y los Z-scores extremos (|Z|>2) — el uso clínico principal — son válidos.
      fam_candidates <- c("NO")
      yy             <- yy_orig
      transform_used <- "none"
      cat(sprintf("[kurt_ex=%.2f skew=%.2f SW_p=%.3f -> NO(neg,spike-slab,%d/3)] ",
                  kurt_excess, skew_val, sw_p, n_crit))
    } else if (n_crit == 1 && all_positive_orig) {
      fam_candidates <- c("BCPE", "NO")
      yy             <- yy_orig
      transform_used <- "none"
      cat(sprintf("[kurt_ex=%.2f skew=%.2f SW_p=%.3f -> BCPE(raw,%d/3)] ",
                  kurt_excess, skew_val, sw_p, n_crit))
    } else {
      # Normal: aplicar transform si fue solicitado
      fam_candidates <- c("NO")
      transform_used <- transform
      yy <- tryCatch(
        apply_transform(yy_orig, transform),
        error = function(e) {
          cat(sprintf("  ERROR transform '%s': %s\n", transform, e$message))
          NULL
        }
      )
      if (is.null(yy)) next
      cat(sprintf("[kurt_ex=%.2f skew=%.2f SW_p=%.3f -> NO(transf='%s',%d/3)] ",
                  kurt_excess, skew_val, sw_p, transform_used, n_crit))
    }

    # ── Nombres únicos en GlobalEnv ───────────────────────────────────────
    env_name_fam <- paste0(".gamlss_fam_", prefix, "_", p)
    env_name_df  <- paste0(".gamlss_df_",  prefix, "_", p)
    env_name_pr  <- paste0(".gamlss_pr_",  prefix, "_", p)

    df_tr <- data.frame(yy = yy, aa = aa)
    df_pr <- data.frame(aa = age_grid)

    fit      <- NULL
    fam_used <- NA_character_

    for (fam_try in fam_candidates) {
      fam <- fam_try
      assign(env_name_fam, fam,   envir = .GlobalEnv)
      assign(env_name_df,  df_tr, envir = .GlobalEnv)
      assign(env_name_pr,  df_pr, envir = .GlobalEnv)
      assign("yy",  yy,  envir = .GlobalEnv)
      assign("aa",  aa,  envir = .GlobalEnv)
      assign("fam", fam, envir = .GlobalEnv)

      fit_try <- tryCatch({
        fit_obj <- gamlss(
          yy ~ pb(log(aa)),
          sigma.formula = ~ pb(log(aa)),
          family  = fam,
          data    = df_tr,
          control = gamlss.control(n.cyc = 50, trace = FALSE)
        )
        # Parchar call para que predict() encuentre los datos
        fit_obj$call$data       <- as.name(env_name_df)
        fit_obj$call$family     <- as.name(env_name_fam)
        fit_obj$sigma.call$data <- as.name(env_name_df)
        fit_obj
      },
      error = function(e) {
        cat(sprintf("[%s->ERR:%s] ", fam_try,
                    gsub("\n","",substr(e$message,1,40))))
        NULL
      })

      if (!is.null(fit_try)) {
        fit      <- fit_try
        fam_used <- fam_try
        # Avisar si hubo fallback
        if (fam_try != fam_candidates[1]) {
          cat(sprintf("[fallback->%s] ", fam_try))
        }
        break
      }
    }

    if (is.null(fit)) { cat("FALLO TOTAL\n"); next }
    fam <- fam_used

    # ── Diagnóstico: residuos normalizados (ya en fit$residuals) ─────────
    pct_out <- mean(abs(fit$residuals) > 1.96, na.rm = TRUE) * 100
    results_diagnostics$pct_outside[p]    <- round(pct_out, 2)
    results_diagnostics$family[p]         <- fam
    results_diagnostics$transform_used[p] <- transform_used
    cat(sprintf("OK (%.1f%%)\n", pct_out))
    n_ok <- n_ok + 1

    # ── Predict mu y sigma sobre la rejilla de edad ───────────────────────
    mu_p <- tryCatch(
      predict(fit, newdata=df_pr, type="response", what="mu",    data=df_tr),
      error = function(e) {
        cat(sprintf("  WARN predict mu: %s\n", e$message))
        rep(NA_real_, N_AGE_GRID)
      }
    )
    sig_p <- tryCatch(
      predict(fit, newdata=df_pr, type="response", what="sigma", data=df_tr),
      error = function(e) {
        cat(sprintf("  WARN predict sigma: %s\n", e$message))
        rep(NA_real_, N_AGE_GRID)
      }
    )

    # ── Predecir nu y tau fuera del loop (una sola vez por variable) ────────
    nu_p  <- NULL
    tau_p <- NULL
    if (fam %in% c("BCT","BCPE")) {
      nu_p <- tryCatch(
        predict(fit, newdata=df_pr, type="response", what="nu",  data=df_tr),
        error = function(e) rep(NA_real_, N_AGE_GRID)
      )
    }
    if (fam == "BCT") {
      tau_p <- tryCatch(
        predict(fit, newdata=df_pr, type="response", what="tau", data=df_tr),
        error = function(e) rep(NA_real_, N_AGE_GRID)
      )
    }

    # ── Centiles en espacio TRANSFORMADO, luego back-transform ────────────
    centile_vals_t <- list()
    for (ci in seq_along(centile_levels)) {
      cv_t <- if (fam == "NO") {
        mu_p + qnorm(centile_levels[ci]) * sig_p
      } else if (fam == "BCT") {
        tryCatch(
          qBCT(centile_levels[ci], mu=mu_p, sigma=sig_p, nu=nu_p, tau=tau_p),
          error = function(e) mu_p + qnorm(centile_levels[ci]) * sig_p
        )
      } else if (fam == "BCPE") {
        tryCatch(
          qBCPE(centile_levels[ci], mu=mu_p, sigma=sig_p, nu=nu_p),
          error = function(e) mu_p + qnorm(centile_levels[ci]) * sig_p
        )
      } else {
        mu_p + qnorm(centile_levels[ci]) * sig_p
      }
      centile_vals_t[[centile_names[ci]]] <- cv_t
      # Back-transform al espacio original antes de guardar
      results_centiles[[paste0(vname, "_", centile_names[ci])]] <-
        apply_inv_transform(cv_t, transform_used)
    }

    # mu y sigma se guardan en espacio TRANSFORMADO (necesario para Z-scores)
    results_centiles[[paste0(vname, "_mu")]]    <- mu_p
    results_centiles[[paste0(vname, "_sigma")]] <- sig_p

    # sigma_eff basado en IQR 2.5–97.5 en espacio transformado
    sig_eff <- (centile_vals_t[["p975"]] - centile_vals_t[["p025"]]) / (2 * 1.96)
    coef_export[[paste0(vname, "_mu")]]    <- mu_p
    coef_export[[paste0(vname, "_sigma")]] <- sig_eff
    if (!is.null(nu_p))  coef_export[[paste0(vname, "_nu")]]  <- nu_p
    if (!is.null(tau_p)) coef_export[[paste0(vname, "_tau")]] <- tau_p


    # ── Limpiar GlobalEnv ─────────────────────────────────────────────────
    rm(list = intersect(c(env_name_fam, env_name_df, env_name_pr,
                          "yy", "aa", "fam"),
                        ls(envir = .GlobalEnv)),
       envir = .GlobalEnv)
  }

  # ── Guardar resultados ────────────────────────────────────────────────────
  write.csv(results_centiles,    file.path(DATA_DIR, paste0(prefix, "_gamlss_centiles.csv")),    row.names=FALSE)
  write.csv(results_diagnostics, file.path(DATA_DIR, paste0(prefix, "_gamlss_diagnostics.csv")), row.names=FALSE)
  write.csv(coef_export,         file.path(DATA_DIR, paste0(prefix, "_gamlss_curves.csv")),      row.names=FALSE)
  # Metadatos del modelo: transform aplicado + nombre del prefix
  write.csv(data.frame(field = c("transform","prefix"),
                       value = c(transform, prefix),
                       stringsAsFactors = FALSE),
            file.path(DATA_DIR, paste0(prefix, "_gamlss_model_meta.csv")),
            row.names=FALSE)

  pct_vals <- results_diagnostics$pct_outside[!is.na(results_diagnostics$pct_outside)]
  cat(sprintf("\n  Ajustadas: %d/%d  mediana: %.1f%%\n  Curvas: %s\n",
              n_ok, P,
              ifelse(length(pct_vals) > 0, median(pct_vals), NA),
              paste0(prefix, "_gamlss_curves.csv")))

  # Limpieza final por si quedó algo
  stale <- ls(envir=.GlobalEnv, pattern=paste0("^\\.gamlss_.*_", prefix, "_"))
  if (length(stale) > 0) rm(list=stale, envir=.GlobalEnv)

  invisible(list(centiles    = results_centiles,
                 diagnostics = results_diagnostics,
                 curves      = coef_export))
}

# ── Ejecutar todos los prefijos ───────────────────────────────────────────────
PREFIXES <- c("AIb", "AIe", "Arousal","Valence","CognAf",
              "TB1325R","TB1325RFC",
              "DB1325R","DB1325RFC",
              "DB1325RF3F4Asym","TB1325RF3F4Asym",
              "DB1325RF7F8Asym","TB1325RF7F8Asym",
              "DB1325RT3T4Asym","TB1325RT3T4Asym",
              "ABR","TBR","TAR","EI",
              "DB1325RT5T6Asym","TB1325RT5T6Asym",
              "AsymIdx","IAF")

for (px in PREFIXES) fit_one_prefix(px)
