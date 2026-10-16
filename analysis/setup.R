# Additional plotting/GO packages are isolated from the running simulation environment.
dir.create('analysis/renv/library', recursive = TRUE, showWarnings = FALSE)
renv::restore(project = 'analysis', lockfile = 'analysis/renv.lock',
  library = normalizePath('analysis/renv/library'), prompt = FALSE)
