
stopifnot(
  as.character(getRversion()) == trimws(readLines(".R-version"))
)

options(timeout = 1200)

repos <- getOption("repos")
portable_repo <- repos[grepl("manylinux", repos, fixed = TRUE)]

if (length(portable_repo)) {
  options(repos = c(CRAN = unname(portable_repo[1])))
} else {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Install renv if necessary
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

if (file.exists("renv.lock")) {

  # Read pinned package versions
  lock <- renv::lockfile_read("renv.lock")

  # Build igraph from source using the locked version
  if (!is.null(lock$Packages$igraph)) {

    igraph_version <- lock$Packages$igraph$Version
    message("Building igraph ", igraph_version, " from source")

    # Force source packages only for this restore
    old_options <- options(pkgType = "source")

    tryCatch(
      renv::restore(
        packages = "igraph",
        rebuild = "igraph",
        repos = c(CRAN = "https://cloud.r-project.org"),
        prompt = FALSE
      ),
      finally = options(old_options)
    )
  }

  # Restore all remaining packages normally
  renv::restore(prompt = FALSE)

} else {

  renv::init(bare = TRUE)
  renv::settings$snapshot.type("explicit")

  # betAS requires Bioconductor metadata
  renv::settings$bioconductor.version("3.22")

  renv::install(
    "DiseaseTranscriptomicsLab/betAS@cd37737c75ab17a64652baf2fdb00d1b7c67c294"
  )

  renv::install("bioc::BiocVersion")

  renv::snapshot(prompt = FALSE)
}

# Verify installation
stopifnot(requireNamespace("betAS", quietly = TRUE))

cat(
  "R:", as.character(getRversion()),
  "betAS:", as.character(packageVersion("betAS")),
  "\n"
)
