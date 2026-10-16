stopifnot(as.character(getRversion()) == trimws(readLines('.R-version')))
options(timeout = 1200)
repos <- getOption('repos')
portable_repo <- repos[grepl('manylinux', repos, fixed = TRUE)]
if (length(portable_repo)) {
  options(repos = c(CRAN = unname(portable_repo[1])))
} else {
  options(repos = c(CRAN = 'https://cloud.r-project.org'))
}
# On a non-portable cluster R, RENV_CONFIG_REPOS_OVERRIDE may select a
# platform-appropriate repository for restoring the same locked versions.
if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv')
if (file.exists('renv.lock')) {
  renv::restore(prompt = FALSE)
} else {
  renv::init(bare = TRUE)
  renv::settings$snapshot.type('explicit')
  # renv detects betAS biocViews and requires BiocVersion metadata.
  renv::settings$bioconductor.version('3.22')
  renv::install('DiseaseTranscriptomicsLab/betAS@cd37737c75ab17a64652baf2fdb00d1b7c67c294')
  renv::install('bioc::BiocVersion')
  renv::snapshot(prompt = FALSE)
}
stopifnot(requireNamespace('betAS', quietly = TRUE))
cat('R:', as.character(getRversion()), 'betAS:', as.character(packageVersion('betAS')), '\n')
