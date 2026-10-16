# Synthetic four-event test; generates no committed data and submits no real jobs.
source('scripts/common.R')
root <- normalizePath('.')
lib <- .libPaths()[1]
work <- tempfile('spire-single-test-'); dir.create(work)
for (d in c('scripts','config')) dir.create(file.path(work,d))
invisible(file.copy(list.files('scripts', full.names=TRUE), file.path(work,'scripts')))
invisible(file.copy(c('.R-version','renv.lock'), work))
config$input <- 'fixture.tab'
writeLines(c('config <-', capture.output(dput(config))), file.path(work,'config/analysis.R'))
raw <- data.frame(GENE = paste0('Gene', 1:4),
  EVENT = c('MmuEX0000001','MmuINT0000001','MmuALTA0000001-1/2','MmuEX6000001'),
  COORD = 'chr1:100-199', LENGTH = 100L, FullCO = 'chr1:1,100-199,300',
  COMPLEX = c('S','IR','Alt3','ANN'))
inc <- matrix(c(10,30,60,100, 15,35,65,95, 12,32,62,98,
                45,10,80,100, 40,12,85,98, 50,15,75,95), nrow=4)
for (i in seq_along(config$sample_groups)) {
  s <- names(config$sample_groups)[i]
  raw[[s]] <- inc[,i]
  raw[[paste0(s,'-Q')]] <- sprintf('OK,OK,0=0=0,OK,S@%s,%s',inc[,i],100-inc[,i])
}
write.table(raw,file.path(work,config$input),sep='\t',quote=FALSE,row.names=FALSE)
setwd(work)
run <- function(script,args=character(),expected=0L) {
  status <- system2(file.path(R.home('bin'),'Rscript'),
    c('--vanilla',shQuote(file.path('scripts',script)),shQuote(args)),
    env=paste0('R_LIBS=',shQuote(lib)))
  stopifnot(status == expected)
}
run('01_filter.R')
run('00_init_run.R','results/test')
run('03_collect.R','results/test',expected=1L)
run('02_simulate.R','results/test')
run('03_collect.R','results/test')
whole <- readRDS('results/test/statistics.rds')$results
stopifnot(identical(whole$EVENT,raw$EVENT), file.exists('results/test/COMPLETE'),
          nrow(read.csv('results/test/csv/spire_all_events.csv')) == 4L)
# Prove the upstream FDR result is unchanged by splitting the row set.
x <- readRDS(config$filtered)
reference <- new.env(); data('maxDevSimulationN100',package='betAS',envir=reference)
separate <- lapply(seq_len(nrow(x$PSI)),function(i)
  betAS::prepareTableVolcanoFDR(psitable=x$PSI[i,,drop=FALSE],
    qualtable=x$Qual[i,,drop=FALSE],npoints=500,colsA=7:9,colsB=10:12,
    labA='Cont',labB='DKO',maxDevTable=reference$maxDevSimulationN100,
    nsim=1000,seed=TRUE,CoverageWeight=FALSE))
stopifnot(identical(whole$FDR,do.call(rbind,separate)$FDR))
checkpoints <- file.path('results/test',paste0(c('FDR','Pdiff','Fstat'),'.rds'))
before <- file.info(checkpoints)$mtime
run('02_simulate.R','results/test')
stopifnot(identical(before,file.info(checkpoints)$mtime),
          identical(whole,readRDS('results/test/statistics.rds')$results))
# Submission must invoke sbatch once, without array/dependency options.
dir.create('bin')
writeLines(c('#!/bin/sh', 'printf "%s\\n" "$*" >> "$SPIRE_SBATCH_CALLS"',
             'echo 12345'), 'bin/sbatch'); Sys.chmod('bin/sbatch','0755')
status <- system2('bash',c('scripts/submit_crg.sh','results/submission'),
  env=c(paste0('PATH=',shQuote(paste(normalizePath('bin'),Sys.getenv('PATH'),sep=':'))),
        paste0('R_LIBS=',shQuote(lib)),
        paste0('SPIRE_SBATCH_CALLS=',shQuote(file.path(work,'sbatch_calls.txt')))))
calls <- readLines('sbatch_calls.txt')
stopifnot(status == 0L,length(calls) == 1L,
          !grepl('--array|--dependency',calls),grepl('scripts/betas.sbatch',calls))
z <- readRDS('results/test/FDR.rds');z$results$EVENT[1] <- 'INVALID'
saveRDS(z,'results/test/FDR.rds')
run('02_simulate.R','results/test',expected=1L)
z <- readRDS('results/test/statistics.rds');z$results$EVENT[1] <- 'INVALID'
saveRDS(z,'results/test/statistics.rds')
run('03_collect.R','results/test',expected=1L)
stopifnot(!file.exists('results/test/COMPLETE'))
setwd(root)
cat('PASS: full-table export, split/whole FDR equality, checkpoint reuse, single sbatch submission and corrupt-result rejection.\n')
