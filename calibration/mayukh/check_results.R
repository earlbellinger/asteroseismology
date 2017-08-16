dirs <- list.dirs('models', recursive=F)
for (dir in dirs) {
    DF <- read.table(file.path(dir, 'LOGS_MS', 'history.data'),
        header=1, skip=5)
    DF <- DF[nrow(DF),]
    cat(paste0(dir, '\t', signif(DF$log_g, 4), '\t', 
                          signif(10**DF$log_Teff, 4), '\n'))
}
