

directories <- c('legacy', 'kages') 
DF <- do.call(rbind, Map(function(directory) { 
        fnames <- list.files(file.path('data', directory)) 
        fnames <- fnames[grepl('freqs', fnames)] 
        do.call(rbind, Map(function(fname) { 
            freqs <- read.table(file.path('data', directory, fname), 
                header=1)
            data.frame(N=nrow(freqs), 
                sigma.mean=mean(freqs$dnu),
                sigma.med=median(freqs$dnu),
                sigma.min=min(freqs$dnu))
        }, fname=fnames))
    }, directory=directories)) 

DF$N.rank <- rank(-DF$N, ties.method='min')
DF$med.rank <- rank(DF$sigma.med, ties.method='min')
DF$min.rank <- rank(DF$sigma.min, ties.method='min')


