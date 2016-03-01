#### Determine which models should be passed to ADIPLS 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(lpSolve)

args <- commandArgs(TRUE)
num_points <- if (length(args)>0) as.integer(args[1]) else 200

# get history file and remove PMS 
log_dir <- file.path("LOGS")
DF <- read.table(file.path(log_dir, 'history.data'), header=TRUE, skip=5)
decreasing_L <- which(diff(DF$L) < 0 & DF$age[-1] < 0.25)
if (any(decreasing_L)) {
    goes_back_up <- diff(decreasing_L) > 1
    pms <- ifelse(any(goes_back_up), 
        which(goes_back_up)[1] + 1, 
        max(decreasing_L))
    DF <- DF[-1:-pms,]
}

## enumerate log files 
# check that the profiles also have associated fgongs 
logs <- list.files(log_dir)
if (length(logs) <= 1) return(NULL)
profile_candidates <- logs[grep('profile.+.data$', logs)]
freq_file_candidates <- logs[grep('profile.+.data.FGONG$', logs)]
has_fgong <- sub('.FGONG$', '', freq_file_candidates) %in% profile_candidates
profiles <- profile_candidates[has_fgong]

# get model numbers from each profile
model_candidates <- do.call(rbind, Map(function(filename) {
        prof <- read.table(filename, header=1, nrows=1, skip=1)
        data.frame(filename=filename,
            model_number=prof$model_number,
            center_h1=prof$center_h1)
    }, file.path(log_dir, profiles)))
in_main_sequence <- model_candidates$model_number %in% DF$model_number
model_candidates <- model_candidates[in_main_sequence,]

x <- model_candidates$center_h1

nrow.DF <- length(x)
if (nrow.DF < num_points) return(NULL)
ideal <- seq(max(x), min(x), length=num_points)
cost.mat  <- outer(ideal, x, function(x, y) abs(x-y))
row.signs <- rep("==", num_points)
row.rhs   <- rep(1, num_points)
col.signs <- rep("<=", nrow.DF)
col.rhs   <- rep(1, nrow.DF)
sol <- lp.transport(cost.mat, "min", row.signs, row.rhs,
col.signs, col.rhs)$solution
new.DF <- model_candidates[apply(sol, 1, which.max),]
for (fname in new.DF$filename) cat(paste0(fname, ".FGONG\n"))

