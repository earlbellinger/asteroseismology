for a in $(ls maybe_sub_logs | grep calibrate); do tail maybe_sub_logs/$a/condor.out -n 1; done

