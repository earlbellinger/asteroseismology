#!/bin/bash
#### Set up slurm job for a given (bash) command 
#### Author: Earl Patrick Bellinger ( bellinger@phys.au.dk ) 

slurm_sub() {
    cmd=$*
    
    name=${cmd//[-. \/]/_}
    ## Make directory for script
    dname="slurm_logs/$name"
    mkdir -p $dname
    cd $dname
    
    ## Create and submit a slurm job 
    echo "#!/bin/sh

#SBATCH --job-name=$name
#SBATCH --partition=$MACHINE
#SBATCH --ntasks=$OMP_NUM_THREADS
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=sbatch.out
#SBATCH --time=$TIME:00:00
#SBATCH --mem=$MEMORY

export PYTHONUNBUFFERED=1

echo `date`

cd ../..
$cmd >| $dname/out

echo `date`
" > "sbatch.sh"
    chmod +x "sbatch.sh"
    sbatch "sbatch.sh"
    cd -
}

## Parse command line arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    -h) HELP=1; break;;
    -p) OMP_NUM_THREADS="$2"; shift 2;;
    -c) MACHINE="$2"; shift 2;;
    -m) MEMORY="$2"; shift 2;;
    -t) TIME="$2"; shift 2;;
    
     *) break;;
  esac
done

if [ -z ${HELP+x} ]; then HELP=0; fi
if [ -z ${OMP_NUM_THREADS+x} ]; then OMP_NUM_THREADS=1; fi
if [ -z ${MACHINE+x} ]; then MACHINE="q20"; fi
if [ -z ${MEMORY+x} ]; then MEMORY="16G"; fi
if [ -z ${TIME+x} ]; then TIME="240"; fi

if [ $HELP -gt 0 ]; then
    echo "  slurm_sub.sh: automatic job queuer for the slurm"
    echo "                queueing system"
    echo
    echo "example:"
    echo "  slurm_sub.sh -p 16 echo hi"
    echo
    echo "flags:"
    echo "  -h   : show this helpful message and quit"
    echo "  -p # : set the number of threads to run in parallel"
    echo "  -m # : memory; default: 16G"
    echo "  -t # : time (in hours) default: 240"
    echo "  -c s : set the machine that this job will run on"
    echo
    exit
fi

slurm_sub $*
