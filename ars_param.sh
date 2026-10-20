#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=192
#SBATCH --time=24:00:00
#SBATCH --job-name=ars_param
#SBATCH --output=/scratch/aripoll/athena_out/outputs/ars_param.out
#SBATCH --mail-user=arispollak@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

cd /scratch/aripoll/athena

module purge
module load intel/2023.2.1
module load openmpi/4.1.5
module load hdf5-mpi/1.14.2


python configure.py  --prob=jet_blast --coord=cylindrical -hdf5 -mpi --hdf5_path="$SCRATCH" -s --flux=hlle --nscalars=1
make clean
make -j 192

cd /scratch/aripoll/athena_out/outputs 

# Same E/mc^2 = 5.006e-3 and h_0 as jb_exph; each run is ~4 h, back to back
for run in hirho hiG hitjet lotjet lot0; do
    echo "===== jb_param_$run ====="
    cp /scratch/aripoll/athena/inputs/mhd/athinput.jet_blast_param_$run .
    mpiexec -n 384 /scratch/aripoll/athena/bin/athena -i athinput.jet_blast_param_$run
done
