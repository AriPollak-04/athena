#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=192
#SBATCH --time=6:00:00
#SBATCH --job-name=ars_exp
#SBATCH --output=/scratch/aripoll/athena_out/outputs/ars_exp.out
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

cp /scratch/aripoll/athena/inputs/mhd/athinput.jet_blast_exp .

mpiexec -n 384 /scratch/aripoll/athena/bin/athena -i athinput.jet_blast_exp

