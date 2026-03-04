#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=16:00:00
#SBATCH --job-name=ars
#SBATCH --output=/scratch/aripoll/athena_out/outputs/ars.out
#SBATCH --mail-user=ari.pollak@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

cd /scratch/aripoll/athena

module purge
module load intel/2023.2.1
module load openmpi/4.1.5
module load hdf5-mpi/1.14.2


python configure.py --prob=jet_blast --coord=cylindrical -hdf5 -mpi --hdf5_path="$SCRATCH" -s --flux=hlle
make clean
make -j 192

cd /scratch/aripoll/athena_out/outputs 

cp /scratch/aripoll/athena/inputs/mhd/athinput.jet_blast .

mpiexec -n 192 /scratch/aripoll/athena/bin/athena -i athinput.jet_blast

