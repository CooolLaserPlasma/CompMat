#!/usr/bin/env bash
#SBATCH -p hebbe 
#SBATCH -A SNIC2020-5-51 # Project
#SBATCH -J AIMD # Name of the job
#SBATCH -N 1 # Use 1 node
#SBATCH -n 10 # Use 10 cores on that node
#SBATCH -t 35:00:00 # Maximum time
#SBATCH -o stdout # stdout goes to this file
#SBATCH -e stderr # stderr goes to this file

module purge
module load intel/2019a GPAW ASE


mpirun gpaw-python task1.py


