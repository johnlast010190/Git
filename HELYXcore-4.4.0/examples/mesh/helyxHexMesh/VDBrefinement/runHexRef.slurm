#!/bin/bash

#======================================================
#
# Job script for running Helyx on multiple nodes
#
#======================================================
#======================================================
# Queue to run on
#SBATCH --partition=standard
#
#SBATCH --qos=commercial
#SBATCH --account=i219
#
# No. of Nodes
#SBATCH --nodes=2
#
# No. MPI ranks per node
#SBATCH --ntasks-per-node=36
#
# No. tasks per MPI rank
#SBATCH --cpus-per-task=1
#
# Distribute processes in round-robin fashion
#SBATCH --distribution=cyclic
#
# Ensure the node is not shared with another job
#SBATCH --exclusive
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=2:00:00
#SBATCH  --mem 250000MB
#
# Job name
#SBATCH --job-name=helyxHexMesh
#
# Output file
#SBATCH --output=slurm-%j.out
#======================================================

# Set environment vars

source /work/i219/i219/cetrarg/.bashrc
compileenv
source /work/i219/i219/cetrarg/HELYXcore-dev/platforms/activeBuild.shrc

cd $SLURM_SUBMIT_DIR

export mpioptions="--kill-on-bad-exit --cpu-bind=verbose,cores"
export myoptions="-parallel"

sed -i "s/castellatedMesh 0/castellatedMesh 1/" system/helyxHexMeshDict
sed -i "s/snap 0/snap 1/" system/helyxHexMeshDict
sed -i "s/VDBrefinement 1/VDBrefinement 0/" system/helyxHexMeshDict

blockMesh 2>&1 | tee log.blockMesh.$SLURM_JOB_ID.out
decomposePar 2>&1 |tee log.decomposePar.$SLURM_JOB_ID.out

srun $mpioptions helyxHexMesh $myoptions 2>&1 | tee log.helyxHexMesh_HexRef.$SLURM_JOB_ID.out
srun $mpioptions checkMesh -allTopology -allGeometry $myoptions 2>&1 | tee log.checkMesh_HexRef.$SLURM_JOB_ID.out
