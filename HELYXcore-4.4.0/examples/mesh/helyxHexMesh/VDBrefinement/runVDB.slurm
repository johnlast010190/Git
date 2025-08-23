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
#SBATCH --ntasks-per-node=2
#
# No. tasks per MPI rank
#SBATCH --cpus-per-task=18
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
#SBATCH --job-name=helyxVDBMesh
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
sed -i "s/snap 1/snap 0/" system/helyxHexMeshDict
sed -i "s/VDBrefinement 0/VDBrefinement 1/" system/helyxHexMeshDict

srun $mpioptions helyxHexMesh $myoptions -decomposeParDict system/decomposeParDict_4 2>&1 | tee log.helyxHexMesh_VDB.$SLURM_JOB_ID.out

export mpioptions="--ntasks-per-node=36 --cpus-per-task=1 --kill-on-bad-exit --cpu-bind=verbose,cores"

sed -i "s/castellatedMesh 1/castellatedMesh 0/" system/helyxHexMeshDict
sed -i "s/snap 0/snap 1/" system/helyxHexMeshDict

srun --ntasks=72 $mpioptions redistributePar -overwrite $myoptions 2>&1 | tee log.redistributePar.$SLURM_JOB_ID.out
srun --ntasks=72 $mpioptions helyxHexMesh $myoptions 2>&1 | tee log.helyxHexMesh_snap.$SLURM_JOB_ID.out
srun --ntasks=36 $mpioptions checkMesh -allTopology -allGeometry $myoptions 2>&1 | tee log.checkMesh.$SLURM_JOB_ID.out
