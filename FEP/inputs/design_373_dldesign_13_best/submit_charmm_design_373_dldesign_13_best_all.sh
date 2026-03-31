#!/bin/bash
# Submit all jobs for design_373_dldesign_13_best
# Generated: 2026-03-15 18:25:31

# Step 1: Submit setup job
SETUP_JOB=$(sbatch --parsable /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup/submit_charmm_design_373_dldesign_13_best_setup.slurm)
echo "Setup job submitted: $SETUP_JOB"

# Step 2: Submit 3 replica jobs in parallel (all depend on setup completion)
REPLICA1_JOB=$(sbatch --parsable --dependency=afterok:$SETUP_JOB /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/submit_charmm_design_373_dldesign_13_best_r1.slurm)
echo "Replica 1 job submitted (parallel): $REPLICA1_JOB"
REPLICA2_JOB=$(sbatch --parsable --dependency=afterok:$SETUP_JOB /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/submit_charmm_design_373_dldesign_13_best_r2.slurm)
echo "Replica 2 job submitted (parallel): $REPLICA2_JOB"
REPLICA3_JOB=$(sbatch --parsable --dependency=afterok:$SETUP_JOB /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/submit_charmm_design_373_dldesign_13_best_r3.slurm)
echo "Replica 3 job submitted (parallel): $REPLICA3_JOB"

echo "All jobs submitted for design_373_dldesign_13_best"
echo "3 replicas will run in parallel after setup completes"
echo "Monitor with: squeue -u $USER"
