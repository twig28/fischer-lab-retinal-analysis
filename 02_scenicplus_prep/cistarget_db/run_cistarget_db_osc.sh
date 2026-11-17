#!/bin/bash
#SBATCH --job-name=cistarget_db
#SBATCH --account=PCON0240       
#SBATCH --time=48:00:00            
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8          
#SBATCH --mem=32G                  

module load python/3.10            # or whatever python module OSC has

source $HOME/cistarget_env/bin/activate

export PATH=/fs/ess/PCON0240/cisTarget_Gallus:$PATH  # make cbust visible


# Safety: don't let libraries oversubscribe threads
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd /fs/ess/PCON0240/cisTarget_Gallus

echo "Running cisTarget motif database build on OSC"
date

python create_cisTarget_databases/create_cistarget_motif_databases.py \
  -f consensus_regions.fa \
  -M v10nr_clust_public/singletons \
  -m motifs_chicken_expr_for_db.lst \
  -o gallus_consensus_chicken_expr \
  --bgpadding 0 \
  -t 8

echo "Done."
date