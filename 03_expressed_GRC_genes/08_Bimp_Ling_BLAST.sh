#!/bin/bash -l

#SBATCH --job-name=BLASTn
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --export=ALL
#SBATCH --time=0-90:00:00
#SBATCH --partition=ac3-compute
#SBATCH --mem=32gb
#SBATCH --output=BLASTn.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s2673271@ed.ac.uk

# Exit on error
set -e
hostname

# Create scratch directory
SCRATCH=/scratch/${USER}/BLASTn.${SLURM_JOB_ID}
mkdir -p ${SCRATCH}
cd ${SCRATCH}

# Activate conda environment
source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

#################################################
# Copy input files
#################################################

echo "Copying input files..."

# Protein queries from B. coprophila
rsync -av \
/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/outputs/GRC_genes.fasta \
${SCRATCH}/

# B. impatiens genome
rsync -av \
/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_impatiens/idBraImpa2.1.primary.curated.fa \
${SCRATCH}/

# L. ingenua genome
rsync -av \
/mnt/loki/ross/assemblies/flies/sciaridae/Lycoriella_ingenua/idLycInge5.1.primary.masked.fa \
${SCRATCH}/

#################################################
# Make BLAST databases
#################################################

echo "Making B. impatiens BLAST database..."

makeblastdb \
-in idBraImpa2.1.primary.curated.fa \
-dbtype nucl \
-parse_seqids \
-out B_imp_DB

echo "Making L. ingenua BLAST database..."

makeblastdb \
-in idLycInge5.1.primary.masked.fa \
-dbtype nucl \
-parse_seqids \
-out L_ing_DB

#################################################
# Run BLASTn
#################################################

echo "Running BLASTn against B. impatiens genome..."

blastn \
-query GRC_genes.fasta \
-db B_imp_DB \
-out Bimp_Bcop_GRC_blastn.tsv \
-evalue 1e-5 \
-num_threads 16 \
-max_target_seqs 10 \
-outfmt '6 std qlen slen'

echo "Running BLASTn against L. ingenua genome..."

blastn \
-query GRC_genes.fasta \
-db L_ing_DB \
-out Ling_Bcop_GRC_blastn.tsv \
-evalue 1e-5 \
-num_threads 16 \
-max_target_seqs 10 \
-outfmt '6 std qlen slen'

#################################################
# Copy outputs back
#################################################

echo "Syncing results back..."

rsync -av *.tsv \
/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/outputs/

#################################################
# Cleanup
#################################################

echo "Cleaning scratch directory..."

rm -rf ${SCRATCH}

echo "Done."

exit 0
