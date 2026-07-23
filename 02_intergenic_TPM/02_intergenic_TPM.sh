#!/bin/bash -l

#SBATCH --job-name=intergenic_TPM
#SBATCH --nodes=1
#SBATCH --nodelist=ac3-n6
#SBATCH --ntasks=16
#SBATCH --export=ALL
#SBATCH --time=0-90:00:00
#SBATCH --partition=ac3-compute
#SBATCH --mem=32gb
#SBATCH --output=intergenic_TPM.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s2673271@ed.ac.uk

set -e  # Exit on error

# Define scratch directory
SCRATCH=/scratch/${USER}/intergenic.260539
mkdir -p $SCRATCH
cd $SCRATCH

# Activate conda
source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

# Sync uniquely mapped BAM files and their indexes
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/01_STAR_Relaxed_Nmax2/outputs/*_uniquely_mapped.bam* $SCRATCH
# Sync intergenic BED file
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/02_TPM_cutoff_intergenic/outputs/bcop_core_GRC.intergenic.gtf $SCRATCH

# Run StringTie for TPM estimation of intergenic regions 
echo "Running StringTie for TPM estimation..."
for file in $(ls *_uniquely_mapped.bam)
do
    base=$(basename "$file" "_uniquely_mapped.bam")
    output_gtf="${base}_intergenic.gtf"
    echo "Calculating TPM for $file"
    stringtie "$file" -p 16 -o "$output_gtf" -G bcop_core_GRC.intergenic.gtf
done

# Sync mapped files and results back
echo "Syncing mapped files and results back..."
rsync -av *.gtf /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/02_TPM_cutoff_intergenic/outputs/

# Remove input files
echo "Removing input files"
rm -rf *

# Clear and delete scratch directory
rm -rf ${SCRATCH}