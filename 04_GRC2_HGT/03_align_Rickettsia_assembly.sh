#!/bin/bash -l

#SBATCH --job-name=minimap2
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --export=ALL
#SBATCH --time=0-90:00:00
#SBATCH --partition ac3-compute
#SBATCH --mem=32gb
#SBATCH --output=minimap2.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s2673271@ed.ac.uk

set -e

SCRATCH=/scratch/${USER}/minimap2.${SLURM_JOB_ID}
mkdir -p $SCRATCH
cd $SCRATCH

source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

echo "Syncing files..."
rsync -av /mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/idBraCopr2.1.primary.masked.fa $SCRATCH/
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/pacbio_long_reads/outputs/Rickettsiaceae.finalassembly.fa $SCRATCH/

# MINIMAP2 Preset info: 
# Same species or highly similar genomes    -x asm20
# Diverged genomes (e.g., different species, distant homologs)	-x asm5

echo "Running minimap2 -x asm5..."
minimap2 -ax asm5 idBraCopr2.1.primary.masked.fa Rickettsiaceae.finalassembly.fa > Bcop_Rickettsia.sam

echo "Converting to BAM, sorting, and indexing..."
samtools view -bS Bcop_Rickettsia.sam | samtools sort -o Bcop_Rickettsia.sorted.bam
samtools index Bcop_Rickettsia.sorted.bam

echo "Syncing mapped files and results back..."
rsync -av *.bam* /mnt/loki/ross/flies/sciaridae/GRCs/pacbio_long_reads/outputs

echo "Cleaning up..."
rm -rf *

rm -rf ${SCRATCH}
