#!/bin/bash -l
#SBATCH --job-name=minimap2
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --export=ALL
#SBATCH --time=00-120:00:00
#SBATCH --partition=ac3-compute
#SBATCH --mem=32gb
#SBATCH --output=minimap2.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s2673271@ed.ac.uk

set -e

SCRATCH=/scratch/${USER}/minimap2_${SLURM_JOB_ID}
mkdir -p $SCRATCH && cd $SCRATCH

source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

# Inputs
rsync -av /mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/idBraCopr2.1.primary.masked.fa $SCRATCH
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/pacbio_long_reads/outputs/ERR12736861.fastq.gz $SCRATCH

# Index and align
minimap2 -d idBraCopr2.1.primary.masked.mmi idBraCopr2.1.primary.masked.fa
minimap2 -ax map-pb idBraCopr2.1.primary.masked.mmi ERR12736861.fastq.gz > ERR12736861.sam

# Convert to sorted BAM
samtools view -@ 16 -bS ERR12736861.sam | samtools sort -@ 16 -o ERR12736861.sorted.bam
samtools index ERR12736861.sorted.bam
rm ERR12736861.sam

# Output
rsync -av *.bam* /mnt/loki/ross/flies/sciaridae/GRCs/pacbio_long_reads/outputs/

echo "Job done. Cleaning scratch..."