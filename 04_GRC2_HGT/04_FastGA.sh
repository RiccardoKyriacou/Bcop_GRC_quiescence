#!/bin/bash -l
#SBATCH --job-name=FastGA
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --export=ALL
#SBATCH --time=02-00:00:00
#SBATCH --partition=ac3-compute
#SBATCH --mem=64gb
#SBATCH --output=fastga.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s2673271@ed.ac.uk

set -e

SCRATCH=/scratch/${USER}/FastGA.${SLURM_JOB_ID}
mkdir -p $SCRATCH
cd $SCRATCH

source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

# Inputs
rsync -av /mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/idBraCopr2.1.primary.masked.fa $SCRATCH
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/pacbio_long_reads/outputs/Rickettsiaceae.finalassembly.fa $SCRATCH

# Run FastGA (produce ALN file for plotting)
FastGA \
  Rickettsiaceae.finalassembly.fa \
  idBraCopr2.1.primary.masked.fa \
  -1:grc2_vs_rickettsia.1aln \
  -T16 -vk

# Force 1 to 1 alignment
ALNchain -v -ogrc2_vs_rickettsia_1to1.1aln grc2_vs_rickettsia.1aln
 
# Convert to PAF
ALNtoPAF grc2_vs_rickettsia_1to1.1aln > grc2_vs_rickettsia_1to1.1aln.paf

# Copy results back
rsync -av *.paf /mnt/loki/ross/flies/sciaridae/GRCs/pacbio_long_reads/outputs/FastGA

echo "Cleaning up..."
rm -rf ${SCRATCH}
