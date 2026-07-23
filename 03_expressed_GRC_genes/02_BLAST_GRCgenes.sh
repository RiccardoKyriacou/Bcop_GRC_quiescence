#!/bin/bash -l

#SBATCH --job-name=BLAST
#SBATCH --nodes=1 #number of nodes requested
#SBATCH --ntasks=16 #number of threads per node
#SBATCH --export=ALL
#              d-hh:mm:ss
#SBATCH --time=0-90:00:00 # Upper time limit for job
#SBATCH --partition ac3-compute
#SBATCH --mem=32gb # How much memory you need.
#SBATCH --output=BLAST.%j.log   # Standard output and error log
#SBATCH --mail-type=END,FAIL # Turn on mail notification
#SBATCH --mail-user=s2673271@ed.ac.uk  # Where to send mail

# Exit script on error
set -e

# Define and create a unique scratch directory for this job
SCRATCH=/scratch/${USER}/BLAST.${SLURM_JOB_ID}
mkdir -p $SCRATCH
cd $SCRATCH

# Activate the conda environment for the job
source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

# Sync files in
# Core-genome genes to form BLASTDB
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST/outputs/core_genes.fasta $SCRATCH
# GRC genes to BLAST 
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST/outputs/GRC_genes_to_BLAST.fasta $SCRATCH

### STEP 1) BLAST GRC genes against core genes ###

# Make a blast data base of core-chromosome genes 
echo "making BLASTDB from core chromosome genes"
makeblastdb \
-in core_genes.fasta \
-dbtype nucl \
-parse_seqids \
-out core_genes_DB

# Blast GRC genes against core genes 
echo "BLASTn for GRC genes"
blastn \
-query GRC_genes_to_BLAST.fasta \
-db core_genes_DB \
-out GRC_gene_BLAST_output.tsv \
-outfmt '6 std qlen slen qseq sseq'


### STEP 2) BLAST GRC genes against masked core genome ###

# Sync files in
# GRC genome
rsync -av /mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/idBraCopr2.1.primary.masked.fa $SCRATCH
# Masking script
rsync -av /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST/scripts/mask_genome.py $SCRATCH

# Make masked GRC genome
python3 mask_genome.py -g idBraCopr2.1.primary.masked.fa -o masked_GRC_genome.fasta -m GRC

# Make a blast data base of core-chromosome genes 
echo "making BLASTDB from maksed genome"
makeblastdb \
-in masked_GRC_genome.fasta \
-dbtype nucl \
-parse_seqids \
-out masked_GRC_genome_DB

# Blast GRC genes against whole genome 
echo "BLASTn for GRC genes"
blastn \
-query GRC_genes_to_BLAST.fasta \
-db masked_GRC_genome_DB \
-out GRC_genome_BLAST_output.tsv \
-outfmt '6 std qlen slen qseq sseq'


# Sync outputs out
echo "Syncing results back..."
rsync -av $SCRATCH/*BLAST_output.tsv /mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST/outputs

# Remove input files
echo "Removing input files"
rm -rf ${SCRATCH}/* 

# Clear and delete scratch
rm -rf ${SCRATCH}

# Finish the script
exit 0
