#!/bin/bash -l

#SBATCH --job-name=complete_GRC_AI
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --export=ALL
#SBATCH --time=0-02:00:00
#SBATCH --partition=ac3-compute
#SBATCH --mem=48gb
#SBATCH --output=complete_GRC_AI.%j.log

set -e
hostname

source /home/s2673271/miniforge3/etc/profile.d/conda.sh
conda activate /home/s2673271/miniforge3/envs/genomics

WORKDIR=/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/outputs/cecidomyiid_v_sciarid
mkdir -p ${WORKDIR}
cd ${WORKDIR}

THREADS=${SLURM_CPUS_PER_TASK}

# ── 1) Download outgroup ───────────────────────────────────────────────────────────────
# echo "Downloading Anopheles gambiae (AgamP4) RefSeq protein set..."
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/943/734/735/GCF_943734735.2_idAnoGambNW_F1_1/GCF_943734735.2_idAnoGambNW_F1_1_protein.faa.gz

# gunzip GCF_943734735.2_idAnoGambNW_F1_1_protein.faa.gz

OUTGROUP_AA=/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/outputs/cecidomyiid_v_sciarid/GCF_943734735.2_idAnoGambNW_F1_1_protein.faa
SCI_AA=/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/outputs/cecidomyiid_v_sciarid/Bcop_core_proteins.faa
CEC_AA=/mnt/loki/ross/flies/cecidomyiidae/Aphidoletes_aphidimyza/aphid_midge_inversion/15_braker/outputs/M_producer_ragtag.Aphidoletes_aphidimyza.braker/Aaph_braker.aa

GRC_AA=/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/outputs/cecidomyiid_v_sciarid/Bcop_all_GRC_proteins.faa

# ── 2) Three way BLASTp  ───────────────────────────────────────────────────────────────
SCRATCH=/scratch/${USER}/BLASTp_3way_GRC.${SLURM_JOB_ID}
mkdir -p ${SCRATCH}
cd ${SCRATCH}

echo "Copying input files..."

echo "Building Cecidomyiid (Aphidoletes aphidimyza) DB..."
makeblastdb -in $CEC_AA -dbtype prot -parse_seqids -out Cecidomyiid_DB

echo "Building Sciarid core (GRC-excluded) DB..."
makeblastdb -in $SCI_AA -dbtype prot -parse_seqids -out Sciarid_core_DB

echo "Building Outgroup (Anopheles gambiae) DB..."
makeblastdb -in $OUTGROUP_AA -dbtype prot -parse_seqids -out Outgroup_DB

rsync -av "${GRC_AA}" ${SCRATCH}/      # All 15,000 translated GRC queries
GRC_AA=${SCRATCH}/Bcop_all_GRC_proteins.faa

run_blast () {
  local db=$1
  local out=$2
  echo "Running BLASTp vs ${db}..."
  blastp \
    -query $GRC_AA \
    -db ${db} \
    -out ${out} \
    -evalue 10 \
    -num_threads $THREADS \
    -max_target_seqs 10 \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs'
}

run_blast Cecidomyiid_DB  All_GRC_vs_Cecidomyiid.tsv
run_blast Sciarid_core_DB All_GRC_vs_SciaridCore.tsv
run_blast Outgroup_DB     All_GRC_vs_Outgroup.tsv

echo "Syncing results back..."
rsync -av *.tsv ${WORKDIR}/

echo "Cleaning scratch directory..."
rm -rf ${SCRATCH}

# ── 4) get Similarity Index  ───────────────────────────────────────────────────────────────
cd ${WORKDIR}
SCRIPT=/mnt/loki/ross/flies/sciaridae/GRCs/GRC_expression/Bradysia_coprophila/03_BLAST_v2/scripts/alien_index.py

python $SCRIPT \
        --cecidomyiid All_GRC_vs_Cecidomyiid.tsv \
        --sciarid     All_GRC_vs_SciaridCore.tsv \
        --outgroup    All_GRC_vs_Outgroup.tsv\
        --out         Full_GRC_proteome_AI_summary.csv

echo "Done."
