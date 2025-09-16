# Load required modules
module purge
module load star

# Define paths to mouse genome fasta and annotation files
fa_dir="/your/reference/genomes/mouse_fasta/mm10.fa"
gtf_dir="/your/reference/genomes/gencode.vM10.annotation.gtf"
genome_dir="/your/reference/genomes/star_index_mm10"

# Create genome directory if it doesn't exist
mkdir -p "${genome_dir}"

# Run STAR genomeGenerate for mouse mm10
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir "${genome_dir}" \
     --genomeFastaFiles "${fa_dir}" \
     --sjdbGTFfile "${gtf_dir}" \
     --sjdbOverhang 150  # RNA-seq read-length from Novogene is 151bp