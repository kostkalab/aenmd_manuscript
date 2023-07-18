#- Running VEP to predict variant impact;
#  to be used to annotate NMD- variants

set -e
set -u

VEP_CACHE=~/.vep/
VEP_FASTA=~/.vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

INPUT=$1
OUTPUT=$INPUT.vep_out.txt

vep     -i $INPUT               \
        -o $OUTPUT              \
        --offline               \
	--cache_version 105	\
        --dir_cache $VEP_CACHE  \
        --fasta $VEP_FASTA      \
        --no_stats              \
        --coding_only           \
        --shift_hgvs 0          \
        --hgvs                  \
        --tsl                   \
        --symbol                \
        --gencode_basic         \
        --plugin Downstream     \
        --plugin NMD            && gzip $OUTPUT
