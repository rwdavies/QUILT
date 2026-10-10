#!/usr/bin/env bash
#
# Build a pre-prepared QUILT2 reference panel (SNPs only) from per-chromosome
# phased VCF/BCF files, e.g. the HGDP+1KG panel from gnomAD.
#
# For every chromosome it
#   1. keeps biallelic SNPs and only the GT field (bcftools)
#   2. splits the chromosome into chunks with QUILT::quilt_chunk_map
#   3. runs QUILT2_prepare_reference.R for the chunks of all chromosomes
#      in one parallel pool, longest chunks first
#   4. writes chunks.tsv, the table used by QUILT2_run.sh
#
# Output layout
#   <outdir>/chunks.tsv                       one row per chunk, all chromosomes
#   <outdir>/<chr>/QUILT_prepared_reference.<chr>.<start>.<end>.RData
#   <outdir>/QUILT2_run.sh                    copy of the imputation runner
#   <outdir>/bcf/<chr>.snps.bcf               SNP-only panel (not needed for imputation)
#
# The script is restartable: finished chunks are skipped.
#
# Example (HGDP+1KG, GRCh38):
#   scripts/build_quilt2_panel.sh \
#     -i 'phased_haplotypes_v2/hgdp1kgp_{CHR}.filtered.SNV_INDEL.phased.shapeit5.bcf' \
#     -m 'maps/hg38/CEU-{CHR}-final.b38.txt.gz' \
#     -c "$(echo chr{1..22})" \
#     -o hgdp1kg_quilt2_snps -t 16

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -i VCF_PATTERN -m MAP_PATTERN -o OUTDIR [options] [-- extra QUILT2_prepare_reference.R args]

Required:
  -i  reference VCF/BCF path, with {CHR} as placeholder for the chromosome
  -m  genetic map path, with {CHR} as placeholder for the chromosome
  -o  output directory

Options:
  -c  chromosomes, space or comma separated  [chr1 ... chr22]
  -t  number of parallel jobs (cores)       [1]
  -b  buffer in bp                           [500000]
  -g  nGen                                   [100]
  -h  show this help

Anything after -- is passed to QUILT2_prepare_reference.R as is,
e.g. -- --reference_exclude_samplelist_file=exclude.txt
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
vcf_pattern=""
map_pattern=""
outdir=""
chrs="$(echo chr{1..22})"
threads=1
buffer=500000
nGen=100

while getopts "i:m:o:c:t:b:g:h" opt; do
    case $opt in
        i) vcf_pattern=$OPTARG ;;
        m) map_pattern=$OPTARG ;;
        o) outdir=$OPTARG ;;
        c) chrs=${OPTARG//,/ } ;;
        t) threads=$OPTARG ;;
        b) buffer=$OPTARG ;;
        g) nGen=$OPTARG ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
[[ "${1:-}" == "--" ]] && shift
EXTRA_ARGS="$*"

if [[ -z "$vcf_pattern" || -z "$map_pattern" || -z "$outdir" ]]; then
    usage
    exit 1
fi

for tool in bcftools Rscript; do
    command -v $tool >/dev/null || { echo "ERROR: $tool not found in PATH" >&2; exit 1; }
done
QUILT2_PREPARE=${QUILT2_PREPARE:-$(command -v QUILT2_prepare_reference.R || echo "$SCRIPT_DIR/../QUILT2_prepare_reference.R")}
[[ -f "$QUILT2_PREPARE" ]] || { echo "ERROR: QUILT2_prepare_reference.R not found, set QUILT2_PREPARE" >&2; exit 1; }

log() { echo "[$(date '+%F %T')] $*" >&2; }

mkdir -p "$outdir/bcf" "$outdir/logs"
outdir=$(cd "$outdir" && pwd)

prepare_chunk() {
    # args: chr chunk start end vcf map
    local chr=$1 chunk=$2 start=$3 end=$4 vcf=$5 map=$6
    local rdata="$outdir/$chr/QUILT_prepared_reference.$chr.$start.$end.RData"
    [[ -s "$rdata" ]] && return 0
    local tmp="$rdata.tmp.RData"
    if Rscript "$QUILT2_PREPARE" \
        --reference_vcf_file="$vcf" \
        --genetic_map_file="$map" \
        --chr="$chr" \
        --regionStart="$start" \
        --regionEnd="$end" \
        --buffer="$buffer" \
        --nGen="$nGen" \
        --outputdir="$outdir/$chr" \
        --output_file="$tmp" \
        $EXTRA_ARGS > "$outdir/logs/prepare.$chr.chunk$chunk.log" 2>&1; then
        mv "$tmp" "$rdata"
        log "done $chr chunk $chunk ($start-$end)"
    else
        log "FAILED $chr chunk $chunk, see $outdir/logs/prepare.$chr.chunk$chunk.log"
        return 1
    fi
}
export -f prepare_chunk log
export outdir buffer nGen QUILT2_PREPARE EXTRA_ARGS

abs_map() { local map=${map_pattern//\{CHR\}/$1}; echo "$(cd "$(dirname "$map")" && pwd)/$(basename "$map")"; }

prepare_chr() {
    # extract biallelic SNPs and split the chromosome into chunks
    local chr=$1
    local in_vcf=${vcf_pattern//\{CHR\}/$chr}
    local map
    map=$(abs_map "$chr")
    local snp_bcf="$outdir/bcf/$chr.snps.bcf"
    mkdir -p "$outdir/$chr"
    if [[ ! -s "$snp_bcf.csi" ]]; then
        log "extracting biallelic SNPs for $chr"
        bcftools view -v snps -m2 -M2 --threads 4 -Ou "$in_vcf" 2>/dev/null |
            bcftools annotate -x INFO,^FORMAT/GT --threads 4 -Ob -o "$snp_bcf.tmp.bcf" &&
            mv "$snp_bcf.tmp.bcf" "$snp_bcf" &&
            bcftools index -f "$snp_bcf" || { log "FAILED extracting SNPs for $chr"; return 1; }
    fi
    local chunks_chr="$outdir/$chr/chunks.txt"
    if [[ ! -s "$chunks_chr" ]]; then
        Rscript -e "d <- QUILT::quilt_chunk_map('$chr', '$map'); \
            r <- do.call(rbind, strsplit(sub('.*:', '', d[['region']]), '-')); \
            write.table(data.frame(d[['chunk']], r), '$chunks_chr.tmp', \
                        quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)" &&
            mv "$chunks_chr.tmp" "$chunks_chr" || { log "FAILED chunking $chr"; return 1; }
    fi
    log "$chr has $(wc -l < "$chunks_chr") chunks"
}
export -f abs_map prepare_chr
export vcf_pattern map_pattern

for chr in $chrs; do
    [[ -f "${vcf_pattern//\{CHR\}/$chr}" ]] || { echo "ERROR: ${vcf_pattern//\{CHR\}/$chr} not found" >&2; exit 1; }
    [[ -f "${map_pattern//\{CHR\}/$chr}" ]] || { echo "ERROR: ${map_pattern//\{CHR\}/$chr} not found" >&2; exit 1; }
done

# step 1: SNP extraction and chunking, several chromosomes at a time (8 threads each)
echo $chrs | tr ' ' '\n' | xargs -P $(( threads / 8 > 1 ? threads / 8 : 1 )) -I{} bash -c 'prepare_chr {}'

# step 2: prepare the chunks of all chromosomes in one pool, longest regions first
for chr in $chrs; do
    awk -v chr="$chr" -v vcf="$outdir/bcf/$chr.snps.bcf" -v map="$(abs_map "$chr")" \
        '{print $3 - $2, chr, $1, $2, $3, vcf, map}' "$outdir/$chr/chunks.txt"
done | sort -k1,1nr | cut -d' ' -f2- |
    xargs -P "$threads" -L 1 bash -c 'prepare_chunk "$@"' _ ||
    { log "some chunks failed, rerun the same command to retry them"; exit 1; }

chunks_all="$outdir/chunks.tsv"
tmp_chunks=$(mktemp)
echo -e "chr\tchunk\tregionStart\tregionEnd\tbuffer\tnGen\tprepared_reference" > "$tmp_chunks"
for chr in $chrs; do
    awk -v chr="$chr" -v b="$buffer" -v g="$nGen" 'BEGIN{OFS="\t"}
        {print chr, $1, $2, $3, b, g, chr"/QUILT_prepared_reference."chr"."$2"."$3".RData"}' \
        "$outdir/$chr/chunks.txt" >> "$tmp_chunks"
done
# keep rows of chromosomes built in earlier runs
if [[ -s "$chunks_all" ]]; then
    tail -n +2 "$chunks_all" | awk -v chrs=" $chrs " 'index(chrs, " "$1" ") == 0' >> "$tmp_chunks"
fi
{ head -n 1 "$tmp_chunks"; tail -n +2 "$tmp_chunks" | sort -k1,1V -k2,2n; } > "$chunks_all"
rm -f "$tmp_chunks"
rmdir "$outdir"/*/RData 2>/dev/null || true  # empty dirs left by QUILT2_prepare_reference.R

cp "$SCRIPT_DIR/../QUILT2_run.sh" "$outdir/"
cat > "$outdir/README.txt" <<EOF
QUILT2 prepared reference panel (biallelic SNPs only)
built $(date '+%F') by build_quilt2_panel.sh from $vcf_pattern
buffer=$buffer nGen=$nGen ${EXTRA_ARGS:+extra QUILT2_prepare_reference.R args: $EXTRA_ARGS}

Impute (BAMs aligned to the same genome build and contig names as the panel):
  ./QUILT2_run.sh -p . -b bamlist.txt -o results -c chr20 -t 8

See https://github.com/rwdavies/QUILT/blob/master/README_QUILT2.org
EOF
log "panel ready in $outdir"
