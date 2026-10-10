#!/usr/bin/env bash
#
# Impute a list of BAM/CRAM files with QUILT2 using a pre-prepared reference
# panel (e.g. the HGDP+1KG panel built by scripts/build_quilt2_panel.sh).
#
# Runs QUILT2.R for every chunk of the requested chromosomes (in parallel),
# then ligates the chunks into one VCF per chromosome:
#   <outdir>/quilt2.<chr>.vcf.gz
#
# The script is restartable: finished chunks are skipped.
#
# Example:
#   ./QUILT2_run.sh -p hgdp1kg_quilt2_snps -b bamlist.txt -o results -c chr20 -t 8

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -p PANEL_DIR -b BAMLIST -o OUTDIR [options] [-- extra QUILT2.R args]

Required:
  -p  directory of the prepared panel (contains chunks.tsv)
  -b  file with one BAM path per line (use -- --cramlist=... --reference=... for CRAM)
  -o  output directory

Options:
  -c  chromosomes, space or comma separated  [all chromosomes in the panel]
  -t  number of chunks imputed in parallel   [1]
  -k  keep per-chunk VCFs after ligation
  -h  show this help

Anything after -- is passed to QUILT2.R as is, e.g.
  -- --nGibbsSamples=7 --n_seek_its=3
  -- --method=nipt --fflist=fflist.txt
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
panel=""
bamlist=""
outdir=""
chrs=""
threads=1
keep=0

while getopts "p:b:o:c:t:kh" opt; do
    case $opt in
        p) panel=$OPTARG ;;
        b) bamlist=$OPTARG ;;
        o) outdir=$OPTARG ;;
        c) chrs=${OPTARG//,/ } ;;
        t) threads=$OPTARG ;;
        k) keep=1 ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
[[ "${1:-}" == "--" ]] && shift
EXTRA_ARGS="$*"

if [[ -z "$panel" || -z "$outdir" ]] || [[ -z "$bamlist" && "$EXTRA_ARGS" != *--cramlist* ]]; then
    usage
    exit 1
fi
[[ -s "$panel/chunks.tsv" ]] || { echo "ERROR: $panel/chunks.tsv not found" >&2; exit 1; }
[[ -z "$bamlist" || -s "$bamlist" ]] || { echo "ERROR: $bamlist not found" >&2; exit 1; }

for tool in bcftools Rscript; do
    command -v $tool >/dev/null || { echo "ERROR: $tool not found in PATH" >&2; exit 1; }
done
QUILT2=${QUILT2:-$(command -v QUILT2.R || echo "$SCRIPT_DIR/QUILT2.R")}
[[ -f "$QUILT2" ]] || { echo "ERROR: QUILT2.R not found, put it in PATH or set QUILT2=/path/to/QUILT2.R" >&2; exit 1; }

log() { echo "[$(date '+%F %T')] $*" >&2; }

panel=$(cd "$panel" && pwd)
mkdir -p "$outdir/chunks" "$outdir/logs"
outdir=$(cd "$outdir" && pwd)
[[ -n "$bamlist" ]] && bamlist=$(cd "$(dirname "$bamlist")" && pwd)/$(basename "$bamlist")
[[ -z "$chrs" ]] && chrs=$(tail -n +2 "$panel/chunks.tsv" | cut -f1 | uniq | tr '\n' ' ')

# make sure the prepared reference of every requested chromosome is there before starting
for chr in $chrs; do
    missing=$(awk -v chr="$chr" 'NR > 1 && $1 == chr {print $7}' "$panel/chunks.tsv" |
        while read -r f; do [[ -s "$panel/$f" ]] || echo "$f"; done | wc -l)
    if [[ $missing -gt 0 ]]; then
        echo "ERROR: $missing prepared reference files of $chr are missing in $panel/$chr/." >&2
        echo "       Download that chromosome, or choose the chromosomes you have with -c" >&2
        exit 1
    fi
done

impute_chunk() {
    # args: chr chunk start end buffer nGen rdata
    local chr=$1 chunk=$2 start=$3 end=$4 buffer=$5 nGen=$6 rdata=$7
    local vcf="$outdir/chunks/quilt2.$chr.chunk$chunk.vcf.gz"
    [[ -s "$vcf.csi" ]] && return 0
    local tmp="$outdir/chunks/quilt2.$chr.chunk$chunk.tmp.vcf.gz"
    local work="$outdir/chunks/work.$chr.chunk$chunk"
    if Rscript "$QUILT2" \
        --prepared_reference_filename="$panel/$rdata" \
        ${bamlist:+--bamlist="$bamlist"} \
        --chr="$chr" \
        --regionStart="$start" \
        --regionEnd="$end" \
        --buffer="$buffer" \
        --nGen="$nGen" \
        --outputdir="$work" \
        --output_filename="$tmp" \
        $EXTRA_ARGS > "$outdir/logs/quilt2.$chr.chunk$chunk.log" 2>&1; then
        mv "$tmp" "$vcf"
        rm -f "$tmp.tbi" "$tmp.csi"
        bcftools index -f "$vcf"
        rm -rf "$work"
        log "done $chr chunk $chunk ($start-$end)"
    else
        log "FAILED $chr chunk $chunk, see $outdir/logs/quilt2.$chr.chunk$chunk.log"
        return 1
    fi
}
export -f impute_chunk log
export panel bamlist outdir QUILT2 EXTRA_ARGS

for chr in $chrs; do
    rows=$(awk -v chr="$chr" 'NR > 1 && $1 == chr' "$panel/chunks.tsv")
    [[ -n "$rows" ]] || { echo "ERROR: $chr not found in $panel/chunks.tsv" >&2; exit 1; }
    final="$outdir/quilt2.$chr.vcf.gz"
    if [[ -s "$final.csi" ]]; then
        log "$final exists, skipping $chr"
        continue
    fi
    log "imputing $chr ($(echo "$rows" | wc -l) chunks)"
    echo "$rows" | xargs -P "$threads" -L 1 bash -c 'impute_chunk "$@"' _

    vcfs=$(echo "$rows" | awk -v d="$outdir/chunks" '{print d"/quilt2."$1".chunk"$2".vcf.gz"}')
    if [[ $(echo "$vcfs" | wc -l) -eq 1 ]]; then
        cp $vcfs "$final"
    elif [[ "$EXTRA_ARGS" == *nipt* ]]; then
        # NIPT output has 3 haplotypes per sample which bcftools can not ligate,
        # so merge-sort the overlapping chunks and drop the duplicated sites
        bcftools concat --allow-overlaps --remove-duplicates -Oz -o "$final" $vcfs
    else
        bcftools concat --ligate -Oz -o "$final" $vcfs
    fi
    bcftools index -f "$final"
    [[ $keep -eq 1 ]] || rm -f $vcfs $(echo "$vcfs" | sed 's/$/.csi/')
    log "written $final"
done
