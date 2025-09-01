#!/bin/bash

set -euo pipefail

wd=$(dirname $0)

seqtk subseq \
    /project/logsdon_shared/projects/Keith/rs-nucflag/core/test/ending_scaffold/input/aln_1.fa \
    /project/logsdon_shared/projects/Keith/rs-nucflag/core/test/ending_scaffold/input/aln_1.bed > $wd/out_1.fa
seqtk subseq \
    /project/logsdon_shared/projects/Keith/rs-nucflag/core/test/hsat/input/aln_1.fa \
    /project/logsdon_shared/projects/Keith/rs-nucflag/core/test/hsat/input/aln_1.bed > $wd/out_2.fa

# 2.30-r1287
minimap2 \
    -t 12 \
    -x asm20 \
    $wd/out_1.fa \
    $wd/out_2.fa \
    --eqx -c --cs > $wd/out.paf

seqtk subseq $wd/out_1.fa <(printf "HG00171_chrX_haplotype2-0000138:55480685-60218583\t4552385\t4552675") > $wd/out_1_subset.fa
seqtk subseq $wd/out_2.fa <(printf "NA18534_chr4_haplotype2-0000075:49478053-53940451\t4168003\t4168292") > $wd/out_2_subset.fa
