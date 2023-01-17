#!/bin/sh

set -ex

## run EVM

singularity exec ../Docker/EVidenceModeler.latest.simg EVidenceModeler --sample_id smalltest.singularity \
                   --genome genome.fasta \
                   --weights ./weights.txt \
                   --gene_predictions gene_predictions.gff3 \
                   --protein_alignments protein_alignments.gff3 \
                   --transcript_alignments transcript_alignments.gff3 \
                   --segmentSize 100000 \
                   --overlapSize 10000

      
echo "Done. See smalltest.singularity.EVM.* outputs"

