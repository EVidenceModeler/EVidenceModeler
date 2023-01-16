#!/bin/sh

set -ex

## run EVM

docker run --rm  -v "$(pwd)":/data brianjohnhaas/evidencemodeler:latest \
       bash -c 'cd /data && EVidenceModeler --sample_id smalltest.docker \
                   --genome genome.fasta \
                   --weights ./weights.txt \
                   --gene_predictions gene_predictions.gff3 \
                   --protein_alignments protein_alignments.gff3 \
                   --transcript_alignments transcript_alignments.gff3 \
                   --segmentSize 100000 \
                   --overlapSize 10000 '

      
echo "Done. See smalltest.EVM.* outputs"

