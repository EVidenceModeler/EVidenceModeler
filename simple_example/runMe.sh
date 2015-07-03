#!/bin/sh


## run EVM

../evidence_modeler.pl --genome genome.fasta \
                       --weights ./weights.txt \
                       --gene_predictions gene_predictions.gff3 \
                       --protein_alignments protein_alignments.gff3 \
                       --transcript_alignments transcript_alignments.gff3 \
                     > evm.out 

echo
echo
echo "*** Created EVM output file: evm.out ***"


## convert output to GFF3 format
./../EvmUtils/EVM_to_GFF3.pl evm.out.orig Contig1 > evm.out.gff3

echo
echo
echo "*** Converted EVM output to GFF3 format: evm.out.gff3 ***"

echo
echo "Done."




