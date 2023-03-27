#!/usr/bin/env python

# contributed by: Kresimir.Krizanovic@fer.hr

# Reformating the miniprot GFF output to match EVM specification
import sys

if len(sys.argv) != 2:
    sys.stderr.write("USAGE: miniprot_GFF_2_EVM_GFF3.py <miniprot GFF file>\n")
    sys.exit()

gff_filename = sys.argv[1]
evmfile = sys.stdout

with open(gff_filename, "r") as csvfile:
    # Read the first line, should be a comment
    line = csvfile.readline()
    i = 1

    # Structures for storing a data for one gene/transcript
    # miniprot currently contains lines for mrna, CDS and stop_codon records
    gene_line = []
    mrna_line = []
    cds_lines = []  # A list of CDS lines
    sc_line = []

    start = True
    while line:
        # This is just in case
        if line == "":
            break

        if line.startswith(
            "##PAF"
        ):  # transform all previously collected data and write it into the output file
            if (
                not start and len(mrna_line) > 0
            ):  # Just not at the start, because data has not been collected
                gene_line = mrna_line[:]  # copy the mrna data
                gene_line[2] = "gene"
                mrna_params = mrna_line[8]
                pos = mrna_params.find(";")
                mrna_id = mrna_params[3:pos]
                gene_id = "G_" + mrna_id
                gene_line[8] = "ID={0};Name={1} model {2}\n".format(
                    gene_id, "miniprot", mrna_id
                )
                mrna_line[8] = (
                    mrna_params[: pos + 1] + "Parent=" + gene_id + mrna_params[pos:]
                )

                evmfile.write("\t".join(gene_line))
                evmfile.write("\t".join(mrna_line))

                # CDS records do not have an ID, have to make one
                cnt = 0
                for cds_line in cds_lines:
                    cnt += 1
                    exon_line = cds_line[:]  # copy all elements
                    exon_line[2] = "exon"
                    exon_line[7] = "."  # set phase field for exons to '.'
                    exon_line[8] = (
                        "ID={0}.{1};".format(mrna_id, cnt) + exon_line[8]
                    )  # set exon id
                    cds_line[8] = (
                        "ID=CDS_{0}.{1};".format(mrna_id, cnt) + cds_line[8]
                    )  # set cds id

                    evmfile.write("\t".join(exon_line))
                    evmfile.write("\t".join(cds_line))

                if len(sc_line) > 0:
                    # Skip stop codon line
                    # evmfile.write('##' + '\t'.join(sc_line))
                    pass

                gene_line = []
                mrna_line = []
                cds_lines = []  # A list of CDS lines
                sc_line = []

        if line.startswith("#"):  # Comment line, just skip
            # evmfile.write(line)
            pass
        else:
            elements = line.split("\t")
            if elements[2].upper() == "MRNA":
                start = False  # We have collected some data - not start any more
                mrna_line = elements
            if elements[2].upper() == "CDS":
                cds_lines.append(elements)

            # Stop_codon usually marks the end of the transcript, but not always

            if elements[2].upper() == "STOP_CODON":
                sc_line = elements

        line = csvfile.readline()
        i += 1

    # Write the last collected set of data
    gene_line = mrna_line[:]  # copy the mrna data
    gene_line[2] = "gene"
    mrna_params = mrna_line[8]
    pos = mrna_params.find(";")
    mrna_id = mrna_params[3:pos]
    gene_id = "G_" + mrna_id
    gene_line[8] = "ID={0};Name={1} model {2}\n".format(gene_id, "miniprot", mrna_id)
    mrna_line[8] = mrna_params[: pos + 1] + "Parent=" + gene_id + mrna_params[pos:]

    evmfile.write("\t".join(gene_line))
    evmfile.write("\t".join(mrna_line))

    # CDS records do not have an ID, have to make one
    cnt = 0
    for cds_line in cds_lines:
        cnt += 1
        exon_line = cds_line[:]  # copy all elements
        exon_line[2] = "exon"
        exon_line[7] = "."  # set phase field for exons to '.'
        exon_line[8] = "ID={0}.{1};".format(mrna_id, cnt) + exon_line[8]  # set exon id
        cds_line[8] = "ID=CDS_{0}.{1};".format(mrna_id, cnt) + cds_line[8]  # set cds id

        evmfile.write("\t".join(exon_line))
        evmfile.write("\t".join(cds_line))

    if len(sc_line) > 0:
        # skip stop codon line
        # evmfile.write('##' + '\t'.join(sc_line))
        pass

sys.stderr.write("Done! Read {0} lines".format(i))
