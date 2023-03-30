#!/usr/bin/env python


# contributed by: Kresimir.Krizanovic@fer.hr


# Reformating the miniprot GFF output to match EVM SPLICED ALIGNMENT specification
import sys

if len(sys.argv) != 2:
    sys.stderr.write(
        "USAGE: miniprot_GFF_2_EVM_alignment_GFF3.py <miniprot GFF file>\n"
    )
    sys.exit()

gff_filename = sys.argv[1]
evmfile = sys.stdout

with open(gff_filename, "r") as csvfile:
    # Read the first line
    line = csvfile.readline()
    i = 1

    # Structures for storing a data for one gene/transcript
    # miniprot currently contains lines for mrna, CDS and stop_codon records
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
                # no gene line for spliced alignments
                mrna_params = mrna_line[8]
                pos = mrna_params.find(";")
                mrna_id = mrna_params[3:pos]

                # CDS records do not have an ID, for spliced alignments, all CDS records have the same ID
                # Using mRNa ID
                cnt = 0
                for cds_line in cds_lines:
                    cnt += 1
                    # no exon lines for spliced alignments
                    cds_line[1] = "miniprot_protAln"
                    cds_line[2] = "nucleotide_to_protein_match"
                    # modifying CDS line params to set ID and setting score fiels as identity percentage
                    cds_elems = cds_line[8].split(";")
                    mrna_elems = mrna_line[8].split(";")
                    identity = float(cds_elems[2][9:])  # Extracting identity parameter
                    cds_line[5] = "{0:.2f}".format(
                        identity * 100
                    )  # Setting identitty percentage as score
                    cds_elems[0] = mrna_elems[0]  # Using MRNA ID for all CDS alignments
                    cds_line[8] = ";".join(cds_elems)

                    evmfile.write("\t".join(cds_line))

                mrna_line = []
                cds_lines = []  # A list of CDS lines
                sc_line = []

        if line.startswith("#"):  # Comment line, just write it in the output file
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
    mrna_params = mrna_line[8]
    pos = mrna_params.find(";")
    mrna_id = mrna_params[3:pos]

    # Using mRNa ID
    cnt = 0
    for cds_line in cds_lines:
        cnt += 1
        # no exon lines for spliced alignments
        cds_line[1] = "miniprot_protAln"
        cds_line[2] = "nucleotide_to_protein_match"
        # modifying CDS line params to set ID and setting score fiels as identity percentage
        cds_elems = cds_line[8].split(";")
        mrna_elems = mrna_line[8].split(";")
        identity = float(cds_elems[2][9:])  # Extracting identity parameter
        cds_line[5] = "{0:.2f}".format(
            identity * 100
        )  # Setting identitty percentage as score
        cds_elems[0] = mrna_elems[0]  # Using MRNA ID for all CDS alignments
        cds_line[8] = ";".join(cds_elems)

        evmfile.write("\t".join(cds_line))

sys.stderr.write("Done! Read {0} lines".format(i))
