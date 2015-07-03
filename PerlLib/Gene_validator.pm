#!/usr/local/bin/perl

package Gene_validator;

## There's no reason to create an instance of this class.  Class method will be used to validate
## a gene model.

## Requirements for annot validator:
#   -gene components must have same orientation as the TU.
#   -gene component coordinates must behave as expected; ie. CDS must reside within the exon.
#   -protein translations must exist for non-pseudogenes and must be complete ORFs.
#   -consensus splice sites must exist at internal exon boundaries.

require Exporter;
use Gene_obj;
use strict;
use DBI;

use vars qw ($SEE $DEBUG @ISA @EXPORT); ## set in script using this module for verbose output.
@ISA = qw(Exporter);
@EXPORT = qw($SEE $DEBUG validate_gene);


## validate gene 
sub validate_gene {
    my ($gene_obj, $asmbl_seq_reference) = @_;
    
    $gene_obj->create_all_sequence_types($asmbl_seq_reference);
    my ($error_text); #hold errors for gene analysis
    my $err_ref = \$error_text; #pass this ref to subs.
    &validate_gene_structure ($gene_obj, $err_ref, $asmbl_seq_reference); ## do all the serious work here.
    
    foreach my $isoform ($gene_obj->get_additional_isoforms()) {
	my $errors = &validate_gene($isoform, $asmbl_seq_reference);
	if ($errors) {
	    $error_text .= $errors;
	}
    }
    
    if ($error_text) {
	return (" errors for $gene_obj->{Model_feat_name} . $error_text\n");
    } else {
	return (undef());
    }
}



################################################################################################################


####
sub validate_gene_structure {
    my ($gene_obj, $err_ref, $asmbl_seq_ref) = @_;
    ## here, we analyse the gene structure, splice site data, and ORF composition.
   
    my @exons = $gene_obj->get_exons();
    unless (@exons) { #no gene structure here
	
	&add_error ($err_ref, "\tNo Gene Structure available for gene\n");
	return;
    }
    
    ## Traverse data structure, identify problems.
    my ($TU_end5, $TU_end3) = $gene_obj->get_coords();
    my $gene_orient = &get_orient ($TU_end5, $TU_end3);
    my $is_pseudogene = $gene_obj->{is_pseudogene};
    my $TU_feat_name = $gene_obj->{TU_feat_name};
    print "TU: $TU_feat_name\t$gene_orient\t$TU_end5\t$TU_end3\n" if $SEE;
    my $model_feat_name = $gene_obj->{Model_feat_name};
    my ($model_end5, $model_end3) = $gene_obj->get_model_span();
    my $model_orient = &confirm_orientation ($gene_orient, $model_end5, $model_end3, "model", $err_ref);
    &confirm_encapsulation ($TU_end5, $TU_end3, $model_end5, $model_end3, "TU/model", $err_ref);
    print "Model: $model_feat_name\t$model_orient\t$model_end5\t$model_end3\n" if $SEE;
    my $num_exons = $#exons;
    print "Model has " . ($num_exons+1) . " number of exons.\n" if $SEE;
    for (my $i = 0; $i <= $num_exons; $i++) {
	my $exon = $exons[$i];
	my $exon_feat_name = $exon->{feat_name};
	my ($exon_end5, $exon_end3) = $exon->get_coords();
	my $exon_orient = &confirm_orientation($gene_orient, $exon_end5, $exon_end3, "exon", $err_ref);
	print "Exon: $exon_feat_name\t$exon_orient\t$exon_end5\t$exon_end3\n" if $SEE;
	&confirm_encapsulation($TU_end5, $TU_end3, $exon_end5, $exon_end3, "TU/exon", $err_ref);
	my $exon_type = &get_exon_type ($i, $num_exons);
	&confirm_splice_sites ($exon_end5, $exon_end3, $gene_orient, $exon_type, 
			       $exon_feat_name, $asmbl_seq_ref, $err_ref) if (!$is_pseudogene); 
	my $cds = $exon->get_CDS_obj();
	if ($cds) {
	    my $CDS_feat_name =  $cds->{feat_name};
	    my ($CDS_end5, $CDS_end3) = $cds->get_coords();
	    my $CDS_orient = &confirm_orientation($gene_orient, $CDS_end5, $CDS_end3, "CDS", $err_ref);
	    &confirm_encapsulation($exon_end5, $exon_end3, $CDS_end5, $CDS_end3, "exon:$exon_feat_name/CDS:$CDS_feat_name", $err_ref);
	    print "CDS: $CDS_feat_name\t$CDS_orient\t$CDS_end5\t$CDS_end3\n\n" if $SEE;
	}
    }
    &examine_protein_sequence ($gene_obj, $err_ref) if (!$is_pseudogene);
    &examine_CDS_sequence ($gene_obj, $err_ref, $is_pseudogene);
}



####
sub confirm_splice_sites {
    my ($exon_end5, $exon_end3, $exon_orient, $exon_type, $exon_feat_name, $asmbl_seq_ref, $err_ref) = @_;
    print "exon_type: $exon_type\n" if $SEE;
    if ($exon_type eq "gene") { return;} ## no splicing here.
    my ($coord1, $coord2) = sort {$a<=>$b} ($exon_end5, $exon_end3);
    ## get two coordinate sets corresponding to potential splice sites
    my $splice_1_start = $coord1-2-1;
    my $splice_2_start = $coord2-1+1;
    print "confirming splice sites at "  . ($splice_1_start +1) . " and " . ($splice_2_start + 1) . "\n"if $SEE;
    my $splice_1 = substr ($$asmbl_seq_ref, $splice_1_start, 2);
    my $splice_2 = substr ($$asmbl_seq_ref, $splice_2_start, 2);
    my ($acceptor, $donor) = ($exon_orient eq '+') ? ($splice_1, $splice_2) : (&reverse_complement($splice_2), &reverse_complement($splice_1)); 
    my $check_acceptor = ($acceptor =~ /ag/i);
    my $check_donor = ($donor =~ /gt|gc/i);
    ## associate results of checks with exon type.
    if ($exon_type eq "initial" || $exon_type eq "internal") {
	unless ($check_donor) {
	    &add_error ($err_ref, "\tnon-consensus $donor donor splice site at $exon_end3\n");
	}
    }
    if ($exon_type eq "internal" || $exon_type eq "terminal") {
	unless ($check_acceptor) {
	    &add_error ($err_ref, "\tnon-consensus $acceptor acceptor splice site at $exon_end5\n");
	}
    }
    return;
}



####
sub get_exon_type {
    my ($curr_exon_num, $num_exons) = @_;

    ## count starts at 0; array based.

    ## types defined
    # gene: single exon gene
    # initial :first exon of a multi-exon gene
    # internal: internal exon of a multi-exon gene
    # terminal: last exon of a multi-exon gene


    if ($curr_exon_num == 0) {
	if ($num_exons == 0) {
	    return "gene";
	} else {
	    return "initial";
	}
    } elsif ($curr_exon_num == $num_exons) {
	return "terminal";
    } else {
	return "internal";
    }
}



####
sub confirm_orientation {
    my ($gene_orient, $coord1, $coord2, $type, $err_ref) = @_;
    my $orient = &get_orient($coord1, $coord2);
    if (($orient) && ($gene_orient ne $orient)) {
	&add_error($err_ref, "\t$type has opposite orientation to gene ($gene_orient)\n");
	return ($gene_orient);
    }
    return ($orient);
}


####
sub confirm_encapsulation {
    my ($coord_big1, $coord_big2, $coord_sm1, $coord_sm2, $type, $err_ref) = @_;
    ($coord_big1, $coord_big2) = sort {$a<=>$b} ($coord_big1, $coord_big2);
    ($coord_sm1, $coord_sm2) = sort {$a<=>$b} ($coord_sm1, $coord_sm2);
    unless ( ($coord_sm1 >= $coord_big1) && ($coord_sm2 <= $coord_big2) ) {
	print "Not Encapsulated: ($coord_sm1, $coord_sm2, $type) within ($coord_big1, $coord_big2)\n" if $SEE;
	&add_error($err_ref, "\t$type, coords not encapsulated\n");
    }
}

####
sub get_orient {
    my ($coord1, $coord2) = @_;
    my $orient;
    if ($coord1 < $coord2) {
	$orient = '+';
    } elsif ($coord1 > $coord2) {
	$orient = '-';
    }
    return ($orient);
}

####
sub get_pseudogene_status {
    my ($dbproc, $TU_feat_name, $err_ref) = @_;
    my $query = "select is_pseudogene from ident where feat_name = \"$TU_feat_name\"\n";
    my $status = &first_result_sql ($dbproc, $query);
    if ($status == 1) {
	return (1);
    } else {
	return(0);
    }
}

####
sub examine_protein_sequence {
    my ($gene_obj, $err_ref) = @_;
    my $prot_seq = $gene_obj->get_protein_sequence();
    unless ($prot_seq) {
	&add_error ($err_ref, "\tNo protein sequence.\n");
	return;
    }

    my $num_stops = &get_number_stops ($prot_seq);
    if ($num_stops > 1) {
	&add_error ($err_ref,  "\tcorrupt protein sequence: [$num_stops] stops in protein sequence\n");
	return;
    }
    my $first_char = substr ($prot_seq, 0, 1);
    unless ($first_char =~ /M/i) {
	&add_error ($err_ref,  "\tProtein seq doesn't start with Methionine.....[$first_char] instead.\n");
    }
    my $prot_length = length ($prot_seq);
    my $last_char = substr ($prot_seq, ($prot_length - 1), 1);
    unless ($last_char =~ /\*/) {
	&add_error ($err_ref,  "\tProtein seq doesn't terminate with a stop codon....[$last_char] instead.\n");
    }
}



####
sub get_number_stops {
    my ($prot_seq) = @_;
    my $stop_num = 0;
    while ($prot_seq =~ /\*/g) {
	$stop_num++;
    } 
    return ($stop_num);
}

####
sub add_error {
    my ($err_ref, $comment) = @_;
    $comment = "\tERROR:\t$comment";
    $$err_ref .= $comment;
}

####
sub examine_CDS_sequence {
    my ($gene_obj, $err_ref, $is_pseudogene) = @_;
    my $sequence = $gene_obj->get_CDS_sequence();
    if (!$sequence) {
	&add_error ($err_ref, "No CDS sequence could be found for this gene. Possibly no CDS-exons?\n");
    } elsif (!$is_pseudogene) {
	chomp $sequence; #shouldn't need to do this, but just to be safe.
	my $cds_length = length ($sequence);
	my $start_codon = substr ($sequence, 0, 3);
	my $stop_codon = substr ($sequence, $cds_length - 3, 3);
	if ($start_codon !~ /atg/i) {
	    &add_error ($err_ref, "CDS sequence does not begin with ATG; instead [$start_codon]\n");
	}
	if ($stop_codon !~ /tag|tga|taa/i) {
	    &add_error ($err_ref, "CDS sequence does not include a STOP codon; instead [$stop_codon]\n");
	}
    }
}


sub reverse_complement { #alternative to revComp, this one makes more sense in the context of the added subroutines
     my($s) = @_;
     my ($rc);
     $rc = reverse ($s);
     $rc =~tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/;
     return($rc);
 }










1; #end of Annot_validator.pm

