#!/usr/bin/env perl

package TIGR_XML_parser;

use strict;
use Gene_obj;
#use lib ("/usr/local/devel/linux_perlib");
use XML::Simple;
use Data::Dumper;
use Assembly::Evidence;

sub new {
    my $self = {gene_objs=>[],
		assembly_evidence => {}, #holds assembly evidence object
		assembly_seq => 0,
		clone_name=>0,
		asmbl_id=>0,
		chromosome=>0,
		clone_id=>0,
		gb_acc=>0,
		seq_group=>0
		};
    bless $self;
    return $self;
}


sub capture_genes_from_assembly_xml {
    my $self = shift;
    my $xml_file = shift;
    print STDERR "Parsing xml file: $xml_file\n";
    my $ref = XMLin($xml_file, forcearray=>1, keeproot=>1, forcecontent=>1, keyattr=>{});
    #print Dumper($ref);

    my $asmbl_ref;
    if (exists ($ref->{TIGR}->[0]->{ASSEMBLY})) {
	$asmbl_ref = $ref->{TIGR}->[0]->{ASSEMBLY};
    } elsif (exists ($ref->{TIGR}->[0]->{PSEUDOCHROMOSOME}->[0]->{ASSEMBLY})) {
	$asmbl_ref = $ref->{TIGR}->[0]->{PSEUDOCHROMOSOME}->[0]->{ASSEMBLY};
    } else {
	die "No assemblies to parse\n";
    }
    


    $self->{clone_name} = &content($asmbl_ref->[0]->{HEADER}->[0]->{CLONE_NAME});
    $self->{asmbl_id} = &content($asmbl_ref->[0]->{ASMBL_ID});
    $self->{chromosome} =  $asmbl_ref->[0]->{CHROMOSOME};
    $self->{clone_id} = $asmbl_ref->[0]->{CLONE_ID};
    $self->{gb_acc} = &content($asmbl_ref->[0]->{HEADER}->[0]->{GB_ACCESSION});
    $self->{seq_group} = &content($asmbl_ref->[0]->{HEADER}->[0]->{SEQ_GROUP});
 

    my $genelist_arrayref = $asmbl_ref->[0]->{GENE_LIST}->[0]->{PROTEIN_CODING}->[0]->{TU};
    my $x = 1;
    print STDERR "Creating Gene objects\n";
    foreach my $gene_ref_xml (@$genelist_arrayref) {
	my $gene_obj = &process_gene ($self, $gene_ref_xml);
	$gene_obj->refine_gene_object();
	$self->add_gene_obj($gene_obj);
	$x++;
    }
     
    # get assembly sequence.
    $self->{assembly_seq} = &content($asmbl_ref->[0]->{ASSEMBLY_SEQUENCE});
    $self->{assembly_seq} =~ s/\s+//g; #rid whitespace if any.
}

####
sub process_gene {
    my ($self, $gene_ref_xml) = @_;
    my $gene_obj = new Gene_obj();
    &process_gene_info ($gene_obj, $gene_ref_xml->{GENE_INFO}->[0]); #send hash ref
    &process_gene_models ($gene_obj, $gene_ref_xml->{MODEL}); #send array ref
    $gene_obj->{TU_feat_name} = &content($gene_ref_xml->{FEAT_NAME});
    &process_gene_evidence ($self, $gene_obj->{TU_feat_name}, $gene_ref_xml->{GENE_EVIDENCE}->[0]); #send hash ref.
    return ($gene_obj);
}

####
sub process_gene_info {
    my ($gene_obj, $gene_info_ref) = @_;
    $gene_obj->{locus} =  &content($gene_info_ref->{LOCUS});
    $gene_obj->{pub_locus} =  &content($gene_info_ref->{PUB_LOCUS});
    $gene_obj->{com_name} =  &content($gene_info_ref->{COM_NAME});
    $gene_obj->{pub_comment} = &content($gene_info_ref->{PUB_COMMENT});
    $gene_obj->{is_pseudogene} =  &content($gene_info_ref->{IS_PSEUDOGENE});
    $gene_obj->{ec_num} =  &content($gene_info_ref->{EC_NUM});
    $gene_obj->{gene_sym} =  &content($gene_info_ref->{gene_sym});
}
								   


####
sub content {
    my ($arrayref) = @_;
    return ($arrayref->[0]->{content});
}

####
sub process_gene_models {
    my ($gene_obj, $gene_models_ref) = @_;
    
    foreach my $gene_model_ref (@$gene_models_ref) {
	my $Model_feat_name = &content($gene_model_ref->{FEAT_NAME});
	$gene_obj->{Model_feat_name} = $Model_feat_name;
	my @exon_refs = @{$gene_model_ref->{EXON}};
	foreach my $exon_ref (@exon_refs) {
	    &get_exon_info($gene_obj, $exon_ref);
	}
	last; #lets not complicate this by including alternative splicing isoforms just yet.
    }

}


####
sub get_exon_info {
    my ($gene_obj, $exon_ref) = @_;
    my ($exon_end5, $exon_end3) = &process_coordset ($exon_ref);
    my $exon_obj = new mRNA_exon_obj($exon_end5, $exon_end3);
    if (exists ($exon_ref->{CDS})) {
	my ($cds_end5, $cds_end3) = &process_coordset ($exon_ref->{CDS}->[0]);
	$exon_obj->add_CDS_exon_obj($cds_end5, $cds_end3);
    }
    $gene_obj->add_mRNA_exon_obj($exon_obj);
}


####
sub process_coordset {
    my ($coords_ref) = @_;
    my $end5 = &content($coords_ref->{COORDSET}->[0]->{END5});
    my $end3 = &content($coords_ref->{COORDSET}->[0]->{END3});
    return ($end5, $end3);
}




sub add_gene_obj {
    my $self = shift;
    my $gene_obj = shift;
    my $index = $#{$self->{gene_objs}};
    $index++;
    $self->{gene_objs}->[$index] = $gene_obj;
}

sub get_genes {
    my $self = shift;
    return (@{$self->{gene_objs}});
}



sub get_assembly_sequence {
    my $self = shift;
    return ($self->{assembly_seq});
}


# if no CDSs are specified, want to create CDSs based on mRNA exons.
sub force_CDS_occupancy {
    my $self = shift;
    my @genes = $self->get_genes();
    foreach my $gene (@genes) {
	my @exons = $gene->get_exons();
	foreach my $exon (@exons) {
	    my ($end5, $end3) = $exon->get_mRNA_exon_end5_end3();
	    $exon->add_CDS_exon_obj($end5, $end3);
	}
	$gene->refine_gene_object();
    }
}




sub toString {
    my $self = shift;
    my $text = "Clone_name: $self->{clone_name}\n"
	. "Clone_id: $self->{clone_id}\n"
	    . "Chromosome: $self->{chromosome}\n"
	. "Seq_group: $self->{seq_group}\n"
	. "Genbank Accession: $self->{gb_acc}\n"
	    . "TIGR Asmbl_id: $self->{asmbl_id}\n";
    $text .= "\n\nGenes:\n";
    my @genes = $self->get_genes();
    my $x = 1;
    foreach my $gene (@genes) {
	$text .= "gene: $x\n" . $gene->toString();
	$x++;
    }
    $text .= "\n\nAssembly Sequence:\n" . $self->get_assembly_sequence();
    return ($text);
}

sub process_gene_evidence {
    my ($self, $TU_feat_name, $xml_ref) = @_;
   
    if (ref $xml_ref) {
	my $evidence_type = $xml_ref->{EVIDENCE_TYPE}->[0];
	## get search data.
	if (exists $evidence_type->{COMPUT_PREDICTION}) {
	    my $prediction_data = $evidence_type->{COMPUT_PREDICTION}->[0];
	    &process_gene_prediction_evidence ($self, $TU_feat_name, $prediction_data);
	}
	## get gene predictions.
	if (exists $evidence_type->{SEQUENCE_DB_MATCH}) {
	    my $search_data = $evidence_type->{SEQUENCE_DB_MATCH}->[0];
	    &process_search_evidence($self, $TU_feat_name, $search_data);
	}
    }
}

sub process_search_evidence {
    my ($self, $TU_feat_name, $search_data) = @_;
    #print Dumper ($search_data);
    foreach my $db_ref (@{$search_data->{SEARCH_DB}}) {
	#print Dumper ($db_ref);
	#exit;
	my $db_type = $db_ref->{DB_TYPE};
	my $db_name = $db_ref->{DB_NAME};
	my $seq_elements = $db_ref->{SEQ_ELEMENT};
	foreach my $seq_element (@$seq_elements) {
	    #print Dumper ($seq_element);
	    my $method = $seq_element->{METHOD};
	    my $score = $seq_element->{SCORE};
	    my $accession = $seq_element->{ACCESSION};
	    my $description = $seq_element->{DESCRIPTION};
	    my ($asmbl_end5, $asmbl_end3) = &process_coordset($seq_element->{ASMBL_COORDS}->[0]);
	    my ($feat_end5, $feat_end3) = &process_coordset ($seq_element->{MATCH_COORDS}->[0]);
	    print "$db_type\t$method\t$score\t$asmbl_end5\t$asmbl_end3\t$feat_end5\t$feat_end3\n";
	    &add_evidence_element($self, $db_name, $db_type, $method, $accession, $TU_feat_name, $score, $asmbl_end5, $asmbl_end3, $feat_end5, $feat_end3, $description);
	}
	#exit;
    }


#    exit;
}


sub process_gene_prediction_evidence {
    my ($self, $TU_feat_name, $prediction_data) = @_;
    #print Dumper ($search_data);
    foreach my $db_ref (@{$prediction_data->{PREDICTION_SET}}) {
	#print Dumper ($db_ref);
	#exit;
	my $db_type = $db_ref->{PREDICTION_TYPE};
	my $db_name = $db_ref->{PREDICTION_TOOL};
	my $seq_elements = $db_ref->{SEQ_ELEMENT};
	foreach my $seq_element (@$seq_elements) {
	    #print Dumper ($seq_element);
	    my $accession = $TU_feat_name . "_${db_name}_" . $seq_element->{ELEMENT_NUM};
	    my ($asmbl_end5, $asmbl_end3) = &process_coordset($seq_element->{ASMBL_COORDS}->[0]);
	    my ($feat_end5, $feat_end3) = &process_coordset ($seq_element->{MATCH_COORDS}->[0]);
	    print "PREDICTION: $accession\t$asmbl_end5\t$asmbl_end3\t$feat_end5\t$feat_end3\n";
	    &add_evidence_element($self, $db_name, $db_type, "", $accession, $TU_feat_name, "", $asmbl_end5, $asmbl_end3, $feat_end5, $feat_end3, "");
	}
	#exit;
    }
}


sub add_evidence_element {
    my $self = shift;
    my ($db_name, $db_type, $method, $accession, $TU_feat_name, $score, $asmbl_end5, $asmbl_end3, $feat_end5, $feat_end3, $description) = @_;
    my $asmbl_ev_key = $db_type . "-" . $method;
    my $assembly_evidence = $self->{assembly_evidence};
    if (!exists ($assembly_evidence->{$asmbl_ev_key})) {
	$assembly_evidence->{$asmbl_ev_key} = new Assembly::Evidence();
    }
    my $evidence_container = $assembly_evidence->{$asmbl_ev_key};
    $evidence_container->add_evidence($db_type, $db_name, $accession, $TU_feat_name, $feat_end5, $feat_end3, $asmbl_end5, $asmbl_end3, $description);
}




1;# end of TIGR_db_parser



