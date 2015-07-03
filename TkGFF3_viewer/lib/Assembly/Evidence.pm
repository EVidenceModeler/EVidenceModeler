#!/usr/bin/env perl

package main;
use strict;
my $UNDEF = undef();
our $SEE;


package Assembly::Evidence;
use strict;
use Data::Dumper;

sub new {
    my $self = {
        elements => [],    ## array of all evidence elements
        evidence_table => {}  ## table for evidence keyed by chain ID.
    };
    bless $self;
    return ($self);
}

sub add_evidence {
    my $self = shift;
    my ($element_type, $db_name, $accession, $match_chain_identifier, 
        $feat_end5, $feat_end3, $asmbl_end5, $asmbl_end3, $description) = @_;
    
    print "Adding evidence: " . join ("\n", @_) . "\n------\n\n" if $SEE;

    my $element = $self->{evidence_table}->{$match_chain_identifier};
    
    my $segment = "";
    ## instantiate new one if doesn't already exist
    if ($element) {
        $segment = Assembly::Evidence::EvidenceElement->new();
    }
    else {
        
        $element = Assembly::Evidence::EvidenceElement->new(); #container for segments.
        $element->{segments} = [];
        $self->{evidence_table}->{$match_chain_identifier} = $element;
        push (@{$self->{elements}}, $element);
      
        $segment = Assembly::Evidence::EvidenceElement->new();
    }
    
    push (@{$element->{segments}}, $segment); # add segment to list
       

    $segment->{element_type} = $element_type;
    $segment->set_type_and_identifier($element_type, $accession, $db_name, $description);
    $segment->set_asmbl_end5_end3($asmbl_end5, $asmbl_end3);
    $segment->set_feat_end5_end3($feat_end5, $feat_end3);
    
    
    return;
}

####
sub refine_evidence_coords {
    my $self = shift;
    foreach my $evidence_ele (@{$self->{elements}}) {
        $evidence_ele->refine_element();
    }

    return;
}



sub get_all_evidence {
    my $self = shift;
    return (@{$self->{elements}});
}

sub get_evidence_via_gene_association {
    my $self = shift;
    my $gene_association = shift;
    my $evidence_table = $self->{evidence_table};
    my $keyed_evidence = $evidence_table->{$gene_association};
    return (@{$keyed_evidence});
}
#private methods to this class
sub within_distance {
    my ($distance, $coord1, $coord2) = @_;
    if (abs ($coord2 - $coord1) <= $distance) {
        return (1); #true
    } else {
        return (0); #false
    }
}

sub same_orientation {
    my ($dir1, $dir2) = @_;
    if ($dir1 eq $dir2) {
	return (1);
    } else {
	return (0);
    }
}

sub same_type_and_identifier {
    my ($element1, $element2) = @_;
    my $type1 = $element1->{type};
    my $identifier1 = $element1->{identifier};
    my $type2 = $element2->{type};
    my $identifier2 = $element2->{identifier};
    if ($type1 eq $type2 && $identifier1 eq $identifier2) {
        return (1);
    } else {
        return (0);
    }
}


###############################################################################################
## Private class EvidenceElement.  Should access only indirectly via Assembly::Evidence object.
package Assembly::Evidence::EvidenceElement;
use strict;
sub new {
    my $packagename = shift;

    my $self = {
	element_type=>0, #descriptive term. Gene prediction program name or search database.
	identifier=>0, #accession or model number
	db_name =>0, #name of the search database.
	description => $UNDEF, #header info for search db match
	feat_end5 => $UNDEF,
	feat_end3 => $UNDEF,
	asmbl_end5 => $UNDEF,
	asmbl_end3 => $UNDEF,
	score => $UNDEF,
	per_id => $UNDEF,
	strand => 0,
	classification => 0, #descriptive term built based on combo of element_type and db_name.
	gene_association => 0,  #holds feat_name of Transcriptional Unit for gene association.
	segments => 0  #an element can be composed of several segments. Contains array ref if assigned.
	};
    bless ($self, $packagename);
    return ($self);
}



sub refine_element {
    my $self = shift;
    print "Refining evidence element\n " if $SEE;
    my @segments = $self->get_segments();
    my (@asmbl_end5s, @asmbl_end3s, $strand, $type, $identifier, $db_name, $description);
    foreach my $segment (@segments) {
        my ($end5, $end3) = $segment->get_coords();
        push (@asmbl_end5s, $end5);
        push (@asmbl_end3s, $end3);
        $strand = $segment->{strand};
        $type = $segment->{element_type};
        $identifier = $segment->{identifier};
        $db_name = $segment->{db_name};
        $description = $segment->{description};
    }
    @asmbl_end5s = sort {$a<=>$b} @asmbl_end5s;
    @asmbl_end3s = sort {$a<=>$b} @asmbl_end3s;
    my ($set_end5, $set_end3);
    if ($strand eq "+") {
        $set_end5 = shift (@asmbl_end5s);
        $set_end3 = pop (@asmbl_end3s);
    } else {
        $set_end5 = pop (@asmbl_end5s);
        $set_end3 = shift (@asmbl_end3s);
    }
    print "\t$type\t$identifier\t$set_end5\t$set_end3\n" if $SEE;
    $self->set_type_and_identifier($type, $identifier, $db_name, $description);
    $self->set_asmbl_end5_end3($set_end5, $set_end3);
    
}


## element set methods.

sub set_type_and_identifier {
    my $self = shift;
    my ($element_type, $identifier, $db_name, $description) = @_;
    $self->{element_type} = $element_type;
    $self->{identifier} = $identifier;
    $self->{db_name} = $db_name;
    $self->{classification} = $element_type . "_" . $db_name;
    $self->{description} = $description;
}

sub set_asmbl_end5_end3 {
    my $self = shift;
    my ($end5, $end3) = @_;
    $self->{asmbl_end5} = $end5;
    $self->{asmbl_end3} = $end3;
    if ($end5 < $end3) {
        $self->{strand} = '+';
    } elsif ($end5 > $end3) {
        $self->{strand} = '-';
    }
}

sub set_feat_end5_end3 {
    my $self = shift;
    my ($end5, $end3) = @_;
    $self->{feat_end5} = $end5;
    $self->{feat_end3} = $end3;
}

sub set_gene_association {
    my $self = shift;
    my $gene_association = shift;
    $self->{gene_association} = $gene_association;
}


## element get methods

sub get_asmbl_end5_end3 {
    my $self = shift;
    return ($self->{asmbl_end5}, $self->{asmbl_end3});
}

sub get_seq_span {
    my $self = shift;
    return ($self->get_asmbl_end5_end3());
}


sub get_feat_end5_end3 {
    my $self = shift;
    return ($self->{feat_end5}, $self->{feat_end3});
}

sub get_coords {
    my $self = shift;
    return ($self->get_asmbl_end5_end3());
}

sub get_segments {
    #not being used.  Currently, an evidence element is a single segment
    my $self = shift;
    if (ref $self->{segments} eq "ARRAY") {
        return (@{$self->{segments}});
    } else {
        return ($self);
    }
}

sub get_annot_text {
    my $self = shift;
    my $element_type = $self->{element_type};
    my $identifier = $self->{identifier};
    my $asmbl_end5 = $self->{asmbl_end5};
    my $asmbl_end3 = $self->{asmbl_end3};
    my $description = $self->{description};
    return ("$element_type $identifier ($asmbl_end5-$asmbl_end3) $description");
}



1;







