package EvidenceFeature;

## an evidence feature is a set of one or more genome segments
## that can be ascribed a single identifier.
## An EvidenceFeature can have several child features to become
## part of a larger evidence feature graph.

use strict;
use warnings;
use Carp;

####
sub new {
    my $packagename = shift;
    my ($feature_ID, $feature_type, $feature_source) = @_;

    unless ($feature_ID && $feature_type && $feature_source) {
        confess "error, need feature_type and feature_source in constructor";
    }


    my $self = { 
        children => undef, #can be made an array ref
        feature_type => $feature_type,
        feature_source => $feature_source,
        segments => [], # holds [end5, end3] values
        orient => undef, # +|-
        ID => $feature_ID,
        parent_ID => undef, #if children are added, children have this set to the parent identifier
        
        ## coordinates maintained as min/max of segments
        lend => undef,
        rend => undef,
        span_length => undef,

    };
    
    bless ($self, $packagename);
    
    return ($self);
}

####
sub add_segment {
    my $self = shift;
    my ($lend, $rend, $orient) = @_;
    
    ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
    
    unless ($self->{orient}) {
        $self->{orient} = $orient if $orient;
    }
    
    push (@{$self->{segments}}, [$lend, $rend]);

    $self->_update_coordinate_span();
    return;
}


####
sub get_segments {
    my $self = shift;
    return (@{$self->{segments}});
}

####
sub get_coords {
    my $self = shift;
    return ($self->{lend}, $self->{rend});
}

####
sub add_child {
    my $self = shift;
    my ($child_feature) = @_;
    $child_feature->{parent_ID} = $self->{ID};
    unless (ref $self->{children}) {
        $self->{children} = []; #make an array ref
    }
    # add the child
    push (@{$self->{children}}, $child_feature);
    return;
}


####
sub overlaps_coords {
    my $self = shift;
    my ($lend, $rend) = @_;
    
    ($lend, $rend) = sort {$a<=>$b} ($lend, $rend); #just to be sure
    
    if ($self->{lend} <= $rend && $self->{rend} >= $lend) {
        ## yes overlaps
        return (1);
    }
    else {
        return (0); 
    }
}

####
sub toString {
    my $self = shift;
    
    my $text = "ID: " . $self->{ID} 
    . " type: " . $self->{feature_type} 
    . " source: " . $self->{feature_source}
    . " orient(" . $self->{orient} . ")";

    if (my $parent = $self->{parent_ID}) {
        $text .= " parent: $parent";
    }
    foreach my $segment ($self->get_segments()) {
        my ($lend, $rend) = @$segment;
        $text .= " $lend-$rend";
    }
   
    $text .= "\n";
    
    return ($text);
}


####
sub retrieve_children {
    my $self = shift;
    if (ref $self->{children}) {
        return (@{$self->{children}});
    }
    else {
        return ();
    }
}


####
sub retrieve_descendants {
    my $self = shift;
    my ($feature_type) = @_;
    
    my @descendants_wanted;
    foreach my $child ($self->retrieve_children()) {
        if ($child->{feature_type} eq $feature_type) {
            push (@descendants_wanted, $child);
        }
        else {
            my @search_child = $child->retrieve_descendants($feature_type);
            if (@search_child) {
                push (@descendants_wanted, @search_child);
            }
        }
    }
    return (@descendants_wanted);
}
        

#### Private methods:

####
sub _update_coordinate_span {
    my $self = shift;
    my @segments = $self->get_segments();
    
    my @coords;
    foreach my $segment (@segments) {
        push (@coords, @$segment);
    }
    @coords = sort {$a<=>$b} @coords;

    my $lend = shift @coords;
    my $rend = pop @coords;
    my $length = $rend - $lend + 1;

    $self->{lend} = $lend;
    $self->{rend} = $rend;
    $self->{span_length} = $length;

    return;
}





1; #EOM

        
