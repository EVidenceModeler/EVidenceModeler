package Trellis;

use strict;
use warnings;
use Carp;

####
sub new {
    my $packagename = shift;
    
    my $self = { 
        ID_to_exon => {},
        exon_list => [],
        compatible_connections => {},
        best_connections => {},
    };
    
    bless ($self, $packagename);
    
    return ($self);
}

####
sub add_exon {
    my $self = shift;
    my ($lend, $rend, $type, $orient, $phase) = @_;
    
    ## reorder lend,rend just to be sure:
    ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
    
    my $id = "$lend,$rend,$type,$orient,$phase";
    
    if (my $exon_struct = $self->{ID_to_exon}->{$id}) {
        return ($id);
    }

    else {
        
        my $exon_struct = { 
            ID => $id,
            lend => $lend,
            rend => $rend,
            type => $type,
            orient => $orient,
            phase => $phase, };
        
        $self->{ID_to_exon}->{$id} = $exon_struct;

        push (@{$self->{exon_list}}, $exon_struct);
        
        return ($id);
    }
}

####
sub add_compatible_connection {
    my $self = shift;
    my ($id_A, $id_B) = @_;
    $self->{compatible_connections}->{$id_A}->{$id_B} = 1;
    return;
}

####
sub add_best_connection {
    my $self = shift;
    my ($id_A, $id_B) = @_;
    
    # first, delete an older one that might have been inserted:
    delete $self->{best_connections}->{$id_A};
    
    # establish current connection:
    $self->{best_connections}->{$id_A}->{$id_B} = 1;

    return;
}

####
sub get_exon_list {
    my $self = shift;
    return (@{$self->{exon_list}});
}

####
sub get_exon_by_ID {
    my $self = shift;
    my ($id) = @_;
    return ($self->{ID_to_exon}->{$id});
}

####
sub get_compatible_connections {
    my $self = shift;
    return (%{$self->{compatible_connections}});
}

####
sub get_best_connections {
    my $self = shift;
    return (%{$self->{best_connections}});
}


1; #EOM
