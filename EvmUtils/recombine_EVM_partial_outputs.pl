#!/usr/bin/env perl

## script uses dynamic programming to find the set of non-overlapping complete predictions, maximizing for gene lengths

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;

my ($SEE, $help, $partitions_file, $output_file_name);

my $usage = "usage: $0 --partitions \$partitions_file --output_file_name \$output_file_name\n\n";

&GetOptions ("partitions=s" => \$partitions_file,
             "output_file_name|O=s" => \$output_file_name,
             "verbose|S" => \$SEE,
             "help|h" => \$help,
             );


if ($help || ! ($partitions_file && $output_file_name) ) {
    die $usage;
}

$output_file_name = basename($output_file_name); #just to be sure...


my %base_directories_to_partitions;
open (my $fh, $partitions_file) or die "Error, cannot open $partitions_file";
while (<$fh>) {
    chomp;
    my ($accession, $base_dir, $partitioned, $partition_dir) = split (/\t/);
    if ($partitioned eq 'Y') {
        $base_directories_to_partitions{$base_dir}->{$partition_dir} = 1;
    }
}
close $fh;

foreach my $base_dir (keys %base_directories_to_partitions) {
    my $partition_dirs_href = $base_directories_to_partitions{$base_dir};

    my @predictions;
    foreach my $partition_dir (keys %$partition_dirs_href) {
        $partition_dir =~ /(\d+)-(\d+)$/ or die "Error, cannot extract coords from $partition_dir";
        my ($range_lend, $range_rend) = ($1, $2);
        
        my $output_file = $partition_dir . "/$output_file_name";
        &parse_and_add_predictions($output_file, $range_lend, \@predictions);
    }
    
    open (my $out_fh, ">$base_dir/$output_file_name") or die $!;
    print STDERR "\twriting output $base_dir/$output_file_name\n\n";
    
    my @final_predictions = &combine_predictions(\@predictions);
    
    ## write predictions:
    foreach my $prediction (@final_predictions) {
        my $pred_text = $prediction->{text};
        print $out_fh "$pred_text\n";

        if (my $aref = $prediction->{intronic_preds}) {
            foreach my $intronic_pred (@$aref) {
                print $out_fh "!! Intron-containing prediction\n";
                print $out_fh $intronic_pred->{text} . "\n";
            }
        }
    }
    
    close $out_fh;
}

exit(0);


####
sub parse_and_add_predictions {
    my ($output_file, $partition_lend, $predictions_aref) = @_;
    
    print STDERR "Parsing $output_file\n";
    
    my $current_prediction_text = "";
    
    my @preds;
    
    open (my $fh, $output_file) or die $!;
    while (<$fh>) {
        if (/^\!/) { # evm comment line. '^!!'
            next; 
        }
        
        if (/^\#/ && ! /EVM/) { next; }
        
        if (/^(\d+|\#)/) {
            $current_prediction_text .= $_;
        } else {
            if ($current_prediction_text) {
                &process_prediction($current_prediction_text, $partition_lend, \@preds);
                $current_prediction_text = "";
            }
        }
    }
    # get last one
    if ($current_prediction_text) {
        &process_prediction($current_prediction_text, $partition_lend, \@preds);
    }
    close $fh;
    
    
    #print Dumper (\@preds);
    


    ## join preds within introns of other preds
    my @final_preds = &join_intronic_preds(\@preds);

    push (@$predictions_aref, @final_preds);
    

}


####
sub join_intronic_preds {
    ## do N^2 comparison among preds.  If a pred is found encapsulated by another pred, it must be because it's within an intron.
    ## not doing any checks now to verify this...

    my ($preds_aref) = @_;
    
    my @preds = sort {$a->{lend}<=>$b->{lend}} @$preds_aref;

    ## add encaps flag to struct if joined to another one.

    for (my $i = 0; $i < $#preds; $i++) {
        
        # since i always comes before j, j can be within an intron of i but not vice-versa

        my $pred_i = $preds[$i];
        my ($pred_i_lend, $pred_i_rend) = ($pred_i->{lend}, $pred_i->{rend});
        
        for (my $j = $i+1; $j <= $#preds; $j++) {
            
            my $pred_j = $preds[$j];
            my ($pred_j_lend, $pred_j_rend) = ($pred_j->{lend}, $pred_j->{rend});
            
            if ($pred_j_lend > $pred_i_lend && $pred_j_rend < $pred_i_rend) {
                ## j encaps in i
                
                $pred_j->{encaps} = 1;
                
                #print "Encaps: $pred_j_lend, $pred_j_rend   to $pred_i_lend $pred_j_rend\n";
                
                my $encaps_list_aref = $pred_i->{intronic_preds};
                unless (ref $encaps_list_aref) {
                    $encaps_list_aref = $pred_i->{intronic_preds} = [];
                }
                push (@$encaps_list_aref, $pred_j);
                ## add its length to i's base score
                $pred_i->{length} += $pred_j->{length};
                $pred_i->{path_score} += $pred_j->{path_score};
                
            }
        }
    }
    
    my @final_preds;
    foreach my $pred (@preds) {
        unless ($pred->{encaps}) {
            push (@final_preds, $pred);
        }
    }


    return (@final_preds);
}



####
sub process_prediction {
    my ($prediction_text, $partition_lend, $predictions_aref) = @_;
    
    my $new_prediction_text = "";
    my @exon_types;
    my @coords;
    
    my @prediction_lines = split (/\n/, $prediction_text);
    my $prediction_header = shift @prediction_lines;
    
    my @header_comps = split (/\s+/, $prediction_header);
    my $coordspan = $header_comps[6];
    my ($lend, $rend) = split (/-/, $coordspan);
    $lend += $partition_lend - 1;
    $rend += $partition_lend -1;
    
    $header_comps[6] = "$lend-$rend";
    
    $prediction_header = join (" ", @header_comps);
    
    $new_prediction_text .= "$prediction_header\n";
    
    foreach my $prediction_line (@prediction_lines) {
        my @x = split (/\t/, $prediction_line);
        $x[0] += $partition_lend -1;
        $x[1] += $partition_lend -1;
        
        push (@coords, $x[0], $x[1]);
        push (@exon_types, $x[2]);
        
        my $new_pred_line = join ("\t", @x) . "\n";
        $new_prediction_text .= $new_pred_line;
    }
    
    my $exon_type_string = join (",", @exon_types);
    
    my $prediction_class = "";
    
    if ( 
         ($exon_type_string =~ /initial/ && $exon_type_string =~ /terminal/) 
         ||
         ($exon_type_string =~ /single/)
         
         ) {
        $prediction_class = "complete";
    } else {
        $prediction_class = "partial";
    }
    
    
    @coords = sort {$a<=>$b} @coords;
    my $gene_lend = shift @coords;
    my $gene_rend = pop @coords;
    
    my $length = $gene_rend - $gene_lend + 1;
    
    
    my $prediction_struct = { lend => $gene_lend,
                              rend => $gene_rend,
                              class => $prediction_class,
                              text => $new_prediction_text,
                              length => $length,
                              path_score => $length,
                              prev_link => 0,
                              intronic_preds => undef,
                              
                              };
    
    push (@$predictions_aref, $prediction_struct);
}


####
sub combine_predictions {
    my ($predictions_aref) = @_;
    @$predictions_aref = sort {$a->{lend}<=>$b->{lend}} @$predictions_aref;
    
    ## Find the maximal set of complete non-overlapping genes.
    ## Allow partials only at beginning or end.
    
    ## doing full dynamic programming.
    
    for (my $i = 1; $i <= $#$predictions_aref; $i++) {
        
        my $struct_i = $predictions_aref->[$i];
        
       # print Dumper ($struct_i);
       # next;
        

        my $class_i = $struct_i->{class};
        my $path_score_i = $struct_i->{path_score};
        my $lend_i = $struct_i->{lend};
        my $length_i = $struct_i->{length};
        
        for (my $j = $i-1; $j >= 0; $j--) {
            
            
            print "Comparing $i vs. $j\n" if $SEE;
            
            my $struct_j = $predictions_aref->[$j];
            
            my $class_j = $struct_j->{class};
            my $path_score_j = $struct_j->{path_score};
            my $prev_link_j = $struct_j->{prev_link};
            my $rend_j = $struct_j->{rend};
            
            if ( 
                 ($rend_j < $lend_i)  
                 #&&
                 #($class_j eq "complete" || $prev_link_j == 0)
                 ) {
                
                my $new_path_score = $path_score_j + $length_i;
                print "Compatible. ($new_path_score)\n" if $SEE;
                if ($new_path_score > $path_score_i) {
                    $struct_i->{path_score} = $path_score_i = $new_path_score;
                    $struct_i->{prev_link} = $struct_j;
                    print "\tmaking best.\n" if $SEE;
                }
            }
        }
    }
    
    ## find the highest path score
    
    my $highest_path_score = 0;
    my $highest_path_struct = 0;
    foreach my $prediction (@$predictions_aref) {
        if ($prediction->{path_score} > $highest_path_score) {
            $highest_path_score = $prediction->{path_score};
            $highest_path_struct = $prediction;
        }
    }
    
    my @final_predictions;
    while ($highest_path_struct != 0) {
        push (@final_predictions, $highest_path_struct);
        $highest_path_struct = $highest_path_struct->{prev_link};
    }
    
    @final_predictions = reverse @final_predictions; #want in ascending order
    return (@final_predictions);
}






