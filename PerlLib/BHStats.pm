package BHStats;

use strict;

sub binomial_probability_sum_k_to_n {
    my ($n,$k,$p) = @_;
    my $sum = 0;
    for (my $i = $k; $i <= $n; $i++) {
	my $binProb = binomial_probability($n,$i,$p);
	$sum += $binProb;
    }
    return ($sum);
}

sub binomial_probability_sum_k_to_0 {
    my ($n,$k,$p) = @_;
    my $sum = 0;
    for (my $i = $k; $i >= 0; $i--) {
	my $binProb = binomial_probability($n,$i,$p);
	$sum += $binProb;
    }
    return ($sum);
}




sub binomial_probability {
    my ($n_observations, $k_successes, $p_probability) = @_;
    
    my ($n, $k, $p) = ($n_observations, $k_successes, $p_probability);
    
    ### Given B(n,p), find P(X=k)
    
    my $binomial_prob = binomial_coefficient($n,$k) * ($p**$k) * (1-$p)**($n-$k);

    return ($binomial_prob);
    
}


sub binomial_coefficient {
    my ($n_things, $k_at_a_time) = @_;
    
    my $number_of_k_arrangements = (factorial($n_things)) / ( factorial($k_at_a_time) * factorial($n_things-$k_at_a_time) );

    return ($number_of_k_arrangements);
}


sub factorial {
    my $x = shift;
    $x = int($x);
    my $factorial = 1;
    while ($x > 1) {
	$factorial *= $x;
	$x--;
    }
    return ($factorial);
}


sub stDev {
    # standard deviation calculation
    my @nums = @_;
    @nums = sort {$a<=>$b} @nums;


    my $avg = avg(@nums);
    my $total = $#nums + 1;
    
    ## sum up the sqr of diff from avg
    my $sum_avg_diffs_sqr = 0;
    foreach my $num (@nums) {
        my $diff = $num - $avg;
        my $sqr = $diff**2;
        $sum_avg_diffs_sqr += $sqr;
    }
    my $stdev = sqrt ( (1/($total-1)) * $sum_avg_diffs_sqr);
    return ($stdev);
}


sub stDev_from_median {
    # standard deviation calculation
    my @nums = @_;
    @nums = sort {$a<=>$b} @nums;
    
    my $total = $#nums + 1;
    
    my $median = median(@nums);
    

    ## sum up the sqr of diff from avg
    my $sum_avg_diffs_sqr = 0;
    foreach my $num (@nums) {
        my $diff = $num - $median;
        my $sqr = $diff**2;
        $sum_avg_diffs_sqr += $sqr;
    }
    my $stdev = sqrt ( (1/($total-1)) * $sum_avg_diffs_sqr);
    return ($stdev);
    
}

sub median {
    my @nums = @_;
    
	@nums = sort {$a<=>$b} @nums;

    my $count = scalar (@nums);
    if ($count %2 == 0) {
        ## even number:
        my $half = $count / 2;
        return (avg ($nums[$half-1], $nums[$half]));
    }
    else {
        ## odd number. Return middle value
        my $middle_index = int($count/2);
        return ($nums[$middle_index]);
    }
}

sub avg {
    my @nums = @_;
    my $total = $#nums + 1;
    my $sum = 0;
    foreach my $num (@nums) {
	$sum += $num;
    }
    my $avg = $sum/$total;
    return ($avg);
}


sub CorrelationCoeff {
    my ($x_aref, $y_aref) = @_;
    my @x = @$x_aref;
    my @y = @$y_aref;
    
    my $total = $#x + 1;
    my $avg_x = avg(@x);
    my $avg_y = avg(@y);
    
    my $stdev_x = stDev(@x);
    my $stdev_y = stDev(@y);
    
    # sum part of equation
    my $summation = 0;
    for (my $i = 0; $i < $total; $i++) {
	my $x_val = $x[$i];
	my $y_val = $y[$i];
	
	my $x_part = ($x_val - $avg_x)/$stdev_x;
	my $y_part = ($y_val - $avg_y)/$stdev_y;
       
	$summation += ($x_part * $y_part);
    }
    
    my $cor = (1/($total-1)) * $summation;
    
    return ($cor);
}


####
sub geometric_mean {
    my @entries = @_;
    
    my $num_entries = scalar (@entries);
    unless ($num_entries) {
	return (undef);
    }
    
    ## All entries must be > 0
    my $logsum = 0;
    foreach my $entry (@entries) {
	unless ($entry > 0) {
	    return (undef);
	}
	$logsum += log ($entry);
    }
    
    my $geo_mean = exp ( (1/$num_entries) * $logsum);
    
    return ($geo_mean);
}


####
sub min {
	my @vals = @_;

	@vals = sort {$a<=>$b} @vals;

	my $min_val = shift @vals;

	return ($min_val);
}

####
sub max {
	my @vals = @_;
	
	@vals = sort {$a<=>$b} @vals;

	my $max_val = pop @vals;

	return ($max_val);
}


####
sub sum {
	my @vals = @_;

	my $x = 0;
	foreach my $val (@vals) {
		$x += $val;
	}
	
	return ($x);
}

1; #EOM
