#!/usr/bin/perl
# 
# Process a TRF dat file to find potential candidates for higher order structure (HOS)
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict; use warnings;

my $usage = "
trf_hos_finder.pl
-----------------

Usage: trf_hos_finder.pl <trf dat file>

Processes a TRF *.dat output file to look for candidate higher order structure (HOS)
Looks for following criteria:

1) More than 1 tandem repeat in a sequence
2) multiple tandem repeats occupy approximately same range in sequence
3) One repeat is approximately double the size of a shorter repeat
4) TRF score of longer repeat is at least 10% higher than score of shorter repeat
5) average %identity of repeat units within longer repeat is >2% higher than in shorter repeat

Reports details of shorter and longer repeats that satisfy criteria 1–3. Additionally adds
'hos' to penultimae column of output if either criterion 4 or 5 is met. Adds 'HOS' to output
if conditions 4 AND 5 are met (these are the most likely contenders for HOS). Use sequence
IDs to go back to raw sequence files to investigate further.

Sample output (tab delimited)
-----------------------------

In this example, sequence 2 is a strong case for HOS, and sequence 1 is a weaker case. The
'hos' or 'HOS' information is only included in the line corresponding to the longer of the
two repeats.

ID	LEVEL	START	END	LENGTH	COPIES	SCORE	%IDENT	HOS?	SEQ_ID
1	2		337		844	122		4.2		299		68				gnl|ti|2250106555 GGZG8286.b1
1	2		337		844	244		2.1		330		69		hos		gnl|ti|2250106555 GGZG8286.b1
2	2		51		948	164		5.5		590		70				gnl|ti|2250104470 GGZG7125.g1
2	2		51		951	328		2.8		735		86		HOS		gnl|ti|2250104470 GGZG7125.g1
3	2		54		641	125		4.7		467		83				gnl|ti|2250102961 GGZG6288.g1
3	2		54		641	250		2.3		479		83				gnl|ti|2250102961 GGZG6288.g1

";

die "$usage" unless (@ARGV == 1);

# keep track of how many sequences in TRF dat file
my $seq_counter = 0;

# keep track of how many repeats are found in any one sequence
# will store this data in an array
my $seq_repeat_counter = 0; # count repeats belonging to each read
my @seq_repeat_data;

# want to capture the (slightly processed) FASTA headers in the trf output
my $header;

# how different can start and end coordinates of 2 diff tandem repeats
# be in order to be considered occupying the same span? Default = 10 nt.
my $offset = 10;

# master hash to store all of the data for potential HOS repeats
my %hos;

# how much higher should the score of the longer repeat be in relation to shorter
# repeat. Try 15%
my $score_threshold = 1.15;

# and what about the absolute difference in %identity
# start with +5% in longer repeat
my $identity_threshold = 5;

################################
# MAIN LOOP OVER DAT FILE
################################

REPEAT: while(<>){

	# skip blank likes
	next if (m/^$/);

	if (m/^Sequence: (.*)/) {
		$header = $1;
		$seq_counter++;
		
		# reset certain counters and data
		$seq_repeat_counter = 0;
		@seq_repeat_data = ();

		# and now move to next line in file
		next REPEAT;
	}


	# the main output that we are interested in will all be on one line which starts with various
	# numerical details of the repeat
	if (m/^\d+ \d+ \d+ \d+\.\d /){

		# capture repeat data into various variables (most will be unsused)
		my ($start,$end,$period,$copies,$consensus,$matches,$indels,$score,$a,$c,$g,$t,$entropy,$repeat_seq) = split(/\s+/);

		# if we get this far then we will capture some of the repeat data in a hash 
		# this is to potentially compare to other repeats in the same sequence
		$seq_repeat_counter++;
		$seq_repeat_data[$seq_repeat_counter]{start}     = "$start";
		$seq_repeat_data[$seq_repeat_counter]{end}       = "$end";
		$seq_repeat_data[$seq_repeat_counter]{consensus} = "$consensus";
		
		# add repeat info to %HOS hash, but only if it is first repeat in sequence
		my $key = "$header,$seq_counter";
		my $info = "$start,$end,$consensus,$copies,1,$score,$matches";
		push(@{$hos{$key}}, $info) if ($seq_repeat_counter == 1);
					
		# no point going any further until we have seen at least 2 repeats within the current sequence
		next REPEAT unless ($seq_repeat_counter > 1);
		
		# loop through previous repeats. Want to see if one is a multiple of another
		for (my $i = 1; $i < $seq_repeat_counter; $i++){

			# want to see if current start/end coordinates are identical or within $offset nt of previous repeat
			next unless (abs($start - $seq_repeat_data[$i]{start}) <= $offset);
			next unless (abs($end   - $seq_repeat_data[$i]{end})   <= $offset);

			# $ratio is the ratio of the longest repeat to the shorter one, need longer repeat to be at least 
			# 1.85x length of shorter repeat. How we calculate this depends on whether current repeat length is 
			# longer than previous one or not
			my $ratio;
			
			if ($consensus < $seq_repeat_data[$i]{consensus}){
				$ratio = $seq_repeat_data[$i]{'consensus'}/ $consensus;
			} else{
				$ratio = $consensus / $seq_repeat_data[$i]{'consensus'};
			}
			my $processed_ratio = $ratio;
			$processed_ratio =~ s/\d+\.(\d+)/0\.$1/;

			next unless (($processed_ratio < 0.15 or $processed_ratio > 0.85) && ($processed_ratio > 1.5));
			
			# we should now be looking at a multiple repeat so add it to hash
			my $info = "$start,$end,$consensus,$copies,$ratio,$score,$matches";
			push(@{$hos{$key}}, $info);
							
			# if we have seen a match, don't need to look any further
			next REPEAT;
		}
	}
}



print "ID\tLEVEL\tSTART\tEND\tLENGTH\tCOPIES\tSCORE\t%IDENT\tHOS?\tSEQ_ID\n";
my $counter = 1;
HOS: foreach my $key (keys %hos){
	my ($seq_id, $seq_counter) = split(/,/,$key);

	# only want to look at sequences with multiple levels of (potential) HOS
	my $levels = @{$hos{$key}};
	next if ($levels == 1);

		
	# now loop over the different levels of HOS present
	LEVELS: for (my $i = 0; $i < @{$hos{$key}}; $i++){
		my ($start, $end, $repeat_length, $copies, $ratio, $score, $identity) = split(/,/, ${$hos{$key}}[$i]);
		print "$counter\t$levels\t$start\t$end\t$repeat_length\t$copies\t$score\t$identity\t";
		

		# can't easily analyze higher level structures (level 3 or greater) for HOS
		# so will just flag these with '???'
		if ($levels > 2){
			print "???\t";
		}
		# compare to previous repeat (but not when we have only seen 1 repeat)
		elsif ($i > 0){
			my (undef, undef, $prev_length, undef, undef, $prev_score, $prev_identity) = split(/,/, ${$hos{$key}}[$i-1]);

			# is score of longer repeat 10% greater compared to shorter repeat?
			# is average percentage identity of longer repeat 2% greater compared to shorter repeat?
			# print 'hos' or 'HOS' in final output to signify weak or high confidence that this is a HOS repeat

			# also need to double check whether first repeat in pair of repeats is shorter (sometimes the longer
			# repeat is reported first). If this happens, reverse scores and identities			
			if($prev_length > $repeat_length){
				($score, $prev_score)       = ($prev_score, $score);
				($identity, $prev_identity) = ($prev_identity, $identity);
			} 
			
			if (($score / $prev_score > $score_threshold) and ($identity - $prev_identity > $identity_threshold)){
				print "HOS\t";
			}
			elsif (($score / $prev_score > $score_threshold) or  ($identity - $prev_identity > $identity_threshold)){
				print "hos\t";
			} else{
				print "\t";
			}
		} else{
			print "\t";
		}
		print "$seq_id\n";			
	}
	$counter++;

}

print "\n";
