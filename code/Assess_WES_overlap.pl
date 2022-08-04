use strict;
use warnings;

### Save UKB (xgen) capture bases - NB chrs represented with the chr
my $capture_file = "xgen_plus_spikein.GRCh38.bed";
my %in_capture;

print "Saving capture regions and 50bp buffer regions either side\n";

open(XGEN, $capture_file)||die "Cannot open $capture_file";
while(<XGEN>)
{
	chomp;
	my ($chr, $start, $end) = split /\t/, $_;
	for(my $n=$start; $n<=$end; $n++)
	{
		my $pos = $chr.":".$n;
		$in_capture{$pos}="Captured";
	}

	for(my $m=$start-50; $m<$start; $m++)
	{
		my $pos = $chr.":".$m;
		$in_capture{$pos}="Buffer";
	}

	for(my $p=($end+1); $p<=$end+50; $p++)
	{
		my $pos = $chr.":".$p;
		$in_capture{$pos}="Buffer";
	}
}

### Look at gff regions

print "Looking at gff data\n";

my %ncRNA_trans;
my %protein_coding_trans;

my $file = "../gff_files/Homo_sapiens.GRCh38.107.gff3.gz";
my %chrs;

print "Saving ncRNA transcript IDs\n";

open(FILE, "gunzip -c $file |")||die "Cannot open $file";
while(<FILE>)
{
	chomp;
	if ($_ =~ /^#/){}
	else
	{
		my @line = split /\t/, $_;
		if ($line[8] =~ /^ID=transcript:(.+?);/)
		{
			my $trans_ID = $1;
			my $biotype = $line[2];

			if (($line[0] =~ /^KZ/)||($line[0] =~ /^KN/)||($line[0] =~ /Unmapped_Scaffold/)||($line[0] =~ /^GL/)||($line[0] =~ /^JH/)||($line[0] =~ /^KI/)){}
			else
			{
				$chrs{$line[0]}=1;

				if ($biotype eq "mRNA")
				{
					$protein_coding_trans{$trans_ID}=$biotype;
				}
				elsif ($biotype =~ /RNA/)
				{
					$ncRNA_trans{$trans_ID}=$biotype;
				}
			}
		}
	}
}
close(FILE);

print "Saving region annotations\n";

my @CDS_regions;
my @fiveUTR_regions;
my @threeUTR_regions;
my @ncRNA_regions;
my @other_regions;

open(REFILE, "gunzip -c $file |")||die "Cannot open $file";
while(<REFILE>)
{
	chomp;
	if ($_ =~ /^#/){}
	else
	{
		my @line = split /\t/, $_;
		my $reg = $line[0].":".$line[3].":".$line[4];
		if ($line[2] eq "CDS")
		{
			push(@CDS_regions, $reg);
		}
		elsif ($line[2] eq "five_prime_UTR")
		{
			push(@fiveUTR_regions, $reg);
		}
		elsif ($line[2] eq "three_prime_UTR")
		{
			push(@threeUTR_regions, $reg);
		}
		elsif ($line[2] eq "exon")
		{
			if ($line[8] =~ /^Parent=transcript:(.+?);/)
			{
				my $this_trans=$1;
				if (exists $ncRNA_trans{$this_trans})
				{
					push(@ncRNA_regions, $reg);
				}
				elsif (exists $protein_coding_trans{$this_trans}){}
				else
				{
					push(@other_regions, $reg);
				}
			}
		}
	}
}
close(REFILE);

print "Calculating overlap\n";

my $total_cds=0;
my $total_rna=0;
my $total_5utr=0;
my $total_3utr=0;
my $total_other=0;

my $capture_cds=0;
my $capture_rna=0;
my $capture_5utr=0;
my $capture_3utr=0;
my $capture_other=0;

my $buffer_cds=0;
my $buffer_rna=0;
my $buffer_5utr=0;
my $buffer_3utr=0;
my $buffer_other=0;

my @chromos = keys %chrs;


foreach my $chromosome (@chromos)
{
	my %found_base;
	print "$chromosome\n";

	print "Looking at CDS\n";

	foreach my $cds_reg (@CDS_regions)
	{
		my ($chr, $from, $to) = split /:/, $cds_reg;
		if ($chr eq $chromosome)
		{
			for (my $n=$from; $n<=$to; $n++)
			{
				my $to_find="chr".$chr.":".$n;
				if (exists $found_base{$to_find}){}
				else
				{
					if (exists $in_capture{$to_find})
					{
						if ($in_capture{$to_find} eq "Captured")
						{
							$capture_cds++;
						}
						elsif ($in_capture{$to_find} eq "Buffer")
						{
							$buffer_cds++;
						}
					}
					$total_cds++;
					$found_base{$to_find}=1;
				}
			}
		}
	}

	print "Looking at 5primeUTR\n";

	foreach my $fiveUTR_reg (@fiveUTR_regions)
	{
		my ($chr, $from, $to) = split /:/, $fiveUTR_reg;
		if ($chr eq $chromosome)
		{
			for (my $n=$from; $n<=$to; $n++)
			{
				my $to_find="chr".$chr.":".$n;
				if (exists $found_base{$to_find}){}
				else
				{
					if (exists $in_capture{$to_find})
					{
						if ($in_capture{$to_find} eq "Captured")
						{
							$capture_5utr++;
						}
						elsif ($in_capture{$to_find} eq "Buffer")
						{
							$buffer_5utr++;
						}
					}
					$total_5utr++;
					$found_base{$to_find}=1;
				}
			}
		}
	}

	print "Looking at 3primeUTR\n";

	foreach my $threeUTR_reg (@threeUTR_regions)
	{
		my ($chr, $from, $to) = split /:/, $threeUTR_reg;
		if ($chr eq $chromosome)
		{
			for (my $n=$from; $n<=$to; $n++)
			{
				my $to_find="chr".$chr.":".$n;
				if (exists $found_base{$to_find}){}
				else
				{
					if (exists $in_capture{$to_find})
					{
						if ($in_capture{$to_find} eq "Captured")
						{
							$capture_3utr++;
						}
						elsif ($in_capture{$to_find} eq "Buffer")
						{
							$buffer_3utr++;
						}
					}
					$total_3utr++;
					$found_base{$to_find}=1;
				}
			}
		}
	}

	print "Looking at ncRNA\n";

	foreach my $ncRNA_reg (@ncRNA_regions)
	{
		my ($chr, $from, $to) = split /:/, $ncRNA_reg;
		if ($chr eq $chromosome)
		{
			for (my $n=$from; $n<=$to; $n++)
			{
				my $to_find="chr".$chr.":".$n;
				if (exists $found_base{$to_find}){}
				else
				{
					if (exists $in_capture{$to_find})
					{
						if ($in_capture{$to_find} eq "Captured")
						{
							$capture_rna++;
						}
						elsif ($in_capture{$to_find} eq "Buffer")
						{
							$buffer_rna++;
						}
					}
					$total_rna++;
					$found_base{$to_find}=1;
				}
			}
		}
	}

	print "Looking at other exonic regions\n";

	foreach my $pseudo_reg (@other_regions)
	{
		my ($chr, $from, $to) = split /:/, $pseudo_reg;
		if ($chr eq $chromosome)
		{
			for (my $n=$from; $n<=$to; $n++)
			{
				my $to_find="chr".$chr.":".$n;
				if (exists $found_base{$to_find}){}
				else
				{
					if (exists $in_capture{$to_find})
					{
						if ($in_capture{$to_find} eq "Captured")
						{
							$capture_other++;
						}
						elsif ($in_capture{$to_find} eq "Buffer")
						{
							$buffer_other++;
						}
					}
					$total_other++;
					$found_base{$to_find}=1;
				}
			}
		}
	}
	undef %found_base;
}

open(OUT, ">WES_captureOverlap_new.txt");
print OUT "ExonType\tTotalSize\tInCapture\tInBuffer\n";
print OUT "CDS\t$total_cds\t$capture_cds\t$buffer_cds\n";
print OUT "NonCodingRNA\t$total_rna\t$capture_rna\t$buffer_rna\n";
print OUT "UTR_5prime\t$total_5utr\t$capture_5utr\t$buffer_5utr\n";
print OUT "UTR_3prime\t$total_3utr\t$capture_3utr\t$buffer_3utr\n";
print OUT "Other\t$total_other\t$capture_other\t$buffer_other\n";
