use strict;
use warnings;

my $full_file = "Saccharomyces_cerevisiae.R64-1-1.107.gff3.gz";

open(OUT, ">>ExonSizes_perOrganism.txt");

# save IDs of ncRNA transcripts

my %ncRNA_trans;
my %protein_coding_trans;
my $count_ncRNA_trans=0;
my $count_protein_coding_trans=0;
my $count_other_trans=0;

my %chrs;

open(FILE, "gunzip -c $full_file |")||die "Cannot open $full_file";
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
					$count_protein_coding_trans++;
					$protein_coding_trans{$trans_ID}=$biotype;
				}
				elsif ($biotype =~ /RNA/)
				{
					$count_ncRNA_trans++;
					$ncRNA_trans{$trans_ID}=$biotype;
				}
				else
				{
					$count_other_trans++;
				}
			}
		}
	}
}
close(FILE);

print "Found $count_protein_coding_trans protein-coding, $count_ncRNA_trans ncRNA, and $count_other_trans other transcripts\n";

# save regions in each category

my @CDS_regions;
my @fiveUTR_regions;
my @threeUTR_regions;
my @ncRNA_regions;
my @other_regions;

open(REFILE, "gunzip -c $full_file |")||die "Cannot open $full_file";
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

# Supplement with additional UTR regions

print "Adding in additional UTR regions\n";

my $count_added_five=0;
my $count_added_three=0;

my @utr_files = ("longest_five_prime_UTRs_nointrons.gff","abundant_five_prime_UTRs_nointrons.gff","abundant_three_prime_UTRs.gff","longest_three_prime_UTRs.gff");

foreach my $utr_file (@utr_files)
{
	open(UTR, $utr_file)||die "Cannot open $utr_file";
	while(<UTR>)
	{
		if ($_ =~ /^#/){}
		else
		{
			my @utr_line = split /\t/, $_;
			my $chr=$utr_line[0];
			if ($utr_line[0] =~ /chr([A-Z]+)/)
			{
				$chr=$1;
			}

			my $reg = $chr.":".$utr_line[3].":".$utr_line[4];
			if ($utr_line[2] eq "three_prime_UTR")
			{
				push(@threeUTR_regions, $reg);
				$count_added_three++;
			}
			elsif ($utr_line[2] eq "five_prime_UTR")
			{
				push(@fiveUTR_regions, $reg);
				$count_added_five++;
			}
		}
	}
	close(UTR);
}

print "Added $count_added_five 5prime and $count_added_three 3prime UTR regions\n";
print "Saved all regions. Now checking overlap and counting.\n";

my $count_CDS=0;
my $count_fiveUTR=0;
my $count_threeUTR=0;
my $count_nc=0;
my $count_other=0;

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
				my $pos=$chr.":".$n;
				if (exists $found_base{$pos}){}
				else
				{
					$count_CDS++;
					$found_base{$pos}=1;
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
				my $pos=$chr.":".$n;
				if (exists $found_base{$pos}){}
				else
				{
					$count_fiveUTR++;
					$found_base{$pos}=1;
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
				my $pos=$chr.":".$n;
				if (exists $found_base{$pos}){}
				else
				{
					$count_threeUTR++;
					$found_base{$pos}=1;
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
				my $pos=$chr.":".$n;
				if (exists $found_base{$pos}){}
				else
				{
					$count_nc++;
					$found_base{$pos}=1;
				}
			}
		}
	}

	print "Looking at all other exons\n";

	foreach my $other_reg (@other_regions)
	{
		my ($chr, $from, $to) = split /:/, $other_reg;
		if ($chr eq $chromosome)
		{
			for (my $n=$from; $n<=$to; $n++)
			{
				my $pos=$chr.":".$n;
				if (exists $found_base{$pos}){}
				else
				{
					$count_other++;
					$found_base{$pos}=1;
				}
			}
		}
	}

	undef %found_base;
}

print OUT "Saccharomyces_cerevisiae\t$count_CDS\t$count_nc\t$count_fiveUTR\t$count_threeUTR\t$count_other\n";


print "Done!!!\n";
