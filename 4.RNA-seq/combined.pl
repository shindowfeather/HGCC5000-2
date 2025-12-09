use strict;

my $input = shift; 
my $input2 = shift;
my %hash;

open IN,$input or die $!;
while(<IN>){
    chomp;
    my @all = split(/\t/,$_);
    my @ids = split(/\|/,$all[3]);
 #   print "$ids[0]\n";
    if (exists $hash{$ids[0]}){
        print "replicates $ids[0]\n";
    }else{
        $hash{$ids[0]} = $all[0];
    }
}
close IN;

open IN, $input2 or die $!;
while(<IN>){
    chomp;
    my @all = split(",",$_);
    my $name = shift @all;
    my @nstrg = split(/\|/,$name);
    my $string = join("\t",@all);
    if (exists $hash{$nstrg[0]}){
        print "$nstrg[0]\t$hash{$nstrg[0]}\t$string\n";
    }else{
        if ( $nstrg[0] !~ /^MSTRG/ ){
            print "$nstrg[0]\t$nstrg[0]\t$string\n";
        }else{
            print "unknown $name\n";
        }
    }
}
close IN;
