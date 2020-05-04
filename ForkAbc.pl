#!/usr/bin/perl
#
# run abc jobs
#

use Parallel::ForkManager;
my $max = 20;
my $pm = Parallel::ForkManager->new($max);

foreach $i (0..9){
    system "sleep 2\n";
    $pm->start and next; ## fork
    $out = "out_abc_set"."$i".".txt";
    system "./abcfs -g example_data/sim_p0.txt -e example_data/sim_env.txt -f example_data/sim_ne.txt -t example_data/sim_trait.txt -n 100000 -o $out -s 0 -a -0.1 -c 0.1 -b -0.1 -d 0.1\n";
    #system "./abcfs -g example_data/sim_p0.txt -e example_data/sim_env.txt -f example_data/sim_ne.txt -t example_data/sim_trait.txt -n 100000 -o $out -s 0 -a -0.25 -c 0.25 -b -0.25 -d .25\n";
    $pm->finish;
}

$pm->wait_all_children;


