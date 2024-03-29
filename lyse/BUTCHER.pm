#!/usr/bin/perl

package Lyse::BUTCHER;

use strict;
use warnings;
use Utility qw(:all);
use Data::Dumper;
use Ilmm;

# Not actually exporting anything
require Exporter;
our @ISA = qw(Exporter);
our %EXPORT_TAGS = ('all' => [qw()]);
our @EXPORT_OK = (@{$EXPORT_TAGS{'all'}});
our @EXPORT = qw();

# CVS VERSION INFORMATION (AUTOMAGICALLY GENERATED BY CVS)
'$Revision: 1.3 $'              =~ /^.Revision: (.*) \$/; our $REVISION   = $1;
'$Date: 2009/01/08 23:59:51 $'  =~ /^.Date: (.*) \$/;     our $CHECKED_IN = $1;
'$Author: scouras $'            =~ /^.Author: (.*) \$/;   our $AUTHOR     = $1;


############################################## DEFAULT CONFIGURATION

our %CONFIG_ANAL = (

  revision        => $REVISION,
  checked_in      => $CHECKED_IN,
  author          => $AUTHOR,

  subtype         => '',
  resolution      => 100, 
  include_solvent => 0,
  
  dir_anal        => 'filet',
  overwrite_anal  => 0,
  overwrite_plots => 0,
  start_emphasis  => 20,
  final_emphasis  => 20,
  width           => 1024,
  height          => 768,
  
  is_plotted      => 1,
  wait_to_finish  => 1,
  queue           => 'opterons',
  plot_command    => "/users/scouras/code/scouras/renumber_pdbs",

);




######################################################### INITIALIZE

sub New {
  return bless {
    name => 'BUTCHER',
    %CONFIG_ANAL,
  };
}

sub Initialize {
  my $config_anal = $_[0];
  
  my %config_global = $_[1] ? %{$_[1]} : ();
  my $ANAL = { 
    %$config_anal, 
    %config_global,
  };

  $ANAL->{'title'} = 'BUTCHER';
  if ( $ANAL->{'subtype'} ) { 
    $ANAL->{'title'} .= " $ANAL->{'subtype'}";
  }
  return bless $ANAL;
}




########################################################### GET JOBS

sub Get_Jobs {

  my $ANAL    = $_[0];
  my $ILMM    = $_[1];
  my $REGION  = $_[2];

  # BUTCHER only executes on whole molecules and the entire simulation
  if ( not $REGION->{'system'} ) { return }
  if ( $REGION->{'period'}     ) { return }

  # SET UP JOB
  my $job = Lyse::Job::new ( $ANAL, $ILMM, $REGION, $ANAL->{'title'} );
  $job->{'targetdir'  } = $job->{'path_anal'};
  $job->{'resolution' } = 1;
  $job->{'line'       } = "pdb_writer 0 $ANAL->{'resolution'} $ANAL->{'include_solvent'}";
  $job->{'data_files' } = [];

  $job->Validate_Initial();

  if ( wantarray ) { return ($job) }
  else             { return [$job] }
}



###################################################### VALIDATE JOBS
# TODO: Make this more efficient.
sub Validate_Job {

  my $ANAL = $_[0];
  my $job  = $_[1];

  my $dir = $job->{'path_anal'};
  my $ILMM = $job->{'ilmm'};

  # Verify that dir and data files exist
  if ( $ANAL->{'overwrite_anal'}  ) { return (0, "forced overwrite" ) }

  # Filet directory exists

  if ( not -e $job->{'path_anal'} ) { return (0, "no filet dir") }
  $job->{'status'}{'existing'} = 1;

  # PDB Files exist at sufficient resolution
  
  my @pdbs = map { /([\d\.]*)ps.pdb$/; $1 } glob ( "$job->{'path_anal'}/*ps.pdb" );
  #print "Found this many pdbs: " . scalar @pdbs . "\n". (join ',', @pdbs) . "\n";
  my $start = $ILMM->{'start_time'};
  my $finish = $ILMM->{'finish_time'};
  my $leniency = $ANAL->{'leniency'};
  my $resolution = $ANAL->{'resolution'};
  my @times = Utility::Range ( $start, $finish, $resolution );
  my @times_sample = grep { defined } @times[0..10,-10..-2];

  my @missed_inc = Utility::Difference(\@times_sample, \@pdbs, 1);
  if ( grep { $_ > ($finish + $leniency) }  @missed_inc ) { return (0, "missing timepoint $missed_inc[0]") }

  #$job->{'data_files'} = [ grep {-e} map { sprintf "%s/%010.2fps.pdb", $job->{'path_anal'}, $_ } @times ];
  #print "Trying to grab these pdbs files:\n" . 
  #  ( join "\n", @{$job->{'data_files'}} )
  #  . "\n";
  #exit;

  $job->{'status'}{'valid'} = 1;
  return 1;
}


###################################################### VALIDATE PLOT
sub Validate_Plot { 

  my $ANAL = $_[0];
  my $job  = $_[1];
  my $region = $job->{'region'};

  my $start_residue = $region->{'start_residue'};

#  print "Starting residue is $start_residue\n";

  # We'll sample pdbs at 1 ns to see if they're properly numbered
  my @pdbs = sort glob ( "$job->{'path_anal'}/*000.00ps.pdb" );

  foreach my $pdb ( @pdbs ) { 
    my $line = `grep ATOM $pdb | head -n 1`;
    #print "$line\n";
    if ( $line =~ /^([\w ]{5,6})([\d ]{5,6})([\w\d ]{6})(\w{3})..([\d ]{4})....([\d.\- ]{8})([\d.\- ]{8})([\d.\- ]{8})(.*)$/ ) {
      my ($type, $a_num, $a_type, $r_type, $r_num, $x, $y, $z, $extra) 
        = ($1, $2, $3, $4, $5, $6, $7, $8, $9);
      if ( $r_num != $start_residue ) { 
        $pdb =~ /(\d+ps.pdb)$/;
        return (0, "pdb needs renumbering ($1)")
      }
    }
  }
  return 1;
}


########################################################### PLOT JOB
sub Plot_Job { 

  my $ANAL = $_[0];
  my $job = $_[1];
  my $region = $job->{'region'};

#  print "Trying to renumber pdbs\n";

  my $plot_command = 
    "$ANAL->{'plot_command'} "
    . "$job->{'path_anal'} "
    . "$region->{'start_residue'} "
    . "&> $job->{'path_anal'}/renumber_pdbs.log";

  my $result = system( $plot_command );
  if ( $result ) { 
    $job->{'status'}{'error'} = "pdb_renumber failed";
    return ( 0, $job->{'status'}{'error'} );
  }
  $job->{'status'}{'plotted'} = 1;
  return 1;
}












1

