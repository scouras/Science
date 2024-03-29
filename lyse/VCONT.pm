#!/usr/bin/perl

package Lyse::VCONT;

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
'$Revision: 1.2 $'              =~ /^.Revision: (.*) \$/; our $REVISION   = $1;
'$Date: 2009/01/09 00:05:53 $'  =~ /^.Date: (.*) \$/;     our $CHECKED_IN = $1;
'$Author: scouras $'            =~ /^.Author: (.*) \$/;   our $AUTHOR     = $1;


############################################## DEFAULT CONFIGURATION

our %CONFIG_ANAL = (

  revision        => $REVISION,
  checked_in      => $CHECKED_IN,
  author          => $AUTHOR,

  subtype         => '',
  resolution      => 10,
  
  priority        => 0,
  threads         => 2,

  dir_anal        => 'vcont',
  file_data       => 'contact_t.dat',
  reference       => 'min.pdb',
  overwrite_anal  => 0,
  overwrite_plots => 0,
  start_emphasis  => 20,
  final_emphasis  => 20,
  width           => 1024,
  height          => 768,
  
  is_plotted      => 0,
  wait_to_finish  => 1,
  queue           => 'xs',
  #queue           => 'analysis',
  #plot            => "/users/scouras/code/scouras/plot_vcont",

);




######################################################### INITIALIZE

sub New {
  return bless {
    name => 'VCONT',
    %CONFIG_ANAL,
  };
}

sub Initialize {
  my $CONFIG_ANAL = $_[0];
  
  my %CONFIG_GLOBAL = $_[1] ? %{$_[1]} : ();
  my $ANAL = { 
    %$CONFIG_ANAL, 
    %CONFIG_GLOBAL,
  };

  $ANAL->{'title'} = $CONFIG_ANAL->{'name'};
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

  # Currently only do VCONT on whole molecules
  if ( not defined $REGION->{'molecule'} ) { return }
  if ( not defined $REGION->{'index'   } ) { return }
  if ( $REGION->{'period'}               ) { return }

  # SET UP JOB
  my $job = Lyse::Job::new ( $ANAL, $ILMM, $REGION, $ANAL->{'title'} );
  $job->{'resolution'} = $ANAL->{'resolution'};
  $job->{'line'      } = "vcont $ANAL->{'reference'}";
  #$job->{'line'      } = "vcont $REGION->{'index'} $ANAL->{'reference'}";
  $job->{'dir_anal'  } = "$job->{'dir_anal'}/$REGION->{'index'}_$REGION->{'name'}";
  $job->{'path_anal' } = "$ILMM->{'dir'}/$job->{'dir_anal'}";
  $job->Add_Image($ANAL->{'title'}, $ANAL->{'file_data'});
  $job->Validate_Initial();

#  use Data::Dumper;
#  $DATA::Dumper::Maxdepth=1;
#  die Dumper($job);

  if ( wantarray ) { return ($job) }
  else             { return [$job] }
}



###################################################### VALIDATE JOBS
sub Validate_Job {

  my $ANAL = $_[0];
  my $job  = $_[1];

  my $dir = $job->{'path_anal'};
  my $ILMM = $job->{'ilmm'};

  # Verify that dir and data files exist
  if ( $ANAL->{'overwrite_anal'}  ) { return (0, "forced overwrite" ) }
  my ($valid, $message);
  ($valid, $message) = 
    Lyse::Job::Basic_Data_Check($job->{'path_anal'}, $job->{'path_data'});
  return ($valid, $message) if not $valid;
  $job->{'status'}{'existing'} = 1;
  ($valid, $message) = 
    Lyse::Job::Basic_Data_Time_Check( $job->{'path_data'}, 
                                $ILMM->{'start_time'},
                                $ILMM->{'finish_time'},
                                $ANAL->{'resolution'},
                                $ANAL->{'leniency'},
                                $ANAL->{'rewind'},
                                 );
  return ($valid, $message) if not $valid;

  $job->{'status'}{'valid'} = 1;
  return 1;
}





1

