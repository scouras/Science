#!/usr/bin/perl

package Lyse::PHIPSI;

use strict;
use warnings;
use Utility qw(:all);
use Ilmm;

# Not actually exporting anything
require Exporter;
our @ISA = qw(Exporter);
our %EXPORT_TAGS = ('all' => [qw()]);
our @EXPORT_OK = (@{$EXPORT_TAGS{'all'}});
our @EXPORT = qw();

# CVS VERSION INFORMATION (AUTOMAGICALLY GENERATED BY CVS)
'$Revision: 1.4 $'              =~ /^.Revision: (.*) \$/; our $REVISION   = $1;
'$Date: 2009/01/29 22:08:06 $'  =~ /^.Date: (.*) \$/;     our $CHECKED_IN = $1;
'$Author: scouras $'            =~ /^.Author: (.*) \$/;   our $AUTHOR     = $1;



################################################# DEFAULT CONFIGURATION

our %CONFIG_ANAL = (

  revision        => $REVISION,
  checked_in      => $CHECKED_IN,
  author          => $AUTHOR,

  subtype         => '',
  resolution      => 10,
  definitions     => '/net/programs/ilmm/lib/phipsi/ss.def',
  file_data       => 'conf.dat',
  dir_anal        => 'phipsi',

  overwrite_anal  => 0,
  overwrite_plots => 0,
  start_emphasis  => 20,
  final_emphasis  => 20,
  width           => 1024,
  height          =>  768,

  is_plotted      => 1,
  wait_to_finish    => 1,
  queue           => 'opterons',
  #queue           => 'analysis',
  plot_command    => '/users/scouras/code/scouras/plot_phipsi',

);


######################################################### INITIALIZE

sub New {
  return bless { 
    name => 'PHIPSI',
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

  # PHIPSI only executes on whole molecules
  if ( not $REGION->{'molecule'} ) { return }
  if ( $REGION->{'period'}       ) { return }


  
  # SET UP JOB
  my $job = Lyse::Job::new ( $ANAL, $ILMM, $REGION, $ANAL->{'title'} );

  $job->{'resolution' } = $ANAL->{'resolution'};
  $job->{'line'       } = "phipsi $ANAL->{'definitions'} $REGION->{'index'}";
  $job->{'dir_anal'  } = "$job->{'dir_anal'}/$REGION->{'index'}_$REGION->{'name'}";
  $job->{'path_anal' } = "$ILMM->{'dir'}/$job->{'dir_anal'}";
  $job->Add_Image("$ANAL->{'title'} - $REGION->{'title_full'}", $ANAL->{'file_data'});
  $job->Validate_Initial();

  if ( wantarray ) { return ($job) }
  else             { return [$job] }
}
   

####################################################### VALIDATE JOB

sub Validate_Job { 

  my $ANAL = $_[0];
  my $job  = $_[1];

  my $path_anal = $job->{'path_anal'};
  my $path_data = $job->{'path_data'};
  my $ILMM      = $job->{'ilmm'};

  # Verify that dir and data files exist
  if ( $ANAL->{'overwrite_anal'}  ) { return (0, "forced overwrite" ) }
  my ($valid, $message);
  ($valid, $message) = 
    Lyse::Job::Basic_Data_Check($job->{'path_anal'}, $job->{'path_data'});
  return ($valid, $message) if not $valid;
  $job->{'status'}{'existing'} = 1;

  # Verify that config file includes every time point
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


###################################################### VALIDATE PLOT

sub Validate_Plot {

  my $ANAL = $_[0];
  my $job  = $_[1];

  my $path_image = $job->{'path_image'};
  my $path_data  = $job->{'path_data'};

  if ( not -e $path_image             ) { return ( 0, "missing plot"            ) }
  if ( -z $path_image                 ) { return ( 0, "plot is empty "          ) }
  if ( $ANAL->{'overwrite_plots'}     ) { return ( 0, "forced overwrite"        ) }
  if ( -M $path_image > -M $path_data ) { return ( 0, "plot is older than data" ) }

  return 1;

}



########################################################### PLOT JOB

sub Plot_Job { 

  my $ANAL = $_[0];
  my $job = $_[1];



  my $plot_command = 
        "$ANAL->{'plot_command'} "
          . "$job->{'path_anal'} "
          . "res=$job->{'region'}{'start_residue'} "
          . "ps=$ANAL->{'resolution'} "
          . "start=$ANAL->{'start_emphasis'} "
          . "final=$ANAL->{'final_emphasis'} "
          . "title=\"$job->{'title_image'}\" "
          . "height=$ANAL->{'height'} "
          . "width=$ANAL->{'width'} "
          . "&> $job->{'path_anal'}/$job->{'region'}{'file_name'}.log"
          ;


  my $result = system ( $plot_command );
  if ( $result ) { 
    $job->{'status'}{'error'} = "plot_dssp failed";
    return 0;
  } else {
    $job->{'status'}{'plotted'} = 1;
    return 1;
  }

}







