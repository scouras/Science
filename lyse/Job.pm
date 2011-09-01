#!/usr/bin/perl


use strict;
use warnings;
use Utility qw(:all);
use Data::Dumper;

# Not actually exporting anything
require Exporter;
our @ISA = qw(Exporter);
our %EXPORT_TAGS = ('all' => [qw()]);
our @EXPORT_OK = (@{$EXPORT_TAGS{'all'}});
our @EXPORT = qw();

# CVS VERSION INFORMATION (AUTOMAGICALLY GENERATED BY CVS)
'$Revision: 1.4 $'      =~ /^.Revision: (.*) \$/; our $REVISION   = $1;
'$Date: 2009/02/06 02:36:42 $'  =~ /^.Date: (.*) \$/;     our $CHECKED_IN = $1;
'$Author: scouras $'            =~ /^.Author: (.*) \$/;   our $AUTHOR     = $1;




####################################################################
#                                                                JOB
#-------------------------------------------------------------------
# Jobs needs to know:
#   1) ilmm module (dir, time, etc.)
#   2) analysis type
#   3) analysis subtype
#   4) completed
#   5) queueable
#   6) command
#
# Job statuses:
#   created     Virtual status, as we typically prep a job as soon as we know of it
#   prepped     Prepped and ready to be queued
#   executed    Analysis now finished
#   valid       Analysis validated
#   plotted     Analysis plotted
#   complete    Analysis complete (for now, same as plotted)
#   error       Something went wrong with the analysis
#   skipped     Intentionally ignored for some reason
####################################################################

package Lyse::Job;

sub new { 

  my $anal        = $_[0];
  my $ILMM        = $_[1];
  my $region      = $_[2];
  my $type        = $_[3];
  my $regional    = $_[4] || 0;
  
  if ( not $type ) { die "Lyse::Job::new - You forgot to add a title." }

  my $job = {
    anal      => $anal,
    ilmm      => $ILMM,
    region    => $region,
    type      => $type,
    protein   => $ILMM->{'protein'},
    period    => $region->{'period'},

    # various status indicators
    status => {
      created   => 1,
      prepped   => 0,
      queued    => 0,
      executed  => 0,
      valid     => 0,
      plotted   => 0,
      complete  => 0,
      skipped   => 0,
      initial   => '',
      error     => '',
    },
  };

  # add standard titles
  $job->{'region_title'     } = $region->{'title'};
  $job->{'region_title_full'} = $region->{'title_full'};
  if ( $regional ) { 
    $job->{'job_title'      } = "$type $region->{'title'}";
    $job->{'job_title_full' } = "$type $region->{'title_full'}";
    $job->{'ilmm_title'     } = lc Utility::Clean_Filename ( 
                                ($anal->{'subtype'} ? "$anal->{'subtype'} " : "") 
                                . $region->{'title'} );
  } else { 
    $job->{'job_title'      } = $type;
    $job->{'job_title_full' } = $type;
    $job->{'ilmm_title'     } = lc Utility::Clean_Filename ( 
                                ($anal->{'subtype'} ? "$anal->{'subtype'} " : "") );
  }

  $job->{'file_name'        } = lc Utility::Clean_Filename($job->{'job_title'});
      
  # add paths
  $job->{'dir'              } = $ILMM->{'dir'};
  if ( $job->{'ilmm_title'} ) { $job->{'dir_anal' } = "$anal->{'dir_anal'}_$job->{'ilmm_title'}"; }
  else                        { $job->{'dir_anal' } = $anal->{'dir_anal'}; }

  if ( $job->{'dir_anal'}   ) { $job->{'path_anal'} = "$job->{'dir'}/$job->{'dir_anal'}" }
  else                        { $job->{'path_anal'} = $job->{'dir'} }

  # add queue information from anal
  if ( $anal->{'queue'} ) { 
    $job->{'runtype'        } = 'queue';
    $job->{'queue'          } = $anal->{'queue'};
  }

  return bless $job;
}


####################################################################
sub Add_Image { 

  my $job   = $_[0];
  my $title = $_[1];
  my $data  = $_[2];
  my $type  = $_[3];

  my $key_suffix = $type ? "_$type" : "";

  #print "ADDING IMAGE: $title\n";

  $job->{'images'                 } = [];
  $job->{"path_data$key_suffix"   } = "$job->{'path_anal'}/$data";
  $job->{"title_image$key_suffix" } = $title;
  $job->{"file_image$key_suffix"  } = lc Utility::Clean_Filename("$title.png");
  $job->{"path_image$key_suffix"  } = $job->{'path_anal'} . "/" 
                                   . $job->{"file_image$key_suffix"};

  #use Data::Dumper;
  #die Dumper($job);

}

sub Add_Data { 
  my $job   = $_[0];
  my $data  = $_[1];
  my $type  = $_[2];

  my $key_suffix = $type ? "_$type" : "";
  $job->{'data'                   } = [];
  $job->{"path_data$key_suffix"   } = "$job->{'path_anal'}/$data";


}



####################################################################
sub Validate_Initial { 
  my $job = $_[0];
  my $anal = $job->{'anal'};

  my ( $valid, $error ) = $anal->Validate_Job ($job);

  if ( $valid ) { 
    $job->{'status'}{'prepped' } = 1;
    $job->{'status'}{'executed'} = 1;
    $job->{'status'}{'valid'   } = 1;
  } else { 
    $job->{'status'}{'initial' } = $error;
  }
}






####################################################################
####################################################################
sub Start { 
  my $job = $_[0];
  return &{$job->{'start'}}($job);
}

####################################################################
sub Refresh { 
  my $job = $_[0];
  if ( not $job->{'refresh'} ) {  
    $job->{'status'}{'executed'} = 1;
    $job->{'status'}{'error'} = "missing refresh operation?";
    return;
  }
  my ($status, $error) =  &{$job->{'refresh'}}($job);
  if    ( $status eq 'running'  ) { $job->{'status'}{'running' } = 1 }
  elsif ( $status eq 'finished' ) { $job->{'status'}{'executed'} = 1 }
  elsif ( $status eq 'error'    ) { 
    $job->{'status'}{'executed'} = $error;
    $job->{'status'}{'error'   } = $error;
  }
}


####################################################################
sub Validate {
  my $job = $_[0];
  my $anal = $job->{'anal'};
  return $anal->Validate_Job($job);
}



####################################################################
sub Basic_Image_Check {

  my $image = $_[0];
  my $data  = $_[1];

  if ( not -e $image        ) { return ( 0, "missing plot"            ) }
  if ( -z $image            ) { return ( 0, "plot is empty"           ) }
  if ( -M $image > -M $data ) { return ( 0, "plot is older than data" ) }

  return 1;
}


sub Basic_Data_Check {

  my $dir  = $_[0];
  my $data = $_[1];

  if ( not $dir  ) { die "Job::Basic_Data_Check - no dir passed" }
  if ( not $data ) { die "Job::Basic_Data_Check - no data passed" }

  # Verify that dir and data files exist
  if ( not $dir     ) { return (0, "no dir") }
  if ( not -e $dir  ) { return (0, "no dir") }
  if ( not -d $dir  ) { return (0, "no dir") }
  if ( not -e $data ) { return (0, "no data file") }
  if ( -z $data     ) { return (0, "empty data file") }

  return 1;

}


sub Basic_Data_Time_Check {

  my $file        = $_[0];
  my $start       = $_[1];
  my $finish      = $_[2];
  my $resolution  = $_[3];
  my $leniency    = $_[4] || 1000;
  my $rewind      = $_[5] || 0;

  #print "Time check on $file, $start, $finish, $resolution\n";
  #print "Testing " . (join ',', @sample_times) . "\n";
  my $sample_size = 10;

  if ( not -e $file ) { return (0, "no data file") }
  if ( -z $file     ) { return (0, "empty data file") }
  my $data_head = `head -n $sample_size $file`;
  my $data_tail = `tail -n $sample_size $file`;
  #my $data = `cat $file`;

  my @head = sort { $a<=>$b } ( $data_head =~ /^\s*(\d+)\./mg );
#  my @all  = sort { $a<=>$b } ( $data      =~ /^\s*(\d+)\./mg );
  my @tail = sort { $a<=>$b } ( $data_tail =~ /^\s*(\d+)\./mg );

  my $final_time = $tail[-1];

  return (0, "short ($final_time)") 
    if ( ( $final_time + $leniency + $resolution ) < $finish );

  my $head_res = abs($head[-2]-$head[-1]);
  if ( $head_res > $resolution ) { 
    return (0, "insufficient resolution in head ($head_res)");
  }
  my $tail_res = abs($tail[-2]-$tail[-1]);
  if ( $tail_res > $resolution ) { 
    return (0, "insufficient resolution in tail ($tail_res)");
  }

  return (0, "long ($final_time)") 
    if ( ( ($final_time - $leniency - $resolution) > $finish ) and $rewind );

  return 1;

#  my $missing;
#  my @times  = Utility::Range ( $start, $finish, $resolution);
#  my @sample = @times[0..$sample_size,-$sample_size..-2];
#  foreach my $time ( @sample ) { 
#    next if $time > $final_time;
#    if ( ( not grep { $time == $_ } @head ) 
#     and ( not grep { $time == $_ } @tail ) 
#     and ( not grep { $time == $_ } @all  ) ) { 
#      $missing = $time;
#      last;
#    }
#  }
#  return 1 if not defined $missing;
#  
#  @times  = Utility::Range ( $start+$resolution-1, $finish, $resolution);
#  @sample = @times[0..$sample_size,-$sample_size..-2];
#  foreach my $time ( @sample ) { 
#    next if $time > $final_time;
#    if ( ( not grep { $time == $_ } @head ) 
#     and ( not grep { $time == $_ } @tail ) 
#     and ( not grep { $time == $_ } @all  ) ) { 
#      $missing = $time;
#      last;
#    }
#  }
#  
#  return 1 if not defined $missing;
#  return (0, "missing timepoint $missing" );
}



####################################################################
#                                                             REGION
####################################################################

package Lyse::Region;

sub new {
  
  my $sim_group = $_[0];
  my $sim_run   = $_[1];
  my $title     = $_[2];
  my $selection = $_[3];
  my $molecule  = $_[4] || 0;
  my $system    = $_[5] || 0;
  my $period    = $_[6] || [];
  my %extra     = $_[7] ? %{$_[7]} : ();
  my $jobs      = [];

  foreach my $n ( 'jobs', 'title', 'whole_molecule', 'selection' ) {
    delete $extra{$n};
  }
  
  my $region = {
    sim_group       => $sim_group,
    sim_run         => $sim_run,
    title           => $title,
    title_full      => "$sim_group $sim_run $title",
    full_name       => $title,
    file_name       => lc Utility::Clean_Filename($title),
    selection       => $selection,
    molecule        => $molecule,
    system          => $system,
    period          => (defined $period and defined $period->[0] ? $period : undef),
    jobs            => $jobs,
    start_residue   => 1,
    %extra,
  };

  return bless $region;
}



