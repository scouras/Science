package RotLib::Library;

#use 5.008006;
use strict;
use warnings;
BEGIN { push @INC, "/net/programs/perl/lib64/perl5/site_perl/5.8.5/x86_64-linux-thread-multi/" }
use Data::Dumper;
use Utility qw(:all);
use PerlIO::gzip;
use Carp;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use RotLib ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);


# CVS VERSION INFORMATION (AUTOMAGICALLY GENERATED BY CVS)
'$Revision: 1.1.1.1 $'              =~ /^.Revision: (.*) \$/; our $VERSION    = $1;
'$Date: 2009/03/25 18:32:30 $'  =~ /^.Date: (.*) \$/;     our $CHECKED_IN = $1;
'$Author: scouras $'            =~ /^.Author: (.*) \$/;   our $AUTHOR     = $1;



my %CONFIG_LIB = (

  files => {

    assignment_dir        => 'assignment',

    # OUTPUT
    library_bbdep         => 'library.bbdep.dat',
    library_bbind         => 'library.bbind.dat',
    library_loose         => 'library.assign_loose.dat',
    library_strict        => 'library.assign_strict.dat',
    dh_hist_all           => 'dihedral_histograms.all.dat',
    dh_hist_assign_strict => 'dihedral_histograms.assign_strict.dat',
    dh_hist_assign_loose  => 'dihedral_histograms.assign_loose.dat',
  },

  bin_width         => 360,
  strict_assignment => 0.5,
  resolution        => 'bbd',
  metadata          => [qw(total_rotamers rotamer_names residue res_name)],

);


###################################### INITIALIZE LIB

sub Initialize {

  my $class     = $_[0];
  my $CONFIG    = $_[1];
  my $AA        = $_[2];
  my $PROTEINS  = $_[3];

  my $LIB = { %CONFIG_LIB };
  
  bless $LIB, $class;

  $LIB->{'AA'} = $AA;
  $LIB->{'main'}                  = $CONFIG;
  $LIB->{'bin_count'}             = POSIX::ceil ( 360 / $LIB->{'bin_width'} );

  my $ROT = $CONFIG->{'ROTAMER_OBJECT'};
  $LIB->{'highest_rotamer_count'} = $ROT->{'highest_rotamer_count'} || 1;

  return $LIB;

};

############################################## INITIALIZE AGGREGATOR

sub Initialize_Aggregator {
  my $self = $_[0];
  my $agg  = $_[1] || {};

  my $resolution = $self->{'resolution'};

  if ( $resolution eq 'bbd' ) { 
    $agg->{'bbd'} = $self->Initialize_Aggregator_BBD ( $agg->{'bbd'} );
  }
  
  if ( $resolution eq 'bbd' or $resolution eq 'bbi' ) {
    $agg->{'bbi'} = $self->Initialize_Aggregator_BBI ( $agg->{'bbi'} );
  } 

  $agg->{'global'} = $self->Initialize_Aggregator_Global ( $agg->{'global'} );
  $self->Initialize_Aggregator_Metadata ( $agg );

  if ( ( not exists $agg->{'for'} ) or ( not $agg->{'for'} =~ /assign/ )) {
    $agg->{'assignment'}{'for'} = 'container_assign';
    $agg->{'assignment'} = $self->Initialize_Aggregator ( $agg->{'assignment'} );
    delete $agg->{'assignment'}{'bbd'};
  }
 
  if ( not exists $agg->{'for'} ) {
    $agg->{'for'} = 'container';
  }

  return $agg;
}

sub Initialize_Aggregator_Metadata {
  my $self = $_[0];
  my $agg  = $_[1] || {};

  $agg->{'freq'       } = 0;
  $agg->{'count'      } = 0;
  $agg->{'normalized' } = 0;
  $agg->{'res_name'   } = '';
  $agg->{'generation' } = 0;
  $agg->{'initialized'} = 1;
  $agg->{'populated'  } = 0;

  return $agg;
}


sub Initialize_Aggregator_BBD {

  my $self = $_[0];
  my $agg  = $_[1] || {};

  if ( defined $agg->{'for'} and $agg->{'for'} ne 'bbd' ) 
    { confess "Trying to init $agg->{'for'} as BBD" . (Dumper ( $agg )) . "\n"; }

  my $bin_count = $self->{'bin_count'};
  foreach my $phi ( 0..$bin_count-1 ) {
    foreach my $psi ( 0..$bin_count-1 ) {
      $agg->{'lib'}[$phi][$psi] = $self->Initialize_Aggregator_BBI ( $agg->{'lib'}[$phi][$psi] );
    }
  }
  $agg->{'for'} = 'bbd';
  $self->Initialize_Aggregator_Metadata ( $agg );
  return $agg;
}


sub Initialize_Aggregator_BBI {
  my $self  = $_[0];
  my $agg   = $_[1] || {};

  if ( defined $agg->{'for'} and $agg->{'for'} ne 'bbi' ) 
    { confess "Trying to init $agg->{'for'} as BBI" . (Dumper ( $agg )) . "\n"; }
  
  my $total_rotamers = $self->{'highest_rotamer_count'};

  foreach my $rot ( 0..$total_rotamers-1) {
    $agg->{'lib'}[$rot] = $self->Initialize_Aggregator_Global ( $agg->{'lib'}[$rot] );
  }
  
  $agg->{'for'} = 'bbi';
  $self->Initialize_Aggregator_Metadata ( $agg );
  return $agg;
}


sub Initialize_Aggregator_Global { 
  my $self = $_[0];
  my $agg  = $_[1] || {};

  if ( defined $agg->{'for'} and $agg->{'for'} ne 'global' ) 
    { confess "Trying to init $agg->{'for'} as global" . (Dumper ( $agg )) . "\n"; }

  $agg->{'degrees'   } = {};
  $agg->{'modes'      } = {};
  $agg->{'histograms' } = {};
#  if ( not defined ( $agg->{'histograms'} ) ) {
#    $agg->{'histograms' } = {};
#  } else {
#    for my $angle ( keys %{$agg->{'histograms'}} ) {
#      for my $i ( 0..$#{$agg->{'histograms'}{$angle}} ) {
#        if ( defined ( $agg->{'histograms'}{$angle}[$i] ) ) {
#          $agg->{'histograms'}{$angle}[$i] = 0;
#        }
#      }
#    }
#  }
  $agg->{'for'        } = 'global';
  $self->Initialize_Aggregator_Metadata ( $agg );
  
  return $agg;
}



############################################################ ANALYZE

sub Analyze {

  my $self = $_[0];
  my $data = $_[1];
  my $RES  = $_[2];

  if ( $data->{'for'} =~ /assign/ ) { confess "Library::Analyze: Analyzing assignment container." }
  if ( $data->{'for'} ne 'container' ) {
    confess "Library::Analyze: Using $data->{'for'} agg for container function.";
  }

  if ( $self->{'resolution'} eq 'bbd' ) { 
    $self->Analyze_BBD              ( $data->{'bbd'}, $RES );
    $self->Copy_Metadata            ( $data, $data->{'bbd'} );
    $self->Normalize                ( $data );
  }

  elsif ( $self->{'resolution'} eq 'bbi' ) {
    $self->Analyze_BBI              ( $data->{'bbi'}, $RES ); 
    $self->Copy_Metadata            ( $data, $data->{'bbi'} );
    $self->Normalize                ( $data );
  }

  else {
    confess "Library::Analyze: Unknown resolution $self->{'resolution'}."; 
  }

  $self->Assign_Conformation ( $data );
  #print Dumper($data->{'bbi'});
  #exit;


}


sub Get_Symmetry_Angle { 

  my $self      = $_[0];
  my $AA        = $_[1];
  my $res_name  = $_[3];
  my $angle     = $_[3];

}


sub Analyze_BBD {
  my $self = $_[0];
  my $data = $_[1];
  my $RES  = $_[2];

  if ( $data->{'for'} ne 'bbd' ) 
    { confess "Library::Analyze_BBD: Using $data->{'for'} agg in BBD function."; }

  my $timepoints  = $RES->{'timepoints'};
  my $rotamer_names  = $RES->{'rotamer_names'};
  my $total_rotamers = $RES->{'total_rotamers'};
  
  my $bin_width   = $self->{'bin_width'};
  my $bin_count   = $self->{'bin_count'};
  
  my %headers     = %{$RES->{'headers'}};
  my $i_time      = $headers{'time'};
  my $i_phi       = $headers{'phi'};
  my $i_psi       = $headers{'psi'};
  my $i_major     = $headers{'rotamer_major'};
  my ($phi_bin, $psi_bin, $major, @minors, $dh_angle);
  my $line;

  my $residue = $RES->{'residue'};
  
  my $bbi;
  my $glo; # corresponding Individual rotamer
  my $degrees;
  
  # Aggregate all timepoint data into a backbone dependent distribution
  foreach my $t ( 0..$timepoints-1 ) {
    $data->{'count'}++;
    $line    = $RES->{'data'}[$t];
    $major   = $line->[$i_major];
    @minors  = @{$residue->{'major_to_minor'}[$major]};
    $phi_bin = (sprintf "%.0f", $line->[$i_phi] / $bin_width) % $bin_count;
    $psi_bin = (sprintf "%.0f", $line->[$i_psi] / $bin_width) % $bin_count;
    $glo = $data->{'lib'}[$phi_bin][$psi_bin]{'lib'}[$line->[$i_major]];
    $glo->{'count'}++;
    foreach my $i ( 0..$#$rotamer_names ) {
      my $angle = $rotamer_names->[$i];
      $degrees = $line->[$headers{$angle}];
      $dh_angle = $residue->{'dihedrals'}{$angle};
      my $minor = $dh_angle->{'primes'}[$minors[$i]];
      my $orig_degrees = $degrees;

      if ( $dh_angle->{'symmetry'} > 1 ) { 
        my $mean = ($minor->{'min_start'} + $minor->{'max_start'}) / 2;
        if ( $degrees < $minor->{'min_start'} 
          or $degrees > $minor->{'max_start'} ) { 
      
          my $step = 360/$dh_angle->{'symmetry'};
          my $diff = sprintf "%.0f", (($mean-$degrees)/$step);
          $degrees += $diff * $step;
          #$degrees = ($degrees - 360/$dh_angle->{'symmetry'}) % (360/$dh_angle->{'symmetry'});
        }
        
        if ( $degrees < $minor->{'min_start'} or $degrees > $minor->{'max_start'} ) { 
          print "$orig_degrees -> $degrees\tMajor: $major\tMinor: $minors[$i]\t$minor->{'min_start'}\t$minor->{'max_start'}\n";
        }
      }

      #print "$line->[$i_phi]\t$phi_bin\t$line->[$i_psi]\t$psi_bin\t$degrees\t$line->[$i_major]\t$glo->{'count'}\n";
      $glo->{'degrees'  }{$angle} += $degrees;
      $glo->{'histograms'}{$angle}[(sprintf "%.0f", $degrees)%360]++;
    }
  }
  

  # Normalize all the leaf data
  foreach my $phi ( 0..$bin_count-1 ) {
    foreach my $psi ( 0..$bin_count-1 ) {
      $bbi = $data->{'lib'}[$phi][$psi];
      foreach my $rot ( 0..$total_rotamers-1 ) {
        $glo = $bbi->{'lib'}[$rot];
        if ( not exists $glo->{'count'} ) {
          confess "Library::Analyze_BBD: Leaf node uninitialized.";
        }
        
        if ( $glo->{'count'} == 0 ) { 
          $glo->{'freq'} = 0;
          next;
        }

        $glo->{'freq'} = $glo->{'count'} / $data->{'count'};
        #$bbi->{'freq'} += $glo->{'freq'};
        foreach my $angle ( keys %{$glo->{'degrees'}} ) {
          $glo->{'degrees'}{$angle} /= $glo->{'count'};
          Normalize_Histogram ( $glo->{'histograms'}{$angle} );
          #Zero_Blanks_In_Histogram ( $glo->{'histograms'}{$angle} );
          $glo->{'normalized'} = 1;
        }
      }
      $bbi->{'normalized'} = 1;
    }
  }

  # Set up basic parameters for the backbone dependent parts.
  $self->Copy_Metadata ( $data, $RES );
  $data->{'freq'          } = 1;
  $data->{'normalized'    } = 1;
  $data->{'generation'    } = 1;
  
}



sub Analyze_BBI {
  my $self = $_[0];
  my $data = $_[1];
  my $RES  = $_[2];

  confess "No Analyze BBI yet.\n";

}

################################################ ASSIGN CONFORMATION

sub Assign_Conformation {
  my $self = $_[0];
  my $data = $_[1];

  if ( $data->{'for'} eq 'container_assign' ) { return; }
  if ( $data->{'for'} ne 'container' ) { confess "Library::Assign_Conformation: Using $data->{'for'} agg for a container function." }

  if ( not $data->{'normalized'} ) { confess "Library::Assign_Conformation: Trying to assign conformation from unnormalized data." }
  if ( not $data->{'assignment'} ) { confess "Library::Assign_Conformation: Assignment not initialized for data." }

  my $rotamer_names  = $data->{'rotamer_names'};
  my $total_rotamers = $data->{'total_rotamers'};

  # Assign Major (and thus minors)
  my $major = 0;
  my $freq  = 0;
  my $bbi = $data->{'bbi'};
  foreach my $rot ( 0..$total_rotamers-1 ) {
    if ( $bbi->{'lib'}[$rot]{'freq'} > $freq ) {
      $major = $rot;
      $freq = $bbi->{'lib'}[$rot]{'freq'};
    }
  }

  my $assign = $data->{'assignment'};

  $assign->{'major'} = $major;
  $assign->{'major_freq' } = $freq;
  $assign->{'freq'} = 1;
  $assign->{'bbi'}{'lib'}[$major]{'freq'} = 1;

  # Assign Dihedral for each Angle
  foreach my $angle ( keys %{$bbi->{'lib'}[$major]{'histograms'}} ) {
    my $max_mode;
    foreach my $mode ( @{$bbi->{'lib'}[$major]{'modes'}{$angle}} ) {
      if ( (not defined $max_mode) or ($mode->{'height'} > $max_mode->{'height'} ) ) { 
        $max_mode = $mode;
      }
    }

    my $m;
    if ( not $max_mode ) { 
      #confess "Library::Assign_Conformation: No mode assigned for major $major angle $angle"; 
      $m = 0;
    } else {
      $m = $max_mode->{'mode'};
    }
    $assign->{'bbi'}{'lib'}[$major]{'degrees'}{$angle} = $m;
    $assign->{'bbi'}{'lib'}[$major]{'histograms'}{$angle}[$m] = 1;
  }
  $self->Copy_Metadata ( $assign->{'bbi'}, $data );
  $assign->{'bbi'}{'normalized'} = 1;
  $self->Copy_Metadata ( $assign, $data );
  $assign->{'normalized'} = 1;

  $self->Aggregate_BBI_to_Global ( $assign );
}




######################################## COPY METADATA TO AGGREGATES

sub Copy_Metadata {
  my $self  = $_[0];
  my $agg   = $_[1];
  my $data  = $_[2];

  if ( not $data->{'populated'} ) {
    confess "Library::Copy_Metadata: Trying to copy metadata from unpopulated source ($data->{'for'} => $agg->{'for'})." . Dumper ( $data );
  }


  if      ( $agg->{'for'} =~ /container/ ) {
    $self->Copy_Metadata_Global ( $agg, $data );
  } elsif ( $agg->{'for'} eq 'bbd' ) {
    $self->Copy_Metadata_BBD ( $agg, $data );
  } elsif ( $agg->{'for'} eq 'bbi' ) { 
    $self->Copy_Metadata_BBI ( $agg, $data );
  } else {
    $self->Copy_Metadata_Global ( $agg, $data );
  }
}


sub Copy_Metadata_BBD {
  my $self  = $_[0];
  my $agg   = $_[1];
  my $data  = $_[2];

  if ( $agg->{'for'} ne 'bbd' ) { confess "Library::Copy_Metadata_BBD: Using $agg->{'for'} agg for BBD function."; }
  $self->Copy_Metadata_Global ( $agg, $data );
  
  my $bin_count = $self->{'bin_count'};
  foreach my $phi ( 0..$bin_count-1 ) {
    foreach my $psi ( 0..$bin_count-1 ) {
      $self->Copy_Metadata_BBI ( $agg->{'lib'}[$phi][$psi], $data );
    }
  }
}


sub Copy_Metadata_BBI {
  my $self  = $_[0];
  my $agg   = $_[1];
  my $data  = $_[2];

  if ( $agg->{'for'} ne 'bbi' ) { confess "Library::Copy_Metadata_BBI: Using $agg->{'for'} agg for BBI function."; }
  $self->Copy_Metadata_Global ( $agg, $data );
}   


sub Copy_Metadata_Global { 
  my $self  = $_[0];
  my $agg   = $_[1];
  my $data  = $_[2];
  
  foreach my $key ( @{$self->{'metadata'}} ) {
    if ( not exists $data->{$key} ) {
      confess "Library::Copy_Metadata_Global: Data is missing key $key ($data->{'for'} => $agg->{'for'}).";
    }
    $agg->{$key} = $data->{$key};
  }
  $agg->{'populated'} = 1;

}
 
#========================================================= NORMALIZE

sub Normalize {
  my $self = $_[0];
  my $data = $_[1];

  return if $data->{'normalized'};

  if    ( $data->{'for'} eq 'bbd' ) { $self->Normalize_BBD ( $data ); }
  elsif ( $data->{'for'} eq 'bbi' ) { $self->Normalize_BBI ( $data ); }

  elsif ( $data->{'for'} eq 'container' ) {
    if    ( $self->{'resolution'} eq 'bbd') { 
      $self->Normalize_BBD ( $data->{'bbd'} );
      $self->Aggregate_BBD_to_BBI     ( $data );
      $self->Normalize_BBI ( $data->{'bbi'} );
      $self->Aggregate_BBI_to_Global  ( $data );
    }
    elsif ( $self->{'resolution'} eq 'bbi') { 
      $self->Normalize_BBI ( $data->{'bbi'} );
      $self->Aggregate_BBI_to_Global  ( $data );    
    }
    
    if ( $data->{'assignment'}{'populated'} ) {
      $self->Normalize_BBI ( $data->{'assignment'}{'bbi'} );
      $self->Aggregate_BBI_to_Global ( $data->{'assignment'} );
    }

  } else {
    confess "Library::Normalize unknown for $data->{'for'} or resolution $self->{'resolution'}.";
  }

  $data->{'normalized'} = 1;

}


sub Normalize_BBD {

  my $self = $_[0];
  my $data = $_[1];

  if ( $data->{'for'} ne 'bbd' ) { confess "Library::Normalize_BBD: Using $data->{'for'} agg for BBD function."; }

  if ( $data->{'normalized'} ) { return; }
  
  my $freq = $data->{'freq'};
  if ( not $freq ) { 
    $data->{'normalized'} = 1;
    return;
  }
  
  my $bin_count = $self->{'bin_count'};
  
  foreach my $phi ( 0..$bin_count-1 ) {
    foreach my $psi ( 0..$bin_count-1 ) {
      $self->Normalize_BBI ( $data->{'lib'}[$phi][$psi], $freq );
    }
  }
  $data->{'normalized'} = 1;
}

sub Normalize_BBI {
  
  my $self  = $_[0];
  my $data  = $_[1];
  my $freq  = $_[2];

  if ( $data->{'for'} ne 'bbi' ) { confess "Library::Normalize_BBI: Using $data->{'for'} agg for BBI function."; }

  if ( $data->{'normalized'} ) { return; }

  if ( not defined $freq ) { 
    if ( not exists $data->{'freq'} or not defined $data->{'freq'} ) {
      confess "Library::Normalize_BBI: No frequency data for the given BBI Rotamer.";
    }
    $freq = $data->{'freq'}; 
  }
  
  foreach my $rot ( 0..$data->{'total_rotamers'}-1 ) {
    $self->Normalize_Global ( $data->{'lib'}[$rot], $freq );
  }
  $data->{'normalized'} = 1;
}

sub Normalize_Global {
  my $self = $_[0];
  my $data = $_[1];
  my $freq = $_[2] || 1;

  return if ( $data->{'normalized'} );
  if ( defined $freq and not $freq ) {
    $data->{'normalized'} = 1;
    return;
  }

  if ( not $data->{'freq'} ) { 
    $data->{'freq'} = 0;
    $data->{'normalized'} = 1;
    return;
  }

  foreach my $angle ( keys %{$data->{'degrees'}} ) {
    $data->{'degrees'}{$angle} /= $data->{'freq'};
    Normalize_Histogram ( $data->{'histograms'}{$angle} );
    #Zero_Blanks_In_Histogram ( $data->{'histograms'}{$angle} );
    $data->{'modes'}{$angle} = 
            Find_Modes ( $data->{'histograms'}{$angle}, 
                         { base=>0, binwidth=>1, modemin=>0.0000001, circular=>1, window=>30 } );
  }
  $data->{'freq'} /= $freq;
  $data->{'normalized'} = 1;
}


############################################## AGGREGATE SIMULATIONS
 
sub Aggregate {

  my $self    = $_[0];
  my $agg     = $_[1];
  my $data    = $_[2];

  if ( $agg->{'for'} ne 'container' ) { confess "Library::Aggregate: Using $agg->{'for'} agg for container function."; }

  if      ( $self->{'resolution'} eq 'bbd' ) { 
    $self->Aggregate_BBD ( $agg->{'bbd'}, $data->{'bbd'} ) ;
    $self->Copy_Metadata ( $agg, $agg->{'bbd'} );
  } elsif ( $self->{'resolution'} eq 'bbi' ) { 
    $self->Aggregate_BBI ( $agg->{'bbi'}, $data->{'bbi'} );
    $self->Copy_Metadata ( $agg, $agg->{'bbi'} );
  } else { 
    confess "Library::Aggregate: Unknown resolution $self->{'resolution'}.";
  }

  if ( not $agg->{'assignment'}{'bbi'} ) { confess Dumper ( $agg ) };
  $self->Aggregate_BBI ( $agg->{'assignment'}{'bbi'}, $data->{'assignment'}{'bbi'} );
  $self->Copy_Metadata ( $agg->{'assignment'}, $data->{'assignment'}{'bbi'} );

}

sub Aggregate_BBD {

  my $self    = $_[0];
  my $agg     = $_[1];
  my $data    = $_[2];

  if ( $agg-> {'for'} ne 'bbd' ) { confess "Library::Aggregate_BBD: Using $agg->{'for'} agg for BBD function."; }
  if ( $data->{'for'} ne 'bbd' ) { confess "Library::Aggregate_BBD: Using $data->{'for'} data for BBD function."; }

  if ( $agg->{'normalized'} ) {
    confess "Rotamer: The aggregator has already been normalized, you can't add more data to it.";}

  if ( $agg->{'res_name'} and ($agg->{'res_name'} ne $data->{'res_name'}) ) {
      confess "State Counts: Aggregator ($agg->{'res_name'}) and Data ($data->{'res_name'}) do not address the same residue.";
  } 

  my $total_rotamers = $data->{'total_rotamers'};
  my $rotamer_names = $data->{'rotamer_names'};
  
  if ( not $data->{'normalized'} ) { $self->Normalize_BBD($data) };

  my $bin_width = $self->{'bin_width'};
  my $bin_count = $self->{'bin_count'};
  my $pps;
  foreach my $phi ( 0..$bin_count-1 ) {
    foreach my $psi ( 0..$bin_count-1 ) {
      $self->Aggregate_BBI ( $agg->{'lib'}[$phi][$psi], $data->{'lib'}[$phi][$psi] );
      $agg->{'freq'} += $agg->{'lib'}[$phi][$psi]{'freq'};
    }
  }

  $agg->{'normalized'    } = 0;
  $agg->{'generation'    } = $data->{'generation'} + 1;
  $self->Copy_Metadata ( $agg, $data );


}

sub Aggregate_BBI {

  my $self    = $_[0];
  my $agg     = $_[1];
  my $data    = $_[2];
  
  if ( $agg->{'for'} ne 'bbi' ) { confess "Library::Aggregate_BBI: Using $agg->{'for'} agg for BBI function." . Dumper ($agg); }

  if ( not $data->{'normalized'} ) { $self->Normalize_BBI($data) }

  my $total_rotamers = $data->{'total_rotamers'};
  
  foreach my $rot ( 0..$total_rotamers-1 ) {
    next if not $data->{'lib'}[$rot]{'freq'};
    $agg->{'lib'}[$rot] = {} 
      if not defined $agg->{'lib'}[$rot];
    $self->Aggregate_Global ( $agg->{'lib'}[$rot], $data->{'lib'}[$rot] );
    $agg->{'freq'} += $data->{'lib'}[$rot]{'freq'};
  }
  $agg->{'normalized'    } = 0;
  $agg->{'generation'    } = $data->{'generation'} + 1;
  $self->Copy_Metadata ( $agg, $data );

}


sub Aggregate_Global {

  my $self    = $_[0];
  my $agg     = $_[1];
  my $data    = $_[2];

  # Add frequency
  if ( not $data->{'freq'} ) { 
    $data->{'freq'} = 0; 
    return;
  }

  $agg->{'freq'} += $data->{'freq'};
  
  foreach my $angle ( keys %{$data->{'degrees'}} ) {

    
    # Avereage contribues as much as it's frequency.
    $agg->{'degrees'}{$angle} += $data->{'degrees'}{$angle} * $data->{'freq'};
    
    # Histogram also contributes proportional to freq.
    if ( not defined $agg->{'histograms'}{$angle} ) {
      $agg->{'histograms'}{$angle} = [];
    }
    Combine_Histograms( $agg->{'histograms'}{$angle},
                        $data->{'histograms'}{$angle}, 
                        1, 
                        $data->{'freq'} );
  }
}





################################################### AGGREGATE X TO Y
sub Aggregate_BBD_to_BBI {
  my $self = $_[0];
  my $agg  = $_[1];

  if ( $agg->{'for'} ne 'container' ) { confess "Library::Aggregate_BBD_to_BBI: Using $agg->{'for'} agg for container function."; }
  
  if ( not $agg->{'bbd'}{'normalized'} ) { $self->Normalize_BBD ( $agg->{'bbd'} ) }

  my $bin_count   = $self->{'bin_count'};
  my $bbd = $agg->{'bbd'};
  my $bbi = $agg->{'bbi'};

  foreach my $phi ( 0..$bin_count-1 ) {
    foreach my $psi ( 0..$bin_count-1 ) {
      $self->Aggregate_BBI ( $bbi, $bbd->{'lib'}[$phi][$psi] );
    }
  }
  $self->Copy_Metadata ( $bbi, $bbd );
  $self->Normalize_BBI ( $bbi, $bbi->{'freq'} );
}


sub Aggregate_BBI_to_Global {
  my $self = $_[0];
  my $agg  = $_[1];
  
  if ( $agg->{'for'} !~ /container/ ) { confess "Library::Aggregate_BBI_to_Global: Using $agg->{'for'} agg for container function."; }

  if ( not $agg->{'bbi'}{'normalized'} ) { $self->Normalize_BBI ( $agg->{'bbi'} ) }
  
  my $total_rotamers = $agg->{'total_rotamers'};
  my $rotamer_names = $agg->{'rotamer_names'};

  my $bbi = $agg->{'bbi'};
  my $global = $agg->{'global'};
  my $freq;

  foreach my $rot ( 0..$total_rotamers-1 ) {
    $freq = $bbi->{'lib'}[$rot]{'freq'};
    $self->Aggregate_Global ( $global, $bbi->{'lib'}[$rot], 1, $freq );
  }
  $self->Normalize_Global ( $global, $global->{'freq'} );
  $self->Copy_Metadata ( $global, $bbi );

}

################################################ GENERATE HISTOGRAMS

sub Generate_Histograms {
  my $self    = $_[0];
  my $RES     = $_[1];

  $self->Generate_Dihedral_Histograms ( $RES );
  $self->Generate_Rotamer_Histograms  ( $RES );

}






####################################################################
#                                                             OUTPUT
####################################################################

sub Output {

  my $self  = $_[0];
  my $agg   = $_[1];
  my $dir   = $_[2];

  #my $agg = $self->Initialize_Aggregator();
  #$self->Aggregate( $pre_agg, $agg );
  #$self->Aggregate_BBD_to_BBI ( $agg );

  if (not $agg->{'res_name'}) { 
    confess "Library::Output_BBI: Misformed data structure.";
  }

  $self->Normalize ( $agg );

  my $assign_dir = "$dir/$self->{'files'}{'assignment_dir'}";
  Make_Directories ($assign_dir);

  # Library Data
  if ( $self->{'resolution'} eq 'bbd' ) { 
    $self->Output_BBD ( $agg->{'bbd'}, $dir );
    $self->Output_BBI ( $agg->{'bbi'}, $dir );
  }
  if ( $self->{'resolution'} eq 'bbi' ) {
    $self->Output_BBI ( $agg->{'bbi'}, $dir );
  }

  # Assigned Library Data
  $self->Output_BBI ( $agg->{'assignment'}{'bbi'}, $assign_dir );

  # Dihedral Histogram Data
  $self->Output_Dihedral_Histograms ( $agg, $dir );
  $self->Output_Dihedral_Histograms ( $agg->{'assignment'}, $assign_dir );
}

sub Output_BBD {

  my $self  = $_[0];
  my $agg   = $_[1];
  my $dir   = $_[2];
  
  if ( $agg->{'for'} ne 'bbd' ) { confess "Library::Output_BBD: Using $agg->{'for'} agg in BBD function."; }
}


sub Output_BBI {

  my $self  = $_[0];
  my $agg   = $_[1];
  my $dir   = $_[2];

  if ( $agg->{'for'} ne 'bbi' ) { confess "Library::Output_BBI: Using $agg->{'for'} agg in BBI function."; }
  if ( not $agg->{'normalized'} ) { $self->Normalize_BBI ( $agg ); }

  my $file = "$dir/$self->{'files'}{'library_bbind'}";
  my $res_name        =   $agg->{'res_name'};
  my @rotamer_names   = @{$agg->{'rotamer_names' }};
  my $total_rotamers  = $agg->{'total_rotamers'};
  my $AA              = $agg->{'residue'};
  #my $AA              = $self->{'AA'};

  if (not $res_name) { 
    confess "Library::Output_BBI: Misformed data structure. " . Dumper ( $agg );
  }

  open OUT, ">:gzip", "$file.gz" or confess "Couldn't open backbone independent library file '$file' for writing. $!";
      
  my $header = "#RES\tMajor\t" 
             . (join "\t", map { "Minor$_" } 1..@rotamer_names )
             . "\tp(r1)\tp(r1234)\tp(r234|r1)\t" . (join "\t", @rotamer_names) . "\n";
  print OUT $header;
  
  my $format = "%s\t%i"                                    # residue name
             . ( join '', map { "\t%i" } 1..@rotamer_names )   # rotamer bins
             . "\t%.6f\t%.6f\t%.6f"                    # probabilities
             . ( join '', map { "\t%.1f" } 1..@rotamer_names )  # average angles
             . "\n"
             ;
             
  my $minors;
  my @minor_confs;
  my $rot;

  my @p_r1 = ();
  foreach my $major ( 0..$total_rotamers-1 ) { 
    $minors = $AA->{'major_to_minor'}[$major];
    $rot = $agg->{'lib'}[$major];
    $p_r1[$minors->[0]] += $rot->{'freq'};
  }

  
  foreach my $major ( 0..$total_rotamers-1 ) {
    $minors = $AA->{'major_to_minor'}[$major];
    @minor_confs = map { $AA->{'conformations'}{ $rotamer_names[$_] }[ $minors->[$_] ] } 0..$#$minors;
    $rot = $agg->{'lib'}[$major];
    my $p_r1 = $p_r1[$minors->[0]];
    printf OUT $format,
                $res_name,          # residue name
                $major+1,           # major
                @$minors,           # minors
                $p_r1,              # probabilities
                $rot->{'freq'},  
                ($p_r1 > 0 ? ($rot->{'freq'}/$p_r1) : 0),
                ( map { exists $rot->{'degrees'} and exists $rot->{'degrees'}{$_} ? $rot->{'degrees'}{$_} : 0 } @rotamer_names),
                ;
  }

  close OUT or confess "Couldn't close backbone indepdent library file '$file' after writing. $!";
    
}



sub Output_Dihedral_Histograms {

  my $self  = $_[0];
  my $agg   = $_[1];
  my $dir   = $_[2];
  
  if ( $agg->{'for'} !~ /container/ ) { confess "Library::Output_Dihedral_Histograms: Using $agg->{'for'} agg in container function." }
 
  my $bbi = $agg->{'bbi'}; 
  $self->Normalize_BBI ( $bbi );

  my $glo = $agg->{'global'};
  $self->Normalize_Global ( $glo );

  my $file = "$dir/$self->{'files'}{'dh_hist_all'}";

  my $res_name        = $agg->{'res_name'};
  my $rotamer_names   = $agg->{'rotamer_names' };
  my $total_rotamers  = $agg->{'total_rotamers'};
  #my $AA              = $self->{'AA'};
  my $AA              = $agg->{'residue'};

  my @MODE_PROPERTIES = (
        [ "MODE",         "%.6f" ],
        [ "HEIGHT",       "%.6f" ],
        [ "WIDTH",        "%.6f" ],
        [ "AREA",         "%.6f" ],
        [ "MEAN",         "%.6f" ],
        [ "STDDEV",       "%.6f" ],
        [ "SKEW",         "%.6f" ],
        [ "KURT",         "%.6f" ],
        [ "HALF_WIDTH",   "%.6f" ],
        [ "HALF_AREA",    "%.6f" ],
        [ "SKEW_RAT",     "%.6f" ],
        [ "KURT_RAT",     "%.6f" ],
        [ "WINDOW_WIDTH", "%.6f" ],
        [ "WINDOW_AREA",  "%.6f" ],
        [ "WINDOW_MEAN",  "%.6f" ],
        [ "WINDOW_STDDEV","%.6f" ],
        [ "WINDOW_SKEW",  "%.6f" ],
        [ "WINDOW_KURT",  "%.6f" ],
        [ "MAX_HEIGHT",   "%s" ],
        [ "MAX_AREA",     "%s" ],
      );


  open OUT, ">:gzip", "$file.gz" or confess "Couldn't open rotameric dihedral distribution file '$file' for writing. $!";

  # PRINT HEADERS  
  my $header_raw = "#RAW\tRES\tMAJOR\t" 
             . ( join "\t", map { "MINOR_$_" } 1..@$rotamer_names )
             . "\tFREQUENCY\tANGLE\t" . (join "\t", 0..359) . "\n";
  print OUT $header_raw;

  my $header_modes = "#MODES\tRES\tMAJOR\t"
              . ( join "\t", map { "MINOR_$_" } 1..@$rotamer_names )
              . "\tFREQUENCY\tANGLE\t"
              . (join "\t", map { $_->[0] } @MODE_PROPERTIES)
              . "\n"
              ;
#             . "\tFREQUENCY\tANGLE\tMODE\tHEIGHT\tWIDTH\tAREA\tHALF_WIDTH\tHALF_AREA\tKURT_RAT\tMEAN\tSTDDEV\tSKEW\tKURT\tW_MEAN\tW_STDDEV\tW_SKEW\tW_KURT\n";
  print OUT $header_modes;


  # LINE FORMATS
  my $format_raw = "RAW\t%s\t%i"
             . ( join '', map { "\t%i" } 1..@$rotamer_names )   # rotamer bins
             . "\t%.6f\t%s\t"
             . ( join '', map { "\t%.6f" } 0..359 )
             . "\n"
             ;

  my $format_modes = "MODES\t%s\t%i"
             . ( join '', map { "\t%i" } 1..@$rotamer_names )   # rotamer bins
             . "\t%.6f\t%s\t"
             . ( join "\t", map { $_->[1] } @MODE_PROPERTIES )
             . "\n"
#             . "\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n"
             ;

  
  my $hist;
  my $minors;

  ##### Print Global Dihedral Angles

  foreach my $angle ( @$rotamer_names ) {

    $hist = $glo->{'histograms'}{$angle} || [];

    printf OUT $format_raw,
             $res_name,
             0, 
             (map { 0 } 1..@$rotamer_names),
             $glo->{'freq'},
             $angle,
             map { defined $hist->[$_] ? $hist->[$_] : 0 } 0..359,
             ;
   
    foreach my $mode ( @{$glo->{'modes'}{$angle}} ) {

      printf OUT $format_modes,
        $res_name,
        0, 
        (map { 0 } 1..@$rotamer_names),
        $glo->{'freq'},
        $angle,
        (map { $mode->{lc $_->[0]} } @MODE_PROPERTIES),
        ;
    }
  }


  ##### Print Dihedral Angles for each Major
  foreach my $major ( 0..$total_rotamers-1 ) {
    $minors = $AA->{'major_to_minor'}[$major];
    #$minors = $self->Get_Rotamer_Minors ( $res_name, $major );
    foreach my $angle ( @$rotamer_names ) {
      
      $hist = $bbi->{'lib'}[$major]{'histograms'}{$angle} || [];

      printf OUT $format_raw,
               $res_name,
               $major+1,
               @$minors,
               $bbi->{'lib'}[$major]{'freq'},
               $angle,
               (map { defined $hist and defined $hist->[$_] ? $hist->[$_] : 0 } 0..359),
               ;
               
      foreach my $mode ( @{$bbi->{'lib'}[$major]{'modes'}{$angle}} ) {
        printf OUT $format_modes,
          $res_name,
          $major+1,
          @$minors,
          $bbi->{'lib'}[$major]{'freq'},
          $angle,
          (map { $mode->{lc $_->[0]} } @MODE_PROPERTIES),
          ;
      }
    }
  }


  close OUT or confess "Couldn't close rotameric dihedral distribution file '$file'. $!";


}


1;
__END__


