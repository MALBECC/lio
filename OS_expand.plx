#!/usr/bin/perl -w

#########################################################################################
#                                 HELPER FUNCTIONS
#########################################################################################

sub print_integral {
  my @integral = @_;

  my $s_int = 1;
  my $m = -1;
  for my $int (@integral) {
    $m = $int if $int == $integral[$#integral];
    $s_int = 0 if $m < 0 and ${$int}[0] != 0;
  }
  if ($s_int) {
    if ($INDICES == 2) {
      return "F_mU[${m}]";
    } else {
      return "F_mT[${m}]";
    }
  }

  my $str = "";
  for my $ref (@integral) {
    if ($ref == $integral[$#integral]) {
      $str = "${str}_$ref";
      last;
    }
    my @inds = @$ref;
    if ($inds[0] == 0) {
      $str = "${str}s";
    } elsif ($inds[0] == 1) {
      $str = "${str}p";
    } elsif ($inds[0] == 2) {
      $str = "${str}d";
    } elsif ($inds[0] == 3) {
      $str = "${str}f";
    }
    for (1..$#inds) {
      $str = "$str$inds[$_]";
    }
  }
  return $str;
}

sub add_integral {
  my ($level_reqs_ref,$integral_ref) = @_;
  my @level_reqs = @$level_reqs_ref;
  my @integral = @$integral_ref;
  my $not_s = 0;
  for (@integral) {
    if ($_ != $integral[$#integral] and ${$_}[0] != 0) {
      $not_s = 1;
      last;
    }
  }
  my $not_req = 1;
  for (@level_reqs) {
    if (&same_integral($integral_ref,\@$_)) {
      $not_req = 0;
      last;
    }
  }
  push(@$level_reqs_ref,$integral_ref) if $not_s and $not_req;
}

sub same_integral {
  my ($integral1_ref,$integral2_ref) = @_;
  my @integral1 = @$integral1_ref;
  my @integral2 = @$integral2_ref;

  return 0 if ((scalar @integral1) != (scalar @integral2));

  for (0..$#integral1) {
    if ($_ == $#integral1) {
      return 0 if $integral1[$_] != $integral2[$_];
    } else {
      my @ind1 = @{$integral1[$_]};
      my @ind2 = @{$integral2[$_]};
      return 0 if ((scalar @ind1) != (scalar @ind2));

      for (0..$#ind1) {
        return 0 if $ind1[$_] != $ind2[$_];
      }
    }
  }

  return 1;
}

sub reduce_integral {
  my ($skip_ref, $integral_ref) = @_;
  my @skip = @$skip_ref;
  my @integral = @$integral_ref;
  my @new_integral = ();

  for my $ind_ref (@integral) {
    if ($ind_ref == $integral[$#integral]) {
      push(@new_integral,$ind_ref);
      last;
    }
    my @inds = @$ind_ref;
    my @new_inds = ();
    my $skip_count = 0;
    my @keep = ();
    for my $ind (1..$#inds) {
      my $skip_this = 0;
      for my $skip_ind (@skip) {
        if ($skip_ind == $inds[$ind]) {
          $skip_this = 1;
          $skip_count += 1;
          last;
        }
      }
      if (not $skip_this) {
        push @keep, $inds[$ind];
      }
    }
    push @new_inds, $inds[0] - $skip_count;
    push @new_inds, $_ for (@keep);
    push @new_integral, \@new_inds;
  }

  return @new_integral;
}

sub inc_m {
  my @integral = @_;
  my @new_integral = ();
  push @new_integral, $integral[$_] for (0..$#integral-1);
  push @new_integral, $integral[$#integral]+1;
  return @new_integral;
}

#########################################################################################
#                               PRINT FUNCTIONS
#########################################################################################

sub print_one_ind_rule {
  my ($integral,$level,$ind_val) = @_;

  my $center_num = $INDEX_MAP[$level-1];
  my $str_int = &print_integral(@$integral);

  my @skip = ($level);
  my @red1_integral = &reduce_integral(\@skip,$integral);
  my $red1_str = &print_integral(@red1_integral);
  my @red1_m1 = &inc_m(@red1_integral);
  my $m1_str = &print_integral(@red1_m1);

  my $str = "scalar_type $str_int =";

  if ($center_num == 1) {
    #$str = "$str PmA[i${level}] * $red1_str";
    $str = "$str PmA[$ind_val] * $red1_str";
  } elsif ($center_num == 2) {
    #$str = "$str PmB[i${level}] * $red1_str";
    $str = "$str PmB[$ind_val] * $red1_str";
  } elsif ($center_num == 3 && $INDICES != 3) {
    #$str = "$str QmC[i${level}] * $red1_str";
    $str = "$str QmC[$ind_val] * $red1_str";
  } elsif ($center_num == 4) {
    #$str = "$str QmD[i${level}] * $red1_str";
    $str = "$str QmD[$ind_val] * $red1_str";
  }

  if ($INDICES == 2) {
    #$str = "$str - PmC[i$level] * $m1_str;";
    $str = "$str - PmC[$ind_val] * $m1_str;";
  } elsif ($INDICES == 3 && $center_num == 3) {
    #$str = "$str WmQ[i$level] * $m1_str;";
    $str = "$str WmQ[$ind_val] * $m1_str;";
  } else { 
    if ($center_num <= 2) {
      #$str = "$str + WmP[i$level] * $m1_str;";
      $str = "$str + WmP[$ind_val] * $m1_str;";
    } else {
      #$str = "$str + WmQ[i$level] * $m1_str;";
      $str = "$str + WmQ[i$ind_val] * $m1_str;";
    }
  }

  return $str;
}

sub print_two_ind_rule {
  my ($integral,$ind1,$ind2) = @_;

  my $center1 = $INDEX_MAP[$ind1 - 1];
  my $center2 = $INDEX_MAP[$ind2 - 1];
  my $str_int = &print_integral(@$integral);

  my @skip = ($ind1,$ind2);
  my @red2_integral = &reduce_integral(\@skip,$integral);
  my $red2_str = &print_integral(@red2_integral);
  my @red2_m1 = &inc_m(@red2_integral);
  my $m1_str = &print_integral(@red2_m1);

  my $str = "$str_int +=";# del_$ind2$ind1 *";
  if ($center1 <= 2 and $center2 <= 2) {
    $str = "$str inv_two_zeta * ($red2_str -";
    if ($INDICES > 2) {
      $str = "$str rho_zeta *";
    }
    $str = "$str $m1_str);";
  } elsif ($center1 > 2 and $center2 > 2) {
    $str = "$str inv_two_eta * ($red2_str - rho_eta * $m1_str);";
  } else {
    $str = "$str inv_two_zeta_eta * $m1_str;";
  }

  return $str;
}

sub print_preterm {
  my ($level,$indent) = @_;

  my $str = "${indent}scalar_type preterm =";
  my $first = 1;
  for (1..$level) {
    if ($NORM_INDICES[$_-1]) {
      if ($first) {
        $str = "$str norm$_";
        $first = 0;
      } else {
        $str = "$str * norm$_";
      }
    }
#    if ($SAME_FUNC_INDICES[$_-1]) {
#      if ($first) {
#        $str = "$str !skip$_";
#        $first = 0;
#      } else {
#        $str = "$str * !skip$_";
#      }
#    }
  }
  $str = "$str;\n";

  return $str;
}

#sub print_skip {
#  my ($level,@ind_vals) = @_;

#  my $str = "bool skip$level =";
  #$str = "$str same_func &&";
#  if ($IND_TO_L[$level-1] == 1) {
    #my $lm1 = $level - 1;
#    if ($ind_vals[$level-1] > $ind_vals[$level-2]) {
      #$str = "$str (i$level > i$lm1)";
#      $str = "$str same_func";
#    } else {
#      $str = "$str false";
#    }
#  } elsif ($IND_TO_L[$level-1] == 2) {
    #my $lm1 = $level - 1;
    #my $lm2 = $level - 2;
    #my $lm3 = $level - 3;
#    if ($ind_vals[$level-2] > $ind_vals[$level-4] or ($ind_vals[$level-2] == $ind_vals[$level-4] and $ind_vals[$level-1] > $ind_vals[$level-3])) {
      #$str = "$str (i$lm1 > i$lm3 || (i$lm1 == i$lm3 && i$level > i$lm2))";
#      $str = "$str same_func";
#    } else {
#      $str = "$str false";
#    }
#  }

#  return $str;
#}

sub print_gradient_inner {
  my ($l,$indent,$ref,$str) = @_;
  my @first_req = @$ref;
  my $grad_l = $l + 1;
  $str = "$str$indent  //START INDEX i$grad_l, GRADIENT\n";
  $str = "$str$indent  for (int i$grad_l = 0; i$grad_l < 3; i$grad_l++)\n";
  $str = "$str$indent  \{\n";
  if ($INDICES == 2) {
    my $int_str = &print_integral(@first_req);
    my @m_p1 = &inc_m(@first_req);
    my $m1_str = &print_integral(@m_p1);
    $str = "$str$indent    scalar_type C_force_term = PmC[i$grad_l] * $m1_str;\n";
    for (1..$l) {
      my @skip = ($_);
      my @red_integral = &reduce_integral(\@skip,\@first_req);
      if (not &same_integral(\@red_integral,\@first_req)) {
        my @red_m_p1 = &inc_m(@red_integral);
        my $red1_str = &print_integral(@red_m_p1);
        $str = "$str$indent    C_force_term += (i$_ == i$grad_l) * inv_two_zeta * $red1_str;\n";
      }
    }
    $str = "$str$indent    scalar_type A_force_term = -C_force_term;\n";
    for (1..$l) {
      my @skip = ($_);
      my @red_integral = &reduce_integral(\@skip,\@first_req);
      if (not &same_integral(\@red_integral,\@first_req)) {
        my $red_str = &print_integral(@red_integral);
        $str = "$str$indent    A_force_term += (i$_ == i$grad_l) * inv_two_zeta * $red_str;\n";
      }
    }
    $str = "$str$indent    scalar_type B_force_term = PmB[i$grad_l] * $int_str + A_force_term;\n";
    $str = "$str$indent    A_force_term += PmA[i$grad_l] * $int_str;\n";
    $str = "$str$indent    A_force_term *= 2.0f * ai;\n";
    $str = "$str$indent    B_force_term *= 2.0f * aj;\n";
    for (1..$l) {
      my @skip = ($_);
      my @red_integral = &reduce_integral(\@skip,\@first_req);
      if (not &same_integral(\@red_integral,\@first_req)) {
        my $red_str = &print_integral(@red_integral);
        if ($INDEX_MAP[$_-1] == 1) {
          $str = "$str$indent    A_force_term -= (i$_ == i$grad_l) * $red_str;\n";
        } else {
          $str = "$str$indent    B_force_term -= (i$_ == i$grad_l) * $red_str;\n";
        }
      }
    }
    $str = "$str$indent    A_force[i$grad_l] += preterm * clatom_sh[j] * dens1[dens1_ind] * A_force_term;\n";
    $str = "$str$indent    B_force[i$grad_l] += preterm * clatom_sh[j] * dens1[dens1_ind] * B_force_term;\n";
    $str = "$str$indent    C_force[i$grad_l][tid] += preterm * clatom_sh[j] * dens1[dens1_ind] * C_force_term;\n";
  } else {
  }
  $str = "$str$indent  \}\n";
}
sub print_energy_inner {
  my ($str,$indent,$ref,$dens1_ind,$dens2_ind) = @_;
  my @first_req = @$ref;
  #$str = "$str$indent  my_fock[dens1_ind] +=";
  $str = "$str${indent}my_fock[$dens1_ind] += (double)(";
  $str = "$str preterm *" if $NORM_SAME_COUNT > 0;
  if ($INDICES == 2) {
    $str = "$str clatom_charge_sh[j] *";
  } else {
    #$str = "$str fit_dens_sh[j+dens2_ind] * prefactor_dens *";
    $str = "$str fit_dens_sh[j+$dens2_ind] * prefactor_dens * ";
  }
  my $int_str = &print_integral(@first_req);
  $str = "$str $int_str );\n";
}

#########################################################################################
#                               O-S EVALUATION
#########################################################################################

sub OS_level {
  #my ($level,$lower_str,$vals_ref,$dens1_ind,$dens2_ind,$reqs_ref) = @_;
  my ($level,$vals_ref,$dens1_ind,$dens2_ind,$reqs_ref) = @_;
  my @level_reqs = @$reqs_ref;
  my $indent = "  ";
  $indent = "$indent  " for (1..$level);
  my $str = "";
  my $lower_str = "";
  #return $lower_str if $level == 0;
  if ($level-1 == $l) {
    my $count = 0;
    my @first_req = ();
    for my $index (@BASE_INTEGRAL) {
      my @tmp = ();
      push @tmp, $index;
      for (0..$index-1) {
        $count++;
        push @tmp, $count;
      }
      push @first_req,\@tmp;
    }
    push @first_req, 0;
    push @$reqs_ref, \@first_req;
    if ($NORM_SAME_COUNT > 0) {
      my $preterm_str = &print_preterm($l,$indent);
      $str = "$str$preterm_str";
    }
    if ($GRADIENT) {
      my @m_p1 = &inc_m(@first_req);
      &add_integral($reqs_ref,\@m_p1);
      for (1..$l) {
        my @skip = ($_);
        my @red_integral = &reduce_integral(\@skip,\@first_req);
        my @red_m_p1 = &inc_m(@red_integral);
        &add_integral($reqs_ref,\@red_integral);
        &add_integral($reqs_ref,\@red_m_p1);
      }
      my $lower_str = &print_gradient_inner($l,$indent,\@first_req,$str);
      return "$str$lower_str";
    } else {
      my $lower_str = &print_energy_inner("",$indent,\@first_req,$$dens1_ind,$$dens2_ind);
      return "$str$lower_str";
    }
  }

  my @ind_vals = @$vals_ref;

  my $left_index = $INDEX_MAP[$level-1] <= 2;

  my $imax = 2;
  if ($LOCAL_LEVEL[$level-1]!=0) {
    #$str = "$str${indent}for (int i$level = 0; i$level < 3; i$level++)\n";
    $imax = $ind_vals[$level-2];
  }# else {
  #  my $last_level = $level-1;
  #  $str = "$str${indent}for (int i$level = 0; i$level <= i$last_level; i$level++)\n";
  #}
  if ($level==1 and $INDEX_MAP[$level-1] >= 3) {
    $$dens2_ind = 0;
  }
  for my $this_ind (0..$imax) {

    $ind_vals[$level-1] = $this_ind;
    #$str = &OS_level($level+1,$str,\@ind_vals,$pass_dens1_ind,$pass_dens2_ind,\@new_level_reqs);
    my @new_level_reqs = ();
    my @calc = ();
    if (($INDEX_MAP[$level-1] <= 2 and ($level < $l and $INDEX_MAP[$level] >= 3))) {
      #$str = "$str${indent}uint dens2_ind = 0;\n" if $INDICES > 2;
      $$dens2_ind = 0;
    }
    $lower_str = &OS_level($level+1,\@ind_vals,$dens1_ind,$dens2_ind,\@new_level_reqs);
    $$dens1_ind++ if $OUTER_DENS_UPDATE[$level-1];
    $$dens2_ind++ if $INNER_DENS_UPDATE[$level-1];

    for my $ref (@new_level_reqs) {
      my @integral = @$ref;
      my @skip = ($level);
      my @red1_integral = &reduce_integral(\@skip,\@integral);
      if (not &same_integral(\@red1_integral,\@integral)) {
        push @calc, \@integral;
        if ($INDICES != 3 or $left_index) {
          &add_integral($reqs_ref,\@red1_integral);
        }
        my @m_plus1 = &inc_m(@red1_integral);
        &add_integral($reqs_ref,\@m_plus1);
        for (1..$level-1) {
          if ($ind_vals[$level-1] == $ind_vals[$_-1]) {
            my $this_left = $INDEX_MAP[$_-1] <= 2;
            @skip = ($_,$level);
            my @red2_integral = &reduce_integral(\@skip,\@integral);
            if ($this_left == $left_index) {
              &add_integral($reqs_ref,\@red2_integral);
            }
            my @red2_m_plus1 = &inc_m(@red2_integral);
            &add_integral($reqs_ref,\@red2_m_plus1);
          }
        }
      } else {
        &add_integral($reqs_ref,$ref);
      }
    }

    $str = "$str${indent}//START INDEX i$level=$this_ind, CENTER $INDEX_MAP[$level-1]\n";
    $str = "$str${indent}\{\n";
    for my $integral (@calc) {
      my $one_ind_rule = &print_one_ind_rule($integral,$level,$this_ind);
      $str = "$str$indent  $one_ind_rule\n";
    }
    for my $ind2 (1..$level-1) {
      $str = "$str$indent  scalar_type norm$level = 1.0f;\n" if $NORM_INDICES[$level-1] and $ind2 == $level-1;
      #$str = "$str$indent    bool del_${ind2}${level} = i$ind2 == i$level;\n";
      if ($ind_vals[$level-1] == $ind_vals[$ind2-1]) {
        #$str = "$str$indent  \{\n";
        for my $integral(@calc) {
          my @skip = ($level);
          my @red1_integral = &reduce_integral(\@skip,$integral);
          @skip = ($ind2,$level);
          my @red2_integral = &reduce_integral(\@skip,$integral);
          if (not &same_integral(\@red1_integral,\@red2_integral)) {
            my $two_ind_rule = &print_two_ind_rule($integral,$level,$ind2);
            $str = "$str$indent  $two_ind_rule\n";
          }
        }
        #$str = "$str$indent    norm$level = del_${ind2}${level} * gpu_normalization_factor + !del_${ind2}${level} * 1.0f;\n" if $NORM_INDICES[$level-1] and $ind2 == $level-1;
        $str = "$str$indent  norm$level = gpu_normalization_factor;\n" if $NORM_INDICES[$level-1] and $ind2 == $level-1;
        #$str = "$str$indent  \}\n";
      }
    }
    $str = "$str$lower_str";
    $str = "$str${indent}\}\n";
  }
#  if ($SAME_FUNC_INDICES[$level-1]) {
#    my $skip_str = &print_skip($level,@ind_vals);
#    $str = "$str$indent  $skip_str;\n";
#  }
  #$str = "$str$lower_str";
#  if ($OUTER_DENS_UPDATE[$level-1]) {
#    if ($SAME_FUNC_INDICES[$level-1]) {
#      $str = "$str$indent  dens1_ind += !skip$level;\n";
#    } else {
      #$str = "$str$indent  dens1_ind++;\n";
#    }
#  }
  #$str = "$str$indent  dens2_ind++;\n" if $INNER_DENS_UPDATE[$level-1];

  return $str;
}

#########################################################################################
#                               MAIN PROGRAM ENTRY
#########################################################################################

($GRADIENT,@BASE_INTEGRAL) = @ARGV;

$INDICES = @BASE_INTEGRAL;
@INDEX_MAP = ();
@NORM_INDICES = ();

my $symm1 = ($BASE_INTEGRAL[0] > 0 and $BASE_INTEGRAL[0] == $BASE_INTEGRAL[1]);
my $symm2 = ($INDICES == 4 and $BASE_INTEGRAL[2] > 0 and $BASE_INTEGRAL[2] == $BASE_INTEGRAL[3]);
#@SAME_FUNC_INDICES = ();
$NORM_SAME_COUNT = 0;

@IND_TO_L = ();

@OUTER_DENS_UPDATE = ();
@INNER_DENS_UPDATE = ();

$l = 0;
$l += $_ for (@BASE_INTEGRAL);
@ind_vals_init = ();
push @ind_vals_init, 0 for (1..$l);
$F_l = $l + 1;
$gamma_l = $l;
if ($GRADIENT) {
    $F_l++;
    $gamma_l++;
}
if ($INDICES>2) {
print <<"END";
{
  scalar_type F_mT[$F_l];
  {
    scalar_type PmQ[3];
    PmQ[0] = P[0] - nuc_pos_dens_sh[j].x;
    PmQ[1] = P[1] - nuc_pos_dens_sh[j].y;
    PmQ[2] = P[2] - nuc_pos_dens_sh[j].z;
    scalar_type T = (PmQ[0] * PmQ[0] + PmQ[1] * PmQ[1] + PmQ[2] * PmQ[2]) * rho;
    lio_gamma<scalar_type,$gamma_l>(F_mT,T);
  }
  {
END
#    uint dens1_ind = 0;
  if (($INDICES == 3 and $BASE_INTEGRAL[2] == 0) or ($INDICES == 4 and ($BASE_INTEGRAL[2] == 0 and $BASE_INTEGRAL[3] == 0))) {
    #print "    uint dens2_ind = 0;\n";
  }
} else {
print <<"END";
{
  scalar_type F_mU[$F_l];
  {
    scalar_type U = (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
    // TODO (maybe): test out storing F(m,U) values in texture and doing a texture fetch here rather than the function calculation
    lio_gamma<scalar_type,$gamma_l>(F_mU,U);
  }
  {
END
#    uint dens1_ind = 0;
}

my @req = ();
my $count = 0;
my @first_req = ();
my $ind_count = 0;
@LOCAL_LEVEL = ();
for my $index (@BASE_INTEGRAL) {
  $ind_count += 1;
  #my @tmp = ();
  #push @tmp, $index;
  for (0..$index-1) {
    $count += 1;
    #push @tmp, $count;
    push @LOCAL_LEVEL, $_;
    push @INDEX_MAP, $ind_count;
    push @NORM_INDICES, 0;
#    push @SAME_FUNC_INDICES, 0;
    push @OUTER_DENS_UPDATE, 0;
    push @INNER_DENS_UPDATE, 0;
    push @IND_TO_L, $index;
  }
  $NORM_INDICES[$#NORM_INDICES] = 1 if $index == 2;
#  $SAME_FUNC_INDICES[$#SAME_FUNC_INDICES] = 1 if ($symm1 and $ind_count == 2) or ($symm2 and $ind_count == 4);
  $NORM_SAME_COUNT += 1 if $index == 2;# or ($symm1 and $ind_count == 2) or ($symm2 and $ind_count == 4);
  if ($index > 0) {
    $OUTER_DENS_UPDATE[$#OUTER_DENS_UPDATE] = 1 if $ind_count == 2;
    $INNER_DENS_UPDATE[$#INNER_DENS_UPDATE] = 1 if ($INDICES == 3 and $ind_count == 3) or $ind_count == 4;
  } elsif ($ind_count == 2 and $BASE_INTEGRAL[0] > 0) {
    $OUTER_DENS_UPDATE[$#OUTER_DENS_UPDATE] = 1;
  } elsif ($ind_count == 4 and $BASE_INTEGRAL[2] > 0) {
    $INNER_DENS_UPDATE[$#INNER_DENS_UPDATE] = 1;
  }
  #push @first_req, \@tmp;
}
#push @first_req, 0;
#push @req, \@first_req;

my $indent = "  ";
$indent = "$indent  " for (1..$l);
my $str = "";
#if ($NORM_SAME_COUNT > 0) {
#  my $preterm_str = &print_preterm($l,$indent);
#  $str = "$str$preterm_str";
#}
#if ($gradient) {
#  my @m_p1 = &inc_m(@first_req);
#  &add_BASE_INTEGRAL(\@req,\@m_p1);
#  for (1..$l) {
#    my @skip = ($_);
#    my @red_BASE_INTEGRAL = &reduce_BASE_INTEGRAL(\@skip,\@first_req);
#    my @red_m_p1 = &inc_m(@red_BASE_INTEGRAL);
#    &add_BASE_INTEGRAL(\@req,\@red_BASE_INTEGRAL);
#    &add_BASE_INTEGRAL(\@req,\@red_m_p1);
#  }
  #$str = &print_gradient_inner($l,$indent,\@first_req,$str);
#}
#if (not $gradient) {
  #$str = &print_energy_inner($str,$indent,\@first_req);
#}

#print &OS_level(1,$str,\@ind_vals_init,0,0,\@req);
my $dens1_ind = 0;
my $dens2_ind = 0;
print &OS_level(1,\@ind_vals_init,\$dens1_ind,\$dens2_ind,\@req);

print "  }\n";
print "}\n";
