#!/bin/bash
#
# Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
# Ensemble Kalman Filter for SHYFEM. 
# 
# 2016 first version
# 2018 important updates
# 2024 other changes
#
# See the README file
# 
#----------------------------------------------------------

# This finds the path of the current script
SCRIPT=$(realpath $0)
SCRIPTPATH=$(dirname $SCRIPT)

FEMDIR=$SCRIPTPATH/..	# fem directory
SIMDIR=$(pwd)		# current dir


#----------------------------------------------------------

Usage()
{
  echo "Usage: enKF.sh [method] [localisation] [bas-file] [nlv] [n] [out]"
  echo
  echo "method = 11|12|13|21|22|23. See: analysis.F90"
  echo "localisation = 0 (disable), 1 (enable)"
  echo "bas-file = name of the basin bas file"
  echo "nlv = number of vertical levels used in the simulations"
  echo "n = number of threads"
  echo "out = 0 saves only mean and std restarts, 1 save all the restarts (needed by enKS)"
  echo
  exit 0
}

#----------------------------------------------------------

Check_file()
{
  if [ ! -s $1 ]; then
     echo "File $1 does not exist or has zero size."
     exit 1
  fi
}

#----------------------------------------------------------

Check_num()
{
  nint='^[0-9]+$'
  nreal='^-?[0-9]+([.][0-9]+)?$'

  if [ $3 = 'int' ]; then
     if ! [[  $4 =~ $nint ]]; then
        echo "$4 is not an integer number"
        exit 1
     fi
  elif [ $3 = 'real' ]; then
     if ! [[ $4 =~ $nreal ]]; then
        echo "$4 is not a real number"
        exit 1
     fi
  fi
  if [[ $(echo "$4 < $1" | bc) = 1 ]] || [[ $(echo "$4 > $2" | bc) = 1 ]]; then
     echo "Number $4 out of range"
     exit 1
  fi 
}

#----------------------------------------------------------

Check_files(){
  echo "Check the exec programs"
  command -v parallel > /dev/null 2>&1 || { echo "parallel it's not installed.  Aborting." >&2; exit 1; }
  [ ! -s $FEMDIR/fem3d/shyfem ] && echo "shyfem exec does not exist. Compile the model first." && exit 1
  # Make here the mod_dimensions and compile main
  
  [ ! -s $FEMDIR/enKF/main ] && echo "main exec does not exist. Compile the enKF first." && exit 1

  echo "Check the input files"
  ens_file_list='ens_list.txt'
  Check_file $ens_file_list  
  obs_file_list='obs_list.txt'
  Check_file $obs_file_list
  an_time_list='antime_list.txt'
  Check_file $an_time_list
}

#----------------------------------------------------------

Read_ens_list(){
# Reads the list of skel and restart files of the ensemble
  echo "Reading the ensemble list"

  rm -f an00001_en*b.rst

  rst1=$(head -1 $ens_file_list | awk '{print $2}')
  rst2=$(head -2 $ens_file_list | tail -1 | awk '{print $2}')

  if [ "$rst1" = "$rst2" ];  then
	  is_new_ens=1
	  echo "Only one initial state..."
  else
	  is_new_ens=0
	  echo "Many initial states..."
  fi

  nrow=0
  while read line
  do
     skelf=$(echo $line | awk '{print $1}')
     rstf=$(echo $line | awk '{print $2}')
     Check_file $skelf
     Check_file $rstf

     skel_file[$nrow]=$skelf

     if [ "$is_new_ens" -eq "0" ] || [ "$nrow" -eq "0" ]; then
        nel=$(printf "%05d" $nrow)
        ln -fs $rstf an00001_en${nel}b.rst
     fi

     nrow=$((nrow + 1))

  done < $ens_file_list

  nrens=$nrow
  echo ""; echo "Number of ensemble members: $nrens"; echo ""
}

#----------------------------------------------------------

Read_an_time_list(){
# Reads the list of the analysis times and determines
# the number of analysis steps (nran)
  echo "Read the list of analysis steps"
  
  # timeo starts from 1 as the analysis steps
  nrow=0
  nran=0
  while read line
  do
    nrow=$((nrow + 1))
    iseven=$((nrow%2))
    if [ "$iseven" -eq "1" ]; then
	    nran=$((nran+1))
	    timeo[$nran]=$(echo $line | awk '{print $1}')
	    nfile[$nran]=$(echo $line | awk '{print $2}')
    else
	    isfile[$nran]=$line
    fi
  done < $an_time_list
  echo ""; echo "Number of analysis steps: $nran"
}


#----------------------------------------------------------

SkelStr(){
# Makes a str file from a skel file
inamesim=$1; iitanf=$2; iitend=$3; irestrt=$4; iskelname=$5; istrname=$6

if [ ! -s $iskelname ]; then
        echo "File $iskelname does not exist"
        exit 1
fi

cat $iskelname | sed -e "s/NAMESIM/$inamesim/g" |  sed -e "s/ITANF/$iitanf/g" \
               | sed -e "s/ITEND/$iitend/g" | sed -e "s/RESTRT/$irestrt/" \
               | sed -e "s/IDTRST/-1/" \
               >  $istrname
}

#----------------------------------------------------------

Make_sim()
{
  basen=$(basename $1 .str)
  $2/shyfem $1 > $basen.log 
}

#----------------------------------------------------------

Write_obs_file(){
# Write a tmp file with the observations used in the current time step

  na=$1

  IFS=' ' read -r -a nisfile <<< "${isfile[$na]}"
  nnfile=${nfile[$na]}

  rm -f obs_list_tmp.txt
  k=0
  while read line; do
    if [ "${nisfile[$k]}" = "1" ]; then
	    echo $line >> obs_list_tmp.txt
    fi
    k=$((k+1))
  done < $obs_file_list

  if [ "$nnfile" -ne "$k" ]; then
	  echo "Error in the length of the obs file list: $nnfile $k"
	  exit 1
  fi
}

#----------------------------------------------------------

Write_info_file(){
# Write a file with informations for the fortran analysis program

  na=$1

  echo $nnlv > analysis.info		# nr of vertical levels
  echo $nrens >> analysis.info		# nr of ens members
  echo $na >> analysis.info		# analysis step
  echo $bas_file >> analysis.info	# name of the basin
  echo ${timeo[$na]} >> analysis.info	# current time
  echo obs_list_tmp.txt >> analysis.info	# obs file list
  echo $is_new_ens >> analysis.info	# if to make a new ens of states
  echo $rmode >> analysis.info          # analysis method
  echo $islocal >> analysis.info        # local analysis
}

#----------------------------------------------------------

Run_ensemble_analysis()
{
nanl=$(printf "%05d" $1)

cd $SIMDIR
$FEMDIR/enKF/main
if [ "$?" -ne "0" ]; then
          echo "Errors while running main."
          exit 1
fi

# Check restart files
for (( ne = 0; ne < $nrens; ne++ )); do
	nensl=$(printf "%05d" $nens)
	filename="an${nanl}_en${nensl}a.rst"
	Check_file $filename
done

#filename="an${nanl}_mean_state_b.rst"	#Average
#Check_file $filename
#filename="an${nanl}_mean_state_a.rst"	#Average
#Check_file $filename
}

#----------------------------------------------------------

Make_ens_str(){
strfiles=""
for (( ne = 0; ne < $nrens; ne++ )); do

   if [ "$is_new_ens" -eq "0" ]; then
      ens_skel_file=${skel_file[$ne]}
   else
      ens_skel_file=${skel_file[0]}
   fi
   Check_file $ens_skel_file

   nel=$(printf "%05d" $ne); nal=$(printf "%05d" $na)
   naa=$((na + 1)); naal=$(printf "%05d" $naa)

   itanf=${timeo[$na]}

   if [ "$na" -ne "$nran" ]; then
        name_sim="an${naal}_en${nel}b"
	itend=${timeo[$naa]}
        rstfile="an${nal}_en${nel}a.rst" 
        strnew="${name_sim}.str"

        SkelStr $name_sim $itanf $itend $rstfile $ens_skel_file $strnew
        strfiles="$strfiles $strnew"
   fi

done
}


#----------------------------------------------------------
#----------------------------------------------------------
#	MAIN
#----------------------------------------------------------
#----------------------------------------------------------

if [ $6 ]; then
   rmode=$1
   islocal=$2
   bas_file=$3
   nnlv=$4
   nthreads=$5
   out_verb=$6
else
   Usage
fi

# Export num of threads for main. Warning! Use the max num of cpu, not threads.
export OMP_NUM_THREADS=$nthreads

# Checking the executable programs
Check_files

# Reading skel file list
Read_ens_list

# Reading obs file list
Read_an_time_list

# Assimilation cycle for every analysis time step
rm -f X5*.uf backKF_*.rst analKF_*.rst 	# old files
for (( na = 1; na <= $nran; na++ )); do

   # make the analysis
   echo; echo "			ANALYSIS STEP $na OF $nran"; echo
   Write_obs_file $na

   Write_info_file $na

   ######################### ANALYSIS ######################
   Run_ensemble_analysis $na
   #########################################################

   # Makes nrens str files for the simulations
   Make_ens_str

   if [ "$na" -ne "$nran" ]; then # not the last one

      # run nrens sims before the obs
      echo; echo "       running $nrens ensemble simulations..."

   ############### MODEL RUN ###############################
      # with nthreads=0 uses the maximum number
      nthsim=$nthreads
      [[ "$nthsim" -gt "$nrens" ]] && nthsim=$nrens
      export -f Make_sim
      parallel --no-notice -P $nthsim Make_sim ::: $strfiles ::: $FEMDIR/fem3d
   #########################################################

   fi

   # merge the rst files
   nanl=$(printf "%05d" $na)
   for (( ne = 0; ne < $nrens; ne++ )); do
        nel=$(printf "%05d" $ne)
        filename1="an${nanl}_en${nel}b.rst"
        filename2="an${nanl}_en${nel}a.rst"
        Check_file $filename1
        Check_file $filename2
	# If not verbose, saves only the last rst for each member
	if [ "$out_verb" -eq "1" ]; then
		[[ "$na" -gt "1" ]] && cat $filename1 >> backKF_en$nel.rst
		cat $filename2 >> analKF_en$nel.rst
	else
	        [[ "$na" -eq "$nran" ]] && mv -f $filename2 analKF_en$nel.rst
	fi
	rm -f $filename1 $filename2
   done
   filename1="an${nanl}_mean_a.rst"
   filename2="an${nanl}_std_a.rst"
   Check_file $filename1
   Check_file $filename2
   cat $filename1 >> analKF_mean.rst
   cat $filename2 >> analKF_std.rst
   rm -f $filename1 $filename2
   rm -f an*_en*b.inf an*_en*.log an*_en*b.str 

done
rm -f X5col.dat X5row.dat X5.uf make.log

exit 0
