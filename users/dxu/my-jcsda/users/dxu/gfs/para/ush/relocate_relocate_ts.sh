
# This script attempts to perform tropical cyclone relocation
#
# It is executed by the script prepobs_makeprepbufr.sh
# ----------------------------------------------------

set -aux

# Positional parameters passed in:
#   1 - path to input tcvitals file from current time - this is the tcvitals
#       file normally in the appropriate /com directory which was generated
#       by the SYNDAT_QCTROPCY program
#   2 - center date/time for the PREPBUFR processing (YYYYMMDDHH)
#   3 - path to input tcvitals file 6-hours prior to the current time - this
#       is the tcvitals file normally in the appropriate /com directory which
#       was generated by the SYNDAT_QCTROPCY program 6-hours ago
#   4 - path to input tcvitals file 12-hours prior to the current time - this
#       is the tcvitals file normally in the appropriate /com directory which
#       was generated by the SYNDAT_QCTROPCY program 12-hours ago


# Imported variables that must be passed in:
#   DATA     - path to working directory
#   RUN      - model run
#   COMOUT   - output /com file location
#   pgmout   - string indicating path to for standard output file 
#   USHRELO  - path to RELOCATE ush files
#   FIXSYND  - path to synthethic data fixed field files
#   POE_OPTS - string indicating options to use with poe command

# Imported variables that can be passed in:
#   SENDCOM - if "YES" copies output tcvitals file to $COMOUT
#              defaults to "YES"
#   RELX  - path to RELOCATE_mv_nvortex program 
#            defaults to
#            "$EXECRELO/relocate_mv_nvortex" if not
#            passed in
#   EXECRELO - path to RELOCATE executables (skipped over by this script if
#              $RELX is passed in)
#   jlogfile - string indicating path to joblog file (skipped over by this
#              script if not passed in)

# The following files are expected to be in the working directory;
#  this script executes programs that read them:
#   $DATA/sgm3prep - global sigma GUESS valid for 3 hrs prior to the current
#                     time
#   $DATA/sgesprep - global sigma GUESS valid for the current time
#   $DATA/sgp3prep - global sigma GUESS valid for 3 hrs after the current time

# Files generated by this script:
#   $DATA/tcvitals - updated tcvitals file based on the relocated tropical
#                     cyclone (this is input to subsequent SYNDAT_SYNDATA
#                     program - if this file is empty then SYNDAT_SYNDATA
#                     will not run)

cd $DATA

SENDCOM=${SENDCOM:-YES}

tcvitals_now=$1
CDATE10=$2
tcvitals_m6=$3
tcvitals_m12=$4

if [ "$tcvitals_m12" != "" ]
then
   cat $tcvitals_m12 > VITL
fi
if [ "$tcvitals_m6" != "" ]
then
   cat $tcvitals_m6 >> VITL
fi
if [ "$tcvitals_now" != "" ]
then
   cat $tcvitals_now >> VITL
fi

NDATE=$NWPROD/util/exec/ndate
MP_PULSE=0
MP_TIMEOUT=600
pdy=$(echo $CDATE10|cut -c1-8)
cyc=$(echo $CDATE10|cut -c9-10)
cycle=t$(echo $CDATE10|cut -c9-10)z
GDATE10=$($NDATE -06 $CDATE10)
TIMEIT=""
[ -s $DATA/timex ] && TIMEIT=$DATA/timex

#  make unique combined tcvitals file for t-12, t-6 and t+0 -- 
#  if tcvitals does not contains record from current time,
#  skip relocation processing
#  ---------------------------------------------------------------

grep "$pdy $cyc" VITL
errgrep=$?
> tcvitals
if [ $errgrep -ne 0 ] ; then
   msg="NO TCVITAL RECORDS FOUND FOR $CDATE10 - EXIT TROPICAL CYCLONE \
RELOCATION PROCESSING"
   set +u
   [ -n "$jlogfile" ] && postmsg "$jlogfile" "$msg"
   set -u
else

   cat VITL >>tcvitals
   grep "$pdy $cyc" VITL > tcvitals.now1 


#  create model forecast track location file
# ------------------------------------------
#               $DATA/$RUN.$cycle.relocate.model_track.tm00

   $TIMEIT $USHRELO/relocate_extrkr.sh
   err=$?
   if [ $err -ne 0 ]; then

#  problem: script relocate_extrkr.sh failed
#  -----------------------------------------

      set +x
      echo
      echo "$USHRELO/relocate_extrkr.sh failed"
      echo "ABNORMAL EXIT!!!!!!!!!!!"
      echo
      set -x
      if [ -s err_exit ]; then
         err_exit "Script $USHRELO/relocate_extrkr.sh failed"
      else
#########kill -9 ${qid}
         exit 555
      fi
      exit 9
   fi

#  relocate model tropical storm vortices in ges sigma files
# ----------------------------------------------------------

#  READ DATA fort.11 fort.60
#  DATA OUT only from _t09 fort.55

##VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV##
##                          HEREFILE RELOCATE_GES                            ##
##VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV##

cat <<\EOF >RELOCATE_GES

{ echo

set -aux


DATA=$1
gesfhr=$2
multi=$3

status=$DATA/mstatus_rel${multi} ; > $status
mp_pgmout=$DATA/mp_rel_pgmout${multi}  ; > $mp_pgmout

{ echo
set +x
echo
echo "********************************************************************"
echo This is task $multi executing on node  `hostname -s`
echo Guess forecast hour is $gesfhr
echo Starting time: `date`
echo "********************************************************************"
echo
set -x
} >> $mp_pgmout

if [ $gesfhr = 03 ] ; then sges=sgm3prep ;fi
if [ $gesfhr = 06 ] ; then sges=sgesprep ;fi
if [ $gesfhr = 09 ] ; then sges=sgp3prep ;fi

cd $DATA

JCAP=`$NWPROD/exec/global_sighdr $sges jcap`
LEVS=`$NWPROD/exec/global_sighdr $sges levs`

RELX=${RELX:-$EXECRELO/relocate_mv_nvortex}

echo "$RELX for gesfhr=${gesfhr}" > relocate_exec_name${multi}

echo "$gesfhr" > gesfhr${multi}

pgm=`basename  $RELX`
#-----mimics prep_step-----
set +x
echo $pgm > pgmname
set +u
[ -z "$mp_pgmout" ] && echo "Variable mp_pgmout not set"
set -u
[ -s $DATA/break ] && paste pgmname $DATA/break >> $mp_pgmout
rm pgmname
[ -f errfile ] && rm errfile
export XLFUNITS=0
unset `env | grep XLFUNIT | awk -F= '{print $1}'`

set +u
if [ -z "$XLFRTEOPTS" ]; then
  export XLFRTEOPTS="unit_vars=yes"
else
  export XLFRTEOPTS="${XLFRTEOPTS}:unit_vars=yes"
fi
set -u

[ -s $DATA/tracer ] && cat $DATA/tracer > errfile
set -x
#--------------------------

rm -f fort.12
ln -s $FIXSYND/global_slmask.t126.grb fort.12
export XLFUNIT_11=tcvitals.now1
export XLFUNIT_20=sgm3prep
export XLFUNIT_21=sgesprep
export XLFUNIT_22=sgp3prep
export XLFUNIT_30=$DATA/model_track.all

export XLFUNIT_52=rel_inform
export XLFUNIT_53=sgm3prep.relocate
export XLFUNIT_55=tcvitals.relocate
export XLFUNIT_56=sgesprep.relocate
export XLFUNIT_59=sgp3prep.relocate

echo $gesfhr $LONA $LATA | $TIMEIT $RELX > relocate.stdout.$gesfhr.$cycle 2>&1
cat relocate.stdout.$gesfhr.$cycle >> $mp_pgmout
errcat=$?
set +x
echo
#echo "The foreground exit status for $RELX is " $errcat
echo "The foreground exit status for $pgm is " $errcat
echo
set -x

##echo "$multi finished -- err$RELX = $errcat" > $status
echo "$multi finished -- err$pgm = $errcat" > $status

{ echo
set +x
echo
echo "********************************************************************"
echo Finished executing on node  `hostname -s`
echo Guess forecast hour is $gesfhr
echo Ending time  : `date`
echo "********************************************************************"
echo
set -x
} >> $mp_pgmout

} 2> $DATA/mp_rel_poe${3}.errfile

exit 0
EOF
chmod 775 RELOCATE_GES

##AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA##
##                       end of HEREFILE RELOCATE_GES                        ##
##AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA##

   nprocs=$(echo $LOADL_PROCESSOR_LIST|wc -w)
   >cmd
   cmd[0]="RELOCATE_GES $DATA 03 \$n"
   cmd[1]="RELOCATE_GES $DATA 06 \$n"
   cmd[2]="RELOCATE_GES $DATA 09 \$n"
   m=-1
   n=-1
   while [ $((n+=1)) -le $nprocs ] ;do
      while [ $((m+=1)) -le 2 ] ;do
         eval echo ${cmd[m]} | tee -a cmd
         ((n+=1))
         echo "echo do-nothing" >>cmd
         ((n+=1))
      done
     echo "echo do-nothing" >>cmd
   done

# Effective with transition to frost/snow, DO NOT execute a time command
#  with poe command!!!
#########/usr/bin/timex /usr/bin/poe -cmdfile cmd $POE_OPTS
   /usr/bin/poe -cmdfile cmd $POE_OPTS

   errSTATUS=0
   n=-1
   while [ $((n+=1)) -le $nprocs ] ;do
      if [ -s $DATA/mp_rel_pgmout${n} ]; then
         cat $DATA/mp_rel_pgmout${n} >> relocate.out
         cat $DATA/mp_rel_pgmout${n} >> $pgmout

#  check for success
# ------------------

         status=$DATA/mstatus_rel${n}
         gesfhr=`cat gesfhr${n}`
         if [ ! -s $status ]; then
            set +x
            echo
   echo "********************************************************************"
   echo "                   P  R  O  B  L  E  M   !   !   !                  "
   echo "********************************************************************"
   echo " ###> RELOCATE_GES Stream (Task) $n FAILED - Cycle date: $CDATE10"
   echo "       Guess forecast hour is $gesfhr "
   echo "       Current working directory: $DATA                             "
   echo "********************************************************************"
            echo
            set -x
            errSTATUS=99
         fi
         set +x
         echo
   echo "********************************************************************"
   echo "    ++  Script trace from RELOCATE_GES for Stream (Task) $n        ++"
   echo "                 Guess forecast hour is $gesfhr "
   echo "********************************************************************"
         echo

         cat mp_rel_poe${n}.errfile

         echo
   echo "********************************************************************"
   echo " ++  End of Script trace from RELOCATE_GES for Stream (Task) $n    ++"
   echo "                 Guess forecast hour is $gesfhr "
   echo "********************************************************************"
         echo
         set -x
         if [ "$errSTATUS" -gt '0' ]; then
            if [ -s err_exit ]; then
               err_exit "Script RELOCATE_GES (herefile) failed"
            else
###############kill -9 ${qid}
               exit 555
            fi
            exit 9
         fi
         err=`cut -f 2 -d = $status`
         RELX=`cat relocate_exec_name${n}`
         pgm=`basename  "$RELX"`
         touch errfile
         if [ -s err_chk ]; then
            err_chk
         else
            if test "$err" -gt '0'
            then
###############kill -9 ${qid}
               exit 555
            fi
         fi
         [ "$err" -gt '0' ]  &&  exit 9
      fi
   done

#  further check for success
#  -------------------------

   for sges in sgm3prep sgesprep sgp3prep; do
      if [ -s $sges.relocate ] ; then
         mv $sges.relocate $sges
      else

#  problem: $sges.relocate does not exist
#  --------------------------------------

         if [ -s err_exit ]; then
            err_exit "The file $sges.relocate does not exist"
         else
############kill -9 ${qid}
            exit 555
         fi
         exit 9
      fi
   done

   if [ -s tcvitals.relocate ]; then
      mv tcvitals.relocate tcvitals
   else
      >tcvitals
   fi
   rm RELOCATE_GES cmd

   if [ "$SENDCOM" = "YES" ]
   then
      cp rel_inform ${COMOUT}/${RUN}.${cycle}.inform.relocate
      cp tcvitals ${COMOUT}/${RUN}.${cycle}.tcvitals.relocate
   fi

   msg="TROPICAL CYCLONE RELOCATION PROCESSING SUCCESSFULLY COMPLETED FOR \
$CDATE10"
   set +u
   [ -n "$jlogfile" ] && postmsg "$jlogfile" "$msg"
   set -u

# end GFDL ges manipulation
#------------------------------------------------------------

fi

exit 0
