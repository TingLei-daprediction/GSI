#! /bin/sh
##SBATCH --account=fv3-cam
#SBATCH --account=fv3-cpu
#SBATCH --qos=debug
##SBATCH --qos=batch
#SBATCH --ntasks=600
#SBATCH -t 00:15:00
#SBATCH --job-name=RRFS_GSI
#SBATCH -o RRFS_gsi.log
##SBATCH --cpus-per-task 2 --exclusive

set -x 

THREAD=yes
export MPICH_ALLTOALL_THROTTLE=0
export MP_COLLECTIVE_OFFLOAD=no
export KMP_STACKSIZE=1024m

#export MP_TASK_AFFINITY=core:2
export OMP_NUM_THREADS=1
. /apps/lmod/lmod/init/sh
#GSISRCDIR="/scratch2/NCEPDEV/fv3-cam/Ting.Lei/dr-misha/dr-emc-gsi/GSI/"
#GSISRCDIR="/scratch1/NCEPDEV/da/Miodrag.Rancic/RRFS-GSI/GSI"

#moddir=${GSISRCDIR}/modulefiles
#modfile=${moddir}/gsi_hera.intel.lua
#modfile=${moddir}/gsi_hera.intel.lua
#gsiexec=$GSISRCDIR/build/src/gsi/gsi.x

module use /scratch1/NCEPDEV/da/Miodrag.Rancic/Modules/modulefiles
module load load gsi_hera.intel

gsiexec=/scratch1/NCEPDEV/da/Miodrag.Rancic/Ting/GSI/build/src/gsi/gsi.x

#module load $modfile

#source /scratch2/NCEPDEV/fv3-cam/Ting.Lei/dr-misha/dr-MGBF-from-wcoss/multigridwork/rtma_gsi_fd/ProdGSI/modulefiles/modulefile.ProdGSI.hera
#source /scratch1/NCEPDEV/da/Miodrag.Rancic/RRFS-GSI/GSI/modulefiles/modulefile.ProdGSI.hera

###################### Start Module load ##################################################
##module load udunits/2.2.25
##module load HDF5-serial/1.10.1
######################   End Module load ##################################################
###################### START MPONDECA EDITS ###############################################
CP=/bin/cp
RUN=RRFS-3km_conus
MCDATE=202202020400
YYYYMMDDHHMU=202202020400
workdir=/scratch1/NCEPDEV/stmp4/${USER}/RRFS_MGBF_orig/tmpreg/dr-$MCDATE
#workdir=/scratch1/NCEPDEV/stmp4/${USER}/RRFS_MGBF/tmpreg/dr-$MCDATE
rm -fr $workdir
mkdir -p /scratch1/NCEPDEV/stmp4/${USER}/RRFS_MGBF_orig/tmpreg/dr-$MCDATE
#mkdir -p /scratch1/NCEPDEV/stmp4/${USER}/RRFS_MGBF/tmpreg/dr-$MCDATE
#bg_and_obs_dir=/scratch2/NCEPDEV/stmp3/Ting.Lei/ufs-weather-app-RRFS_B_exp0_det_dir/stmp/tmpnwprd/RRFS_conus_3km_for_MGBF/2022020204/anal_conv_gsi_spinup_test
bg_and_obs_dir=/scratch1/NCEPDEV/da/Miodrag.Rancic/RRFS-GSI/RRFS_conus_3km_for_MGBF/2022020204/anal_conv_gsi_spinup_test
cd $workdir
cp $bg_and_obs_dir/* .
cp fv3_dynvars.bg fv3_dynvars &
cp fv3_tracer.bg fv3_tracer & 
cp fv3_sfcdata.bg fv3_sfcdata
wait 
cp $gsiexec gsi.x
pgmout=stdout_gsi
#
# Copy here 568-572
#
#vdef_namelist=/scratch1/NCEPDEV/da/Miodrag.Rancic/Ting/GSI/aux/gsiparm.anl_mgbf
vdef_namelist=/scratch1/NCEPDEV/da/Miodrag.Rancic/Ting/GSI/aux_orig/gsiparm.anl
source $vdef_namelist
cat << EOF > gsiparm.anl
$gsi_namelist
EOF
#
#
#
srun --label gsi.x >$pgmout  2> errfile
export err=$?
export error=$err


export endianness=Big_Endian

#===========================================================#
# error checking used in GSD script
#===========================================================#
export pgmout_stdout="stdout"
cat ${pgmout} > ${pgmout_stdout}
if [ -f errfile ] ; then
    cat errfile >> ${pgmout_stdout}
fi

if [ ${error} -ne 0 ]; then
  ${ECHO} "ERROR: ${GSI} crashed  Exit status=${error}"
  cp -p ${pgmout_stdout}  ../.
  exit ${error}
fi

ls -l > GSI_workdir_list

# Look for successful completion messages in rsl files
nsuccess=`tail -200 ${pgmout_stdout} | awk '/PROGRAM GSI_ANL HAS ENDED/' | wc -l`
ntotal=1 
echo  "Found ${nsuccess} of ${ntotal} completion messages"
if [ ${nsuccess} -ne ${ntotal} ]; then
   ${ECHO} "ERROR: ${GSI} did not complete sucessfully  Exit status=${error}"
   cp -p ${pgmout_stdout}  ../.
   cp GSI_workdir_list ../.
   if [ ${error} -ne 0 ]; then
     exit ${error}
   else
     exit 1
   fi
fi

# Loop over first and last outer loops to generate innovation
# diagnostic files for indicated observation types (groups)
#
# NOTE:  Since we set miter=2 in GSI namelist SETUP, outer
#        loop 03 will contain innovations with respect to 
#        the analysis.  Creation of o-a innovation files
#        is triggered by write_diag(3)=.true.  The setting
#        write_diag(1)=.true. turns on creation of o-g
#        innovation files.
#

loops="01 02"
for loop in $loops; do

case $loop in
  01) string=ges;;
  02) string=anl;;
   *) string=$loop;;
esac

#  Collect diagnostic files for obs types (groups) below
#  listall="hirs2_n14 msu_n14 sndr_g08 sndr_g11 sndr_g11 sndr_g12 sndr_g13 sndr_g08_prep sndr_g11_prep sndr_g12_prep sndr_g13_prep sndrd1_g11 sndrd2_g11 sndrd3_g11 sndrd4_g11 sndrd1_g12 sndrd2_g12 sndrd3_g12 sndrd4_g12 sndrd1_g13 sndrd2_g13 sndrd3_g13 sndrd4_g13 hirs3_n15 hirs3_n16 hirs3_n17 amsua_n15 amsua_n16 amsua_n17 amsub_n15 amsub_n16 amsub_n17 hsb_aqua airs_aqua amsua_aqua imgr_g08 imgr_g11 imgr_g12 pcp_ssmi_dmsp pcp_tmi_trmm conv sbuv2_n16 sbuv2_n17 sbuv2_n18 omi_aura ssmi_f13 ssmi_f14 ssmi_f15 hirs4_n18 hirs4_metop-a amsua_n18 amsua_metop-a mhs_n18 mhs_metop-a amsre_low_aqua amsre_mid_aqua amsre_hig_aqua ssmis_las_f16 ssmis_uas_f16 ssmis_img_f16 ssmis_env_f16 iasi_metop-a"
   listall="conv"
   for type in $listall; do
      count=`ls pe*.${type}_${loop}* | wc -l`
      if [[ $count -gt 0 ]]; then
         `cat pe*.${type}_${loop}* > diag_${type}_${string}.${YYYYMMDDHHMU}`
      fi
   done
done

# save results from 1st run
${CP} fort.201    fit_p1.${YYYYMMDDHHMU}
${CP} fort.202    fit_w1.${YYYYMMDDHHMU}
${CP} fort.203    fit_t1.${YYYYMMDDHHMU}
${CP} fort.204    fit_q1.${YYYYMMDDHHMU}
${CP} fort.207    fit_rad1.${YYYYMMDDHHMU}
cat   fort.* >    fits_${YYYYMMDDHHMU}.txt
${CP} -p fort.220 minimization_fort220.${YYYYMMDDHHMU}
# cat fort.* > ${COMOUT}/fits_${YYYYMMDDHH}.txt


 #Comment out / MPondeca
# Saving ANALYSIS, DIAG, Obs-Fitting files TO COM2 DIRECTORY AS PRODUCT for archive
#${CP} -p ${DATA}/wrf_inout                  ${COMOUTgsi_rtma3d}/${ANLrtma3d_FNAME}
#${CP} -p ${pgmout_stdout}                   ${COMOUTgsi_rtma3d}/${pgmout_stdout}_gsianl.${YYYYMMDDHHMU}
#${CP} -p fits_${YYYYMMDDHHMU}.txt             ${COMOUTgsi_rtma3d}/fits_${YYYYMMDDHHMU}.txt
#${CP} -p minimization_fort220.${YYYYMMDDHHMU} ${COMOUTgsi_rtma3d}/minimization_fort220.${YYYYMMDDHHMU}
#${CP} -p gsiparm.anl                        ${COMOUTgsi_rtma3d}/gsiparm.anl.${YYYYMMDDHHMU}
#${CP} -p diag_*                             ${COMOUTgsi_rtma3d}/

 #Comment out / MPondeca
#tar -zcvf obsfit_fort220.tgz  ./fort.* ./fit_*
#${CP} -p  obsfit_fort220.tgz                 ${COMOUTgsi_rtma3d}
#tar -zcvf misc_info.tgz       ./*info ./errtable ./prepobs_prep.bufrtable  ./*bias*  ./current_bad_aircraft ./gsd_sfcobs_uselist.txt ./gsd_sfcobs_provider.txt ./GSI_workdir_list
#${CP} -p  misc_info.tgz                      ${COMOUTgsi_rtma3d}
#gzip ${COMOUTgsi_rtma3d}/diag_*

# extra backup (NOT necessary)
#${LN} -sf ${COMOUTgsi_rtma3d}/${ANLrtma3d_FNAME} ${COMOUT}/${ANLrtma3d_FNAME}
#${CP} -p  ${pgmout_stdout}                       ${COMOUT}/${pgmout_stdout}_gsianl.${YYYYMMDDHHMU}
#${CP} -p  fits_${YYYYMMDDHHMU}.txt                 ${COMOUT}/fits_${YYYYMMDDHHMU}.txt

#/bin/rm -f ${DATA}/wrf_inout
#/bin/rm -f ${DATA}/sig*
#/bin/rm -f ${DATA}/obs*
/bin/rm -f ${DATA}/pe*

exit 0
