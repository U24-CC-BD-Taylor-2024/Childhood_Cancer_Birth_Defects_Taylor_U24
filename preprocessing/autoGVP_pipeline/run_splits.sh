
OPEN_PEDCAN_DIR="/mnt/isilon/opentargets/OpenPedCan_Data"

GNINT_DIR="${OPEN_PEDCAN_DIR}/CHD_KF_phs001846/"
CHD_DIR="${OPEN_PEDCAN_DIR}/CHD_KF_phs001138/vcf/vep_105/"
NBL_DIR="${OPEN_PEDCAN_DIR}/CHD_KF_phs001436/vcf/"
RSBD_DIR="${OPEN_PEDCAN_DIR}/CHD_KF_phs002590/vcf/"
MMC_DIR="${OPEN_PEDCAN_DIR}/CHD_KF_phs002187/"

chd_probands="${OPEN_PEDCAN_DIR}/single_vcfs/proband_ids/chd_probands.txt"
nbl_probands="${OPEN_PEDCAN_DIR}/single_vcfs/proband_ids/nbl_probands.txt"
gnint_probands="${OPEN_PEDCAN_DIR}/single_vcfs/proband_ids/gnint_probands.txt"
rsbd_probands="${OPEN_PEDCAN_DIR}/single_vcfs/proband_ids/rsbd_probands.txt"
mmc_probands="${OPEN_PEDCAN_DIR}/single_vcfs/proband_ids/mmc_probands.txt"


chd_out_dir="${OPEN_PEDCAN_DIR}/single_vcfs/split_vcfs/chd/"
gnint_out_dir="${OPEN_PEDCAN_DIR}/single_vcfs/split_vcfs/gnint/"
nbl_out_dir="${OPEN_PEDCAN_DIR}/single_vcfs/split_vcfs/nbl/"
rsbd_out_dir="${OPEN_PEDCAN_DIR}/single_vcfs/split_vcfs/rsbd/"
mmc_out_dir="${OPEN_PEDCAN_DIR}/single_vcfs/split_vcfs/mmc/"

#echo 'Running split vcf script in GNINT directory'
#sleep 5s
#bash $OPEN_PEDCAN_DIR/single_vcfs/split_multivcf.sh $GNINT_DIR $gnint_probands $gnint_out_dir

#echo 'Running split vcf script in CHD directory'
#sleep 5s
#bash $OPEN_PEDCAN_DIR/single_vcfs/split_multivcf.sh $CHD_DIR $chd_probands $chd_out_dir

#echo 'Running split vcf script in NBL directory'
#sleep 5s
#bash $OPEN_PEDCAN_DIR/single_vcfs/split_multivcf.sh $NBL_DIR $nbl_probands $nbl_out_dir

#echo ''
#echo 'Running split vcf script in RSBD directory'
#sleep 5s

#bash $OPEN_PEDCAN_DIR/single_vcfs/split_multivcf.sh $RSBD_DIR $rsbd_probands $rsbd_out_dir


echo 'Running split vcf script in MMC directory'
sleep 5s

echo $OPEN_PEDCAN_DIR
echo $MMC_DIR
echo $mmc_probands
echo $mmc_out_dir

bash $OPEN_PEDCAN_DIR/single_vcfs/split_multivcf.sh $MMC_DIR $mmc_probands $mmc_out_dir



