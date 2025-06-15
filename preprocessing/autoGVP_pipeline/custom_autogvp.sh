#!/bin/bash
#SBATCH --mem=50G
#SBATCH --time=04:00:00
#SBATCH --array=1-78
#SBATCH --output=/scr1/users/stearb/U24/autoGVP/AutoGVP-main/results/%a_output.o

# Define cohort level output directories. These live under the results/ 
# directory of the respective prereq steps (InterVar, ANNOVAR and AutoPVS1)
# The files that contain the filepaths of the vcfs are in the data subdir of the Intervar dir.


#COHORT="CHD_NBL"
COHORT="TALL"
#COHORT="RSBD"

if [ "$COHORT" = "CHD_NBL" ]
then
        INTERVAR_COHORT_OUT_DIR="results_chd_nbl"
        ANNOVAR_COHORT_OUT_DIR="results_chd_nbl"
        AUTOPVS1_COHORT_OUT_DIR="results_chd_nbl"
        AUTOGVP_COHORT_OUT_DIR="results_chd_nbl"
        VCF_COHORT_FILES_LIST="chd_nbl_vcfs.txt"
elif [ "$COHORT" = "GNINT_MMC" ]
then
        INTERVAR_COHORT_OUT_DIR="results_gnint_mmc"
        ANNOVAR_COHORT_OUT_DIR="results_gnint_mmc"
        AUTOPVS1_COHORT_OUT_DIR="results_gnint_mmc"
        AUTOGVP_COHORT_OUT_DIR="results_gnint_mmc"
        VCF_COHORT_FILES_LIST="gnint_mmc_need_to_run_563.txt"
elif [ "$COHORT" = "RSBD" ]
then
        INTERVAR_COHORT_OUT_DIR="results_rsbd"
        ANNOVAR_COHORT_OUT_DIR="results_rsbd"
        AUTOPVS1_COHORT_OUT_DIR="results_rsbd"
        AUTOGVP_COHORT_OUT_DIR="results_rsbd"
        VCF_COHORT_FILES_LIST="rsbd_vcfs.txt"
elif [ "$COHORT" = "TALL" ]
then
        INTERVAR_COHORT_OUT_DIR="results_tall"
        ANNOVAR_COHORT_OUT_DIR="results_tall"
        AUTOPVS1_COHORT_OUT_DIR="results_tall"
        AUTOGVP_COHORT_OUT_DIR="results_tall"
        VCF_COHORT_FILES_LIST="tall_need_to_run_78.txt"  #"tall_vcfs.txt"
else
        echo 'No matching Cohort was found.';
        exit;
fi

# Generate the list of files using ls command
# create array of file paths to input VCFs
# the files that contain the filepaths for each cohort currently live at this directory:
readarray -t files < "/scr1/users/stearb/U24/InterVar/data/${VCF_COHORT_FILES_LIST}"

# Get the file for this array task
current_file=${files[$SLURM_ARRAY_TASK_ID-1]}
# drop path and just get file name
basefilename=$(basename $current_file)
echo "Running pipeline for cohort ${COHORT} on file: ${basefilename}"

# get just the uuid part of the filename to append the hg38_multianno.txt suffix to.
# also to append .vcf.gz to (the annotated vcf).
uuid_basename="$( cut -d '.' -f 1 <<< "$basefilename" )"

AUTO_GVP_DIR="/scr1/users/stearb/U24/autoGVP/AutoGVP-main"
AUTOGVP_OUTDIR="${AUTO_GVP_DIR}/results/${AUTOGVP_COHORT_OUT_DIR}"

# check that this VCF file hasnt been run through the pipeline already.
echo "checking if ${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-full.tsv.gz exists.."
if [ -f "${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-full.tsv.gz" ]; then
        echo "Input VCF ${basefilename} has already been run through this pipeline!"; echo "Exiting..."; exit;
fi

#echo "current_file = ${current_file}"; echo "basefilename = ${basefilename}";  echo "uuid_basename = ${uuid_basename}"
# check that $current_file exists
if [ ! -f "${current_file}" ]; then echo "Input VCF ${current_file} does not exist.";  echo "Exiting..."; exit; fi

#################################
######## Run InterVar ###########
#################################

INTERVAR_OUT_DIR="/scr1/users/stearb/U24/InterVar/results/${INTERVAR_COHORT_OUT_DIR}/"
cd /scr1/users/stearb/U24/InterVar/

#echo "Running InterVar..."
python /scr1/users/stearb/U24/InterVar/Intervar.py -b hg38 -i $current_file --input_type=VCF -o $INTERVAR_OUT_DIR$basefilename

# check if intervar file already exists (.hg38_multianno.txt.intervar)
### Remove INTERVAR output files that we dont need (we only need the 3rd file below) ###
# .hg38_multianno.txt     .hg38_multianno.txt.grl_p    .hg38_multianno.txt.intervar
intervar_text_file="${INTERVAR_OUT_DIR}${basefilename}.hg38_multianno.txt"
intervar_grl_file="${INTERVAR_OUT_DIR}${basefilename}.hg38_multianno.txt.grl_p"
if [[ -f $intervar_text_file ]]; then echo "Removing INTERVAR .txt output file ${intervar_text_file}"; rm $intervar_text_file; fi
if [[ -f $intervar_grl_file ]]; then echo "Removing INTERVAR .grl_p output file ${intervar_grl_file}"; rm $intervar_grl_file; fi


#################################
#########  Run ANNOVAR  #########
#################################

# Input here is the same as for InterVar, the original VCF ($current_file)
ANNOVAR_OUT_DIR="/scr1/users/stearb/U24/annovar/results/${ANNOVAR_COHORT_OUT_DIR}/"

echo 'Running ANNOVAR...'
cd /scr1/users/stearb/U24/annovar/

perl /scr1/users/stearb/U24/annovar/table_annovar.pl $current_file  hg38 --buildver hg38 --out "${ANNOVAR_OUT_DIR}${uuid_basename}" --remove --protocol gnomad211_exome,gnomad211_genome --operation f,f --vcfinput

### Remove ANNOVAR output files we dont need (only need the 3rd one below)
# $uuid_basename.avinput    $uuid_basename.hg38_multianno.vcf   $uuid_basename.hg38_multianno.txt
annovar_avinput_file="${ANNOVAR_OUT_DIR}${uuid_basename}.avinput"
annovar_vcf_file="${ANNOVAR_OUT_DIR}${uuid_basename}.hg38_multianno.vcf"

if [[ -f $annovar_avinput_file ]]; then echo "Removing ANNOVAR .avinput file ${annovar_avinput_file}"; rm $annovar_avinput_file; fi
if [[ -f $annovar_vcf_file ]]; then echo "Removing ANNOVAR .vcf file ${annovar_vcf_file}"; rm $annovar_vcf_file; fi

ANNOVAR_FILE="${ANNOVAR_OUT_DIR}${uuid_basename}.hg38_multianno.txt"

# check that ANNOVAR output file exists
if [ ! -f $ANNOVAR_FILE ]; then echo "ANNOVAR output file ${ANNOVAR_FILE}  does not exist.";  fi 
# remove intermediary files and then exit if ANNOVAR file was not produced!

##############################
######## Run AutoPVS1  #######
##############################

cd /scr1/users/stearb/U24/autoGVP/AutoGVP-main/prereqs/autopsv1/D3b-autoPVS1/
# Input here is the original VCF ($current_file)

AUTO_PVS1_OUT_DIR="/scr1/users/stearb/U24/autoGVP/AutoGVP-main/prereqs/autopsv1/D3b-autoPVS1/results/${AUTOPVS1_COHORT_OUT_DIR}/"
AUTO_PVS1_SCRIPT="/scr1/users/stearb/U24/autoGVP/AutoGVP-main/prereqs/autopsv1/D3b-autoPVS1/"

echo "----------------------------------------------------------"; echo ''
echo "Running AutoPVS1..."

# Input for AutoPVS1 is the original VCF ( $current_file)
python "${AUTO_PVS1_SCRIPT}autoPVS1_from_VEP_vcf.py" --genome_version hg38 --vep_vcf $current_file  > "${AUTO_PVS1_OUT_DIR}${uuid_basename}_autopvs1.txt"

echo "----------------------------------------------------------"

##############################
########  Run AutoGVP  #######
##############################

cd $AUTO_GVP_DIR

### Define input args: ###
# use original VCF as input ($current_file)
INTERVAR_FILE="${INTERVAR_OUT_DIR}${basefilename}.hg38_multianno.txt.intervar"
ANNOVAR_FILE="${ANNOVAR_OUT_DIR}${uuid_basename}.hg38_multianno.txt"
AUTOPVS1_FILE="${AUTO_PVS1_OUT_DIR}${uuid_basename}_autopvs1.txt"

echo "checking that AutoGVP input files exist..."
files=("${current_file}"  "${INTERVAR_FILE}" "${ANNOVAR_FILE}" "${AUTOPVS1_FILE}")
for file in "${files[@]}"; do
  if [[ ! -f "$file" ]]; then
    echo "AutoGVP input File ${file} does not exist."; echo "Exiting..."; exit
  fi
done

# Define path to where VCFs are so we can bind it to the singularity image
VCF_DIR="/mnt/isilon/opentargets/OpenPedCan_Data/"
# the GNINT and MCC cohort group did not have to be split, so those files are in their original directory.
if [ $COHORT = "CHD_NBL" ];  then VCF_DIR="/mnt/isilon/opentargets/OpenPedCan_Data/single_vcfs/split_vcfs/"; fi
if [ $COHORT = "GNINT_MMC" ] || [ $COHORT = "RSBD" ] || [ $COHORT = "TALL" ]; then VCF_DIR="/mnt/isilon/opentargets/OpenPedCan_Data/"; fi

echo 'pulling AutoGVP docker container...'
singularity pull docker://pgc-images.sbgenomics.com/diskin-lab/autogvp:latest

echo 'Running AutoGVP...'
# HAD TO ADD DOUBLE QUOTES AROUND PASS AND . IN `filter_criteria` arg:   (gnomad_3_1_1_FILTER=PASS|gnomad_3_1_1_FILTER=.)
singularity exec --bind $PWD:$PWD,$INTERVAR_OUT_DIR:$INTERVAR_OUT_DIR,$ANNOVAR_OUT_DIR:$ANNOVAR_OUT_DIR,$AUTO_PVS1_OUT_DIR:$AUTO_PVS1_OUT_DIR,$VCF_DIR:$VCF_DIR \
docker://pgc-images.sbgenomics.com/diskin-lab/autogvp:latest bash run_autogvp.sh \
--workflow="custom" \
--vcf=$current_file \
--clinvar=$AUTO_GVP_DIR/data/clinvar.vcf.gz \
--selected_clinvar_submissions=$AUTO_GVP_DIR/data/ClinVar-selected-submissions.tsv \
--filter_criteria='FORMAT/DP>=10 (FORMAT/AD[0:1-])/(FORMAT/DP)>=0.08 (gnomad_3_1_1_AF_non_cancer_popmax<0.001|gnomad_3_1_1_AF_non_cancer_popmax=".") (gnomad_3_1_1_FILTER="PASS"|gnomad_3_1_1_FILTER=".")' \
--intervar=$INTERVAR_FILE \
--multianno=$ANNOVAR_FILE \
--autopvs1=$AUTOPVS1_FILE \
--variant_summary=$AUTO_GVP_DIR/data/variant_summary.txt.gz \
--submission_summary=$AUTO_GVP_DIR/data/submission_summary.txt.gz \
--conceptIDs=$AUTO_GVP_DIR/data/clinvar_all_disease_concept_ids.txt \
--conflict_res="latest" \
--outdir=$AUTOGVP_OUTDIR \
--out=$uuid_basename

echo '---------------- Done running AutoGVP! -------------------'


##########################
####### Clean up #########
##########################

# Check that final output file (ie uuid-autogvp-annotated-full.tsv) exists and if it does gzip it
if [[ ! -f "${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-full.tsv" ]]; then
        echo "AUTOGVP FILE NOT FOUND!"; 
else
        echo 'Zipping AutoGVP full output file...'; 
	gzip "${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-full.tsv"; 
	echo 'Done zipping.';
fi

# Remove the large ANNOVAR file
if [[ -f $ANNOVAR_FILE ]]; then echo "Removing ${ANNOVAR_FILE}"; rm "${ANNOVAR_FILE}";  fi

# Remove AutoPSV1 file
if [[ -f $AUTOPVS1_FILE ]]; then echo "Removing ${AUTOPVS1_FILE}"; rm "${AUTOPVS1_FILE}";  fi

# Remove abridged file
#echo "Checking if ${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-abridged.tsv or"
#echo "${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-full.tsv exists..."

if [[ -f "${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-abridged.tsv" ]]; then
	#echo 'Removing unzipped autoGVP full-annotation file...';
        #rm "${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-full.tsv"; echo 'Done.'
        echo "Removing ${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-abridged.tsv"
	rm "${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-abridged.tsv"; echo 'Done.'
fi

# If autoGVP or ANNOVAR does not run correctly there will be output files that were not properly removed. 
# Check if they exist and remove them if they do.
files=("${AUTOGVP_OUTDIR}/${uuid_basename}_intervar_filtered.txt"  "${AUTOGVP_OUTDIR}/${uuid_basename}_autopvs1_filtered.tsv" "${AUTOGVP_OUTDIR}/${uuid_basename}_multianno_filtered.txt" "${AUTOGVP_OUTDIR}/${uuid_basename}.filtered.vcf" "${ANNOVAR_OUT_DIR}${uuid_basename}.hg38_gnomad211_exome_*" "${ANNOVAR_OUT_DIR}${uuid_basename}.hg38_gnomad211_genome_*" )

for file in "${files[@]}"; do
  echo "Checking if ${file} exists..." 
  if [[ -f "$file" ]]; then echo "Removing autoGVP/ANNOVAR intermediary Files ${file}."; fi
done

chmod 755 ${AUTOGVP_OUTDIR}/${uuid_basename}-autogvp-annotated-full.tsv.gz
