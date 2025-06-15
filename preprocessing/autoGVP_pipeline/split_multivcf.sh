
# this script gets called by run_splits.sh

# first input is directory with multivcfs
# second input is list of proband ids
# third input is output directory



#TEMP_DIR=$1/multi_to_single
#mkdir $TEMP_DIR

#echo 'Creating directory: ' $TEMP_DIR
echo ''
echo 'Starting multi-VCF to single-VCF conversion...placing outputs in ' $3
echo ''


for file in $1/*.multi.*.vcf.gz; do
    
     filename=$(basename $file)
     
     # Use sed to insert '_single' at the end
     NEW_FILE=$(echo $filename | sed 's/\.vcf\.gz/_single\.vcf\.gz/g'); #echo 'NEW_FILE = ' $NEW_FILE
     
     if [ -f "$3/$NEW_FILE" ]; then
        echo "$3/$NEW_FILE exists. Skipping..."
     else
	echo 'Converting ' $filename '...'
        sleep 3
        # Use bcftools to get just the sample we want. it will search through the list of probands and 
        # filter for just the IDs in there.  -Oz is for output = gz, and -S is for sample list.
        /mnt/isilon/dbhi_bfx/bin/bcftools-1.3 view --force-samples  -Oz -S $2 $file > $3/$NEW_FILE
     fi
done

echo "Finished $3!"
