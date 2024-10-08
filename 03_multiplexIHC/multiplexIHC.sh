#!/bin/bash +x
start_time=$(date)
# create directory
# cd $(dirname $0)
#-------------------------------------------------------------------------------
# definitions
source masterParameters.txt
#-------------------------------------------------------------------------------
# create folders

mkdir Analysis

mkdir -p $PWD/Analysis/1_raw_input  #store all downloaded input files
mkdir -p $PWD/Analysis/2_crop_output
mkdir -p $PWD/Analysis/3_aligned_stacks
mkdir -p $PWD/Analysis/4_artifacts_output
mkdir -p $PWD/Analysis/5_artifacts_removal
mkdir -p $PWD/Analysis/lut
mkdir -p $PWD/Analysis/6_deconvolution
mkdir -p $PWD/Analysis/7_pseudocolor
mkdir -p $PWD/Analysis/8_montages
mkdir -p $PWD/Analysis/9_zStacks

mkdir -p $PWD/Analysis/10_image_cytometry_output
mkdir -p $PWD/Analysis/10_image_cytometry_output/1_images
mkdir -p $PWD/Analysis/10_image_cytometry_output/1_images/ImageAfterMath
mkdir -p $PWD/Analysis/10_image_cytometry_output/1_images/multiplex
#mkdir -p $PWD/Analysis/10_image_cytometry_output/1_images/crop

mkdir -p $PWD/Analysis/10_image_cytometry_output/2_relationships
mkdir -p $PWD/Analysis/10_image_cytometry_output/3_experiments

mkdir -p $PWD/Analysis/10_image_cytometry_output/4_analysis
mkdir -p $PWD/Analysis/10_image_cytometry_output/4_analysis/fibrosis
mkdir -p $PWD/Analysis/10_image_cytometry_output/4_analysis/fibrosisImages
mkdir -p $PWD/Analysis/10_image_cytometry_output/4_analysis/multiplexNuclei
mkdir -p $PWD/Analysis/10_image_cytometry_output/4_analysis/multiplexCells
mkdir -p $PWD/Analysis/10_image_cytometry_output/4_analysis/multiplexCytoplasm

#-------------------------------------------------------------------------------
## convert to bioformat
## ref: https://www.glencoesoftware.com/blog/2019/12/09/converting-whole-slide-images-to-OME-TIFF.html
/Volumes/Multiplexing_G_Drive/Multiplexing/bioformats2raw-0.4.0/bin/bioformats2raw \
/Volumes/Multiplexing_G_Drive/Multiplexing/multiplexing/TMA_Ahmed_new/QuPath/tma-MULTIPLEX-022022/0-H\&E/1.mrxs \
/Volumes/Multiplexing_G_Drive/Multiplexing/multiplexing/TMA_Ahmed_new/QuPath/new/ \
--resolutions 6

#-------------------------------------------------------------------------------
# Step 1: copy data from QuPath project to Analysis folder
echo -e "\033[1;36m...starting raw_data download...\033[0m"
## cp -R $PWD/${folder_name}/* $PWD/1_raw_input
cp -n $PWD/QuPath/**/raw_cores/*.tiff $PWD/Analysis/1_raw_input/
echo -e "\033[1;35mraw_data download completed....\033[0m"
#-------------------------------------------------------------------------------
# Step 2: Autocrop all images
echo -e "\033[1;36m...starting autocroping...\033[0m";{
  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  --headless --heap $memory --console \
  --run $PWD/scripts/autocrop_images.py \
  "input='$PWD/Analysis/1_raw_input',
  output='$PWD/Analysis/2_crop_output',
  suffix='.tiff'"
} 2> /dev/null
echo -e "\033[1;35mAutocrop completed....\033[0m"
#-------------------------------------------------------------------------------
# Step 3: Image registration
echo -e "\033[1;36m...coalesce imges per core...\033[0m"
## coalesce all images from one core to one folder
for file in $(find $PWD/Analysis/2_crop_output -type f -iname "*.tiff");
do
  data_core=${file%_*}
  data_core=${data_core%_*} ## only if second last underscore
  core=${data_core}
  if [ ! -d $core ]; then
    mkdir $core
  fi
  mv $file $core
done
## registration of images using SIFT
echo -e "\033[1;36m...starting registration...\033[0m";{
  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  --headless --heap $memory --console \
  --run $PWD/scripts/stackAlign_SIFT.py \
  "input='$PWD/Analysis/2_crop_output',
  output='$PWD/Analysis/3_aligned_stacks',
  suffix='.tiff',
   nuImages='$nuImages'"
} 2> /dev/null
echo -e "\033[1;35mRegistration completed....\033[0m"
#-------------------------------------------------------------------------------
# Step 4: import all definations
for file in $(find $PWD/Analysis/2_crop_output -type f -iname "*.tiff");
do
  fileNameExt=${file%%.*}
  ## fileName=${fileNameExt##*_} ## for last underscore
  fileName= echo $fileNameExt | rev | cut -f1,2 -d"_" | rev
  echo $fileName
done >> fileName.txt


grep '^.' fileName.txt > fileNameList.txt
rm fileName.txt

## create definations.txt file for all staining
cat -n fileNameList.txt | sort -uk2 | sort -nk2 | cut -f2- | sort -f -V > definitions.txt
rm fileNameList.txt
#-------------------------------------------------------------------------------
# Step 5: refining names from image cytometry
echo -e "\033[1;36m...starting image renaming...\033[0m"
## rename files by arranging name to
## (core number)_(replicate number)_(staining)
## SIFT generate names with digits
loop=1
while IFS='' read -r LINE || [ -n "${LINE}" ];
do
  count="0000$loop"
  ## echo ${count:(-4)}
  loop=$((loop+1))
  ## echo "${LINE}"
  for file in $(find $PWD/Analysis/3_aligned_stacks/* -type f -iname "*.tif");
  do
    if [[ $file =~ ${count:(-4)} ]]; then
      mv "$file" "${file//_${count:(-4)}/_${LINE}}"
    fi
  done
done < definitions.txt
echo -e "\033[1;35mImage renaming completed....\033[0m"
rm definitions.txt
#-------------------------------------------------------------------------------
# Step 6: Artifacts definations
## define Artifacts
echo -e "\033[1;36m...defining artifacts...\033[0m";{

  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  --heap $memory --console \
  --run $PWD/scripts/deconvolution_lut.py \
  "output='$PWD/Analysis/lut'"

  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  --headless --heap $memory --console \
  --run $PWD/scripts/defineArtifacts_boost.ijm \
  "input='$PWD/Analysis/3_aligned_stacks',
  inputLut='$PWD/Analysis/lut',
  output='$PWD/Analysis/4_artifacts_output',
  suffix='.tif',
  mask='$mask'"
} 2> /dev/null
echo -e "\033[1;35mArtifcacts definitions completed....\033[0m"
## coalesce all images from one core to one folder
for file in $(find $PWD/Analysis/3_aligned_stacks -type f -iname "*.tif");
do
  data_core=${file%_*}
  data_core=${data_core%_*} ## only for second underscore
  core=${data_core}
  if [ ! -d $core ]; then
    mkdir $core
  fi
  mv $file $core
done
#-------------------------------------------------------------------------------
# Step 7: Masks substraction
## coalesce all images from one core to one folder
for file in $(find $PWD/Analysis/4_artifacts_output -type f -iname "*.tif");
do
  data_core=${file%_*}
  data_core=${data_core%_*} ## only for second underscore
  core=${data_core}
  if [ ! -d $core ]; then
    mkdir $core
  fi
  mv $file $core
done
# define masks
echo -e "\033[1;36m...defining masks...\033[0m";{
  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  --headless --heap $memory --console \
  --run $PWD/scripts/defineMasks.py \
  "input='$PWD/Analysis/4_artifacts_output',
  output='$PWD/Analysis/3_aligned_stacks',
  suffix='.tif',
  nuImages=$nuImages"
} 2> /dev/null
echo -e "\033[1;35mMasks defination completed....\033[0m"
# clean names
# for file in $(find $PWD/Analysis/3_aligned_stacks/* -maxdepth 0 -type f -iname "*.tif");
# do
#   mv "$file" "${file//_${first_stain}}"
# done
# ## move file to folder of same name
# for file in $PWD/Analysis/3_aligned_stacks/*.tif;
# do
#   [[ -file "$file" ]] || continue
#   dir="${file%.*}"
#   mv "$file" "$dir"
# done
#-------------------------------------------------------------------------------
# Step 8: Artifacts substraction
echo -e "\033[1;36m...substracting artifacts...\033[0m";{
  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  --headless --heap $memory --console \
  --run $PWD/scripts/artifactsSubstraction.ijm \
  "input='$PWD/Analysis/3_aligned_stacks',
  output='$PWD/Analysis/5_artifacts_removal',
  suffix='.tif'"
} 2> /dev/null
echo -e "\033[1;35mArtifacts substraction completed....\033[0m"
## white balance correction
# echo -e "\033[1;36m...Correcting white balance...\033[0m";{
#   /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
#   --headless --heap $memory --console \
#   --run $PWD/scripts/whiteBalance.ijm \
#   "input='$PWD/Analysis/5_artifacts_removal',output='$PWD/Analysis/5_artifacts_removal',suffix='.tif'"
# } &> /dev/null
# echo -e "\033[1;35mWhite Balancing completed....\033[0m"
#-------------------------------------------------------------------------------
# Step 9: Deconvolution
echo -e "\033[1;36m...starting Deconvolution...\033[0m";{

  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  --headless --heap $memory --console \
  --run $PWD/scripts/deconvolution_boost.ijm \
  "input='$PWD/Analysis/5_artifacts_removal',
  inputLut='$PWD/Analysis/lut',
  output='$PWD/Analysis/6_deconvolution',
  suffix='.tif',
  costain='$costain',
  costain1='$costain1',
  costain2='$costain2',
  heCostain='$heCostain'"

  # /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  # --heap $memory --console \
  # --run $PWD/scripts/deconvolution.py \
  # "input='$PWD/Analysis/5_artifacts_removal',
  # output='$PWD/Analysis/6_deconvolution',
  # suffix='.tif',
  # costain='$costain',
  # heCostain='$heCostain'"
} 2> /dev/null
echo -e "\033[1;35mDeconvolution completed....\033[0m"
#-------------------------------------------------------------------------------
# Step 10: Assign psudocolor
echo -e "\033[1;36m...assigning pseudocolor...\033[0m";{
  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  --headless --heap $memory --console \
  --run $PWD/scripts/colorAssignment.ijm \
  "input='$PWD/Analysis/6_deconvolution',
  output='$PWD/Analysis/7_pseudocolor',
  suffix='.tif',
  costain='$costain',
  costain1='$costain1',
  costain2='$costain2',
  stain1costain='$stain1costain',
  stain2costain='$stain2costain',
  stain1costain1='$stain1costain1',
  stain2costain1='$stain2costain1',
  stain1costain2='$stain1costain2',
  stain2costain2='$stain2costain2',
  stain1Namecostain='$stain1Namecostain',
  stain2Namecostain='$stain2Namecostain',
  stain1Namecostain1='$stain1Namecostain1',
  stain2Namecostain1='$stain2Namecostain1',
  stain1Namecostain2='$stain1Namecostain2',
  stain2Namecostain2='$stain2Namecostain2',
  mask='$mask'"
} &> /dev/null
## coalesce all images from one core to one folder
for i in $(find $PWD/Analysis/7_pseudocolor -type f -iname "*.tif");
do
  data_core=${i%_*}
  data_core=${data_core%_*} ## only for second underscore
  core=${data_core}
  if [ ! -d $core ]; then
    mkdir $core
  fi
  mv $i $core
done
echo -e "\033[1;35mPseudocolor assignment completed....\033[0m"
#-------------------------------------------------------------------------------
# Step 11: Montage
echo -e "\033[1;36m...performing montage...\033[0m";{
  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 \
  --headless --heap $memory --console \
  --run $PWD/scripts/montage.ijm \
  "input='$PWD/Analysis/7_pseudocolor',output='$PWD/Analysis/8_montages', suffix='.tif'"
} &> /dev/null
## rename montage files by core names
# for file in $(find $PWD/Analysis/8_montages/* -type f -iname "*.tif");
# do
#   mv "$file" "${file//_${first_stain}}"
# done
echo -e "\033[1;35mMontage completed....\033[0m"
#-------------------------------------------------------------------------------
# Step 12: Z-stacking
echo -e "\033[1;36m...performing Z-stacking...\033[0m";{
  /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx --ij2 --headless \
  --heap $memory --console \
  --run $PWD/scripts/zStack.py \
  "input='$PWD/Analysis/7_pseudocolor',output='$PWD/Analysis/9_zStacks', suffix='.tif', nuPseudo='$nuPseudo'"
} &> /dev/null
## rename montage files by core names
# for file in $(find $PWD/Analysis/9_zStacks/* -type f -iname "*.tif");
# do
#   mv "$file" "${file//_${first_stain}}"
# done
# echo -e "\033[1;35mZ-stack completed....\033[0m"
#-------------------------------------------------------------------------------
# Step 13: Image_cytometry pipeline
for file in $(find $PWD/Analysis/5_artifacts_removal -type f -iname "*.tif");
do
  data_core=${file%_*}
  data_core=${data_core%_*} ## only for second underscore
  core=${data_core}
  if [ ! -d $core ]; then
    mkdir $core
  fi
  mv $file $core
done

echo -e "\033[1;36m...starting image cytometry...\033[0m"
## coalesce all images from one core to one folder
# for i in $(find $PWD/Analysis/6_deconvolution/* -type f -iname "*.tif");
# do
#   data_core=${i%_*}
#   core=${data_core}
#   if [ ! -d $core ]; then
#     mkdir $core
#   fi
#   mv $i $core
# done
## run cell_profiler
echo "\033[1;35mRunning cell cytometry pipeline............\033[0m"; {
  /Volumes/Multiplexing_G_Drive/Multiplexing/programs/CellProfiler.app/Contents/MacOS/cp -r -c \
  -p $PWD/scripts/multiplexIntesityAnalysis_finalv002.cppipe \
  -i $PWD/Analysis/5_artifacts_removal \
  -o $PWD/Analysis/10_image_cytometry_output \
  -L 20
} &> /dev/null
## sort cell_profiler results
### move files to ImageAfterMath
cp -r $PWD/Analysis/10_image_cytometry_output/1_images/5_artifacts_removal/. \
$PWD/Analysis/10_image_cytometry_output/1_images/ImageAfterMath && rm -rf \
$PWD/Analysis/10_image_cytometry_output/1_images/5_artifacts_removal
find $PWD/Analysis/10_image_cytometry_output/1_images/ImageAfterMath -mindepth 2 -type f -exec mv {} $PWD/Analysis/10_image_cytometry_output/1_images/ImageAfterMath \;
find $PWD/Analysis/10_image_cytometry_output/1_images/ImageAfterMath -depth -type d -exec rmdir {} \; 2>/dev/null
### move results
find $PWD/Analysis/10_image_cytometry_output -type f -name '*rel*' \
-exec mv {} $PWD/Analysis/10_image_cytometry_output/2_relationships \;
find $PWD/Analysis/10_image_cytometry_output -type f -iname '*exp*' \
-exec mv  {} $PWD/Analysis/10_image_cytometry_output/3_experiments \;
find $PWD/Analysis/10_image_cytometry_output -type f -iname '*_image*' \
-exec mv  {} $PWD/Analysis/10_image_cytometry_output/4_analysis/fibrosisImages \;
find $PWD/Analysis/10_image_cytometry_output -maxdepth 1 -type f -name 'fibrosis_*' \
-exec mv  {} $PWD/Analysis/10_image_cytometry_output/4_analysis/fibrosis \;
find $PWD/Analysis/10_image_cytometry_output -maxdepth 1 -type f -name '*_NucleiObjects.csv' \
-exec mv {} $PWD/Analysis/10_image_cytometry_output/4_analysis/multiplexNuclei \;
find $PWD/Analysis/10_image_cytometry_output -maxdepth 1 -type f -name '*_CellsObjects.csv' \
-exec mv {} $PWD/Analysis/10_image_cytometry_output/4_analysis/multiplexCells \;
find $PWD/Analysis/10_image_cytometry_output -maxdepth 1 -type f -name '*_CytoplasmObjects.csv' \
-exec mv {} $PWD/Analysis/10_image_cytometry_output/4_analysis/multiplexCytoplasm \;
# find $PWD/Analysis/10_image_cytometry_output -maxdepth 1 -type f -name '*_Cells.csv' \
# -exec mv {} $PWD/Analysis/10_image_cytometry_output/4_analysis/multiplexCells \;
echo -e "\033[1;35mImage cytometry completed....\033[0m"
#-------------------------------------------------------------------------------
end_time=$(date)
echo "processing started on: ${start_time}" > summary.txt
echo "Details of analysis:" >> summary.txt
tree -L 1 --du -h  >> summary.txt
echo "Details of Image cytometry:" >> summary.txt
tree $PWD/7_image_cytometry_output -L 2 --du -h  >> summary.txt
echo "processing ended on: ${end_time}" >> summary.txt
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@

find $PWD/Analysis/5_artifacts_removal/ -mindepth 2 -type f -print -exec mv {} $PWD/Analysis/5_artifacts_removal/ \;
