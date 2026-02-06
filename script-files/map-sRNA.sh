#!/usr/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=5g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --time=30:00
#SBATCH --qos=short

hostname

set -ux

###################################################################################################
###################################################################################################
# Extract variables from input argument (comma-separated, possibly quoted)
VARI=$(echo $1 | sed 's/,/\t/g;s/"//g')
eval $VARI
VARI=$1
echo $VARI | sed 's/,/\n/g'

# Source utility functions (assumes $SCRIPTdir is set)
source ${SCRIPTdir}tools

# Determine number of threads and available memory for the job
CORES=$(( $SLURM_CPUS_PER_TASK * 2 ))
MEM=$(scontrol show job $SLURM_JOBID | grep TRES | awk '{ split($NF, X, /,|=|G/);{print X[5]-5}}' | head -n 1)

###################################################################################################
# Setup phase

# Extract fasta from pass and fail individually and merge to single fasta file

# Extract name of current sensor and related info from a lookup file
libNAME=$(cat ${OPENdir}library-fasta-pairing.txt |  sed -n ${SLURM_ARRAY_TASK_ID}p | tr ' ' '\t' | cut -f 1 )
SENSOR=$(cat ${OPENdir}library-fasta-pairing.txt |  sed -n ${SLURM_ARRAY_TASK_ID}p | tr ' ' '\t' | cut -f 2 )
QPCRnorm=$(cat ${OPENdir}library-fasta-pairing.txt |  sed -n ${SLURM_ARRAY_TASK_ID}p | tr ' ' '\t' | cut -f 3 )

echo $libNAME $SENSOR
locTMP=${TMPdir}TMP_map-data/$libNAME/
mkdir -p $locTMP

###################################################################################################
# Build index for current sensor to determine boundaries

rm -rf ${locTMP}index*
bowtieBuild --quiet ${OPENdir}fasta-files/$SENSOR ${locTMP}index_origPlasmid

###################################################################################################
# Determine start and end of sensor in the plasmid and limit alignment

# Map sensor features to plasmid and convert to BED format
bowtie -m 1 -v 1 -f -S -x ${locTMP}index_origPlasmid ${UTILITYdir}sensor_parts.fa |
  samtools view -bS - |
  bedtools bamtobed -i - >${locTMP}sensor_parts.bed

# Extract start and end coordinates for sensor and UTRs from BED file
START=$(awk -v OFS="\t" '{if($4=="5UTR") print $2-100}' ${locTMP}sensor_parts.bed)
END=$(awk -v OFS="\t" '{if($4=="SV40") print $3+100}' ${locTMP}sensor_parts.bed)

STARTutr=$(awk -v OFS="\t" '{if($4=="GFP") print $3}' ${locTMP}sensor_parts.bed)
ENDutr=$(awk -v OFS="\t" '{if($4=="SV40") print $2}' ${locTMP}sensor_parts.bed)

# Extract sensor sequence from plasmid using coordinates
seqkit subseq -r $START:$END ${OPENdir}fasta-files/$SENSOR > ${locTMP}sensor.fa

# If UTR coordinates are available, extract UTR sequence and analyze U content
if [[ ! -z $STARTutr && ! -z $ENDutr ]]; then
  STARTutr=$(awk -v OFS="\t" '{if($4=="GFP") print $3}' ${locTMP}sensor_parts.bed)
  ENDutr=$(awk -v OFS="\t" '{if($4=="SV40") print $2}' ${locTMP}sensor_parts.bed)

  # Extract UTR sequence from plasmid
  seqkit subseq -r $STARTutr:$ENDutr ${OPENdir}fasta-files/$SENSOR > ${locTMP}UTR.fa

  # Determine number of T nucleotides in the UTR (single, double, triple, total)
  seqkit seq --upper-case ${locTMP}UTR.fa | seqkit fx2tab | 
  mawk -v OFS="\t" '
  BEGIN{
    print "NAME","totalU","U","UU","UUU","UTRlength"
  }
  {
    SEQ=$NF
    gsub(/A|C|G/,"", SEQ)
    totalU+=length(SEQ)

    SEQ=$NF
    gsub(/A|C|G/," ", SEQ)
    n=split(SEQ, splitSEQ, " ")

    T=0
    TT=0
    TTT=0
    for(i=1; i<=n; i++){
    if(splitSEQ[i]=="T"){ T++}
    if(splitSEQ[i]=="TT"){ TT++}
    if(splitSEQ[i]~"TTT"){ TTT++}
    }
    print $1, totalU,T,TT,TTT , length($NF)
  }' > ${locTMP}U-content.UTR.txt
fi

# Determine nucleotide content along the sensor using sliding window
seqkit sliding -s 1 -W 50 ${locTMP}sensor.fa > ${locTMP}sliding.fa
faCount  ${locTMP}sliding.fa | 
  head -n -2  | 
  mawk -v OFS="\t" '
  {
  if(NR==1){
    for(i=3; i<=6;i++){
    $i="NUC_"$i  
   }
    print "POS~"$3,$4,$5,$6
  }else{
   split($1, splitNAME, /_sliding:/)
   split(splitNAME[2], splitPOS, /-/)
   print (splitPOS[1]+splitPOS[2])/2-0.5"~"$3*100/$2,$4*100/$2,$5*100/$2,$6*100/$2
  }
  }' > ${locTMP}nuc-content.txt

###################################################################################################
# Align reads to plasmid and controls

# Build index for current sensor and tj-locus
rm -rf ${locTMP}input.fa
cat ${locTMP}sensor.fa  > ${locTMP}input.fa
cat ${UTILITYdir}tj-locus.fa >> ${locTMP}input.fa

# Build bowtie index
bowtieBuild --quiet ${locTMP}input.fa ${locTMP}index

# Map reads to sensor and tj-locus, normalize shRNA, filter for piRNAs, and sort
gunzip -c ${APpath}individual-libraries/${libNAME}/${libNAME}_annotated.fa.gz |
  # Count normalization shRNA 
  mawk -v OFS="\t" -v TMP=$locTMP '
  {
  # Extract count for shRNA for normalization 
  if( $1~">" && $1 ~ "TTCTTAAAGACCTTCACGGATG"){ 
    split($1, splitNAME, /:|=/)
    COUNT+=splitNAME[4]
  }
  print
  }
  END{    
  print COUNT > TMP "shRNAnorm.txt"
  }' | 
  # Map to sensor
  bowtie -p $CORES -m 1 -v 0 -f -S -x ${locTMP}index - | tee ${locTMP}mapped.sam |
  samtools view -bS - |
  bedtools bamtobed -i - |
  # Filter for piRNAs by length (24-35nt)
  mawk -v OFS="\t" '
  {
  split($4, splitNAME, /:|=/)
  if(length(splitNAME[2])>23 && length(splitNAME[2])<=35){
    for( i=1; i<=splitNAME[4]; i++) {
    print 
    }
  }
  }' |
  sort -k1,1 -k2,2n >${locTMP}mapped.bed

###################################################################################################
# Determine UTR boundaries and quantify UTR reads

# Map sensor features again to highlight in the graph
bowtie -m1 -v 1 -f -S -x ${locTMP}index ${UTILITYdir}sensor_parts.fa |
  samtools view -bS - |
  bedtools bamtobed -i - >${locTMP}sensor_parts.bed

# Extract UTR boundaries on the sensor
STARTutrOnSensor=$(awk -v OFS="\t" '{if($4=="GFP") print $3}' ${locTMP}sensor_parts.bed)
ENDutrOnSensor=$(awk -v OFS="\t" '{if($4=="SV40") print $2}' ${locTMP}sensor_parts.bed)

# Count number of reads mapping to the UTR region
UTRcount=$(mawk -v OFS="\t" -v STARTutr=$STARTutr -v ENDutr=$ENDutr '
  {
  if($2>=STARTutr && $2<=ENDutr && $1 !~ "tj-locus"){
    X+=1
  }
  }
  END{
  print X
  }' ${locTMP}mapped.bed )

###################################################################################################
# Generate bedgraph files

# Create length file for genomecov
seqkit fx2tab --name --length ${locTMP}input.fa | awk -v OFS="\t" '{print $1,$NF}' > ${locTMP}SENSORlength.txt

# Unnormalized sense coverage
bedtools genomecov -d -strand + -g ${locTMP}SENSORlength.txt -i ${locTMP}mapped.bed |
  mawk -v OFS="\t" ' 
  {
  $3=$3
  print 
  }' >${locTMP}${libNAME}_sense.bg

# Unnormalized antisense coverage (negative values)
bedtools genomecov -d -strand - -g ${locTMP}SENSORlength.txt -i ${locTMP}mapped.bed |
  mawk -v OFS="\t" ' 
  {
  $3=-$3
  print 
  }' >${locTMP}${libNAME}_antisense.bg

##################################################################################################
# Extract normalization factors

# miRNA normalization: from file if available, else default to 1
if [[ -s ${APpath}individual-libraries/${libNAME}/normalization.txt ]]; then
  echo "Extracting miRNA normalization factor from file"
  miRNA_NORM=$(tail -n 1 ${APpath}individual-libraries/${libNAME}/normalization.txt | tr ' ' '\t' | cut -f 1)
else
  echo "No miRNA normalization factor found, using 1 as default"
  miRNA_NORM=1
fi

# piRNA normalization: count piRNA reads (24-34nt, not tRNA/miRNA/rRNA/snRNA/snoRNA)
piRNA_NORM=$(seqkit fx2tab ${APpath}individual-libraries/${libNAME}/${libNAME}_annotated.fa.gz | 
  awk -v OFS="\t" '
  {
  if(length($NF)>23 && length($NF)<35){
    split($1, splitNAME, /:|=/)
    if(splitNAME[10] !~"tRNA" && splitNAME[10] !~"miRNA" && splitNAME[10] !~"rRNA" && splitNAME[10] !~"snRNA" && splitNAME[10] !~"snoRNA"){
    X+=splitNAME[4]
    }
  }
  }
  END{
  print X/1000000
  }' 
)

# shRNA normalization: from shRNAnorm.txt, scaled
shRNA_NORM=$(awk -v miRNA_NORM=${miRNA_NORM} '{if(NR==1){print $1/1000}}' ${locTMP}shRNAnorm.txt)

# QPCR normalization: use value if available, else default to 1
if [[ $QPCRnorm == "" ]]; then
  QPCRnorm=1
fi

# Normalize sense strand data and add nucleotide content
mawk -v miRNA_NORM=${miRNA_NORM} -v piRNA_NORM=${piRNA_NORM} -v shRNA_NORM=${shRNA_NORM} -v QPCRnorm=${QPCRnorm} -v OFS="\t" -v NAME=${libNAME} -v NUCcont=${locTMP}nuc-content.txt '
  BEGIN{
  SWITCH="N"
  while((getline LINE < NUCcont)){
    split(LINE, splitLINE, /~/)
    if(SWITCH=="N"){
    HEADER=splitLINE[2]
    SWITCH="Y"
    }else{
    NUC[splitLINE[1]]=splitLINE[2]
    }
  }
  print "SENSOR","POSITION","COUNT_RAW","COUNT_miRNAnorm","COUNT_piRNAnorm","COUNT_shNORM","COUNT_QPCRnorm","libNAME",HEADER
  }
  {
    # Normalize for sequencing depth
    $4=$3/miRNA_NORM
    $5=$3/piRNA_NORM
    $6=$3/shRNA_NORM
    $7=$4/QPCRnorm
    if($2 in NUC){
    currNUC=NUC[$2]
    }else{
    currNUC="0 0 0 0"
    }
    print $0,NAME,currNUC
  }' ${locTMP}${libNAME}_sense.bg | tr ' ' '\t' > ${locTMP}${libNAME}_sense.norm.bg

# Print all normalization factors and other stats to file
echo $libNAME $miRNA_NORM $piRNA_NORM $shRNA_NORM $QPCRnorm $STARTutrOnSensor $ENDutrOnSensor $UTRcount | tr ' ' '\t' > ${locTMP}/stats.tmp

# Plot data using R script
Rscript ${SCRIPTdir}plot_hist.R NAME=$libNAME TMP=$locTMP OPENdir=${OPENdir}/plots/ RAWdir=${OPENdir}raw/  

exit
