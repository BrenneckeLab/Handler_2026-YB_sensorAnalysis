#!/usr/bin/bash

###############################################################################################
set -u

# Argument = -i input -I input_directory -c chunksize -D blast-database -v
usage() {
  cat <<EOF
  usage: $0 options

  This tool calculates read counts on a template of choice. It reports the
  counts for both, the sense and the antisense strand.

  usage: [PATH]/selectUTRs [options] -F

  OPTIONS:
      -i, --input         hyperlink to the AnnotationPipeline run to analyze
                          please copy the full link from the top page of the analysis
                          (exclusive with -I)

      -I, --input-dir     directory containing fasta files to analyze
                          library names will be extracted from filenames
                          (exclusive with -i)

      -N --name           name of the analysis
      -a, --all           set flag for local processing of all libraries
                          by default only the first one will be analyzed if option -C is set
      -f, --force         delete all old data to start a fresh analysis
                          also sample to plasmid reference file will be deleted
      -h, --help          Show this message
      -C, --computing     set flag for local processing (use only if multiple cores available)
      -D, --debug         set debug mode - does not trigger git commit
EOF
}

INPUT=
INPUTDIR=
analysisNAME=
localALL=N
FORCE=N
COMPUTING=C
DEBUG=N

OPTS=$(getopt -o i:I:N:afhCD --long input:,input-dir:,name:,all,force,help,computing,debug -n 'parse-options' -- "$@")
if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

eval set -- "$OPTS"

while true; do
  case "$1" in
    -i | --input ) INPUT="$2"; shift; shift ;;
    -I | --input-dir ) INPUTDIR="$2"; shift; shift ;;
    -N | --name ) analysisNAME="$2"; shift; shift ;;
    -a | --all ) localALL=Y; shift ;;
    -f | --force ) FORCE=Y; shift ;;
    -h | --help ) usage; exit 1; shift ;;
    -C | --computing ) COMPUTING=L; shift ;;
    -D | --debug ) DEBUG=Y; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done
###################################################################################################
#check if inputs are correct

# Check that -i and -I are mutually exclusive
if [[ -n $INPUT ]] && [[ -n $INPUTDIR ]]; then
  printf "\n\nERROR: Options -i and -I are mutually exclusive. Please specify only one.\n\n"
  exit 1
fi

# Check that at least one input option is provided
if [[ -z $INPUT ]] && [[ -z $INPUTDIR ]]; then
  printf "\n\nplease specify either:
  -i: hyperlink to the AnnotationPipeline run to analyze
  OR
  -I: directory containing fasta files to analyze
  \n\n"
  exit 1
fi

if [[ -z $analysisNAME ]]; then
  printf "\n\nplease specify the name of the analysis\n\n"
  exit 1
fi
###################################################################################################

#hard-coded variables

#location of all singularity images
SINGULARITYdir=""

#base-location of AnnotationPipeline
APutil=""

###################################################################################################
#variable-setup
DATE_OF_DAY=$(date +%F)
FULL_DATE=$(date)

# Conditional parsing based on which input option was used
if [[ -n $INPUT ]]; then
  # Original AnnotationPipeline logic
  echo "Using AnnotationPipeline input mode..."
  
  #AnnotationPipeline analysis folder
  APfolder=$(echo $INPUT | awk '{split($1, X, "/tmp/"); print X[2] }')
  echo $APfolder

  APname=$(echo $INPUT | awk '{n=split($1, X, "/"); print X[n] }')

  #path to tmp-storage
  TMP="${APname}/"

  #base-location of AnnotationPipeline
  APpath=""

  ##new style deposits output directly into AP folder
  OPENdir=${APpath}/${APfolder}/sensor-analysis/${analysisNAME}/
  TMPdir=${TMP}/${analysisNAME}/
  APpath=${APpath}/${APfolder}/

  if [[ $FORCE == Y ]]; then
    rm -rf $OPENdir
    rm -rf $TMPdir
  fi

  mkdir -p $OPENdir
  mkdir -p $TMPdir

  #print warning messages and exit if some problem occurs
  if [[ -z $analysisNAME ]]; then
    echo no name for Analysis given!!!
    exit 1
  fi
  if [[ ! -d ${APpath} ]]; then
    echo AP folder non existent!!!!!
    exit 1
  fi

  nORIG=$(find "${APpath}individual-libraries/" -maxdepth 1 -mindepth 1 -type d | wc -l)
  LIBnames=$(find "${APpath}individual-libraries/" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" | tr '\n' '~' | sed 's/.$//')
  nNAME=$(echo "$LIBnames" | tr '~' '\n' | wc -l)

elif [[ -n $INPUTDIR ]]; then
  # New fasta directory input logic
  echo "Using fasta directory input mode..."
  
  # Check that input directory exists
  if [[ ! -d $INPUTDIR ]]; then
    echo "ERROR: Input directory does not exist: $INPUTDIR"
    exit 1
  fi
  
  # Set variables to maintain compatibility with downstream code
  APfolder=$(basename "$INPUTDIR")
  echo "APfolder set to: $APfolder"
  
  APname=$analysisNAME
  
  #path to tmp-storage
  TMP="${APname}/"
  
  # Set APpath to the input directory for compatibility
  APpath="$INPUTDIR"
  
  OPENdir="${TMP}/${analysisNAME}/"
  TMPdir="${TMP}/${analysisNAME}/"
  
  if [[ $FORCE == Y ]]; then
    rm -rf "$OPENdir"
    rm -rf "$TMPdir"
  fi
  
  mkdir -p "$OPENdir"
  mkdir -p "$TMPdir"
  
  #print warning messages and exit if some problem occurs
  if [[ -z $analysisNAME ]]; then
    echo "No name for Analysis given!!!"
    exit 1
  fi
  
  # Extract library names from fasta filenames in the input directory
  # Assuming fasta files have extensions like .fa, .fasta, .fna, etc.
  LIBnames=$(find "$INPUTDIR" -maxdepth 1 -type f \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fna" \) -printf "%f\n" | sed -E 's/\.(fa|fasta|fna)$//' | tr '\n' '~' | sed 's/.$//')
  nORIG=$(echo "$LIBnames" | tr '~' '\n' | wc -l)
  nNAME=$nORIG
  
  echo "Found $nORIG libraries: $(echo $LIBnames | tr '~' ' ')"
fi

#test if ~ was contained in a filename
if [[ ! $nORIG -eq $nNAME ]]; then
  printf " some filenmae contained ~ which causes problems"
  exit
fi


if [[ ! -s ${OPENdir}library-fasta-pairing.txt  ]]; then
  printf "\n\nplease add the name of the fasta file into the file listed below and copy the fasta files into the fasta-files directory in the parent-dir
(samples can be left empty - these will not get processed)

${OPENdir}library-fasta-pairing.txt
  \n\n"
  echo $LIBnames | tr '~' '\n' | sort > ${OPENdir}library-fasta-pairing.txt
  mkdir -p ${OPENdir}fasta-files
  exit
fi


SWITCH="N"
errNAMES=""
nLIBs=0
for longNAME in $(echo $LIBnames | tr '~' '\t'); do
  TEST=""
  
  TEST=$(grep  "$longNAME" ${OPENdir}library-fasta-pairing.txt | tr ' ' '\t' | cut -f 2)
  if [[ -z $TEST ]]; then
    echo $longNAME
    if [[ -z $errNAMES ]]; then
      errNAMES=${longNAME}
    else
      errNAMES=${errNAMES}~${longNAME}
    fi
    SWITCH=Y
  else
    nLIBs=$(( $nLIBs + 1 ))
  fi
done

#ask user if it is OK to not process all files
if [[ $SWITCH == Y ]]; then
  printf "\n\nthe following sequences are missing or misspelled in the sequence file
${OPENdir}sensor-sequences.fa
  \n"
  echo ${errNAMES} | tr '~' '\n'
  printf "\n\n"
  while true; do
    read -r -p "Do you want to continue with available sensor sequences? [y or n]" yn
    case $yn in
      [Yy])
        break
      ;;
      [Nn])
        printf "please add the  sensor sequences into the  following fasta file and restart
(the fasta header has to correspond to the file-names)
${OPENdir}sensor-sequences.fa
        \n\n"
        
        exit
      ;;
      *) echo "Please answer yes [y] or no [n]." ;;
    esac
  done
  
fi


###################################################################################################
#preset scripts

#determine script-location
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done

#generate final SCRIPTdir variables
SCRIPTdirraw="$(cd -P "$(dirname "$SOURCE")" && pwd)"
SCRIPTdirraw="${SCRIPTdirraw}/"
SCRIPTdir="${SCRIPTdirraw}script-files/"
UTILITYdir="${SCRIPTdirraw}utility-files/"

#move scripts to TMP-directory
rm -rf ${TMPdir}script-files
cp -r ${SCRIPTdir} ${TMPdir}
SCRIPTdir=${TMPdir}script-files/
mv ${SCRIPTdir}main.sh ${SCRIPTdir}${analysisNAME}.sh
chmod 777 ${SCRIPTdir}${analysisNAME}.sh
MAIN=${SCRIPTdir}${analysisNAME}.sh

###################################################################################################
#push git version and commit

cd ${SCRIPTdirraw}

if [[ $DEBUG != Y ]]; then
  
  #ask for commit-message
  while true; do
    read -r -p "Plese specify a commit-message: " msg
    case $msg in
      [Nn]) break ;;
      *)
        commitMESSAGE=$msg
        break
      ;;
    esac
  done
  
  #commit all changes
  git add .
  if [[ -z $commitMESSAGE ]]; then
    git commit -m "automatic commit on submission"
  else
    git commit -m "$commitMESSAGE"
  fi
  #git push --all
fi
commitID=$(git log -1 --pretty=format:"%h")

mkdir -p ${OPENdir}settings/
cat ${SCRIPTdirraw}*.sh |
awk -v RS="#@!@#" '{if (NR==1) print }' >${OPENdir}settings/SETTINGS_${commitID}.log

mkdir -p $SINGULARITYdir
cd $SINGULARITYdir
singularity pull --disable-cache --name APmaster.simg shub://singularity.vbc.ac.at/brenneckelab/ap_master:master
singularity pull --disable-cache --name R.simg shub://singularity.vbc.ac.at/brenneckelab/ap_r:master

###################################################################################################
#submit main-run script
LOG=${OPENdir}/LOGs/
mkdir -p ${LOG}
rm -rf ${LOG}*.txt
cd $LOG


COMMAND="${MAIN}"
VARI="analysisNAME=${analysisNAME},commitID=${commitID},COMPUTING=${COMPUTING},localALL=${localALL},OPENdir=${OPENdir},TMPdir=${TMPdir},SCRIPTdir=${SCRIPTdir},UTILITYdir=${UTILITYdir},LOG=${LOG},SINGULARITYdir=${SINGULARITYdir},APpath=${APpath},APutil=${APutil},LIBnames=${LIBnames},nLIBs=${nLIBs}"

if [[ $COMPUTING == C ]]; then
  sbatch $COMMAND ${VARI}
else
  if [[ -z ${SLURM_CPUS_PER_TASK+x} ]]; then
    srun --cpus-per-task=10 --mem-per-cpu=5g --qos=short $COMMAND ${VARI}
  else
    $COMMAND ${VARI}
  fi
  wait
fi

exit
