#!/bin/bash
### Script for gathering RFQAmodel features

### PATHS
SERVER=$(echo $HOSTNAME | cut -d "." -f1)
export SAINT2=/homes/west/SAINT2/SAINT_eleanor/
#export SAINT2=/homes/west/SAINT2/SAINT_laura/
export QA=/data/$SERVER/not-backed-up/west/QualityAssessment/

### DEFAULT VALUES FOR COMMAND LINE ARGUMENTS
RUN_PPV=false
RUN_MAPALIGN=false
RUN_TDB=false
RUN_SAINT2=false
RUN_TMSCORE=false
RUN_LDDT=false
RUN_BEFF=false
PREP_EIGEN=false
RUN_EIGEN=false
RUN_PCONS=false
RUN_PROQ=false
GATHER_EIGEN=false
GATHER_PCONS=false
GATHER_PROQ=false
SCAFFOLD=false
RUN_STMSCORE=false
RUN_SLDDT=false
GATHER_SAMPLED=false
#GATHER_FNAT=false
#GATHER_FLEX=false
#GATHER_END=false
#GATHER_LDDT=false


usage="\n./$(basename "$0") [options] <NAME>\n\n
Gather features for RFQAmodel\n\n
Positional arguments:\n
\t-t | --target <NAME> \tprefix of the input files\n\n
Optional arguments:\n
\t-h | --help\t\tshow this help text\n
\t-d | --dir\t\tdirectory containing models (default: target)\n
\t-c | --chain\t\tchain identifier of native structure (default: final character\n
\t-p | --ppv\t\tcalculate PPV of models\n
\t-m | --mapalign\t\trun mapalign on models\n
\t-pe | --prepeigen\tprepare files for EigenTHREADER\n
\t-s | --saint2\t\tcalculate SAINT2 score for all models (requires config file)\n
\t-e | --eigenthreader\trun EigenTHREADER on models and gather scores\n
\t-p | --pcons\t\trun Pcons for all models and gather scores\n
\t-q | --proq\t\trun ProQ3D for all models and gather scores \n
\t-m | --tmscore\t\tcalculate TM-score for all models\n
\t-l | --lddt\t\tcalculate lddt for all models\n
\t-b | --beff\t\tcalculate Beff for the target (requires TARGET.aln)\n
\t-ge | --gathereigen\tgather EigenTHREADER scores from output files\n
\t-gp | --gatherpcons\tgather Pcons scores from output files\n
\t-gq | --gatherproq\tgather ProQ3D scores from output files\n
\t--scaffold <SEG> <TERM> flag as ScafFold and define segment length and terminal sampled\n
\t--gathersampled\tgather ProQ3D, LDDT and PCons for sampled region only (requires --scaffold)\n\n"

### GET COMMAND LINE ARGUMENTS
while [ ! $# -eq 0 ]; do
  case "$1" in
    --help | -h)
      echo -e $usage
      exit
      ;;
    --target | -t) 
      TARGET=$2
      ## Default decoy directory and chain:
      DECOY_DIR=$2              
      CHAIN_TARGET="${2: -1}"
      shift
      ;;
    --dir | -d)
      DECOY_DIR=$2
      shift
      ;;
    --chain | -c)
      CHAIN_TARGET=$2
      shift
      ;;
    --ppv | -p)
      RUN_PPV=true
      ;;
    --mapalign | -m)
      RUN_MAPALIGN=true
      ;;
    --prepeigen | -pe)
      PREP_EIGEN=true
      ;;
    --saint2 | -s)
      RUN_SAINT2=true
      ;;
    --eigenthreader | -e)
      RUN_EIGEN=true
      GATHER_EIGEN=true
      ;;
    --pcons | -c)
      RUN_PCONS=true
      GATHER_PCONS=true
      ;;
    --proq | -q)
      RUN_PROQ=true
      GATHER_PROQ=true
      ;;
    --tmscore | -m)
      RUN_TMSCORE=true
      ;;
    --stmscore | -sm)
      RUN_STMSCORE=true
      ;;
    --lddt | -l)
      RUN_LDDT=true
      ;;
    --beff | -b)
      RUN_BEFF=true
      ;;
    --gathereigen | -ge)
      GATHER_EIGEN=true
      ;;
    --gatherpcons | -gp)
      GATHER_PCONS=true
      ;;
    --gatherproq | -gq)
      GATHER_PROQ=true
      ;;
    --scaffold)
      SCAFFOLD=true
      SEG=$2
      TERM=$3
      LENGTH=$(tail -n1 $TARGET.fasta.txt | tr -d '\n' | wc -c)
      if [ $TERM = "C" ]; then
        BEGIN=$((SEG+1))
        END=$LENGTH
      elif [ $TERM = "N" ]; then
        BEGIN=1
        END=$((LENGTH-SEG-1))
      else
        echo "Not a valid terminus"
        exit
      fi
      shift
      shift
      ;;
    --gathersampled)
      if [ $SCAFFOLD != true ]; then
        echo "Please provide the segment and terminus information using the --scaffold option."
        exit
      fi
      GATHER_SAMPLED=true
      ;;
    --* | -*)
      echo "Unknown option:" $1
      exit
      ;;
  esac
  shift
done

## Progress bar ##
watcher () {
    PERCENT=$((i/5))
    echo -ne "$i $PERCENT"%"[...]\r"
}


### PPV: Get true contacts and calculate target PPVs ###
if [ $RUN_PPV = true ]; then
  echo Calculating PPVs...
  if [ -f $TARGET.ppvs ] ; then rm $TARGET.ppvs; fi
  ## Get only predicted contacts >= 0.5
  if [ ! -f $TARGET.models.metapsicov.stage1 ];
  then 
    if [ $SCAFFOLD = true ]; then
      awk -v START=$BEGIN -v STOP=$END '{if ($1 >= START && $1 <= STOP && $2>=START && $2 <= STOP && $5 >= 0.5) print $0}' $TARGET.metapsicov.stage1 > $TARGET.models.metapsicov.stage1
    else
      awk '{ if ($5>=0.5) print $0}' $TARGET.metapsicov.stage1 > $TARGET.models.metapsicov.stage1
    fi
  fi
  i=0
  for DECOY in $(cat $TARGET.lst);
  do
    PPV=NA
    if [ ! -f $DECOY_DIR/$DECOY.proxy_fasta ] || [ ! -f $DECOY_DIR/$DECOY.proxy_map ]
    then
      $QA/bin/getcontactsSparse3.py $DECOY_DIR/$DECOY $CHAIN_TARGET 8.0 2>> $TARGET.log
    fi
    PPV=$($QA/bin/convert $DECOY_DIR/$DECOY.proxy_fasta $TARGET.aln $TARGET.models.metapsicov.stage1 $DECOY_DIR/$DECOY.proxy_map 2>&1 | grep "^True Positives" | awk '{print $3}');
    echo $DECOY $PPV >> $TARGET.ppvs
    i=$(($i+1))
    watcher
  done 
echo Done: $TARGET.ppvs
fi


### MAPALIGN ###
if [ $RUN_MAPALIGN = true ]; then
  echo "Calculating map_align..."
  if [ ! -f $TARGET.meta_map ]; then
    $QA/bin/clare_parse_metapsicov $TARGET.metapsicov.stage1 $LENGTH 0.5 > $TARGET.meta_map
  fi
  i=0
  if [ -f $TARGET.mapalign ]; then rm $TARGET.mapalign; fi
  for DECOY in $(cat $TARGET.lst); do
      $QA/bin/parse_proxy_map $DECOY_DIR/$DECOY.proxy_map > $DECOY_DIR/$DECOY.map

      $QA/bin/map_align -a $TARGET.meta_map -b $DECOY_DIR/$DECOY.map > $DECOY_DIR/$DECOY.out_map
      tail -n 1 $DECOY_DIR/$DECOY.out_map > $DECOY_DIR/$DECOY.align

    if [ $(awk '{print NF}' $DECOY_DIR/$DECOY.align) -gt 9 ] 
    then
      LEN=$(awk '{print $8}' $DECOY_DIR/$DECOY.align)
      FIRST=$(awk '{print $9}' $DECOY_DIR/$DECOY.align | awk -v FS=":" '{print $1}')
      FIRST2=$(awk '{print $9}' $DECOY_DIR/$DECOY.align | awk -v FS=":" '{print $2}')
      LAST=$( expr $LEN + $FIRST )
      LAST2=$( expr $LEN + $FIRST2 )
      SCORE=$(cat $DECOY_DIR/$DECOY.align | grep "^MAX" | awk '{print $5,$8}' );
    else
      LEN=0
      FIRST=0
      FIRST2=0
      LAST=0
      LAST2=0
      SCORE="0.0 0"
    fi
    echo $DECOY $SCORE >> $TARGET.mapalign
    i=$((i+1))
    watcher
  done
echo "Done: $TARGET.mapalign"
fi 

### PREPARE EIGENTHREADER ###
if [ $PREP_EIGEN = true ]; then
  echo "Preparing EigenTHREADER files..."
  if [ ! -d $DECOY_DIR/TDB_DIR/ ]
  then 
    mkdir $DECOY_DIR/TDB_DIR
  fi
  TDB_DIR=$DECOY_DIR/TDB_DIR
  i=0
  for DECOY in $(cat $TARGET.lst); do
    if [ ! -f $DECOY_DIR/$DECOY.dssp ]
    then
      $QA/bin/dssp-2.0.4-linux-amd64 -i $DECOY_DIR/$DECOY > $DECOY_DIR/$DECOY.dssp
    fi

    if [ ! -f $DECOY_DIR/$DECOY.pdb ]
    then
      cp  $DECOY_DIR/$DECOY $DECOY_DIR/$DECOY.pdb
    fi

    if [ ! -f $DECOY_DIR/TBD_DIR/$DECOY.eig ]
    then
      $QA/bin/strsum_eigen $DECOY_DIR/$DECOY.pdb $DECOY_DIR/$DECOY.dssp $TDB_DIR/$DECOY.tdb $TDB_DIR/$DECOY.eig
    fi
    i=$((i+1))
    watcher
  done
  echo "Done."
fi

### TM-SCORE ###
if [ $RUN_TMSCORE = true ]; then
  if [ -f $TARGET.tmscores ]; then rm $TARGET.tmscores ; fi
  echo "Calculating TM-scores..."
  i=0
  for DECOY in $(cat $TARGET.lst); do
    if [ -f $TARGET.pdb ]; then
      TM2=$($SAINT2/3rdparty/TMalign $DECOY_DIR/$DECOY $TARGET.pdb | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')
    fi
    if [ -z $TM2 ]; then
      TM2=0.00
    fi
    echo $DECOY $TM2 >> $TARGET.tmscores
    i=$((i+1))
    watcher
  done
 echo "Done: $TARGET.tmscores" 
fi

### sampled TM-SCORE ###
if [ $RUN_STMSCORE = true ]; then
  if [ -f $TARGET.sampledtmscores ]; then rm $TARGET.sampledtmscores ; fi
  echo "Calculating TM-scores for the sampled region..."
  i=0
  for DECOY in $(cat $TARGET.lst); do
    if [ -f sampled_$TARGET.pdb ]; then
      STM=$($SAINT2/3rdparty/TMalign $DECOY_DIR/$DECOY.cut sampled_$TARGET.pdb | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')
    fi
    if [ -z $STM ]; then
      STM=0.00
    fi
    echo $DECOY $STM >> $TARGET.sampledtmscores
    i=$((i+1))
    watcher
  done
 echo "Done: $TARGET.sampledtmscores" 
fi

if [ $RUN_BEFF = true ]; then
  echo "Calculating Beff..."
  $QA/bin/calculate_bf $TARGET.aln | awk '{print $2}' > $TARGET.beff
  echo "Done: $TARGET.beff"
fi


### SAINT2 SCORE ###
# incomplete score files can cause missing values that aren't caught
# so if the file exists (and is not empty), check they contain the right scores
# if they don't, do it again
# If the scores still don't exist, assign "NA"
if [ $RUN_SAINT2 = true ]; then
  echo "Calculating SAINT2 scores..."
  if [ -f $TARGET.saintscores ]; then rm $TARGET.saintscores ; fi
  i=0
  for DECOY in $(cat $TARGET.lst); do 
    SAULO=NA
    COMB=NA
    $SAINT2/bin/saint2 config_"$TARGET"* -- "$DECOY_DIR/$DECOY"  > $DECOY_DIR/"$DECOY"_scores
    if [ -s $DECOY_DIR/"$DECOY"_scores ] 
    then
      SAULO=`cat $DECOY_DIR/"$DECOY"_scores | awk '/^Saulo =/ { print $NF; }'`
      COMB=`cat $DECOY_DIR/"$DECOY"_scores | awk '/^Combined score =/ { print $NF; }'`
    fi

    if [ -z $COMB ]; then
      SAULO="NA"
      COMB="NA"
    fi
    echo $DECOY $SAULO $COMB >> $TARGET.saintscores
    i=$((i+1))
    watcher
  done 
echo "Done: $TARGET.saintscores"
fi 

if [ $RUN_EIGEN = true ]; then
  echo "Running EigenTHREADER..."
  $QA/bin/run_eigen.sh $TARGET ${PWD##*/} > $TARGET.eigenlog
  mkdir -p Eigen_$TARGET 
  mv "$TARGET"_*.model.pdb Eigen_$TARGET/ 
  cp $TARGET.out $TARGET.out.backup;
  $QA/bin/paste_eig_lst.sh $TARGET > $TARGET.out.tmp;
  mv $TARGET.out.tmp $TARGET.out 
  echo "Done."
fi

if [ $RUN_PCONS = true ]; then
  if [ -f $TARGET.pcons_lst ]; then 
    echo "Running Pcons..."
    $QA/bin/run_pcons.sh $TARGET 
    echo "Done."
  else
    echo "No pcons_lst file. Ignoring pcons request."
  fi 
fi

if [ $RUN_PROQ = true ]; then
  echo "Running ProQ3D..."
  $QA/bin/runproq3D.sh $TARGET
  echo "Done." 
fi

if [ $RUN_LDDT = true ]; then
  echo "Calculating LDDT..."
  lddt -c -v0 $(cat $TARGET.pcons_lst) $TARGET.pdb > $TARGET.all_lddts
  grep -e File -e Global $TARGET.all_lddts | sed 's/ //g' | cut -d ":" -f2 | paste -d ' ' - - | cut -d "/" -f2 > $TARGET.lddts
  echo "Done: $TARGET.lddts"
fi

### Parts to include once EigenTHREADER Pcons and PROQ3D have been run separately ###
if [ $GATHER_EIGEN = true ]; then
  echo "Gathering EigenTHREADER scores..."
  if [ -f $TARGET.eigen ]; then rm $TARGET.eigen ; fi
  i=0
  for DECOY in $(cat $TARGET.lst); do
    if [ -f $TARGET.out ]; then
      EIGEN=$(grep -w "$DECOY" $TARGET.out | awk '{print $1}')
    else EIGEN=0.00;
    fi
    echo $DECOY $EIGEN >> $TARGET.eigen
#    i=$((i+1))
#    watcher
  done 
  echo "Done: $TARGET.eigen"
fi
  
if [ $GATHER_PCONS = true ]; then
  echo "Gathering pcons scores..."
  if [ -f $TARGET.pcons ]; then rm $TARGET.pcons ; fi
  i=0
  for DECOY in $(cat $TARGET.lst); do 
    PCONS=$(grep -w "$DECOY" $TARGET.pcons.txt | awk '{print $2}')
    echo $DECOY $PCONS >> $TARGET.pcons
#    i=$((i+1))
#    watcher
  done 
  echo "Done: $TARGET.pcons"
fi

if [ $GATHER_PROQ = true ]; then
  echo "Gathering ProQ3D scores..."
  if [ -f $TARGET.proq3dscores ]; then rm $TARGET.proq3dscores; fi
  i=0
  for DECOY in $(cat $TARGET.lst); do 
    if [ -f ProQ3D_$TARGET/"$DECOY".pdb.proq3.global ]; then
      PROQ=$(tail -n1 ProQ3D_$TARGET/"$DECOY".pdb.proq3.global)
      echo $DECOY $PROQ >> $TARGET.proq3dscores
    else
      echo $DECOY NA >> $TARGET.proq3dscores
    fi
#    i=$((i+1))
#    watcher
  done
  echo "Done: $TARGET.proq3dscores"
fi

if [ $GATHER_SAMPLED = true ]; then 
  echo "Gathering sampled LDDT scores..."
  $QA/bin/get_local_lddts $TARGET.all_lddts $BEGIN $END | cut -d "/" -f2- > $TARGET.local_lddts
  echo "Done: $TARGET.local_lddts"
  echo "Gathering sampled PCons scores..."
  python3 $QA/bin/get_sampled_pcons.py $TARGET $BEGIN $END
  echo "Done: $TARGET.sampledpcons"
  echo "Gathering sampled ProQ3D scores..."
  for DECOY in $(cat $TARGET.lst); do 
    sed "$((BEGIN+1)),$((END+1))!d" ProQ3D_$TARGET/$DECOY.pdb.proq3.local | awk -v DECOY=$DECOY '{PROQ2+=$1; PROQCEN+=$2; PROQFAD+=$3; PROQ3+=$4} END {print DECOY, PROQ2/NR,PROQCEN/NR, PROQFAD/NR, PROQ3/NR}' 
  done > $TARGET.localproq3dscores
  echo "Done: $TARGET.localproq3dscores"
fi

