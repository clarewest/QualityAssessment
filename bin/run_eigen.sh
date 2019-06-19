export QA=/data/icarus/not-backed-up/west/QualityAssessment
export TDB_DIR=$QA/data/$2/$1/TDB_DIR/

$QA/bin/eigenthreader -c12 -C0 -t20 -z1250 -m$1 -F$1.fasta.ss2 $1.fasta.txt $1.metapsicov.stage1 $1.out $1.lst
