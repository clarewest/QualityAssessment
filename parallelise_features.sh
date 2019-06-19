maxjobs=18
parallelize () {
        while [ $# -gt 0 ] ; do
                jobcnt=(`jobs -p`)
                if [ ${#jobcnt[@]} -lt $maxjobs ] ; then
#      ../../bin/run_pcons.sh $1 &  
#                  ../../bin/runproq3D.sh $1 &
              #   bash ../../get_features.sh --target $1 --ppv --mapalign --saint2 --tmscore --gatherpcons --gatherproq --gathereigen > $1.progress &
 #                bash ../../get_features.sh --target $1 --tmscore --lddt --beff  > $1.progress &
                 bash ../../get_features.sh --target $1 --proq > $1.progress &
#                 bash ../../get_features.sh --target $1 --saint2 > $1.progress &
                shift
                fi
        done
        wait
}

LIST=`cat ~/catchupvalidation.txt`
parallelize $LIST

