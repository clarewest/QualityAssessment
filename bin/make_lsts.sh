#!/usr/bin/bash
#for DECOY in $(ls $1 |  cut -d "." -f1-7  | grep -v ".pdb" | sort | uniq); do if [ -f $1/$DECOY ]; then echo $DECOY; fi; done |  head -n500 > $1.lst
sed "s/^/$1\//" $1.lst > $1.pcons_lst
sed 's/$/.pdb/' $1.pcons_lst > $1.proq_lst
