#!/bin/sh
python --version
python ../statphase2slowvec.py  \
         --reflat  67.93590   \
         --reflon  25.83491   \
         --phaselistfile Finland_phaselist.txt   \
         --outfile       Finland_ak135_slovecs.txt
