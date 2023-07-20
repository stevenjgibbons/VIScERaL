#!/bin/sh
for runID in \
   00800 00801  00802 00803  00804 00805  00806 00807  00808 00809 \
   00810 00811  00812 00813  00814 00815  00816 00817  00818 00819 \
   00820 00821  00822 00823  00824 00825  00826 00827  00828 00829 \
   00830 00831  00832 00833  00834 00835  00836 00837  00838 00839 \
   00840 00841  00842 00843  00844 00845  00846 00847  00848 00849 \
   00850 00851  00852 00853  00854 00855  00856 00857  00858 00859 \
   00860 00861  00862 00863  00864 00865  00866 00867  00868 00869 \
   00870 00871  00872 00873  00874 00875  00876 00877  00878 00879 \
   00880 00881  00882 00883  00884 00885  00886 00887  00888 00889 \
   00890 00891  00892 00893  00894 00895  00896 00897  00898 00899
do
  outfile=rellocs_${runID}.txt
  if test ! -r ${outfile}
  then
    ./locate_all_vs_all_variableslovecs_DPRK.sh ${runID} | grep TOTAL
  else
    echo Already calculated file ${outfile}
  fi
done
