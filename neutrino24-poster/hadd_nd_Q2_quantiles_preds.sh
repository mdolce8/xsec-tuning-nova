#!/bin/bash

# hadd_cafana to the ND Enu Reco Numu preds in the quantiles for the residual difference fit.

# to run: 
#   $ cp hadd_nd_Q2_quantiles_preds.sh </path/to/grid/files>
#   $ source /path/to/grid/output/files/hadd_nd_Q2_quantiles_preds.sh <var>"

#NOTE: the FHC and RHC preds are in the same directory, so can do both in this script 

echo "hadd-ing ROOT files"

OUTDIR=$PWD/"hadded-files"
mkdir -v $OUTDIR
# 5 quantiles bc 5 is the total. (even though we only use 4)

for beam in fhc rhc
do
  for quantile in Q1 Q2 Q3 Q4 Q5
  do
      hadd_cafana $OUTDIR/pred_interp_${quantile}_nxp_xsec24_nd_${beam}_trueQ2.root pred_interp_${quantile}_nxp_xsec24_nd_${beam}_trueQ2.*of*.root
  done

  echo "All " $beam " predictions hadded."
done

echo "done."
 

