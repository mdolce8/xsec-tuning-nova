# Investigating effects of cRPA model on NOvA Nue sample.

From cRPA paper: https://journals.aps.org/prd/supplemental/10.1103/PhysRevD.108.L031301

## May 2024

### M. Dolce

--- 

## Order of Operations:

1. `generate_fd_nue.C` macro 

    * generate the `PredictionNoExtrap` object.
      * produced in the same binning scheme from the paper. (So easy to do step 3.)
    * save the signal, signal + beam, and "all" FD Nue predictions.

2. `plot_fd_nue.C` macro

    * read in the ROOT file and make plots of the signal + beam, and "all" FD Nue.
    * save them to ROOT file

3. `apply_crpa_nue.C` macro

    * take the ratio from the cRPA paper and apply it to the ROOT file from 2). 
    * This shows the impact of the cRPA effect on NOvA FD Nue sample.

NOTE: there is a separate macro for the LowE FHC sample too. (that follows the same order/workflow)

NOTE: it is not explicitly stated, but you should make `Clone()`s of all histograms when you open them from one file and attempt to save them into another. Because you will have segfaults when you attempt to read that hist from the new file. This is a known nuissance from ROOT...