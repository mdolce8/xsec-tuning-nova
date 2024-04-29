/*
 *  maketextlistfile_prod51_2020_vs_2024_cuts.C
 *
 *  use MakeTextListFile() to create a text file for
 *  macro to give text files of the events, run, subrun,
 *  etc to debug the 2020 vs 2024 cuts.
 *
 *
 * August 2022
 * M. Dolce
 */

#include <iostream>
#include <3FlavorAna/Cuts/NumuCuts2024.h>


#include "3FlavorAna/Cuts/NumuCuts2020.h"
#include "3FlavorAna/Systs/3FlavorAna2020Systs.h"
#include "3FlavorAna/Vars/NumuVars.h"

#include "CAFAna/Core/EventList.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionXSecTuning.h"
#include "CAFAna/Weights/GenieWeights.h"



using namespace ana;

// NOTE: I have only re-run the Reconstruction Chain for `prodFilePair` = "a"....... FHC 11385, RHC 12536.


using namespace ana;

// =====================================================================================================
void maketextlistfile_prod51_2020_vs_2024_cuts(
											 const std::string& outDir="", 				// two dirs: one for ehadvis, another for its components
                       const std::string& prodFilePair="a" // "a" is the only option at the moment.
)
// =====================================================================================================
{

//  const std::string& outDir = "/nova/ana/users/mdolce/mcmc/plot/validation/maketextlistfile_ehadvis_prod5p1_prod5_data/";

  std::cout << " -- Using MakeTextListFile() for Prod5.1 MC comparing 2020 and 2024 cuts...." << std::endl;
  std::cout << "**********************************************************" << std::endl;
  std::cout << "'outdir' == " << outDir << std::endl;

// -------------------------------------------------------------------------------------------

 // NOTE:
 /*
  * After running Reco-chain, file is here:
  *   -- /nova/ana/users/mdolce/prod5.1-data-revert-CleanUpTrackAlg/fhc-11385-10-file/neardet_r00011385_s10_t00_R20-11-25-prod5.1.revertCleanUpTrackAlg.caf.root
  *   -- /nova/ana/users/mdolce/prod5.1-data-revert-CleanUpTrackAlg/rhc-12536-10-file/neardet_r00012536_s10_t00_R20-11-25-prod5.1.revertCleanUpTrackAlg.caf.root
  */

  // The original definitions:
  // samweb locate-file "def_snapshot <defname> and run_number <runNumber>" (11385 FHC, 12536 RHC)
  // -- ND Prod5.1 Data
  // ----> prod_caf_R20-11-25-prod5.1reco.d_nd_numi_fhc_full_v1_goodruns
  // ----> prod_caf_R20-11-25-prod5.1reco.d_nd_numi_rhc_full_v1_goodruns
  //
  // -- ND Prod5 Data
  // ------> prod_caf_R19-11-18-prod5reco.d.f.h.l_nd_numi_fhc_full_v1_goodruns
  // ------> prod_caf_R19-11-18-prod5reco.g_nd_numi_rhc_full_v1_goodruns




//   final file and Run Nos.
  std::string finalFileP5p1FHC, finalFileP5p1RHC;
  int fhcRunNo, rhcRunNo;

  if (prodFilePair == "a") {
    fhcRunNo = 11385, rhcRunNo = 12536;
    std::cout << "Using FHC Run Number....." << fhcRunNo << std::endl;
    std::string fileP5p1FHC_R11385S10 = "/nova/ana/users/mdolce/prod5.1-data-revert-CleanUpTrackAlg/fhc-11385-10-file/neardet_r00011385_s10_t00_R20-11-25-prod5.1.revertCleanUpTrackAlg.caf.root";
    finalFileP5p1FHC = fileP5p1FHC_R11385S10;

    std::cout << "Using RHC Run Number....." << rhcRunNo << std::endl;
    const std::string &fileP5p1RHC_R12536S10 = "/nova/ana/users/mdolce/prod5.1-data-revert-CleanUpTrackAlg/rhc-12536-10-file/neardet_r00012536_s10_t00_R20-11-25-prod5.1.revertCleanUpTrackAlg.caf.root";
    finalFileP5p1RHC = fileP5p1RHC_R12536S10;
  }


  else if (prodFilePair == "b") {
    fhcRunNo = 13279; rhcRunNo = 13125;
    std::cout << "Using FHC Run Number....." << fhcRunNo << std::endl;
    const std::string &fileP5P1FHC_R13279S18 = "";
    finalFileP5p1FHC = fileP5P1FHC_R13279S18;

    std::cout << "Using RHC Run Number....." << rhcRunNo << std::endl;
    const std::string &fileP5p1RHC_R13125S10 = "";
    finalFileP5p1RHC = fileP5p1RHC_R13125S10;
  }

  else { std::cerr << "Unknown filePair string. " << prodFilePair << std::endl; abort();}


  std::unordered_map<std::string, ana::Cut> map_cuts
          {
                  {"kNumu2020DecafCut", kNumu2020FDDecafCut},
                  {"kNumu2024DecafCut", kNumu2024FDDecafCut},
          };


  std::cout << finalFileP5p1FHC << std::endl;
  std::cout << finalFileP5p1RHC << std::endl;

  // loop through the cuts and add the EHadVis slices...
  for (const auto &edition : {"2020", "2024"}){
    const auto cutName = edition;
    std::cout << "Looping through cut....." << cutName;

    // make the text files for each production and for each topology -- record the EHadVis -- Run -- Subrun -- Event -- Slice

    // July 22. Swap out kNUmuHadVisE for the SR objects that make up kNumuHadVisE = kNumuHadCalE + kNumuHadTrkE
    MakeTextListFile({finalFileP5p1FHC}, {kNumu2020ND && cutName}, {Form("%s/fhc_Prod5p1_revert_CleanUpTrackAlg_%i_%s.txt", outDir.c_str(), fhcRunNo, cutName.c_str())}, {&kNumuHadVisE, &kNumuHadCalE, &kNumuHadTrkE}, {&kEvt}); // {&kRun, &kSubrun, &kEvt, &kSlc}

  } // fhc topologies


  for (const auto &rhcTopoCut : rhcTopologicalCuts){
    const auto topoCut = rhcTopoCut.first;
    const auto topoCutName = rhcTopoCut.second;
    std::cout << "Looping through cut....." << topoCutName;

    // make the text files for each production and for each topology

    // July 22. Swap out kNUmuHadVisE for the SR objects that make up kNumuHadVisE = kNumuHadCalE + kNumuHadTrkE
    MakeTextListFile({finalFileP5p1RHC}, {kNumu2020ND && topoCut}, {Form("%s/rhc_Prod5p1_revert_CleanUpTrackAlg_%i_%s.txt", outDir.c_str(), rhcRunNo, topoCutName.c_str())}, {&kNumuHadVisE, &kNumuHadCalE, &kNumuHadTrkE}, {&kEvt});

  } // rhc topologies


}
