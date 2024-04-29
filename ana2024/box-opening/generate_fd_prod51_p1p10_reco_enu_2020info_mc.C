/*
 * generate_fd_prod5.1_p1p10_reco_enu_2020info_mc.C:
 *    Create Spectrum objects from Prod5.1 MC
 *    in Reco Enu in the FD 2020 Quantiles.
 *    Use the 2020 cuts and vars to compare with the 2020 ones.
 *    No systematics are involved here.
 *
 *    Author: M. Dolce
 *    Date:  April 2024
 *
 */

#include "3FlavorAna/Cuts/QuantileCuts2020.h"
#include "3FlavorAna/Cuts/NumuCuts2020.h"
#include "3FlavorAna/NDFit/Samples/UsefulCutsVars.h"
#include "3FlavorAna/Vars/HistAxes.h"
#include "3FlavorAna/Vars/NumuVars.h"

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Weights/XsecTunes.h"
#include "CAFAna/Weights/PPFXWeights.h"

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"


//  ------ cafe -bq -l <file-number> --stride too?
// NOTE: this uses Prod5.1 !

using namespace ana;

// =====================================================================================================
void generate_fd_prod51_p1p10_reco_enu_2020info_mc(const std::string& beam,        // fhc or rhc
                                                     const std::string& outDir,      // $data/preds+spectra/ana2024/box-opening ("." for grid: -o $scratch/data )
                                                     const bool gridSubmission = false
)
// =====================================================================================================
{

  std::cout << "Ana2024 Box Opening........" << std::endl;
  std::cout << "Plotting Prod5.1 FD Numu Quantile(s) MC with 2020 cuts in Reco Enu........" << std::endl;


  // we are only looking at p1-10. no new periods

  std::string defMC;
  if (beam == "fhc") defMC = "prod_sumdecaf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_numu2020"; // 3F concat -- 150 files
  else if (beam == "rhc") defMC = ""; // 3F concat
  else {std::cerr << "Unknown 'beam'. exit..." << std::endl; exit(1);}

  SpectrumLoader loader(defMC);
  loader.SetSpillCut(kStandardSpillCuts);

  std::vector<Cut> cutQuantiles = GetNumuEhadFracQuantCuts2020(beam != "fhc");

  std::map<std::string, Spectrum*> spectrumMap;

  // quantiles
  for (unsigned int quantileIdx = 0; quantileIdx < cutQuantiles.size(); quantileIdx++) {
    spectrumMap.try_emplace(Form("pred_interp_Q%d", quantileIdx + 1),
                        new Spectrum(loader, kNumuCCOptimisedAxis2020, kNumu2020FD && cutQuantiles[quantileIdx], kNoShift, kXSecCVWgt2024 * kPPFXFluxCVWgt));
  }

  spectrumMap.try_emplace("pred_interp_Q5", new Spectrum(loader, kNumuCCOptimisedAxis2020, kNumu2020FD, kNoShift, kXSecCVWgt2024 * kPPFXFluxCVWgt));

  loader.Go();



  std::string out_dir;
  if (gridSubmission) out_dir = "."; // for grid: " -o /pnfs/nova/scratch/users/mdolce/<path>/<outDir> "
  else {out_dir = outDir;}   // for local
  if ( gSystem->AccessPathName( out_dir.c_str() ) ) gSystem->mkdir( out_dir.c_str(), true );


  // save the spectra to each Quantile ROOT file
  int quantileCount = 1;
  for (const auto &specPair : spectrumMap){
    std::string fileName = Form("spectra_fd_prod5.1_p1p10_reco_enu_2020info_mc_%s_numu_Q%i.root",  beam.c_str(), quantileCount);
    const std::string& finalOutDir = out_dir + "/" + fileName;
    TFile ofile(Form("%s", finalOutDir.c_str()), "recreate");

    specPair.second->SaveTo(&ofile, specPair.first);

    const double pot = specPair.second->POT();
    specPair.second->ToTH1(pot)->Write(Form("h1_%s", specPair.first.c_str()));
    std::cout << "saving TH1 with intrinsic POT: " << pot << std::endl;

    ofile.Close();
    std::cout << "Wrote file: " << finalOutDir << std::endl;

    quantileCount++;
  } // preds

}
