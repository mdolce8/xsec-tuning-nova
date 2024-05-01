/*
 * generate_fd_decaf_cuts_prod5.1_p1p10_reco_enu_2024info_mc.C:
 *    Create Spectrum objects from Prod5.1 MC
 *    in Reco Enu in the FD 2024 Quantiles.
 *    Use the 2020 cuts and vars to compare with the 2024 ones.
 *    No systematics are involved here.
 *
 *    Author: M. Dolce
 *    Date:  April 2024
 *
 */

#include <CAFAna/Prediction/PredictionNoExtrap.h>
#include <CAFAna/Core/Loaders.h>
#include <OscLib/OscCalcPMNSOpt.h>
#include "3FlavorAna/Cuts/QuantileCuts2024.h"
#include "3FlavorAna/Cuts/NumuCuts2024.h"
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
void generate_fd_decaf_cuts_prod51_p1p10_reco_enu_2024info_mc(const std::string& beam,        // fhc or rhc
                                                     const std::string& outDir,      // $data/preds+spectra/ana2024/box-opening ("." for grid: -o $scratch/data )
                                                     const bool gridSubmission = false
)
// =====================================================================================================
{

  // Asimov A. The 2020 best fit.
  auto calc = new osc::OscCalcPMNSOpt();
  calc->SetL(810);
  calc->SetRho(2.84);
  calc->SetDmsq21(7.53e-5);
  calc->SetTh12(asin(sqrt(0.307)));
  calc->SetDmsq32(2.41e-3);
  calc->SetTh23(asin(sqrt(0.57)));
  calc->SetdCP(0.82*M_PI);
  calc->SetTh13(asin(sqrt(2.18e-2)));

  struct Component
  {
      Flavors::Flavors_t flav;
      Current::Current_t curr;
      Sign::Sign_t sign;
  };

  std::map<std::string, Component> flavors = {
          {"numucc", {Flavors::kAllNuMu, Current::kCC, Sign::kBoth}},
  };


  std::cout << "Ana2024 Box Opening........" << std::endl;
  std::cout << "Plotting Prod5.1 FD Numu Quantile(s) MC with 2024 cuts in Reco Enu........" << std::endl;


  // we are only looking at p1-10. no new periods

  std::string defMC_Non, defMC_Flux, defMC_Tau;
  if (beam == "fhc") {
    defMC_Non = "prod_sumdecaf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_numu2020"; // 3F concat -- 150 files
    defMC_Flux = "prod_sumdecaf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_fluxswap_fhc_nova_v08_full_v1_numu2020";
    defMC_Tau = "prod_sumdecaf_R20-11-25-prod5.1reco.j_fd_genie_N1810j0211a_tauswap_fhc_nova_v08_full_v1_numu2020";
  }
  else if (beam == "rhc") {
    defMC_Non = "prod_sumdecaf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_v1_numu2020"; // 3F concat -- 150 files
    defMC_Flux = "prod_sumdecaf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_fluxswap_rhc_nova_v08_full_v1_numu2020";
    defMC_Tau = "prod_sumdecaf_R20-11-25-prod5.1reco.j_fd_genie_N1810j0211a_tauswap_rhc_nova_v08_full_v1_numu2020";
  }
  else {std::cerr << "Unknown 'beam'. exit..." << std::endl; exit(1);}

  Loaders loader;
  loader.SetLoaderPath(defMC_Non, caf::kFARDET, ana::Loaders::kMC, DataSource::kBeam, Loaders::kNonSwap);
  loader.SetLoaderPath(defMC_Flux, caf::kFARDET, ana::Loaders::kMC, DataSource::kBeam, ana::Loaders::kFluxSwap);
//  loader.SetLoaderPath(defMC_Tau, caf::kFARDET, ana::Loaders::kMC, DataSource::kBeam, ana::Loaders::kTauSwap);

  loader.SetSpillCut(kStandardSpillCuts);



  std::map<std::string, const PredictionNoExtrap*> predNxp;


  loader.Go();



  std::string out_dir;
  if (gridSubmission) out_dir = "."; // for grid: " -o /pnfs/nova/scratch/users/mdolce/<path>/<outDir> "
  else {out_dir = outDir;}   // for local
  if ( gSystem->AccessPathName( out_dir.c_str() ) ) gSystem->mkdir( out_dir.c_str(), true );


  // save the spectra to each Quantile ROOT file
  int quantileCount = 1;
  for (const auto &prednxp : predNxp){
    std::string fileName = Form("pred_nxp_fd_prod5.1_p1p10_reco_enu_2024info_mc_%s_numu_Q%i.root",  beam.c_str(), quantileCount);
    const std::string& finalOutDir = out_dir + "/" + fileName;
    TFile ofile(Form("%s", finalOutDir.c_str()), "recreate");

    // save the PredNoExtrap info
    std::cout << "saving PredNxp: " << std::endl;
    prednxp.second->SaveTo(&ofile, prednxp.first);

    // save the Spectrum info.
    std::cout << "saving Spectrum: " << std::endl;
    std::map<std::string, Spectrum> specs_preds;
    for (auto flavor : flavors){
      auto tmp_spec = prednxp.second->PredictComponent(calc, flavor.second.flav, flavor.second.curr, flavor.second.sign);
      std::string key = prednxp.first+"_"+flavor.first+"_all";
      specs_preds.insert({key, tmp_spec});
    } // flavors

    for(auto spec: specs_preds)
      spec.second.SaveTo(&ofile, spec.first.c_str());


    ofile.Close();
    std::cout << "Wrote file: " << finalOutDir << std::endl;

    quantileCount++;
  } // preds

}
