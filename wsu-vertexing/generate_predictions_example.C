/*
 * generate_predictions_examples.C:
 *    Create PredInterp objects from Prod5.1 MC
 *    as an example to get started with CAFAna.
 *
 *    Author: M. Dolce
 *    Date:  April 2025
 *
 */

#include "3FlavorAna/Cuts/QuantileCuts2024.h"
#include "3FlavorAna/Cuts/NumuCuts2024.h"
#include "3FlavorAna/NDFit/Samples/UsefulCutsVars.h"
#include "3FlavorAna/Vars/HistAxes.h"

#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionNoOsc.h"
#include "CAFAna/Systs/RESSysts.h"
#include "CAFAna/Systs/DISSysts.h"
#include "CAFAna/Weights/XsecTunes.h"
#include "CAFAna/Weights/PPFXWeights.h"

#include "OscLib/OscCalcPMNSOpt.h"

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"


//  ------ cafe -bq -l <file-number> --stride <if desired>

using namespace ana;

// =====================================================================================================
void generate_predictions_example(const std::string& beam,        // fhc or rhc
                                  const bool gridSubmission = false
)
// =====================================================================================================
{
  std::string outDir = "/exp/nova/data/users/mdolce/preds+spectra/wsu-vertexer/testing";
  double pot = -5.;

  std::cout << "Predictions will be made and placed into..." << outDir << std::endl;


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

  // the xsec systs
  std::vector<const ISyst*> xsecSystsSubset {
    &kRESLowQ2SuppressionSyst2020,
    &kRESDeltaScaleSyst,
    GetGenieKnobSyst(rwgt::fReweightZNormCCQE),
    &kDISvpCC0pi_2020,
  };

  std::cout << "Generating some example predictions to prep for WSU Vertexing analysis........" << std::endl;


// 		Definitions:
  // only use nonswap for ND
  std::string defNonSwap;
  if (beam == "fhc") {
    std::cout << "Using FHC definitions...." << std::endl;
    defNonSwap = "prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_numu2020";
    pot = kProd5p1NDFHCPOT;
  }
  if (beam == "rhc") {
    std::cout << "Using RHC Definitions...." << std::endl;
    defNonSwap = "prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_v1_numu2020prod5.1";
    pot = kProd5p1NDRHCPOT;
  }

  // make sure the def isn't empty
  if (defNonSwap.empty()) throw std::runtime_error( "MC SAM Definition is empty" );

  // create a Loader object, this will "load" the files from the definition.
  // Make sure the arguments in SetLoaderPath match your definition!
  Loaders loader;
  loader.SetLoaderPath(defNonSwap, caf::kNEARDET,  Loaders::kMC, kBeam, Loaders::kNonSwap);
  // this is a basic quality cut, not strictly essential, but good practice to use.
  loader.SetSpillCut(kStandardSpillCuts);

  // create our HistAxis object. One for True W and one for Reco W.
  // inside the HistAxis constructor, we create a binning scheme: 30 bins, from 0 - 1.5 GeV.
  HistAxis histAxisTrueW("True W (GeV)", Binning::Simple(30, 0.0, 1.5), kTrueW);
  HistAxis histAxisRecoW("Reco W (GeV)", Binning::Simple(30, 0.0, 1.5), kRecoW);

  // create the prediction "Generator" and the actual prediction, "PredictionInterp"
  std::map<std::string, const PredictionInterp*> predInterps;
  std::map<std::string, NoOscPredictionGenerator> predGens;


  // Create the two Prediction Generators we want
  predGens.try_emplace("pred_interp_TrueW",
                       NoOscPredictionGenerator(loader.GetLoader(caf::kNEARDET, Loaders::kMC), histAxisTrueW, kNumu2024ND, kPPFXFluxCVWgt * kXSecCVWgt2024));

  predGens.try_emplace("pred_interp_RecoW",
                       NoOscPredictionGenerator(loader.GetLoader(caf::kNEARDET, Loaders::kMC), histAxisRecoW, kNumu2024ND, kPPFXFluxCVWgt * kXSecCVWgt2024));

  // Create the actual Predictions
  for (const auto &predGen : predGens) {
    predInterps.try_emplace(predGen.first,
                            new PredictionInterp(xsecSystsSubset, calc, predGen.second, loader));
  }

  // the Go() call fills the spectra.
  // this could take A WHILE, depending on: (1) how many files you use from the def, (2) how many bins you use, (3), how many systematics you are filling.
  loader.Go();


  // now let's save the predictions to a ROOT file so we can use them later.
  std::string out_dir;
  if (gridSubmission) out_dir = "."; // for grid: " -o /pnfs/nova/scratch/users/mdolce/<path>/<outDir> "
  else {out_dir = outDir;}   // for local
  if ( gSystem->AccessPathName( out_dir.c_str() ) ) gSystem->mkdir( out_dir.c_str(), true );  // make out_dir if it doesn't exist.


  // save the PredInterps to each separate ROOT file
  for (const std::pair<std::string, const PredictionInterp*> predPair : predInterps){

    // create ROOT file.
    std::string fileName = Form("%s_nxp_nd_systSubset_%s_.root", predPair.first.c_str(), beam.c_str());
    const std::string finalOutDir = out_dir + "/" + fileName;
    TFile ofile(Form("%s", finalOutDir.c_str()), "recreate");

    // save the prediction to the ROOT file
    predPair.second->SaveTo(&ofile, predPair.first);

    // Not necessary, but we can make the TH1 right here, so let's save that too.
    predPair.second->Predict(calc).ToTH1(pot)->Write(Form("h1_%s", predPair.first.c_str()));
    std::cout << "saving TH1 with POT: " << pot << std::endl;

    ofile.Close();
    std::cout << "Wrote file: " << finalOutDir << std::endl;
  }


}
