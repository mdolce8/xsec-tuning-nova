/*
 *  generate_nd_reco_enu_quantiles_predictions.C
 *
 *	Generate ND Reco. Ev quantiles predictions.
 *	These predictions will be used to produce some
 *	validation plots for the NOvA 2024 analysis.
 *	Specifically, the "new" RES/DIS systs and
 *	MEC systs (shape and Double Gaussian).
 *
 *  Saves each Quantile prediction into its own ROOT file.
 *
 *	Oct. 2023
 *	M. Dolce
 *
 *
 */

#include "3FlavorAna/Cuts/NumuCuts2020.h"
#include "3FlavorAna/Cuts/QuantileCuts2020.h"
#include "3FlavorAna/Systs/DummyRockScaleSyst.h"
#include "3FlavorAna/Systs/EnergySysts2020.h"
#include "3FlavorAna/Systs/3FlavorAna2020Systs.h"
#include "3FlavorAna/Systs/3FlavorSystHelper.h"
#include "3FlavorAna/Systs/GeniePCASyst.h"
#include "3FlavorAna/Vars/Binnings.h"
#include "3FlavorAna/Vars/HistAxes.h"

#include "3FlavorAna/NDFit/NDFitHelper.h"
#include "3FlavorAna/NDFit/Samples/UsefulCutsVars.h"
#include "3FlavorAna/NDFit/Systs/LoadSysts.h"

#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoOsc.h"
#include "CAFAna/Systs/DISSysts.h"
#include "CAFAna/Systs/FSISysts.h"
#include "CAFAna/Systs/RESSysts.h"
#include "CAFAna/Systs/RPASysts.h"
#include "CAFAna/Systs/Systs.h"
#include "CAFAna/Weights/PPFXWeights.h"

#include "OscLib/OscCalcPMNSOpt.h"

#include "TSystem.h"
#include "TFile.h"

using namespace ana;


namespace save
{    void WriteTH1ToFile(TH1 * h1, const std::string& histName, TDirectory * dirName)
    {
      h1->SetName(Form("%s", histName.c_str()));
      h1->SetTitle(Form("%s", histName.c_str()));
      h1->SetDirectory(dirName);
      h1->Write(h1->GetName());
    }
}

namespace ndfit {
    std::vector<const ana::ISyst *> GetNDProd51MCMCDolceSysts() {
      std::vector<const ana::ISyst *> systs_ptrs;

      auto mecTuneReduced = ana::getMECtuneSystsCorrected_GSFProd5p1();
      for (
        auto &syst: mecTuneReduced) {
        for (
                unsigned int i = 0; i < mecTuneReduced.size(); i++) {
          if (mecTuneReduced[i]->ShortName() == "MECDoubleGaussEnhSystNorm_2_GSFProd5p1")
            mecTuneReduced.erase(mecTuneReduced.begin() + i); //rm this element
        }
        for (
                unsigned int i = 0; i < mecTuneReduced.size(); i++) {
          if (mecTuneReduced[i]->ShortName() == "MECDoubleGaussEnhSystMeanQ0_2_GSFProd5p1")
            mecTuneReduced.erase(mecTuneReduced.begin() + i);
        }
        for (
                unsigned int i = 0; i < mecTuneReduced.size(); i++) {
          if (mecTuneReduced[i]->ShortName() == "MECDoubleGaussEnhSystCorr_2_GSFProd5p1")
            mecTuneReduced.erase(mecTuneReduced.begin() + i); //rm this element
        }
      }
      for (auto &syst: mecTuneReduced) systs_ptrs.push_back(syst);

      // the MEC Enu Shape, and MEC NP frac (for nu and nubar)
      for (auto &syst: ana::MECsysts()) systs_ptrs.push_back(syst);

      // QE
      systs_ptrs.push_back(ana::GetGenieKnobSyst(rwgt::fReweightZNormCCQE));
      systs_ptrs.push_back(&ana::kZExpEV1Syst2020);
      systs_ptrs.push_back(&ana::kZExpEV2Syst2020);
      systs_ptrs.push_back(&ana::kZExpEV3Syst2020);
      systs_ptrs.push_back(&ana::kRPACCQEEnhSyst2020);
      systs_ptrs.push_back(&ana::kRPACCQESuppSyst2020);

      // RES
      systs_ptrs.push_back(&ana::kRESLowQ2SuppressionSyst2020);
      systs_ptrs.push_back(ana::GetGenieKnobSyst(rwgt::fReweightMaCCRES));
      systs_ptrs.push_back(ana::GetGenieKnobSyst(rwgt::fReweightMvCCRES));
      systs_ptrs.push_back(ana::GetGenieKnobSyst(rwgt::fReweightTheta_Delta2Npi));


      // DIS
      systs_ptrs.push_back(&ana::kFormZone_2020);
      systs_ptrs.push_back(&ana::kDISvbarnCC3pi_2020);
      systs_ptrs.push_back(&ana::kDISvbarpCC1pi_2020);
      systs_ptrs.push_back(&ana::kDISvbarpCC3pi_2020);
      systs_ptrs.push_back(&ana::kDISvnCC1pi_2020);
      systs_ptrs.push_back(&ana::kDISvnCC2pi_2020);
      systs_ptrs.push_back(&ana::kDISvnCC3pi_2020);
      systs_ptrs.push_back(&ana::kDISvpCC0pi_2020);
      systs_ptrs.push_back(&ana::kDISvpCC2pi_2020);
      systs_ptrs.push_back(&ana::kDISvpCC3pi_2020);



      // RES/FSI custom systs (ResScale{Delta,Other}, DISHadro{nu,nubar}, RESvpvn{nu,nubar}ratio
      for (auto &syst: ana::NewRESDISSysts()) systs_ptrs.push_back(syst);

      // ND + FD PCA Flux 2023 systs
      for (int i = 0; i <= 4; ++i) systs_ptrs.push_back(ana::GetFluxPrincipalsNDFD2023(i));


      //  neutronVisEPrimaries2018.
      for (auto &syst: ana::getNeutronSyst_2020()) systs_ptrs.push_back(syst);

      // FSI
      systs_ptrs.push_back(&ana::khNFSISyst2020_EV1);
      systs_ptrs.push_back(&ana::khNFSISyst2020_MFP);

      // radiative corrections, for nue.
      systs_ptrs.push_back(&ana::kRadCorrNue);
      systs_ptrs.push_back(&ana::kRadCorrNuebar);
      systs_ptrs.push_back(&ana::k2ndClassCurrs);

      // correlated Muon energy scale syst
      systs_ptrs.push_back(&ana::kCorrMuEScaleSyst2020);

//      systs_ptrs.push_back(&ana::kAnaCalibShapeSyst);
//      systs_ptrs.push_back(&ana::kAnaCalibrationSyst);
//      systs_ptrs.push_back(&ana::kAnaCherenkovSyst);
//      systs_ptrs.push_back(&ana::kAnaLightlevelNDSyst);

      return systs_ptrs;
    }
}

// for snapshot do: cafe -ss

// NOTE: these are reweightable systematics only!

//todo: maybe create a README.txt which lists the `opt` indicating which systematics were used in the predictions?

using namespace ana;

// =====================================================================================================
void generate_nd_reco_enu_quantiles_predictions(const std::string& beam = "fhc", // or "rhc"
                                               const std::string& outDir="", 				// local --> $ana/preds+spectra/<outDir> , grid --> mkdir outDir, fill arg with "none"
                                               const bool gridSubmission=false, 			// outdir for grid is "."
                                               const bool fillSysts = true            // generate preds w/o systs
)
// =====================================================================================================
{

  ana::BeamType2020 beamType;
  beam == "fhc" ? beamType = BeamType2020::kFHC : beamType = BeamType2020::kRHC;


  // --------------------------- get systs in place -----------------------------------------
  std::vector<const ISyst *> systs_ptrs;

  // will store all fit parameters
  if (fillSysts) {
    std::cout << "Filling preds with systs..." << std::endl;

    // these are the reweightable systematics only! (no detector systs)
    systs_ptrs = ndfit::GetNDProd51MCMCDolceSysts();

    std::cout << "******* We are filling predictions for  " << systs_ptrs.size() << " systematic parameters *******" << std::endl;
    for (auto &syst: systs_ptrs) std::cout << syst->ShortName() << std::endl;
  } // fillSysts = true

  else {
    std::cout << "Generating predictions with NO systematics." << std::endl;
    std::cout << "number of systs: " << systs_ptrs.size() << std::endl;
  }
  std::cout << "**********************************************************" << std::endl;
//  ------------------------------------------------------------------------------

// 		Definitions:
  // only use nonswap for ND
  std::string defNonSwap;
  if (beam == "fhc") {
    std::cout << "Using FHC definitions...." << std::endl;
    defNonSwap = "prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_numu2020";

  }
  if (beam == "rhc") {
    std::cout << "Using RHC Definitions...." << std::endl;
    defNonSwap = "prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_v1_numu2020prod5.1";
  }

  Loaders loaders;
  loaders.SetLoaderPath(defNonSwap, caf::kNEARDET,  Loaders::kMC, kBeam, Loaders::kNonSwap);
  loaders.SetSpillCut(kStandardSpillCuts);

  if (defNonSwap.empty()) throw std::runtime_error( "MC SAM Definition is empty" );


  //important notes: kCCE = kNumuE2020 (Reco Ev), kTrueE (true Ev).
  // kIsNumuCC is a truth-level cut, so we prob don't want that for our MC.

  // Cuts.
  // quantiles (different for FHC or RHC).
  std::vector<Cut> cutQuantiles = GetNumuEhadFracQuantCuts2020(!(beam == "fhc"));


  // calc values don't really matter for a ND prediction that has no oscillations in it...
  osc::OscCalcPMNSOpt calc;
  ndfit::Calculator2020BestFit(calc);



  // ==== Create Spectra ====
  std::map<std::string, NoOscPredictionGenerator> predGens;
  std::map<std::string, const PredictionInterp*> predInterps;

// for consistency in LoadNDTopologicalPreds and its file structure
  for (unsigned int quantileIdx = 0; quantileIdx < cutQuantiles.size(); quantileIdx++) {
    predGens.try_emplace(Form("pred_interp_Q%d", quantileIdx+1),
                         NoOscPredictionGenerator(loaders.GetLoader(caf::kNEARDET, Loaders::kMC), kNumuCCOptimisedAxis2020, kNumu2020ND && cutQuantiles[quantileIdx], kPPFXFluxCVWgt*kXSecCVWgt2020GSFProd51)); // kNumuE -- reco Enu
  }

  // Create the FD "AllNumu" predinterp as well. It is Q5.
  // jeremy says this is typically quantile 4.... (and 1 is 0), but this is my way.
  predGens.try_emplace(Form("pred_interp_Q%d", (int) cutQuantiles.size()+1),
                       NoOscPredictionGenerator(loaders.GetLoader(caf::kNEARDET, Loaders::kMC), kNumuCCOptimisedAxis2020, kNumu2020ND, kPPFXFluxCVWgt*kXSecCVWgt2020GSFProd51));

  for (const auto &predGen : predGens) {
    predInterps.try_emplace(predGen.first,
                            new PredictionInterp(systs_ptrs, &calc, predGen.second, loaders));
  }




  loaders.Go();


  // save to ROOT file.
  std::string out_dir;
  if (gridSubmission) out_dir = "."; // for grid: " -o /pnfs/nova/scratch/users/mdolce/<path>/<outDir> "
  else {out_dir = outDir;}   // for local
  if ( gSystem->AccessPathName( out_dir.c_str() ) ) gSystem->mkdir( out_dir.c_str(), true );



  // save the PredInterps to each Quantile ROOT file
  int quantileCount = 1;
  for (const auto &predPair : predInterps){

    // create ROOT file.
    std::string fileName = Form("pred_interp_nxp_nd_%s_numu_Q%i.root",  beam.c_str(), quantileCount);
    const std::string& finalOutDir = out_dir + "/" + fileName;
    TFile ofile(Form("%s", finalOutDir.c_str()), "recreate");

    predPair.second->SaveTo(&ofile, predPair.first);
    const double pot = predPair.second->Predict(&calc).POT(); //keep this POT, and use `kAna2020{F,R}HCPOT` when plotting elsewhere....
    predPair.second->Predict(&calc).ToTH1(pot)->Write(Form("h1_%s", predPair.first.c_str()));
    std::cout << "saving TH1 with intrinsic POT: " << pot << std::endl;

    ofile.Close();
    std::cout << "Wrote file: " << finalOutDir << std::endl;

    quantileCount++;
  } // preds



  std::cout << "done." << std::endl;

}

