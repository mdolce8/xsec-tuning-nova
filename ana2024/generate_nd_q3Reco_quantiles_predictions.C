/*
 *  generate_nd_q3Reco_quantiles_predictions.C
 *
 *	Generate ND q3 Reco quantiles predictions.
 *	These predictions will be used to produce some
 *	validation plots for the NOvA 2024 analysis.
 *	Specifically:
 *	  -- the "new" RES/DIS systs and
 *	  -- MEC systs
 *	    (shape and Double Gaussian).
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
#include "3FlavorAna/Systs/EnergySysts2020.h"
#include "3FlavorAna/Systs/3FlavorAna2020Systs.h"
#include "3FlavorAna/Systs/3FlavorSystHelper.h"
#include "3FlavorAna/Vars/Binnings.h"
#include "3FlavorAna/Vars/HistAxes.h"

#include "3FlavorAna/NDFit/NDFitHelper.h"
#include "3FlavorAna/NDFit/Samples/UsefulCutsVars.h"
#include "3FlavorAna/NDFit/Systs/LoadSysts.h"

#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoOsc.h"
#include "CAFAna/Systs/MECSysts.h"
#include "CAFAna/Systs/RESSysts.h"
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
    // systString takes either: "resdis", "mecdg", or "mecshape", or "all"
    std::vector<const ana::ISyst *> GetSystematics(const std::string& systString, const bool reducedDG = false) {
      std::vector<const ana::ISyst *> systs_ptrs;
      std::cout << "'systString' is == " << systString << std::endl;

      if (systString == "mecdg"){
        auto mecDGSysts = ana::getMECtuneSystsCorrected_GSFProd5p1();

        if (reducedDG) {
          std::cout << "Removing the MEC DG systs that were omitted from the ND fitting..." << std::endl;
          for (auto &syst: mecDGSysts) {
            for (unsigned int i = 0; i < mecDGSysts.size(); i++) {
              if (mecDGSysts[i]->ShortName() == "MECDoubleGaussEnhSystNorm_2_GSFProd5p1")
                mecDGSysts.erase(mecDGSysts.begin() + i); //rm this element
            }
            for (unsigned int i = 0; i < mecDGSysts.size(); i++) {
              if (mecDGSysts[i]->ShortName() == "MECDoubleGaussEnhSystMeanQ0_2_GSFProd5p1")
                mecDGSysts.erase(mecDGSysts.begin() + i);
            }
            for (unsigned int i = 0; i < mecDGSysts.size(); i++) {
              if (mecDGSysts[i]->ShortName() == "MECDoubleGaussEnhSystCorr_2_GSFProd5p1")
                mecDGSysts.erase(mecDGSysts.begin() + i); //rm this element
            }
          } // Double Gaussian systs
        } // reduced DG

        for (auto &syst: mecDGSysts) systs_ptrs.push_back(syst);
      } // systString == mecdg

      // the MEC Enu Shape, and MEC NP frac (for nu and nubar)
      else if (systString == "mecshape"){
        systs_ptrs.push_back(&kMECInitStateNPFracSyst2020Nu);
        systs_ptrs.push_back(&kMECInitStateNPFracSyst2020AntiNu);
        systs_ptrs.push_back(&kMECShapeSyst2020GSFNu);
        systs_ptrs.push_back(&kMECShapeSyst2020GSFAntiNu);
      }


      // RES/FSI custom systs (ResScale{Delta,Other}, DISHadro{nu,nubar}, RESvpvn{nu,nubar}ratio)
      else if (systString == "resdis") for (auto &syst: ana::NewRESDISSysts()) systs_ptrs.push_back(syst);

      else {std::cerr << "ERROR. Unknown 'systString' = " << systString << std::endl; exit(1);}

      return systs_ptrs;
    }
}

// for snapshot do: cafe -ss

// NOTE: these are reweightable systematics only!

//todo: maybe create a README.txt which lists the `opt` indicating which systematics were used in the predictions?

using namespace ana;

// =====================================================================================================
void generate_nd_q3Reco_quantiles_predictions(const std::string& beam = "fhc", // or "rhc"
                                                const std::string& systString = "",
                                                const bool test = true,
                                                const bool gridSubmission = false, 			// outdir for grid is "."
                                                const bool fillSysts = true            // generate preds w/o systs
)
// =====================================================================================================
{

  std::string outDir;
  if (test)
  outDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/generate_nd_q3Reco_quantiles_predictions_ana2024/test/";
  else {
    outDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/generate_nd_q3Reco_quantiles_predictions_ana2024/";
  }
  std::cout << "Predictions will be made and placed into..." << outDir << std::endl;

  ana::BeamType2020 beamType;
  beam == "fhc" ? beamType = BeamType2020::kFHC : beamType = BeamType2020::kRHC;


  // --------------------------- get systs in place -----------------------------------------
  std::vector<const ISyst *> systs_ptrs;

  // will store all fit parameters
  if (fillSysts) {
    std::cout << "Filling preds with systs..." << std::endl;

    // these are the reweightable systematics only! (no detector systs)
    systs_ptrs = ndfit::GetSystematics(systString);

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
  if (beamType == ana::BeamType2020::kFHC) {
    std::cout << "Using FHC definitions...." << std::endl;
    defNonSwap = "prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_numu2020";

  }
  if (beamType == ana::BeamType2020::kRHC) {
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
  // NOTE: this Ana2024 Cut is in QuantileCuts2020.{cxx,h}. No wonder no one can find anything...
  std::vector<Cut> cutQuantiles = GetNumuEhadFracQuantCuts2024(!(beam == "fhc"));


  // calc values don't really matter for a ND prediction that has no oscillations in it...
  osc::OscCalcPMNSOpt calc;
  ndfit::Calculator2020BestFit(calc);



  // ==== Create Spectra ====
  std::map<std::string, NoOscPredictionGenerator> predGens;
  std::map<std::string, const PredictionInterp*> predInterps;


  // NOTE: using the "standard" NOvA 2020 kRecoQmag.
  HistAxis histaxisKRecoQmag("Reconstructed |#vec{q}| (GeV)", Binning::Simple(40, 0.0, 2.0), kRecoQmag);
  // NOTE: kNumuND2020 should be a fine cut for now. The binning may also be changed too...

// for consistency in LoadNDTopologicalPreds and its file structure
  for (unsigned int quantileIdx = 0; quantileIdx < cutQuantiles.size(); quantileIdx++) {
    predGens.try_emplace(Form("pred_interp_Q%d", quantileIdx+1),
                         NoOscPredictionGenerator(loaders.GetLoader(caf::kNEARDET, Loaders::kMC), histaxisKRecoQmag, kNumu2020ND && cutQuantiles[quantileIdx], kPPFXFluxCVWgt*kXSecCVWgt2020GSFwFSIProd51)); // kNumuE -- reco Enu
  }

  // Create the FD "AllNumu" predinterp as well. It is Q5.
  // jeremy says this is typically quantile 4.... (and 1 is 0), but this is my way.
  predGens.try_emplace(Form("pred_interp_Q%d", (int) cutQuantiles.size()+1),
                       NoOscPredictionGenerator(loaders.GetLoader(caf::kNEARDET, Loaders::kMC), histaxisKRecoQmag, kNumu2020ND, kPPFXFluxCVWgt*kXSecCVWgt2020GSFwFSIProd51));

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
  for (const std::pair<std::string, const PredictionInterp*> predPair : predInterps){

    // create ROOT file.
    std::string fileName = Form("pred_interp_nxp_%s_nd_%s_numu_Q%i.root",  systString.c_str(), beam.c_str(), quantileCount);
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

