/*
 *  validate_new_RES_DIS_systs.C
 *
 *	We want to validate the "new" RES / DIS systematics.
 *	They are now written into NOvARwgt, and so we should
 *	be diligent and check them again.
 *
 *
 *
 *	Feb. 2024
 *	M. Dolce
 *
 *
 */

#include "3FlavorAna/Cuts/NumuCuts2020.h"
#include "3FlavorAna/Cuts/QuantileCuts2020.h"
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
#include "CAFAna/Systs/XSecSystLists.h"
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



// for snapshot do: cafe -ss

// NOTE: these are reweightable systematics only!

// NOTE: These predictions contain ALL Ana2024 xsec systematics.
// -- Moreover, they will be used for a blessing package.

using namespace ana;

// =====================================================================================================
void validate_new_RES_DIS_systs(const std::string& beam = "fhc", // or "rhc"
                                                const bool test = true,
                                                const bool gridSubmission = false 			// outdir for grid is "."
)
// =====================================================================================================
{

  std::string outDir;
  if (test)
    outDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/validate_new_RES_DIS_systs/test/";
  else {
    outDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/validate_new_RES_DIS_systs/";
  }
  std::cout << "Predictions will be made and placed into..." << outDir << std::endl;

  ana::BeamType2020 beamType;
  beam == "fhc" ? beamType = BeamType2020::kFHC : beamType = BeamType2020::kRHC;


  // --------------------------- get systs in place -----------------------------------------
  std::vector<const ISyst *> systs_ptrs;

  // will store all fit parameters
  std::cout << "Filling preds with systs..." << std::endl;

  // these are the reweightable systematics only! (no detector systs)
  systs_ptrs = getAllXsecSysts_2024();

  std::cout << "******* We are filling predictions for  " << systs_ptrs.size() << " systematic parameters *******" << std::endl;
  for (auto &syst: systs_ptrs) std::cout << syst->ShortName() << std::endl;

  std::cout << "number of systs: " << systs_ptrs.size() << std::endl;
  std::cout << "**********************************************************" << std::endl;
//  ------------------------------------------------------------------------------

// 		Definitions:
//   only use nonswap for ND
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



  // calc values don't really matter for a ND prediction that has no oscillations in it...
  osc::OscCalcPMNSOpt calc;
  ndfit::Calculator2020BestFit(calc);

  // ==== Create Spectra ====
  std::map<std::string, NoOscPredictionGenerator> predGens;
  std::map<std::string, const PredictionInterp*> predInterps;


  // NOTE: using the "standard" NOvA 2020 ehadvis.
  HistAxis histaxisEHadVisE("E_{had}^{vis} (GeV)", Binning::Simple(40, 0.0, 0.8), kNumuHadVisE);
  // NOTE: kNumuND2020 should be a fine cut for now. The binning may also be changed too...


  // Create the FD "AllNumu" predinterp as well. It is Q5.
  // jeremy says this is typically quantile 4.... (and 1 is 0), but this is my way.
  predGens.try_emplace("pred_interp_Q5"),
                       NoOscPredictionGenerator(loaders.GetLoader(caf::kNEARDET, Loaders::kMC), histaxisEHadVisE, kNumu2020ND, kPPFXFluxCVWgt*kXSecCVWgt2020GSFProd51);

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
  for (const std::pair<std::string, const PredictionInterp*> predPair : predInterps){

    // create ROOT file.
    std::string fileName = Form("pred_interp_nxp_xsec24_nd_%s_numu_Q5.root", beam.c_str());
    const std::string& finalOutDir = out_dir + "/" + fileName;
    TFile ofile(Form("%s", finalOutDir.c_str()), "recreate");

    predPair.second->SaveTo(&ofile, predPair.first);
    const double pot = predPair.second->Predict(&calc).POT(); //keep this POT, and use `kAna2020{F,R}HCPOT` when plotting elsewhere....
    predPair.second->Predict(&calc).ToTH1(pot)->Write(Form("h1_%s", predPair.first.c_str()));
    std::cout << "saving TH1 with intrinsic POT: " << pot << std::endl;

    ofile.Close();
    std::cout << "Wrote file: " << finalOutDir << std::endl;

  } // preds



  std::cout << "done." << std::endl;

}

