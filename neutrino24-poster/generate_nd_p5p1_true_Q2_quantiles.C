/*
 * generate_nd_p5p1_true_Q2_quantiles.C:
 *    Create PredInterp objects from Prod5.1 MC
 *    in True W and Q^2 ND 2024 Quantiles.
 *    Xsec systematics are involved here.
 *
 *    Author: M. Dolce
 *    Date:  May 2024
 *
 */

#include "3FlavorAna/Cuts/QuantileCuts2024.h"
#include "3FlavorAna/Cuts/NumuCuts2024.h"
#include "3FlavorAna/NDFit/Samples/UsefulCutsVars.h"
#include "3FlavorAna/Vars/HistAxes.h"

#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionNoOsc.h"
#include "CAFAna/Systs/XSecSystLists.h"
#include "CAFAna/Weights/XsecTunes.h"
#include "CAFAna/Weights/PPFXWeights.h"

#include "OscLib/OscCalcPMNSOpt.h"

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"


//  ------ cafe -bq -l <file-number> --stride too?
// NOTE: this uses Prod5.1 !

using namespace ana;

// =====================================================================================================
void generate_nd_p5p1_true_Q2_quantiles(const std::string& beam,        // fhc or rhc
                                          const bool test = true,
                                          const bool gridSubmission = false
)
// =====================================================================================================
{

  std::string outDir;
  if (test)
    outDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/neutrino24-poster/test/";
  else {
    outDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/neutrino24-poster/";
  }
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
  std::vector<const ISyst*> xsecSysts = getAllXsecSysts_2024();


  std::cout << "Ana2024 Box Opening........" << std::endl;
  std::cout << "Plotting Prod5.1 FD Numu Quantile(s) MC with 2024 cuts in Reco Enu........" << std::endl;


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

  Loaders loader;
  loader.SetLoaderPath(defNonSwap, caf::kNEARDET,  Loaders::kMC, kBeam, Loaders::kNonSwap);
  loader.SetSpillCut(kStandardSpillCuts);

  if (defNonSwap.empty()) throw std::runtime_error( "MC SAM Definition is empty" );

  std::vector<Cut> cutQuantiles = GetNumuEhadFracQuantCuts2024(beam != "fhc");

  HistAxis histAxisTrueQ2("True W (GeV)", Binning::Simple(30, 0.0, 1.5), kTrueQ2);


  std::map<std::string, NoOscPredictionGenerator> predGens;
  std::map<std::string, const PredictionInterp*> predInterps;

	if (!test){
		for (unsigned int quantileIdx = 0; quantileIdx < cutQuantiles.size(); quantileIdx++) {
			predGens.try_emplace(Form("pred_interp_Q%d", quantileIdx + 1),
													 NoOscPredictionGenerator(loader.GetLoader(caf::kNEARDET, Loaders::kMC), histAxisTrueQ2,
																										kNumu2024ND && cutQuantiles[quantileIdx],
																										kPPFXFluxCVWgt * kXSecCVWgt2024));
		}
	}

  // Q5 is Inclusive sample.
  predGens.try_emplace(Form("pred_interp_Q%d", (int) cutQuantiles.size()+1),
                       NoOscPredictionGenerator(loader.GetLoader(caf::kNEARDET, Loaders::kMC), histAxisTrueQ2, kNumu2024ND, kPPFXFluxCVWgt * kXSecCVWgt2024));

  for (const auto &predGen : predGens) {
    predInterps.try_emplace(predGen.first,
                            new PredictionInterp(xsecSysts, calc, predGen.second, loader));
  }

  loader.Go();



  std::string out_dir;
  if (gridSubmission) out_dir = "."; // for grid: " -o /pnfs/nova/scratch/users/mdolce/<path>/<outDir> "
  else {out_dir = outDir;}   // for local
  if ( gSystem->AccessPathName( out_dir.c_str() ) ) gSystem->mkdir( out_dir.c_str(), true );


  // save the PredInterps to each Quantile ROOT file
  int quantileCount = 1;
  for (const std::pair<std::string, const PredictionInterp*> predPair : predInterps){

    // create ROOT file.
    std::string fileName = Form("%s_nxp_xsec24_nd_%s_trueQ2.root", predPair.first.c_str(), beam.c_str());
    const std::string& finalOutDir = out_dir + "/" + fileName;
    TFile ofile(Form("%s", finalOutDir.c_str()), "recreate");

    predPair.second->SaveTo(&ofile, predPair.first);
    const double pot = predPair.second->Predict(calc).POT(); //keep this POT, and use `kAna2020{F,R}HCPOT` when plotting elsewhere....
    predPair.second->Predict(calc).ToTH1(pot)->Write(Form("h1_%s", predPair.first.c_str()));
    std::cout << "saving TH1 with intrinsic POT: " << pot << std::endl;

    ofile.Close();
    std::cout << "Wrote file: " << finalOutDir << std::endl;

    quantileCount++;
  } // preds

}
