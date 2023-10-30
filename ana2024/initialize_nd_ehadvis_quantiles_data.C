/*
 * initialize_nd_ehadvis_quantiles_data.C:
 *    Create Spectrum objects from Prod5.1 data
 *    in EHadVis in the ND 2024 Quantiles
 *
 *    Author: M. Dolce
 *    Date:  Oct. 2023
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

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"


//  ------ cafe -bq -l <file-number> --stride too?
// NOTE: this uses Prod5.1 !

using namespace ana;

// =====================================================================================================
void initialize_nd_ehadvis_quantiles_data(const std::string& beam,        // fhc or rhc
                                           const bool test = false,
                                           const bool gridSubmission = false
)
// =====================================================================================================
{

  std::cout << "Producing Prod5.1 ND Data Spectra with FD Numu 2024 Quantile cuts on EHadVis........" << std::endl;

  std::string outDir;
  if (test)
    outDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/generate_nd_ehadvis_quantiles_data_ana2024/test/";
  else {
    outDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/generate_nd_ehadvis_quantiles_data_ana2024/";
  }
  std::cout << "Predictions will be made and placed into..." << outDir << std::endl;

  std::string defData;
  if (beam == "fhc") defData = "prod_sumdecaf_R20-11-25-prod5.1reco.d_nd_numi_fhc_full_v1_goodruns_numu2020"; // 3F concat
  else if (beam == "rhc") defData = "prod_sumdecaf_development_nd_numi_rhc_full_v1_goodruns_numu2020"; // 3F concat
  else {std::cerr << "Unknown 'beam'. exit..." << std::endl; exit(1);}

  SpectrumLoader dataLoader(defData);
  dataLoader.SetSpillCut(kStandardSpillCuts);

  std::vector<Cut> cutQuantiles = GetNumuEhadFracQuantCuts2024(beam != "fhc");

  // NOTE: using the "standard" NOvA 2020 ehadvis.
  HistAxis histaxisEHadVisE("E_{had}^{vis} (GeV)", Binning::Simple(40, 0.0, 0.8), kNumuHadVisE);

  std::map<std::string, Spectrum*> spectrumMap;

  // quantiles
  for (unsigned int quantileIdx = 0; quantileIdx < cutQuantiles.size(); quantileIdx++) {
    spectrumMap.try_emplace(Form("pred_interp_Q%d", quantileIdx + 1),
                            new Spectrum(dataLoader, histaxisEHadVisE, kNumu2020ND && cutQuantiles[quantileIdx]));
  }

  // Create the "AllNumu" data as well. It is Q5.
  spectrumMap.try_emplace(Form("pred_interp_Q%d", (int) cutQuantiles.size() + 1),
                          new Spectrum(dataLoader, histaxisEHadVisE, kNumu2020ND));

  dataLoader.Go();



  std::string out_dir;
  if (gridSubmission) out_dir = "."; // for grid: " -o /pnfs/nova/scratch/users/mdolce/<path>/<outDir> "
  else {out_dir = outDir;}   // for local
  if ( gSystem->AccessPathName( out_dir.c_str() ) ) gSystem->mkdir( out_dir.c_str(), true );


  // save the spectra to each Quantile ROOT file
  int quantileCount = 1;
  for (const auto &specPair : spectrumMap){
    std::string fileName = Form("nd_ehadvis_prod5.1_data_%s_numu_Q%i.root",  beam.c_str(), quantileCount);
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
