/*
 *  plot_fd_fhc_numu.C
 *
 *  macro to plot any FD Numu Quantile spectra.
 *
 *  Apr. 2024
 *  M. Dolce
 */

#include <iostream>
#include <CAFAna/Analysis/Exposures.h>


#include "3FlavorAna/NDFit/NDFitHelper.h"
#include "3FlavorAna/NDFit/InitializeFit.h"
#include "3FlavorAna/NDFit/LoadTopoPreds.h"
#include "3FlavorAna/NDFit/Systs/LoadSysts.h"

#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Core/LoadFromFile.h"

#include "OscLib/OscCalcPMNSOpt.h"

#include "TAxis.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "Utilities/rootlogon.C"

namespace histogram
{
    // correct histogram: bin_content/bin_width
    void NormalizeBinContent(TH1 * h1){
      for (int binIdx = 0; binIdx < h1->GetNbinsX(); binIdx++) {
        auto contentNormalized = h1->GetBinContent(binIdx) / h1->GetBinWidth(binIdx);
        h1->SetBinContent(binIdx, contentNormalized);
      }
    }
}

using namespace ana;

// =====================================================================================================
void plot_fd_fhc_numu(
                      const std::string& outDirSuffix = "" // description of the preds used.
)
// =====================================================================================================
{

  const double POT = 4;

  // dir of the FD Numu Data ROOT files
  const std::string inputDir = "/exp/nova/data/users/mdolce/preds+spectra/ana2024/box-opening";
  const std::string outDir = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/box-opening/plot_fd_fhc_numu";

  const std::vector<std::string> filenames
          {
    "spectra_fd_prod5.1_p1p10_reco_enu_2020info_data_fhc_numu_Q1.root",
    "spectra_fd_prod5.1_p1p10_reco_enu_2020info_data_fhc_numu_Q2.root",
    "spectra_fd_prod5.1_p1p10_reco_enu_2020info_data_fhc_numu_Q3.root",
    "spectra_fd_prod5.1_p1p10_reco_enu_2020info_data_fhc_numu_Q4.root",
    "spectra_fd_prod5.1_p1p10_reco_enu_2020info_data_fhc_numu_Q5.root"
          };


  std::unordered_map<std::string, ana::Spectrum> map_spec;

  int qCount = 1;
  for (const std::string & fname : filenames){
    const std::string filePath = inputDir + "/" + fname;
    const std::string specName = Form("pred_interp_Q%i", qCount);
    TFile * f = TFile::Open(filePath.c_str());
    Spectrum spec = *ana::Spectrum::LoadFrom(f, specName);
    map_spec.try_emplace(specName, spec);


    // do plotting
    TCanvas c;

    spec.ToTH1(spec.POT());

    c.SaveAs(Form("%s/%s.png", outDir.c_str(), )); // todo: fix the filename.

    qCount++;
  }




}