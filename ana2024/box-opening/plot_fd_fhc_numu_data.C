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
void plot_fd_fhc_numu_data(
                      const std::string& outDirSuffix = "" // description of the preds used.
)
// =====================================================================================================
{

  const double POT = 4;

  // dir of the FD Numu Data ROOT files
  const std::string inputDir = "/exp/nova/data/users/mdolce/preds+spectra/ana2024/box-opening";
  const std::string outDir = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/box-opening/plot_fd_fhc_numu";


  for (int qCount = 1; qCount <= 5; qCount++){
    const std::string filePath = inputDir + "/spectra_fd_prod5.1_p1p10_reco_enu_2020info_data_fhc_numu_Q" + std::to_string(qCount) + ".root";
    const std::string specName = Form("pred_interp_Q%i", qCount);
    TFile * f = TFile::Open(filePath.c_str());
    Spectrum spec = *ana::Spectrum::LoadFrom(f, specName);

    // do plotting
    TCanvas c;

    TH1D * h = spec.ToTH1(spec.POT());
    h->Draw("e");
    h->SetMaximum(h->GetMaximum() * 1.3);

    TLatex latex;
    latex.DrawLatexNDC(0.15, 0.8, Form("Neutrino Beam"));
    const std::string qStr = "Q" + std::to_string(qCount);
    TLatex latex2;
    if (qCount < 5) latex2.DrawLatexNDC(.15,.85,(Form("Quantile %i", qCount)));

    h->SetXTitle("Reconstructed Neutrino Energy (GeV)");
    h->SetYTitle("Events / GeV");

    c.SaveAs(Form("%s/plot_fd_fhc_prod5.1_p1p10_reco_enu_2020info_data_fhc_numu_%s_unnormalized.png", outDir.c_str(), qStr.c_str()));

    histogram::NormalizeBinContent(h);
    h->SetMaximum(h->GetMaximum() * 1.3);

    c.SaveAs(Form("%s/plot_fd_fhc_prod5.1_p1p10_reco_enu_2020info_data_fhc_numu_%s_normalized.png", outDir.c_str(), qStr.c_str()));

  }




}