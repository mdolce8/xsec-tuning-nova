/*
 *  plot_fd_fhc_numu_mc_sloshing.C
 *
 *  macro to plot any FD Numu Quantile spectra.
 *  Created to plot 2020 vs 2024 version of Prod5.1 MC.
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
  const std::string outDir = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/box-opening/plot_fd_fhc_numu_mc_sloshing";


  for (int qCount = 1; qCount <= 5; qCount++){
    const std::string filePath20 = inputDir + "/spectra_fd_prod5.1_p1p10_reco_enu_2020info_mc_fhc_numu_Q" + std::to_string(qCount) + ".root";
    const std::string filePath24 = inputDir + "/spectra_fd_prod5.1_p1p10_reco_enu_2024info_mc_fhc_numu_Q" + std::to_string(qCount) + ".root";

    const std::string specName = Form("pred_interp_Q%i", qCount);
    TFile * f20 = TFile::Open(filePath20.c_str());
    TFile * f24 = TFile::Open(filePath24.c_str());

    Spectrum spec20 = *ana::Spectrum::LoadFrom(f20, specName);
    Spectrum spec24 = *ana::Spectrum::LoadFrom(f24, specName);

    // do plotting
    TCanvas c;

    // we are looking at only p1-10 results. So use 2020 POT.
    TH1D * h20 = spec20.ToTH1(kAna2020FHCPOT);
    TH1D * h24 = spec24.ToTH1(kAna2020FHCPOT);

    h20->Draw("e");
    h24->Draw("e");

    h24->SetLineColor(kRed);

    h20->SetMaximum(h20->GetMaximum() * 1.3);
    h24->SetMaximum(h20->GetMaximum() * 1.3);

    TLatex latex;
    latex.DrawLatexNDC(0.15, 0.8, Form("Neutrino Beam"));
    const std::string qStr = "Q" + std::to_string(qCount);
    TLatex latex2;
    if (qCount < 5) latex2.DrawLatexNDC(.15,.85,(Form("Quantile %i", qCount)));

    h20->SetXTitle("Reconstructed Neutrino Energy (GeV)");
    h20->SetYTitle("Events");

    float evts20 = h20->Integral();
    float evts24 = h24->Integral();

    TLegend leg(0.65, 0.7, 0.9, 0.9);
    leg.AddEntry(h20, Form("2020. Events = %.2f", evts20), "l");
    leg.AddEntry(h24, Form("2024. Events = %.2f", evts24), "l");

    c.SaveAs(Form("%s/plot_fd_fhc_prod5.1_p1p10_reco_enu_numu_mc_sloshing_%s_unnormalized.png", outDir.c_str(), qStr.c_str()));

  } // quantiles




}