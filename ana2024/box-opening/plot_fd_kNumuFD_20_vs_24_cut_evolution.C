/*
 *  plot_fd_kNumuFD_20_vs_24_cut_evolution.C
 *
 *  macro to plot inclusive FD Numu spectra
 *  after each round of cuts.
 *  Created to plot 2020 vs 2024 version of Prod5.1 MC.
 *
 *  Apr. 2024
 *  M. Dolce
 */

#include <iostream>
#include <CAFAna/Analysis/Exposures.h>


#include "3FlavorAna/NDFit/NDFitHelper.h"
#include "3FlavorAna/NDFit/InitializeFit.h"
#include "3FlavorAna/Cuts/NumuCuts2020.h"

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
void plot_fd_kNumuFD_20_vs_24_cut_evolution(const std::string& beam)
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
          //{"nuecc", {Flavors::kAllNuE, Current::kCC, Sign::kBoth}},
//          {"beam_nuecc", {Flavors::kNuEToNuE, Current::kCC, Sign::kBoth}},
//          {"app_nuecc", {Flavors::kNuMuToNuE, Current::kCC, Sign::kNu}},
//          {"app_nuebarcc", {Flavors::kNuMuToNuE, Current::kCC, Sign::kAntiNu}},
//          {"numucc", {Flavors::kAllNuMu, Current::kCC, Sign::kBoth}},
          {"nc", {Flavors::kAll, Current::kNC, Sign::kBoth}}
  };

  std::unordered_map<std::string, const ana::Cut> mapCuts
          {
//                  {"kNumuQuality", kNumuQuality}, // made this already with the 2024 macro, they are completely identical.
                  {"kNumuQuality_kNumuContainFD2020", kNumuQuality&&kNumuContainFD2020},
                  {"kNumuQuality_kNumuContainFD2020_kNumu2020PID", kNumuQuality&&kNumuContainFD2020&&kNumu2020PID},
                  {"kNumuQuality_kNumuContainFD2020_kNumu2020PID_kNumu2020CosRej", kNumuQuality&&kNumuContainFD2020&&kNumu2020PID&&kNumu2020CosRej},
                  {"kNumuQuality_kNumuContainFD2020_kNumu2020PID_kNumu2020CosRej_kCosRejVeto", kNumuQuality&&kNumuContainFD2020&&kNumu2020PID&&kNumu2020CosRej&&kCosRejVeto}
          };


  const double POT = 4;

  // dir of the FD Numu Data ROOT files
  const std::string inputDir = "/exp/nova/data/users/mdolce/preds+spectra/ana2024/box-opening";
  const std::string outDir = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/box-opening/plot_fd_kNumuFD_20_vs_24_cut_evolution" + beam;


  for (const auto& cutPair : mapCuts){
    const std::string filename20 = Form("pred_nxp_fd_prod5.1_reco_enu_mc_%s_numu_Q5_cut_%s.root", beam.c_str(), cutPair.first.c_str());
    const std::string filePath20 = inputDir + "/" + filename20 ;

    std::cout << "Quantile " << qCount << "........" << std::endl;

    // these are the "numucc_all" category
    const std::string specName = Form("pred_nxp_Q%i_numucc_all", qCount); // this is a dir.
    TFile * f20 = TFile::Open(filePath20.c_str());
    TFile * f24 = TFile::Open(filePath24.c_str());


//    auto s20 = ana::LoadFrom<ana::Spectrum>(f20, specName);
    Spectrum spec24 = *ana::Spectrum::LoadFrom(f24, specName);
    Spectrum spec20 = *ana::Spectrum::LoadFrom(f20, specName);

    // do plotting
    TCanvas c;

    // we are looking at only p1-10 results. So use 2020 POT.
    TH1D * h20 = spec20.ToTH1(kAna2020FHCPOT);
    TH1D * h24 = spec24.ToTH1(kAna2020FHCPOT);

    h24->Draw("same hist e");
    h20->Draw("same hist e");

    h24->SetLineColor(kRed);

    h20->SetMaximum(h20->GetMaximum() * 1.5);
    h24->SetMaximum(h20->GetMaximum() * 1.5);

    TLatex latex;
    latex.DrawLatexNDC(0.15, 0.8, Form(beam == "fhc" ? "Neutrino Beam" : "AntiNeutrino Beam"));
    const std::string qStr = "Q" + std::to_string(qCount);
    TLatex latex2;
    if (qCount < 5) latex2.DrawLatexNDC(.15,.85,(Form("Quantile %i", qCount)));

    h20->SetXTitle("Reconstructed Neutrino Energy (GeV)");
    h20->SetYTitle("Events");

    float evts20 = h20->Integral();
    float evts24 = h24->Integral();
    std::cout << "2020 version events integral: " << evts20 << std::endl;
    std::cout << "2024 version events integral: " << evts24 << std::endl;

    h20->Scale(0.1, "width");
    h24->Scale(0.1, "width");

    TLegend leg(0.65, 0.7, 0.9, 0.9);
    leg.SetFillStyle(0);
    leg.AddEntry(h20, Form("2020. Events = %.2f", evts20), "l");
    leg.AddEntry(h24, Form("2024. Events = %.2f", evts24), "l");

    leg.Draw("same");

    c.SaveAs(Form("%s/plot_fd_%s_prod5.1_p1p10_reco_enu_numu_mc_sloshing_%s_unnormalized.png", outDir.c_str(),  beam.c_str(), qStr.c_str()));

  } // quantiles




}