/* plot_preds_nd_ehadvis_quantiles_errorbands.C
 *
 * preparation for Ana2024.
 *
 * produce some plots of the ND predictions
 * of the "new" RES and DIS systematics,
 * in the inclusive and quantile distributions
 * in 1 sigma error bands.
 *
 * Large chunk taken from: plot_fd_systematic_shifts_from_numu_quantile_predictions.C
 *
 *
 * Oct. 2023.
 * M. Dolce
*/

#include <iostream>
#include <fstream>

#include "3FlavorAna/MCMC/MCMC3FShared.h"
#include "3FlavorAna/NDFit/LoadTopoPreds.h"
#include "3FlavorAna/NDFit/Systs/LoadSysts.h"
#include "3FlavorAna/NDFit/InitializeFit.h"
#include "3FlavorAna/NDFit/NDFitHelper.h"
#include "3FlavorAna/NDFit/NDFitEnums.h"
#include "3FlavorAna/Systs/EnergySysts2020.h"

#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Systs/XSecSystLists.h"

#include "OscLib/OscCalcAnalytic.h"

#include "Utilities/rootlogon.C"

#include "TCanvas.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TPaveText.h"

#include "boost/algorithm/string.hpp"




namespace files
{
    struct PredIdentifier
    {
        ana::Loaders::FluxType horn;
        std::string topology;

        bool operator==(const PredIdentifier & other) const { return other.horn == horn && other.topology == topology; }
    };
    const std::unordered_map<std::string, std::vector<PredIdentifier>> ND_QUANTILE_PREDS
            {
                    {
                            "nd-quantiles",
                            {
                                    {ana::Loaders::kRHC, "Q1"},
                                    {ana::Loaders::kRHC, "Q2"},
                                    {ana::Loaders::kRHC, "Q3"},
                                    {ana::Loaders::kRHC, "Q4"},
                                    {ana::Loaders::kRHC, "Q5"},

                                    {ana::Loaders::kFHC, "Q1"},
                                    {ana::Loaders::kFHC, "Q2"},
                                    {ana::Loaders::kFHC, "Q3"},
                                    {ana::Loaders::kFHC, "Q4"},
                                    {ana::Loaders::kFHC, "Q5"},
                            }
                    },
            };
    const std::unordered_map<std::string, std::string> FILE_PATTERNS
            {
                    // Prod5.1 ND Reco Enu Quantile Cut pred Files
                    {"xsec24", "pred_interp_nxp_%s_nd_%s_numu_%s.root"},
                    {"resdis", "pred_interp_nxp_%s_nd_%s_numu_%s.root"},
                    {"mecdg", "pred_interp_nxp_%s_nd_%s_numu_%s.root"},
                    {"mecshape", "pred_interp_nxp_%s_nd_%s_numu_%s.root"}
            };
}

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

void NeutrinoLabel(const ndfit::NeutrinoType nu, const bool antiParticle = false){
  std::string nuLtx;
  nu == ndfit::NeutrinoType::kNue ? nuLtx = "#nu_{e}" : nuLtx = "#nu_{#mu}";
  if (antiParticle) nuLtx = "#bar{" + nuLtx + "}";
  auto * prelim = new TLatex(.25, .65, Form("%s", nuLtx.c_str())); // top left: below topology, horn
  prelim->SetTextColor(nu == ndfit::NeutrinoType::kNumu ? kBlue - 3 :  kRed - 3);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

using namespace ana;

// TODO: include data flag for these plots... because we want to be able to make both: set of plots with and WITHOUT data
// TODO: next need to make captions too...

// =====================================================================================================
void plot_preds_nd_ehadvis_quantiles_errorbands(const bool plotData = false,
                                                const bool saveCaptions = false)
// =====================================================================================================
{

  const std::string systString = "xsec24"; // --> all xsec24 systs.

  std::string outDirPlot = "/nova/ana/users/mdolce/xsec-tuning-nova/plots/ana2024/plot_preds_nd_ehadvis_quantiles_errorbands";
  if (plotData) outDirPlot += "/data/";

  //load all systs that exist in the preds ROOT file
  NewRESDISSysts();
  getAllXsecSysts_2020_GSF();






  // Load the ND Quantile predictions.
  const std::string& inputDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/generate_nd_ehadvis_quantiles_predictions_ana2024/";
  // code from LoadFDNumuPreds(), but adapted to read in any input directory.
  auto it_topo = files::ND_QUANTILE_PREDS.find("nd-quantiles");
  std::vector<ana::FitPredictions> preds;
  for (const auto &predId: it_topo->second) {
    std::vector <std::pair<std::string, std::string>> sampleFilenamePairs;
    std::vector <std::string> predObjNames;

    std::string hc = predId.horn == ana::Loaders::kFHC ? "fhc" : "rhc";
    sampleFilenamePairs.emplace_back(hc + "_nd",
                                     inputDir + "/" + Form(files::FILE_PATTERNS.at(systString).c_str(),
                                                           systString.c_str(),
                                                           hc.c_str(),
                                                           predId.topology.c_str()));
    predObjNames.emplace_back("pred_interp_" + predId.topology);

    std::cout << "File location: " << inputDir + "/" + Form("%s", Form(files::FILE_PATTERNS.at(systString).c_str(), systString.c_str(), hc.c_str(), predId.topology.c_str())) << std::endl;
    auto pred = ndfit::LoadPreds(sampleFilenamePairs, predObjNames, caf::kNEARDET);
    preds.emplace_back(std::move(pred[0]));
    std::cout << "Added to preds: " << hc << "_nd_pred_interp_" << predId.topology << std::endl;
  } // load preds from LoadFDNumuPreds()




  // Load the ND Reco Enu Quantile data to get the pot.
  std::cout << "Loading the ND Reco Enu Quantile data for the pot...." << std::endl;
  std::map<std::string, ana::Spectrum> dataSpectra;
  const std::string& dataPath = "/nova/ana/users/mdolce/preds+spectra/ana2024/initialize_nd_ehadvis_quantiles_data_ana2024/";
  for (std::string beam : {"fhc","rhc"}){
    for (const std::string& q : {"Q1", "Q2", "Q3", "Q4", "Q5"}){
      const std::string& dataFilename = Form("nd_ehadvis_prod5.1_data_%s_numu_%s.root", beam.c_str(), q.c_str());
      const std::string& filenameStr = dataPath + "/" + dataFilename;
      TFile* infileData= TFile::Open(filenameStr.c_str(), "read");
      const std::string spectraName = Form("%s_nd_pred_interp_%s", beam.c_str(), q.c_str());
      // I am giving this map an R-value, which is a temporary object and has no name.
      dataSpectra.try_emplace(spectraName, *ana::Spectrum::LoadFrom(infileData, Form("pred_interp_%s", q.c_str()))); // this works, ignore the CLion error
    } // quantiles
  } // beam


  auto calc2020BF = std::make_unique<osc::OscCalcAnalytic>();
  ndfit::Calculator2020BestFit(*calc2020BF);

  // Systematics from the Prediction
  auto systs = dynamic_cast<const ana::PredictionInterp *>(preds[0].pred)->GetAllSysts();
  std::sort(systs.begin(),
            systs.end(),
            [](const auto syst1, const auto syst2) {
                return syst1->ShortName() < syst2->ShortName();
            });





  TCanvas c("c","c", 600,600); // 900, 600
  TPad * p1, * p2; //p1 upper, p2 lower


  // set scale factors here.
  const double scaleFactor = 1e6;

  int sampleType = 0;

  std::cout << "Plotting the ratio and fitted predictions with data" << std::endl;
  std::cout << "Looping through systematics. Total number of systematics: " << systs.size() << std::endl;

  for (const auto &predBundle : preds) {
    std::cout << predBundle.name << "......." << std::endl;

    const ndfit::Quantiles q = ndfit::visuals::GetQuantileEnum(predBundle.name);
    const std::string quantileString = ndfit::visuals::GetQuantileString(q);
    const std::string beamType = ndfit::visuals::GetHornCurrent(predBundle.name);

    double POT;
    if (plotData) {
      if (predBundle.name.find("fhc") != std::string::npos)
        POT = kAna2024SensitivityFHCPOT;
      else {
        POT = kAna2024SensitivityRHCPOT;
      }
      std::cout << "Setting POT to MC....." << POT << std::endl;
    }
    else {
      POT = dataSpectra.at(predBundle.name).POT();
      std::cout << "Setting POT to data....." << POT << std::endl;
    }

    /// create the error bands -- one vector for each prediction.
    std::vector<TH1*> up1Shifts, dn1Shifts;

    // use the systs from each specific topology
    // Use ONLY the systs that were used in the fitting...
    for (const auto &syst : systs) {
      std::cout << "Looping through syst....." << syst->ShortName() << std::endl;


      SystShifts pm1SigmaShift;
      pm1SigmaShift.SetShift(syst, +1.);
      TH1 *hUp1 = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(POT,
                                                                                      EExposureType::kPOT,
                                                                                      kBinDensity);
      std::cout << "Up integral: " << hUp1->Integral() << std::endl;
      up1Shifts.push_back(hUp1);

      pm1SigmaShift.SetShift(syst, -1.);
      TH1 *hDn1 = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(POT,
                                                                                      EExposureType::kPOT,
                                                                                      kBinDensity);
      std::cout << "Down integral: " << hDn1->Integral() << std::endl;
      dn1Shifts.push_back(hDn1);

    } // systs to create error bands


      c.cd();
      c.Clear();




      auto * xAxisEHad = new TGaxis(0.001, 0.5, 0.8, 0.501, 0., 0.8, 10, "");
      xAxisEHad->SetLabelOffset(-0.015); // default is 0.005
      xAxisEHad->SetLabelFont(42);
      xAxisEHad->SetTitle("E_{had}^{vis} (GeV)");
      xAxisEHad->SetTitleOffset(1.2);
      xAxisEHad->CenterTitle();
      xAxisEHad->SetTitleFont(42);



      //pavetext to print out the events for each topology
      TPaveText ptEnuEvents(0.7, 0.60, 0.85, 0.67, "ARC NDC");
      ptEnuEvents.SetFillColor(0);
      ptEnuEvents.SetFillStyle(0);
      ptEnuEvents.SetBorderSize(0);
      ptEnuEvents.SetTextSize(0.032);
      ptEnuEvents.SetTextFont(102);


      // contains all systs
      PlotWithSystErrorBand((IPrediction *&) predBundle.pred, systs, calc2020BF.get(), POT, kGray + 2, kGray);
      c.SaveAs(ndfit::FullFilename(outDirPlot, "profiled_error_bands_plot_" + predBundle.name + ".png").c_str());
      c.Clear();
      // 2D profile




      /// Plot comparison and ratio on save canvas
      SplitCanvas(0.25, p1, p2);
      // EHadVis
      //create the histograms for the PlotWithSystErrorBand() function
      std::cout << "Producing EHadVis plots for " << predBundle.name << "......" << std::endl;
      TH1 * hData = dataSpectra.at(predBundle.name).ToTH1(POT, EExposureType::kPOT, kBinDensity);
      TH1 * hCVPred = predBundle.pred->PredictSyst(calc2020BF.get(), SystShifts::Nominal()).ToTH1(POT,
                                                                                                  EExposureType::kPOT,
                                                                                                  kBinDensity);

      // Rescale
      hCVPred->SetLineColor(kGray + 2);
      hCVPred->SetLineWidth(3);
      hCVPred->Scale(scaleFactor);
      for (TH1 * hist : up1Shifts)
        hist->Scale(scaleFactor);
      for (TH1 * hist : dn1Shifts)
        hist->Scale(scaleFactor);

      hData->SetMarkerColor(kBlack);
      hData->SetMarkerStyle(kFullCircle);
      hData->SetLineWidth(2);
      hData->Scale(scaleFactor);


      p1->cd();
      auto ErrorBand = PlotWithSystErrorBand(hCVPred, up1Shifts, dn1Shifts, kGray + 2, kGray);
      hCVPred->Draw("same hist e");
      if (plotData) hData->Draw("same hist p"); // draw data as points
      hCVPred->GetYaxis()->SetTitle("Events / GeV");
      hCVPred->GetYaxis()->SetTitleSize(0.036);
      hCVPred->GetYaxis()->SetTitleOffset(1.1);
      hCVPred->SetMaximum(hCVPred->GetMaximum() * 2.0);
      hCVPred->GetXaxis()->SetLabelSize(0.0);
      hCVPred->GetXaxis()->SetTitleSize(0.0);

      TLegend leg(0.45, 0.65, 0.9, 0.9);
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      const std::string errorBands = "#pm1#sigma xsec24";
      const std::string cvPred = "NOvA 2024 MC";
      std::string legCVText;
      leg.AddEntry(hCVPred, cvPred.c_str(), "l");
      up1Shifts.at(0)->SetFillColor(kGray);
      up1Shifts.at(0)->SetLineColor(kGray);
      leg.AddEntry(up1Shifts.at(0), Form("%s", errorBands.c_str()), "f");
      if (plotData) leg.AddEntry(hData, "Prod5.1 Data", "p");
      leg.Draw("same");
      TLatex latex;
      latex.SetTextSize(0.04);
      latex.SetTextAlign(13);
      latex.DrawLatexNDC(.15, .85, (Form("%s", beamType.c_str())));
      latex.DrawLatexNDC(.15, .8, Form("%s", quantileString.c_str()));
      latex.Draw("same");
//    ptEnuEvents.Draw("same");
      if (plotData) Preliminary();
      else {Simulation();}
      NeutrinoLabel(ndfit::NeutrinoType::kNumu, beamType == "AntiNeutrino Beam");
      ndfit::visuals::DetectorLabel(caf::kNEARDET);

      /// EHadVis ratio
      p2->cd();
      p2->SetGridy(1);
      TH1 *hUnity = (TH1F *) hCVPred->Clone("hEUnity");
      hUnity->Divide(hCVPred);
      TH1 *hDataRatio = (TH1F *) hData->Clone("hDataRatio");
      hDataRatio->Divide(hCVPred);


      ///create the ratios for the error bands
      std::vector<TH1*> up1ShiftsRatio = up1Shifts;
      std::vector<TH1*> dn1ShiftsRatio = dn1Shifts;
      for (auto &hist: up1ShiftsRatio)
        hist->Divide(hCVPred);
      for (auto &hist: dn1ShiftsRatio)
        hist->Divide(hCVPred);


      PlotWithSystErrorBand(hUnity, up1ShiftsRatio, dn1ShiftsRatio, kGray + 2, kGray);

      hUnity->GetXaxis()->CenterTitle();
      hUnity->GetXaxis()->SetTitleOffset(1.);
      hUnity->GetXaxis()->SetTitleSize(0.045);
      hUnity->SetXTitle(""); // set from the TAxis object
      hUnity->GetYaxis()->CenterTitle();
      hUnity->GetYaxis()->SetRangeUser(0.5, 1.5);
      hUnity->GetYaxis()->SetTitleSize(0.02);
      hUnity->GetYaxis()->SetLabelSize(0.02);
      hUnity->GetYaxis()->SetTitleOffset(1.5);
      hUnity->SetYTitle("#frac{Prod5.1 Data}{NOvA MC}");
      hUnity->GetYaxis()->CenterTitle();
      xAxisEHad->Draw("same");

      if (plotData) hDataRatio->Draw("hist same pe");


      ndfit::visuals::DetectorLabel(predBundle.det);
      for (const auto &ext: {".png", ".pdf"}) // ".root"
        c.SaveAs(ndfit::FullFilename(outDirPlot, "plot_" + predBundle.name + "_EHadVis_xsec24_errorbands" + ext).c_str());

    // write captions here...
    if (saveCaptions) {
      std::ofstream ofile;
      ofile.open(ndfit::FullFilename(outDirPlot, "plot_" + predBundle.name + "_EHadVis_xsec24_errorbands.txt").c_str());
      ofile << "Near Detector " << beamType << " Prod5.1 Ana2024 Monte Carlo prediction in the hadronic energy fraction:  " << quantileString
                << ". The variable in this plot is reconstructed hadronic visible energy (in dark grey)."
                   " The light grey band is the 1 sigma error from all NOvA cross-section uncertainties for the Ana2024 analysis."
                   " This includes the new RES and DIS uncertainties from the ND fitting work."
                   " The top distribution is the number of events, and the bottom is the ratio from the MC." << std::endl;
                   if (plotData) {
                     ofile << "The black points are the Prod5.1 data." << std::endl;
                   }
      ofile.close();
    }


    sampleType++;
  } //predBundle in preds



}
