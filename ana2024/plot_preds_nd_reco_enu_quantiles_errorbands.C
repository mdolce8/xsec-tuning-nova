/* plot_preds_nd_reco_enu_quantiles_errorbands.C
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

#include "3FlavorAna/MCMC/MCMC3FShared.h"
#include "3FlavorAna/NDFit/LoadTopoPreds.h"
#include "3FlavorAna/NDFit/Systs/LoadSysts.h"
#include "3FlavorAna/NDFit/InitializeFit.h"
#include "3FlavorAna/NDFit/NDFitHelper.h"
#include "3FlavorAna/NDFit/NDFitEnums.h"
#include "3FlavorAna/Systs/EnergySysts2020.h"

#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Experiment/MultiExperiment.h"

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

using namespace ana;

// =====================================================================================================
void plot_preds_nd_reco_enu_quantiles_errorbands(const std::string& systString)          // "resdis", "mecdg", "mecshape"
// =====================================================================================================
{

  const std::string& outDirPlot = "/nova/ana/users/mdolce/xsec-tuning-nova/plots/ana2024/plot_preds_nd_reco_enu_quantiles_errorbands";

  //load all systs that exist in the preds ROOT file
  NewRESDISSysts();
  getMECtuneSystsCorrected_GSFProd5p1();
  MECsysts();






  // Load the ND Reco Enu Quantile predictions.
  const std::string& inputDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/generate_nd_reco_enu_quantiles_predictions_ana2024/";
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
  const std::string& dataPath = "/nova/ana/users/mdolce/mcmc/data/nd_reco_enu_quantiles/";
  for (std::string beam : {"fhc","rhc"}){
    for (const std::string& q : {"Q1", "Q2", "Q3", "Q4", "Q5"}){
      const std::string& dataFilename = Form("nd_reco_enu_prod5.1_data_%s_numu_%s.root", beam.c_str(), q.c_str());
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


  ///TODO: could add this in later if we want...
//  // load the expts to get the LL for each topology
//  std::vector<std::pair<const std::string, std::unique_ptr<ana::IExperiment>>> expts = ndfit::BuildExperiments(dataSpectra, preds);
//  std::unique_ptr<ana::MultiExperiment> expt = ndfit::BuildMultiExperiment(ndfit::ExptPtrs(expts), acceptedSystsSet, false);

//  // total ChiSq. Value printed at end of macro -- easier to read.
//  double chiSqNominalTotal = expt->ChiSq(calc2020BF.get(), SystShifts::Nominal());
//  double chiSqShiftedTotal = expt->ChiSq(calcRepSample.get(), *shiftsRepSample); // NO DoF
//  int totalBins = 0; // every single bin.
//  int totalMCBins = 0; // all non-zero MC bins (used in the chiSq calculation).
//  // value to use to check to confirm the total expt->ChiSq() = expts.at(sampleType)->ChiSq() for Rep Sample
//  double chiSqRepSampleTotal_Check = 0.;
//
//  //save the ChiSq for each topology into here
//  std::unordered_map<std::string, double> chiSqMCMCMap {};
//  std::unordered_map<std::string, double> chiSqNominalMap {};


  TCanvas c("c","c", 600,600); // 900, 600
  TPad * p1, * p2; //p1 upper, p2 lower


  // set the max and scale factors here at the start.
  double maxFactor = 1. , scaleFactor = 1.;

  int sampleType = 0;

  std::cout << "Plotting the ratio and fitted predictions with data" << std::endl;
  std::cout << "Looping through systematics. Total number of systematics: " << systs.size() << std::endl;

  for (const auto &predBundle : preds) {
    std::cout << predBundle.name << "......." << std::endl;
    const double dataPOT = dataSpectra.at(predBundle.name).POT();

    /// create the error bands
    std::vector<TH1*> up1ShiftEnuReco, dn1ShiftEnuReco;

    // use the systs from each specific topology
    // Use ONLY the systs that were used in the fitting...
    for (const auto &syst : systs) {
      SystShifts pm1SigmaShift;
      pm1SigmaShift.SetShift(syst, +1.);
      TH1 *hUp1ShiftEnuReco = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(dataPOT,
                                                                                                  EExposureType::kPOT,
                                                                                                  kBinDensity);
      std::cout << "Up integral: " << hUp1ShiftEnuReco->Integral() << std::endl;
      up1ShiftEnuReco.push_back(hUp1ShiftEnuReco);

      pm1SigmaShift.SetShift(syst, -1.);
      TH1 *hDn1ShiftEnuReco = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(dataPOT,
                                                                                                  EExposureType::kPOT,
                                                                                                  kBinDensity);
      std::cout << "Down integral: " << hDn1ShiftEnuReco->Integral() << std::endl;
      dn1ShiftEnuReco.push_back(hDn1ShiftEnuReco);


      // retrieve the ChiSq()s for each topological sample.
//    double tmpChiSqMCMC = expts.at(sampleType).second->ChiSq(calcRepSample.get(), *shiftsRepSample); // this is the REAL ChiSq.
//    double tmpChiSqNom  = expts.at(sampleType).second->ChiSq(calc2020BF.get(), SystShifts::Nominal());
//    chiSqMCMCMap.try_emplace(predBundle.name, tmpChiSqMCMC);
//    chiSqNominalMap.try_emplace(predBundle.name, tmpChiSqNom);
//    chiSqRepSampleTotal_Check += tmpChiSqMCMC;
//
//    // calculate total bins for: chi2/D.o.F.
//    TH1 * hPred = predBundle.pred->PredictSyst(calcRepSample.get(), SystShifts::Nominal()).ToTH1(dataPOT);
//    for (int i = 1; i <= hPred->GetNbinsX(); ++i) {
//      if (hPred->GetBinContent(i) > 0.) totalMCBins++;
//      totalBins++;
//    } // i

      c.cd();
      c.Clear();


      const ndfit::Quantiles q = ndfit::visuals::GetQuantileEnum(predBundle.name);
      const std::string quantileString = ndfit::visuals::GetQuantileString(q);
      const std::string beamType = ndfit::visuals::GetHornCurrent(predBundle.name);


      auto *xAxisENu = new TGaxis(0.001, 0.5, 5.0, 0.501, 0., 5.0, 10, "");
      xAxisENu->SetLabelOffset(-0.015); // //	  std::cout << xAxisENu->GetLabelOffset() --> 0.005
      xAxisENu->SetLabelFont(42);




      //pavetext to print out the events for each topology
      TPaveText ptEnuEvents(0.7, 0.60, 0.85, 0.67, "ARC NDC");
      ptEnuEvents.SetFillColor(0);
      ptEnuEvents.SetFillStyle(0);
      ptEnuEvents.SetBorderSize(0);
      ptEnuEvents.SetTextSize(0.032);
      ptEnuEvents.SetTextFont(102);


      // contains all systs
      PlotWithSystErrorBand((IPrediction *&) predBundle.pred, systs, calc2020BF.get(), dataPOT, kGray + 2, kGray);
      c.SaveAs(ndfit::FullFilename(outDirPlot, "profiled_error_bands_plot_" + predBundle.name + ".png").c_str());
      c.Clear();
      // 2D profile




      /// Plot comparison and ratio on save canvas
      SplitCanvas(0.25, p1, p2);
      // Enu-Reco
      //create the histograms for the PlotWithSystErrorBand() function
      std::cout << "Producing Enu-Reco plots for " << predBundle.name << "......" << std::endl;
      TH1 * hEnuData = dataSpectra.at(predBundle.name).ToTH1(dataPOT, EExposureType::kPOT, kBinDensity);
      TH1 * hEnuCVPred = predBundle.pred->PredictSyst(calc2020BF.get(), SystShifts::Nominal()).ToTH1(dataPOT,
                                                                                                    EExposureType::kPOT,
                                                                                                    kBinDensity);


      //get the events BEFORE the re-scaling
      int dataEnuEvents = hEnuData->Integral();
      int nomEnuEvents = hEnuCVPred->Integral();
      std::cout << "Events (data, CV): " << dataEnuEvents << ", " << nomEnuEvents << std::endl;

      hEnuCVPred->SetLineColor(kGray + 2);
      hEnuCVPred->SetLineWidth(3);
//      hEnuCVPred->Scale(scaleFactor);
      hEnuData->SetMarkerColor(kBlack);
      hEnuData->SetMarkerStyle(kFullCircle);
      hEnuData->SetLineWidth(2);
//      hEnuData->Scale(scaleFactor);
//      for (TH1 * hist: up1ShiftEnuReco)
//        hist->Scale(scaleFactor);
//      for (TH1 * hist: dn1ShiftEnuReco)
//        hist->Scale(scaleFactor);


      p1->cd();
      auto EnuErrorBand = PlotWithSystErrorBand(hEnuCVPred, up1ShiftEnuReco, dn1ShiftEnuReco, kGray + 2, kGray);
      hEnuCVPred->Draw("same hist e");
      hEnuData->Draw("same hist p"); // draw as points (to distinguish with data)

      hEnuCVPred->GetYaxis()->SetTitle("Events / GeV");
      hEnuCVPred->GetYaxis()->SetTitleSize(0.036);
      hEnuCVPred->GetYaxis()->SetTitleOffset(1.1);
      hEnuCVPred->SetMaximum(hEnuCVPred->GetMaximum() * 1.8);
      hEnuCVPred->GetXaxis()->SetLabelSize(0.0);
      hEnuCVPred->GetXaxis()->SetTitleSize(0.0);

      TLegend leg(0.45, 0.65, 0.9, 0.9);
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      const std::string errorBands = Form("#pm1#sigma %s", syst->LatexName().c_str());
      const std::string cvPred = "NOvA 2024 MC"; // 2020 tune (MEC+FSI)";
      std::string legCVText;
      leg.AddEntry(hEnuCVPred, cvPred.c_str(), "l");
      up1ShiftEnuReco.at(0)->SetFillColor(kGray);
      up1ShiftEnuReco.at(0)->SetLineColor(kGray);
      leg.AddEntry(up1ShiftEnuReco.at(0), Form("%s", errorBands.c_str()), "f");
      leg.AddEntry(hEnuData, "Prod5.1 Data", "p");
      leg.Draw("same");
      TLatex latex;
      latex.SetTextSize(0.04);
      latex.SetTextAlign(13);
      latex.DrawLatexNDC(.15, .85, (Form("%s", beamType.c_str())));
      latex.DrawLatexNDC(.15, .8, Form("%s", quantileString.c_str()));
      latex.Draw("same");
//    ptEnuEvents.Draw("same");
      Preliminary();
      ndfit::visuals::NeutrinoLabel(ndfit::NeutrinoType::kNumu);
      ndfit::visuals::DetectorLabel(caf::kNEARDET);

      /// Enu-Reco ratio
      p2->cd();
      p2->SetGridy(1);
      TH1 *hEnuUnity = (TH1F *) hEnuCVPred->Clone("hEnuUnity");
      hEnuUnity->Divide(hEnuCVPred);
      TH1 *hEnuDataRatio = (TH1F *) hEnuData->Clone("hEnuDataRatio");
      hEnuDataRatio->Divide(hEnuCVPred);


      ///create the ratios for the error bands
      std::vector<TH1*> up1ShiftEnuRecoRatio = up1ShiftEnuReco;
      std::vector<TH1*> dn1ShiftEnuRecoRatio = dn1ShiftEnuReco;
      for (auto &hist: up1ShiftEnuRecoRatio)
        hist->Divide(hEnuCVPred);
      for (auto &hist: dn1ShiftEnuRecoRatio)
        hist->Divide(hEnuCVPred);


      PlotWithSystErrorBand(hEnuUnity, up1ShiftEnuRecoRatio, dn1ShiftEnuRecoRatio, kGray + 2, kGray);

      hEnuUnity->GetXaxis()->CenterTitle();
      hEnuUnity->GetXaxis()->SetTitleOffset(1.);
      hEnuUnity->GetXaxis()->SetTitleSize(0.045);
      hEnuUnity->SetXTitle("Reco. E_{#nu} (GeV)");
      hEnuUnity->GetYaxis()->CenterTitle();
      hEnuUnity->GetYaxis()->SetRangeUser(0.5, 1.5);
      hEnuUnity->GetYaxis()->SetTitleSize(0.02);
      hEnuUnity->GetYaxis()->SetLabelSize(0.02);
      hEnuUnity->GetYaxis()->SetTitleOffset(1.5);
      hEnuUnity->SetYTitle("#frac{Prod5.1 Data}{NOvA MC}");
      hEnuUnity->GetYaxis()->CenterTitle();
      xAxisENu->Draw("same");

      hEnuDataRatio->Draw("hist same pe");


      ndfit::visuals::DetectorLabel(predBundle.det);
      for (const auto &ext: {".png", ".pdf"}) // ".root"
        c.SaveAs(ndfit::FullFilename(outDirPlot,
                                     "plot_" + predBundle.name + "_" + syst->ShortName() + "_RecoEnu_errorbands" +
                                     ext).c_str());

    } // systs

    // write captions here...


    sampleType++;
  } //predBundle in preds

//
//  std::cout << "Printing X^2 numbers for the FD spectra only...." << std::endl;
//  std::cout << "========== Nominal X^2 =========== (NO DoF): " << std::endl;
//  for (auto const & pair : chiSqNominalMap) std::cout << pair.first << " --> " << pair.second << std::endl;
//  std::cout << "========== MCMC X^2 =========== (NO DoF): " << chiSqShiftedTotal << ". and the check: " << chiSqRepSampleTotal_Check << std::endl;
//  for (auto const & pair : chiSqMCMCMap) std::cout << pair.first << " --> " << pair.second << std::endl;
//  std::cout << "============================================" << std::endl;
//  std::cout << "============================================" << std::endl;
//  std::cout << "Total X^2 (NO DoF) Nominal is: " << chiSqNominalTotal << std::endl;
//  std::cout << "Total X^2 (NO DoF) Rep. Sample is: " << chiSqShiftedTotal << std::endl;
//  std::cout << "============================================" << std::endl;
//  std::cout << "============================================" << std::endl;
//  std::cout << "Total LL for Nominal is: " << chiSqNominalTotal / -2. << std::endl;
//  std::cout << "Total LL for Rep. Sample is: " << chiSqShiftedTotal / -2. << std::endl;
//  // NOTE: very dangerous not to use the function directly, but this is the correct math....
//  std::cout << "============================================" << std::endl;
//  std::cout << "Total bins = " << totalBins << std::endl;
//  std::cout << "Total Non-zero MC bins = " << totalMCBins << std::endl;
//  std::cout << "Total systs = " << totalSysts << std::endl;
//  std::cout << "Total vars = " << totalfitVars << std::endl;
//  const int totalDoF = totalSysts + totalfitVars;
//  std::cout << "Total DoF = " << totalDoF << std::endl;
//  double chiSqDoFNominal = chiSqNominalTotal;
//  chiSqDoFNominal /= (totalMCBins - (double) totalDoF);
//  double chiSqDoFRepSample = chiSqShiftedTotal;
//  chiSqDoFRepSample /= (totalMCBins - (double) totalDoF);
//  std::cout << "============================================" << std::endl;
//  std::cout << "Total X^2/DoF Nominal is: " << chiSqDoFNominal << std::endl;
//  std::cout << "Total X^2/DoF Rep. Sample is: " << chiSqDoFRepSample << std::endl;
//// printout the X^2 per degree of freedom.


}
