/* plot_preds_nd_reco_enu_quantiles_RES_DIS_shifts.C
 *
 * preparation for Ana2024.
 *
 * produce some plots of the ND predictions
 * of the "new" RES and DIS systematics,
 * in the inclusive and quantile distributions.
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
#include "CAFAna/Fit/MCMCSamples.h"

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
void plot_preds_nd_reco_enu_quantiles_RES_DIS_shifts(const std::string& systString)          // "resdis", "mecdg", "mecshape"
// =====================================================================================================
{

  const std::string& outDirPlot = "/nova/ana/users/mdolce/xsec-tuning-nova/plots/ana2024/plot_preds_nd_reco_enu_quantiles_RES_DIS_shifts";

  //load all systs that exist in the preds ROOT file
  NewRESDISSysts();

  auto calc2020BF = std::make_unique<osc::OscCalcAnalytic>();
  ndfit::Calculator2020BestFit(*calc2020BF);






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
                                                           hc.c_str(),
                                                           predId.topology.c_str()));
    predObjNames.emplace_back("pred_interp_" + predId.topology);

    std::cout << "File location: " << inputDir + "/" + Form("%s", Form(files::FILE_PATTERNS.at(systString).c_str(), hc.c_str(), predId.topology.c_str())) << std::endl;
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



  TCanvas c("c","c", 600,650); // this will be a square plot...
  c.cd();
  c.Clear();

  // Loop through the ND Numu Quantile predictions
  for (const ana::FitPredictions &predBundle : preds){
    std::cout << "Looping through " << predBundle.name << std::endl;

    auto hc = ndfit::visuals::GetHornCurrent(predBundle.name);
    const double POT = dataSpectra.at(predBundle.name).POT();
    std::cout << "horn: " << hc  << ", data pot: " << POT << std::endl;
    std::cout << "MC pot is from predBundle.pot = " << predBundle.pot << std::endl;

    // Systematics from the Prediction
    auto systs = dynamic_cast<const ana::PredictionInterp *>(preds[0].pred)->GetAllSysts();
    std::sort(systs.begin(),
              systs.end(),
              [](const auto syst1, const auto syst2) {
                  return syst1->ShortName() < syst2->ShortName();
              });


    // We are bin-wdith normalizing these plots!!
    TH1 * hnom = predBundle.pred->PredictSyst(calc2020BF.get(), SystShifts::Nominal()).ToTH1(POT, EExposureType::kPOT, EBinType::kBinDensity);

    TH1 * hData = dataSpectra.at(predBundle.name).ToTH1(POT, EExposureType::kPOT, EBinType::kBinDensity);





    std::string predBundleName = predBundle.name;
    std::vector<std::string> result;
    boost::split(result, predBundleName, boost::is_any_of("_"));
    const ndfit::Quantiles topologyEnum = ndfit::visuals::GetQuantileEnum(predBundleName);
    const std::string& topologyName = ndfit::visuals::GetQuantileString(topologyEnum);
    const std::string beamType = ndfit::visuals::GetHornCurrent(predBundle.name);


    /// this is crucial to getting the systs and shifting them
    auto pred = dynamic_cast<const ana::PredictionInterp *>(predBundle.pred);

    /// for recording the up and down error bands for all systs (in Ev)
    std::vector<TH1*> histsUp1, histsDn1;

    int systsCount = 0;
    std::cout << "Now looping over systematics......" << systs.size() << std::endl;
    for (const ISyst* syst : systs) {

      /// Save the information for the 1sigma shifts for error bands
      SystShifts pm1SigmaShift;
      // -------------------

      /// Create the hists
      SystShifts shifts;
      shifts.SetShift(syst, +1.0);
      TH1 * up1 = pred->PredictSyst(calc2020BF.get(), shifts).ToTH1(POT, EExposureType::kPOT, EBinType::kBinDensity);
      up1->SetLineColor(kBlue+1);
      histsUp1.push_back(up1);

      shifts.SetShift(syst, +2);
      TH1 * up2 = pred->PredictSyst(calc2020BF.get(), shifts).ToTH1(POT, EExposureType::kPOT, EBinType::kBinDensity);
      up2->SetLineColor(kBlue); up2->SetLineStyle(kDashed);

      shifts.SetShift(syst, +3);
      TH1 * up3 = pred->PredictSyst(calc2020BF.get(), shifts).ToTH1(POT, EExposureType::kPOT, EBinType::kBinDensity);
      up3->SetLineColor(kBlue-4); up3->SetLineStyle(kDotted);


      shifts.SetShift(syst, -1);
      TH1 * down1 = pred->PredictSyst(calc2020BF.get(), shifts).ToTH1(POT, EExposureType::kPOT, EBinType::kBinDensity);
      down1->SetLineColor(kRed+1);
      histsDn1.push_back(down1);

      shifts.SetShift(syst, -2);
      TH1 * down2 = pred->PredictSyst(calc2020BF.get(), shifts).ToTH1(POT, EExposureType::kPOT, EBinType::kBinDensity);
      down2->SetLineColor(kRed); down2->SetLineStyle(kDashed);

      shifts.SetShift(syst, -3);
      TH1 * down3 = pred->PredictSyst(calc2020BF.get(), shifts).ToTH1(POT, EExposureType::kPOT, EBinType::kBinDensity);
      down3->SetLineColor(kRed-4); down3->SetLineStyle(kDotted);
      // -------------------



      double maxContent= 1.6 * std::max(up3->GetMaximum(), down3->GetMaximum());
      hnom->SetMaximum(maxContent);

      double legx1, legx2, legy1, legy2;
      legx1=.6;legx2=0.9; legy1=0.55; legy2=0.85;
      auto *leg = new TLegend(legx1,legy1,legx2,legy2);
      leg->SetFillStyle(0);
      leg->SetLineColor(0);
      leg->AddEntry(hnom, "NOvA N18_10j_00_000 tune","l");
      leg->AddEntry(up3,"+3 #sigma shift","l");
      leg->AddEntry(up2,"+2 #sigma shift","l");
      leg->AddEntry(up1,"+1 #sigma shift","l");
      leg->AddEntry(down1,"-1 #sigma shift","l");
      leg->AddEntry(down2,"-2 #sigma shift","l");
      leg->AddEntry(down3,"-3 #sigma shift","l");

      TH1* ratio_up1 = (TH1*)up1->Clone("ratio_up1");       ratio_up1->Divide(hnom);
      TH1* ratio_up2 = (TH1*)up2->Clone("ratio_up2");       ratio_up2->Divide(hnom);
      TH1* ratio_up3 = (TH1*)up3->Clone("ratio_up3");       ratio_up3->Divide(hnom);
      TH1* ratio_down1 = (TH1*)down1->Clone("ratio_down1"); ratio_down1->Divide(hnom);
      TH1* ratio_down2 = (TH1*)down2->Clone("ratio_down2"); ratio_down2->Divide(hnom);
      TH1* ratio_down3 = (TH1*)down3->Clone("ratio_down3"); ratio_down3->Divide(hnom);

      ratio_up1->SetMinimum(0.5);
      ratio_up1->SetMaximum(1.5);
      ratio_up1->GetYaxis()->SetTitle("Shifted / CV");
      ratio_up1->GetYaxis()->SetTitleFont(43);
      ratio_up1->GetYaxis()->SetTitleOffset(1.4);
      ratio_up1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      ratio_up1->GetYaxis()->SetLabelSize(15);
      ratio_up1->GetYaxis()->SetTitleSize(15);
      ratio_up1->GetXaxis()->SetTitle("Reco. E_{#nu} (GeV)");
      ratio_up1->GetXaxis()->CenterTitle();
      ratio_up1->GetXaxis()->SetTitleFont(43);
      ratio_up1->GetXaxis()->SetTitleOffset(3.4);
      ratio_up1->GetXaxis()->SetNdivisions(510, "X");
      ratio_up1->GetXaxis()->SetTitleOffset(3.5);
      ratio_up1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      ratio_up1->GetXaxis()->SetLabelSize(15);
      ratio_up1->GetXaxis()->SetTitleSize(21);

      TLatex latex;
      latex.SetTextSize(0.04);
      latex.SetTextAlign(13);  //align at top


      // gPad->SetBottomMargin(0.15);
      c.Draw();
      // Upper plot will be in pad1
      TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
      pad1->SetBottomMargin(0); // Upper and lower plot are joined
      pad1->Draw();             // Draw the upper pad: pad1
      c.cd();          // Go back to the main canvas before defining pad2
      TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
      pad2->SetTopMargin(0);
      pad2->SetBottomMargin(0.3);
      pad2->Draw();
      pad1->cd();

      hnom->GetYaxis()->SetTitle("Events");
      hnom->GetYaxis()->SetTitleOffset(.8);
      hnom->GetYaxis()->CenterTitle();
      hnom->Draw("hist");
      up1->Draw("hist same");
      up2->Draw("hist same");
      up3->Draw("hist same");
      down1->Draw("hist same");
      down2->Draw("hist same");
      down3->Draw("hist same");
      c.Update();
      leg->Draw();
      latex.DrawLatexNDC(.15,.85,("#color[1]{"+ndfit::visuals::GetHornCurrent(predBundle.name)+"}").c_str());
      latex.DrawLatexNDC(.15,.8,("#color[1]{"+topologyName+"}").c_str());
      latex.DrawLatexNDC(.15,.75,("#color[1]{"+syst->LatexName()+"}").c_str());
      ndfit::visuals::NeutrinoLabel(ndfit::NeutrinoType::kNumu);
      Simulation();
      ndfit::visuals::DetectorLabel(caf::kFARDET);
      pad1->Update();
      c.Update();
      c.cd();
      pad2->cd();
      pad2->SetGridy();

      ratio_up1->Draw("hist");
      ratio_up2->Draw("hist same");
      ratio_up3->Draw("hist same");
      ratio_down1->Draw("hist same");
      ratio_down2->Draw("hist same");
      ratio_down3->Draw("hist same");

      pad2->Update();
      pad2->Update();
      for (const auto &ext: {"pdf", "png"})
        c.SaveAs(ndfit::FullFilename(outDirPlot, "plot_" + predBundleName + "_" + syst->ShortName() + "_reco_enu_shifts_ratio." + ext).c_str());
      c.Clear();



      /// Draw with error bands
      // --------------------
      c.cd();
      pad1->cd();

      TPaveText ptTopology(0.7, 0.68, 0.85, 0.75, "ARC NDC");
      ptTopology.SetFillStyle(0);
      ptTopology.SetFillColor(0);
      ptTopology.SetBorderSize(0);
      ptTopology.AddText(beamType.c_str());
      ptTopology.AddText(topologyName.c_str());

      hnom->SetLineColor(kGray + 2);
      hnom->Scale(1e-5);
      hnom->SetLineWidth(3);
      hData->SetMarkerColor(kBlack);
      hData->SetMarkerStyle(kFullCircle);
      hData->Scale(1e-5);
      hData->SetLineWidth(2);
      for (TH1* h : histsUp1) h->Scale(1e-5);
      for (TH1* h : histsDn1) h->Scale(1e-5);

      auto plotErrorBand = PlotWithSystErrorBand(hnom, histsUp1, histsDn1, kGray, kGray + 2);
      hData->Draw("same hist p");

      TLegend leg2(0.6, 0.75, 0.9, 0.9);
      leg2.SetFillColor(0);
      leg2.SetFillStyle(0);
      const std::string systLegName = "CV Prediction & 1#sigma: " + syst->ShortName();
      leg2.AddEntry(hnom, Form("%s", systLegName.c_str()), "l");
      leg2.AddEntry(hData, "Prod5.1 Data", "p");

      leg2.Draw("same");
      ptTopology.Draw("same");
      Preliminary();

      // ratio
      pad2->cd();
      pad2->SetGridy(1);
      TH1 * hUnity = (TH1F*) hnom->Clone("hUnity");
      hUnity->Divide(hnom);
      TH1 * hDataRatio = (TH1F*) hData->Clone("hDataRatio");
      for (TH1* h : histsUp1) h->Divide(hnom);
      for (TH1* h : histsDn1) h->Divide(hnom);


      hUnity->GetXaxis()->CenterTitle();
      hUnity->GetXaxis()->SetTitleOffset(1.);
      hUnity->GetXaxis()->SetTitleSize(0.045);
      hUnity->SetXTitle("#vec{q}_{3} Reco (GeV)");
      hUnity->GetYaxis()->CenterTitle();
      hUnity->GetYaxis()->SetRangeUser(0.5, 1.5);
      hUnity->GetYaxis()->SetTitleSize(0.02);
      hUnity->GetYaxis()->SetLabelSize(0.02);
      hUnity->GetYaxis()->SetTitleOffset(1.5);
      hUnity->SetYTitle("#frac{Prod5.1 Data}{MC}");
      hUnity->GetYaxis()->CenterTitle();
//      xAxisq3->Draw("same");


      hDataRatio->SetMarkerColor(kBlack);
      hDataRatio->SetMarkerStyle(kFullCircle);
      hDataRatio->Draw("hist same pe");

      hUnity->Draw("same hist pe");
      hDataRatio->Draw("same hist p");


      for (const auto & ext : {"png", "pdf"})
        c.SaveAs(ndfit::FullFilename(outDirPlot, "plot_" + predBundle.name + "NOvA_errorband_RecoEnu." + ext).c_str());





      systsCount++;
    } // systs


  } // predBundle


}
