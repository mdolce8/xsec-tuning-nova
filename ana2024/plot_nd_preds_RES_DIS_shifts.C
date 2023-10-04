/* .C
 *
 * preparation for Ana2024.
 *
 * produce some plots of the ND predictions
 * of the "new" RES and DIS systematics,
 * in the inclusive and quantile distributions.
 *
 *
 * Oct. 2023.
 * M. Dolce
*/


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
                    {"nd-quantiles", "pred_interp_nxp_nd_%s_numu_%s.root"}
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
void apply_ndfit_rep_sample_to_nd_reco_enu_quantile_predictions()
// =====================================================================================================
{

  const std::string& outDirROOT = "/nova/ana/users/mdolce/mcmc/residual-difference-fit";
  const std::string& outDirPlot = "/nova/ana/users/mdolce/mcmc/plot/residual-difference-fit/apply_ndfit_rep_sample_to_nd_reco_enu_quantile_predictions/";

  //load all systs that exist in the preds ROOT file
  getNDXsecSysts();
  getNeutronSyst_2020();
  getMECtuneSysts();
  NDFDpca();
  DetectorSystsProd5p1(caf::kFARDET);
  getLightSystsND_2020();
  getMECtuneSystsCorrected_GSFProd5p1();
  NewRESDISSysts();

  auto calc2020BF = std::make_unique<osc::OscCalcAnalytic>();
  ndfit::Calculator2020BestFit(*calc2020BF);


  // Systematics from the MCMC Samples
  std::vector<const ana::ISyst*> systsSamples = shiftsRepSample->ActiveSysts();
  std::sort(systsSamples.begin(),
            systsSamples.end(),
            [](const auto syst1, const auto syst2) {
                return syst1->ShortName() < syst2->ShortName();
            });




  // Load the ND Reco Enu Quantile predictions.
  const std::string& inputDir = "/nova/ana/users/mdolce/preds+spectra/Merged/";
  // code from LoadFDNumuPreds(), but adapted to read in any input directory.
  auto it_topo = files::ND_QUANTILE_PREDS.find("nd-quantiles");
  std::vector<ana::FitPredictions> preds;
  for (const auto &predId: it_topo->second) {
    std::vector <std::pair<std::string, std::string>> sampleFilenamePairs;
    std::vector <std::string> predObjNames;

    std::string hc = predId.horn == ana::Loaders::kFHC ? "fhc" : "rhc";
    sampleFilenamePairs.emplace_back(hc + "_nd",
                                     inputDir + "/" + Form(files::FILE_PATTERNS.at("nd-quantiles").c_str(),
                                                           hc.c_str(),
                                                           predId.topology.c_str()));
    predObjNames.emplace_back("pred_interp_" + predId.topology);

    std::cout << "File location: " << inputDir + "/" + Form("%s", Form(files::FILE_PATTERNS.at("nd-quantiles").c_str(), hc.c_str(), predId.topology.c_str())) << std::endl;
    auto pred = ndfit::LoadPreds(sampleFilenamePairs, predObjNames, caf::kNEARDET);
    preds.emplace_back(std::move(pred[0]));
    std::cout << "Added to preds: " << hc << "_nd_pred_interp_" << predId.topology << std::endl;
  } // load preds from LoadFDNumuPreds()

  // Systematics from the Predictions
  auto systsPreds = dynamic_cast<const ana::PredictionInterp *>(preds[0].pred)->GetAllSysts();
  std::sort(systsPreds.begin(),
            systsPreds.end(),
            [](const auto syst1, const auto syst2) {
                return syst1->ShortName() < syst2->ShortName();
            });


// NOTE: the only difference between systematics is: Light_Level_FD.
// -- the merged preds do, but the ND Fit samples do NOT.
// dated March 24, 2023.



  // Load the ND Reco Enu Quantile data to get the POT.
  std::cout << "Loading the ND Reco Enu Quantile data for the POT...." << std::endl;
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


  std::map<std::string, TH1*> mapHistsRepSample;
  std::map<std::string, TH1*> mapHistsRepSample_NoMuEScaleSyst;
  std::map<std::string, TH1*> mapHistsNominal;

  // Loop through the ND Numu Quantile predictions
  // and apply the Rep. Sample pulls from the ND Topological Fit.
  // and then save to a ROOT file.
  TFile ofile(Form("%s/NDFit_RepSamplePulls_AppliedNDRecoENuQuantilePreds.root", outDirROOT.c_str()), "recreate");
  for (const ana::FitPredictions &predBundle : preds){
    std::cout << "Looping through " << predBundle.name << std::endl;

    auto hc = ndfit::visuals::GetHornCurrent(predBundle.name);
    const double pot = dataSpectra.at(predBundle.name).POT();
    std::cout << "horn: " << hc  << ", data POT: " << pot << std::endl;

    // NOTE: these histograms NOT normalized!!!
    TH1 * hQuantileRepSample = predBundle.pred->PredictSyst(static_cast<osc::IOscCalc*>(nullptr), *shiftsRepSample).ToTH1(pot, EExposureType::kPOT, kBinContent);
    TH1 * hQuantileRepSample_noCorrMuEScaleSyst = predBundle.pred->PredictSyst(static_cast<osc::IOscCalc*>(nullptr), *shiftsRepSample_NoCorrMuEScaleSyst).ToTH1(pot, EExposureType::kPOT, kBinContent);
    TH1 * hNom = predBundle.pred->Predict(static_cast<osc::IOscCalc*>(nullptr)).ToTH1(pot, EExposureType::kPOT, kBinContent);

    mapHistsRepSample.try_emplace(predBundle.name, hQuantileRepSample);
    mapHistsRepSample_NoMuEScaleSyst.try_emplace(predBundle.name, hQuantileRepSample_noCorrMuEScaleSyst);
    mapHistsNominal.try_emplace(predBundle.name, hNom);

    hQuantileRepSample->Write(Form("%s_RepSample", predBundle.name.c_str()));
    hQuantileRepSample_noCorrMuEScaleSyst->Write(Form("%s_RepSample_NoCorrMuEScaleSyst", predBundle.name.c_str()));

  } // predBundle

  // Crucially, save the Rep. Sample pulls to the ROOT file for easy access.
  ana::SaveTo(static_cast<SystShifts>(*shiftsRepSample), &ofile, "NDFit_systRepSample");
  ana::SaveTo(static_cast<SystShifts>(*shiftsRepSample_NoCorrMuEScaleSyst), &ofile, "NDFit_systRepSample_NoCorrMuEScaleSyst");
  ofile.Close();
  std::cout << "ROOT file created: " << ofile.GetName() << std::endl;








  // Next is to create simple plots of the Nominal and Rep. Sample
  // ND Reco Enu Quantile predictions...
  auto * xAxisENu = new TGaxis(0.001, 0.5, 5.0, 0.501, 0., 5.0, 10, "");
  xAxisENu->SetLabelOffset(-0.015);
  xAxisENu->SetLabelFont(42);
  for (const auto &pairHistRepSample : mapHistsRepSample){
    TH1* hRepSample = pairHistRepSample.second;
    TH1* hRepSample_noCorrMuEScaleSyst = mapHistsRepSample_NoMuEScaleSyst.at(pairHistRepSample.first);
    TH1* hNominal = mapHistsNominal.at(pairHistRepSample.first);

    // Normalize all the histograms here -- for drawing ONLY.
    for (const auto &h : {hNominal, hRepSample, hRepSample_noCorrMuEScaleSyst}) histogram::NormalizeBinContent(h);

    const std::string& predName = pairHistRepSample.first;
    
    const ndfit::Quantiles q = ndfit::visuals::GetQuantileEnum(predName);
    const std::string quantileString = ndfit::visuals::GetQuantileString(q);
    const std::string beamType = ndfit::visuals::GetHornCurrent(predName);

    //pavetext to print out the events for each topology
    TPaveText ptEnuEvents(0.7, 0.60, 0.85, 0.67, "ARC NDC");

    TCanvas c("c", "c", 600, 600); // 900, 600
    TPad *p1, *p2; //p1 upper, p2 lower

    // Plot comparison and ratio on save canvas
    SplitCanvas(0.25, p1, p2);

    ndfit::visuals::PredPreDrawAesthetics(hRepSample, 1e-4, true);
    ndfit::visuals::PredPreDrawAesthetics(hRepSample_noCorrMuEScaleSyst, 1e-4, true);
    hRepSample_noCorrMuEScaleSyst->SetLineColor(kGreen - 2);
    ndfit::visuals::PredPreDrawAesthetics(hNominal, 1e-4, false);

    p1->cd();
    hRepSample->Draw("same hist e");
    hRepSample_noCorrMuEScaleSyst->Draw("same hist e");
    hNominal->Draw("same hist e");

    hRepSample->GetYaxis()->SetTitle("10^{4} Events / GeV");
    hRepSample->GetYaxis()->SetTitleSize(0.036);
    hRepSample->GetYaxis()->CenterTitle();
    hRepSample->GetYaxis()->SetTitleOffset(1.1);
    hRepSample->SetMaximum(hRepSample->GetMaximum() * 1.7);
    hRepSample->GetXaxis()->SetLabelSize(0.0);
    hRepSample->GetXaxis()->SetTitleSize(0.0);

    TLegend leg(0.45, 0.7, 0.9, 0.9);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    const std::string strRepSample = "ND Fit MCMC 'Rep. Sample' pulls applied";
    const std::string strNominal = "Prod5.1 Nominal Prediction";
    leg.AddEntry(hNominal, Form("%s", strNominal.c_str()), "l");
    leg.AddEntry(hRepSample, Form("%s", strRepSample.c_str()), "l");
    leg.AddEntry(hRepSample_noCorrMuEScaleSyst, "Rep. Sample (kCorrMuEScaleSyst2020 = 0#sigma", "l");
    leg.Draw("same");
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);
    latex.DrawLatexNDC(.15, .85, (Form("%s", beamType.c_str())));
    latex.DrawLatexNDC(.15, .8, Form("%s", quantileString.c_str()));
    latex.Draw("same");
    Simulation();
    ndfit::visuals::DetectorLabel(caf::kNEARDET);


    /// Enu-Reco ratio
    p2->cd();
    p2->SetGridy(1);
    TH1 *hEnuUnity = (TH1F *) hRepSample->Clone("hEnuUnity");
    hEnuUnity->Divide(hRepSample);

    TH1 *hRatioRepByNom = (TH1F *) hRepSample->Clone("hRatioRepByNom");
    hRatioRepByNom->Divide(hNominal);

    TH1 *hRatioRepByNom_NoMuEScaleSyst = (TH1F *) hRepSample_noCorrMuEScaleSyst->Clone("hRatioRepByNom_NoMuEScaleSyst");
    hRatioRepByNom_NoMuEScaleSyst->Divide(hNominal);

    hEnuUnity->SetLineColor(kBlack);
    hEnuUnity->SetLineWidth(1);
    hEnuUnity->Draw("same hist");
    hRatioRepByNom->Draw("same hist e");
    hRatioRepByNom_NoMuEScaleSyst->Draw("same hist e");

    hEnuUnity->GetXaxis()->CenterTitle();
    hEnuUnity->GetXaxis()->SetTitleOffset(1.);
    hEnuUnity->GetXaxis()->SetTitleSize(0.045);
    hEnuUnity->SetXTitle("Reco. E_{#nu} (GeV)");
    hEnuUnity->SetYTitle("#frac{Rep. Sample}{Nominal}");
    hEnuUnity->GetYaxis()->CenterTitle();
    hEnuUnity->GetYaxis()->SetRangeUser(0.5, 1.5);
    hEnuUnity->GetYaxis()->SetTitleSize(0.02);
    hEnuUnity->GetYaxis()->SetLabelSize(0.02);
    hEnuUnity->GetYaxis()->SetTitleOffset(1.5);
    hEnuUnity->GetYaxis()->CenterTitle();
    xAxisENu->Draw("same");
    
    for (const auto &ext: {"pdf", "png"})
      c.SaveAs(ndfit::FullFilename(outDirPlot, "plot_" + predName + "_RepSample_Nominal_CorrMuEScaleSyst_comparison_ratio." + ext).c_str());

  } // TH1 pairs, plotting



}
