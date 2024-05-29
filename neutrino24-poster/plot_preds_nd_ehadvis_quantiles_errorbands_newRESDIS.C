/* plot_preds_nd_ehadvis_quantiles_errorbands_newRESDIS.C
 *
 * preparation for Neutrino 2024 poster.
 *
 * produce some plots of the ND predictions
 * of two error bands:
 * 1. all xsec + "new" RES and DIS systematics,
 * 2. xsec WITHOUT the "new" RES/DIS
 * in the inclusive and quantile distributions
 * in 1 sigma error bands.
 *
 * Large chunk taken from: plot_fd_systematic_shifts_from_numu_quantile_predictions.C
 *
 *
 * May. 2024.
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
#include "3FlavorAna/Systs/DummySystStorage.h"
#include "3FlavorAna/Systs/EnergySysts2020.h"

#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Systs/RESSysts.h"
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
  std::cout << "antiNeutrino ? " << antiParticle << std::endl;
  nu == ndfit::NeutrinoType::kNue ? nuLtx = "#nu_{e}" : nuLtx = "#nu_{#mu}";
  if (antiParticle) nuLtx = "#bar{" + nuLtx + "}";
  auto * prelim = new TLatex(.25, .7, Form("%s", nuLtx.c_str())); // top left: below topology, horn
  prelim->SetTextColor(nu == ndfit::NeutrinoType::kNumu ? kBlue - 3 :  kRed - 3);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

using namespace ana;


// =====================================================================================================
void plot_preds_nd_ehadvis_quantiles_errorbands_newRESDIS(const bool saveCaptions = false)
// =====================================================================================================
{

  const std::string systString = "xsec24"; // --> all xsec24 systs.

  std::string outDirPlot = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/neutrino24-poster/plot_preds_nd_ehadvis_quantiles_errorbands_newRESDIS";

  //load all systs that exist in the preds ROOT file
  getAllXsecSysts_2024();

	// these are older preds (2023), before we changed the name, so need to instantiate as Dummy systs.
	DummyAnaSyst kDummyDIS0Syst = DummyAnaSyst("DIS_nubar_hadro_Q0_syst", "DIS_nubar_hadro_Q0_syst");
	DummyAnaSyst kDummyDIS1Syst = DummyAnaSyst("DIS_nu_hadro_Q1_syst", "DIS_nu_hadro_Q1_syst");
	DummyAnaSyst kDummyRESvpvnNuSyst = DummyAnaSyst("RES_vpvn_Nu_ratio_xsec_syst", "RES #frac{#nu+p}{#nu+n} Ratio XSec Syst");
	DummyAnaSyst kDummyRESvpvnNuBarSyst = DummyAnaSyst("RES_vpvn_NuBar_ratio_xsec_syst", "RES #frac{#bar{#nu}+p}{#bar{#nu}+n} Ratio XSec Syst");

	DummyAnaSyst khNFSI_MFP_2020GSF = DummyAnaSyst("hNFSI_MFP_2020GSF", "hNFSI_MFP_2020GSF");
	DummyAnaSyst khNFSI_FateFracEV1_2020GSF = DummyAnaSyst("hNFSI_FateFracEV1_2020GSF", "hNFSI_FateFracEV1_2020GSF");
	DummyAnaSyst khNFSI_FateFracEV2_2020GSF = DummyAnaSyst("hNFSI_FateFracEV2_2020GSF", "hNFSI_FateFracEV2_2020GSF");
	DummyAnaSyst khNFSI_FateFracEV3_2020GSF = DummyAnaSyst("hNFSI_FateFracEV3_2020GSF", "hNFSI_FateFracEV3_2020GSF");

	DummyAnaSyst kMECShape2020GSFNu = DummyAnaSyst("MECShape2020GSFNu", "MECShape2020GSFNu");
	DummyAnaSyst kMECShape2020GSFAntiNu = DummyAnaSyst("MECShape2020GSFAntiNu", "MECShape2020GSFAntiNu");

	DummyAnaSyst kFormZone2020GSF = DummyAnaSyst("FormZone2020GSF", "FormZone2020GSF");

	std::vector<const ISyst*> systsRESDIS;
	systsRESDIS.push_back(&kDummyDIS0Syst);
	systsRESDIS.push_back(&kDummyDIS1Syst);
	systsRESDIS.push_back(&kDummyRESvpvnNuSyst);
	systsRESDIS.push_back(&kDummyRESvpvnNuBarSyst);
	systsRESDIS.push_back(&kRESDeltaScaleSyst);
	systsRESDIS.push_back(&kRESOtherScaleSyst);





  // Load the ND Quantile predictions.
  const std::string& inputDir = "/exp/nova/data/users/mdolce/preds+spectra/ana2024/generate_nd_ehadvis_quantiles_predictions_ana2024/";
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





  auto calc2020BF = std::make_unique<osc::OscCalcAnalytic>();
  ndfit::Calculator2020BestFit(*calc2020BF);

  // Systematics from the Prediction
  auto systs = dynamic_cast<const ana::PredictionInterp *>(preds[0].pred)->GetAllSysts();
  std::sort(systs.begin(),
            systs.end(),
            [](const auto syst1, const auto syst2) {
                return syst1->ShortName() < syst2->ShortName();
            });

	auto systswRESDIS = systs;
	std::cout << "systs count: " << systs.size() << std::endl;
	// erase the RES/DIS systs from the systs big vector.
	for (const ISyst* &syst : systsRESDIS) {
		std::cout << "Found: " << syst->ShortName();
		auto it = std::find(systs.begin(), systs.end(), syst);
		systs.erase(it);
		std::cout << ".......erased: " << std::endl;
	}
	std::cout << "systs count: " << systs.size() << std::endl;

  TCanvas c("c","c", 600,600); // 900, 600
  TPad * p1, * p2; //p1 upper, p2 lower


  // set scale factors here.
  const double scaleFactor = 1e-6;

  int sampleType = 0;

  std::cout << "Plotting the ratio and fitted predictions with data" << std::endl;
  std::cout << "Looping through systematics. Total number of systematics: " << systs.size() << std::endl;

  for (const auto &predBundle : preds) {
    std::cout << predBundle.name << "......." << std::endl;

    const ndfit::Quantiles q = ndfit::visuals::GetQuantileEnum(predBundle.name);
    const std::string quantileString = ndfit::visuals::GetQuantileString(q);
    const std::string beamType = ndfit::visuals::GetHornCurrent(predBundle.name);

    double POT = -5.;

		if (predBundle.name.find("fhc") != std::string::npos)
			POT = kAna2024FHCPOT;
		else {
			POT = kAna2024RHCPOT;
		}
		std::cout << "Setting POT to MC....." << POT << std::endl;



    /// create the error bands -- one vector for each prediction.
    std::vector<TH1*> up1Shifts, dn1Shifts; // NO RES/DIS
    std::vector<TH1*> up1Shifts_RESDIS, dn1Shifts_RESDIS; // with RES/DIS

		for (const ISyst* &newSyst : systsRESDIS){
			std::cout << "Found systematic: " << newSyst->ShortName() << std::endl;
			SystShifts pm1SigmaShift;
			pm1SigmaShift.SetShift(newSyst, +1.);

			TH1 *hUp1 = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(POT,
																																												EExposureType::kPOT,
																																												kBinDensity);
	//      std::cout << "Up integral: " << hUp1->Integral() << std::endl;
				up1Shifts_RESDIS.push_back(hUp1);

				pm1SigmaShift.SetShift(newSyst, -1.);
				TH1 *hDn1 = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(POT,
																																												EExposureType::kPOT,
																																												kBinDensity);
	//      std::cout << "Down integral: " << hDn1->Integral() << std::endl;
				dn1Shifts_RESDIS.push_back(hDn1);
		} // RES DIS systs

    // Use ONLY the systs that were used in the fitting...
		// NOTE: these systematics now have the RES/DIS systs removed!
    for (const auto &syst : systs) {
      std::cout << "Looping through syst....." << syst->ShortName() << std::endl;

      SystShifts pm1SigmaShift;
      pm1SigmaShift.SetShift(syst, +1.);
      TH1 *hUp1 = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(POT,
                                                                                      EExposureType::kPOT,
                                                                                      kBinDensity);
//      std::cout << "Up integral: " << hUp1->Integral() << std::endl;
      up1Shifts.push_back(hUp1);

      pm1SigmaShift.SetShift(syst, -1.);
      TH1 *hDn1 = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(POT,
                                                                                      EExposureType::kPOT,
                                                                                      kBinDensity);
//      std::cout << "Down integral: " << hDn1->Integral() << std::endl;
      dn1Shifts.push_back(hDn1);
    } // systs to create error bands


		// organize the systs TH1s properly now.
		auto up1Shifts_total = up1Shifts;
		auto dn1Shifts_total = dn1Shifts;

		for (TH1* &h : up1Shifts_RESDIS)
			up1Shifts_total.push_back(h);
		for (TH1* &h : dn1Shifts_RESDIS)
			dn1Shifts_total.push_back(h);


      c.cd();
      c.Clear();




      auto * xAxisEHad = new TGaxis(0.001, 0.5, 0.8, 0.501, 0., 0.8, 10, "");
      xAxisEHad->SetLabelOffset(-0.015); // default is 0.005
      xAxisEHad->SetLabelFont(42);
      xAxisEHad->SetTitle("E_{had}^{vis} (GeV)");
      xAxisEHad->SetTitleOffset(1.2);
      xAxisEHad->CenterTitle();
      xAxisEHad->SetTitleFont(42);


      /// Plot comparison and ratio on save canvas
      SplitCanvas(0.25, p1, p2);
      // EHadVis
      //create the histograms for the PlotWithSystErrorBand() function
      std::cout << "Producing EHadVis plots for " << predBundle.name << "......" << std::endl;
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
			for (TH1 * hist : up1Shifts_total)
				hist->Scale(scaleFactor);
		for (TH1 * hist : dn1Shifts_total)
				hist->Scale(scaleFactor);


      p1->cd();
			// draw the larger error band first.
		auto ErrorBandRESDIS = PlotWithSystErrorBand(hCVPred, up1Shifts_total, dn1Shifts_total, kGray + 2, kGreen + 2, 1.3, false, 0.5,
																								 false);
		auto ErrorBand = PlotWithSystErrorBand(hCVPred, up1Shifts, dn1Shifts, kGray + 2, kGray, 1.3, false, 0.8, false); // w/o RES DIS
      hCVPred->Draw("same hist e");
      hCVPred->GetYaxis()->SetTitle("10^{6} Events / GeV");
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
			up1Shifts_total.at(0)->SetFillColor(kGreen + 2);
			up1Shifts_total.at(0)->SetLineColor(kGreen + 2);
      leg.AddEntry(up1Shifts.at(0), Form("%s", errorBands.c_str()), "f");
			leg.AddEntry(up1Shifts_total.at(0), "#1pm#sigma #pi^{#pm} unc.", "f");
      leg.Draw("same");
      TLatex latex;
      latex.SetTextSize(0.04);
      latex.SetTextAlign(13);
      latex.DrawLatexNDC(.15, .85, (Form("%s", beamType.c_str())));
      latex.DrawLatexNDC(.15, .8, Form("%s", quantileString.c_str()));
      latex.Draw("same");
      Simulation();
      NeutrinoLabel(ndfit::NeutrinoType::kNumu, beamType == "Antineutrino Beam");
      ndfit::visuals::DetectorLabel(caf::kNEARDET);

      /// EHadVis ratio
      p2->cd();
      p2->SetGridy(1);
      TH1 *hUnity = (TH1F *) hCVPred->Clone("hEUnity");
      hUnity->Divide(hCVPred);


      ///create the ratios for the error bands
      std::vector<TH1*> up1ShiftsRatio = up1Shifts;
      std::vector<TH1*> dn1ShiftsRatio = dn1Shifts;
			std::vector<TH1*> up1ShiftsRatio_total = up1Shifts_total;
      std::vector<TH1*> dn1ShiftsRatio_total = dn1Shifts_total;
      for (auto &hist: up1ShiftsRatio)
        hist->Divide(hCVPred);
      for (auto &hist: dn1ShiftsRatio)
        hist->Divide(hCVPred);
			for (auto &hist: up1ShiftsRatio_total)
        hist->Divide(hCVPred);
      for (auto &hist: dn1ShiftsRatio_total)
        hist->Divide(hCVPred);



			// draw the largest error first.
		PlotWithSystErrorBand(hUnity, up1ShiftsRatio_total, dn1ShiftsRatio_total, kGray + 2, kGreen + 2, 1.3, false, 0.8, false);
		PlotWithSystErrorBand(hUnity, up1ShiftsRatio, dn1ShiftsRatio, kGray + 2, kGray, 1.3, false, 0.5, false);

      hUnity->GetXaxis()->CenterTitle();
      hUnity->GetXaxis()->SetTitleOffset(1.);
      hUnity->GetXaxis()->SetTitleSize(0.045);
      hUnity->SetXTitle(""); // set from the TAxis object
      hUnity->GetYaxis()->CenterTitle();
      hUnity->GetYaxis()->SetRangeUser(0.6, 1.4);
      hUnity->GetYaxis()->SetTitleSize(0.02);
      hUnity->GetYaxis()->SetLabelSize(0.02);
      hUnity->GetYaxis()->SetTitleOffset(1.5);
      hUnity->SetYTitle("NOvA MC Ratio");
      hUnity->GetYaxis()->CenterTitle();
      xAxisEHad->Draw("same");


      ndfit::visuals::DetectorLabel(predBundle.det);
      for (const auto &ext: {".png", ".pdf"}) // ".root"
        c.SaveAs(ndfit::FullFilename(outDirPlot, "plot_" + predBundle.name + "_EHadVis_xsec24_errorbands_newRESDIS" + ext).c_str());

    // write captions here...
    if (saveCaptions) {
      std::ofstream ofile;
      ofile.open(ndfit::FullFilename(outDirPlot, "plot_" + predBundle.name + "_EHadVis_xsec24_errorbands_newRESDIS.txt").c_str());
      ofile << "Near Detector " << beamType << " Prod5.1 Ana2024 Monte Carlo prediction in the hadronic energy fraction:  " << quantileString
                << ". The variable in this plot is reconstructed hadronic visible energy (in dark grey)."
                   " The light grey band is the 1 sigma error from all existing NOvA cross-section uncertainties."
									 " The green band is the 1 sigma error from the same uncertainties plus the new RES and DIS uncertainties for the Ana2024 analysis."
                   " The top distribution is the number of events, and the bottom is the ratio from the MC." << std::endl;
      ofile.close();
    }


    sampleType++;
  } //predBundle in preds



}
