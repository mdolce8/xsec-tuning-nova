/* plot_preds_nd_ehadvis_allnumu_errorbands_charged_vs_neutral_pions.C
 *
 * preparation for Ana2024.
 *
 * produce plots of the ND AllNumu EHadVis pred
 * of the "new" RES and DIS systematics,
 * with the true Npi^+- > 0 component
 * for Neutrino 2024 (talk) and poster.
 * Contains ONLY the RESvpn + DISHadro systs.
 *
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





// map of the two files.
std::map<std::string, std::string> FILE_NAMES
{
				{"pred_interp_Q5", "pred_interp_Q5_nd_allnumu_%s_resvpvn_dishadro_ehadvis.root"},
				{"pred_interp_Q5_chg_pi", "pred_interp_Q5_chg_pi_nd_allnumu_%s_resvpvn_dishadro_ehadvis.root"}
};

void NeutrinoLabel(const ndfit::NeutrinoType nu, const bool antiParticle = false){
  std::string nuLtx;
//  std::cout << "antiNeutrino ? " << antiParticle << std::endl;
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
void plot_preds_nd_ehadvis_allnumu_errorbands_charged_vs_neutral_pions(const std::string beam = "fhc",
				const bool saveCaptions = false)
// =====================================================================================================
{

  const std::string systString = "xsec24"; // --> all xsec24 systs.

  std::string outDirPlot = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/neutrino24-poster/plot_preds_nd_ehadvis_allnumu_errorbands_charged_vs_neutral_pions";

  //load all systs that exist in the preds ROOT file
  NewRESDISSysts();



  // Load the ND predictions.
  const std::string& inputDir = "/exp/nova/data/users/mdolce/preds+spectra/ana2024/neutrino24-poster/generate_nd_allnumu_ehadvis_charged_vs_neutral_hadrons/test";
  std::vector<ana::FitPredictions> preds;
  for (const auto &pairPredName: FILE_NAMES) {

		const std::string path = inputDir + "/" + Form(FILE_NAMES.at(pairPredName.first).c_str(), beam.c_str());
		TFile f(path.c_str());
		auto p = ana::LoadFrom<PredictionInterp>(&f, pairPredName.first);

		if (!p)
		{
			// no prediction loaded from file
			std::cerr << "File not not found. Failed to find file: " << &f << std::endl;
			exit(1);
		}

		// Ensure the pred exists and spline gets initialized, and save some memory.
		p->PredictSyst(static_cast<osc::IOscCalc *>(nullptr), ana::kNoShift);
		p->MinimizeMemory();

		const double pot = beam == "fhc" ? kAna2024FHCPOT : kAna2024RHCPOT;

		ana::FitPredictions pred {
			"fhc_nd_" + pairPredName.first,
			p.release(),  // ignore, this is the right thing to do.
			{nullptr, 0.},
			pot,
			0.,
			caf::kNEARDET,
			ndfit::kNumu
		};

		preds.emplace_back(pred);
  } // taken partially from LoadPreds()



	// not really needed, since ND pred.
  auto calc2020BF = std::make_unique<osc::OscCalcAnalytic>();
  ndfit::Calculator2020BestFit(*calc2020BF);

  // Systematics from the Prediction
  auto systs = dynamic_cast<const ana::PredictionInterp *>(preds[0].pred)->GetAllSysts();
  std::sort(systs.begin(),
            systs.end(),
            [](const auto syst1, const auto syst2) {
                return syst1->ShortName() < syst2->ShortName();
            });

	std::cout << "Total number of systs found: " << systs.size() << std::endl;



  TCanvas c("c","c", 600,600); // 900, 600
  TPad * p1, * p2; //p1 upper, p2 lower


  // set scale factors here.
  const double scaleFactor = 1e-6;


	int predCounter = 0;
  std::cout << "Plotting the ratio and predictions" << std::endl;
  for (const auto &predBundle : preds) {
    std::cout << predBundle.name << "......." << std::endl;

    const ndfit::Quantiles q = ndfit::visuals::GetQuantileEnum(predBundle.name);
    const std::string quantileString = ndfit::visuals::GetQuantileString(q);
    const std::string beamType = ndfit::visuals::GetHornCurrent(predBundle.name);

    /// create the error bands -- one vector for each prediction.
    std::vector<TH1*> up1Shifts, dn1Shifts;

    // use the systs to create error bands
    // Use ONLY the systs that were used in the fitting...
		std::cout << "Looping through systematics. Total number of systematics: " << systs.size() << std::endl;
		for (const auto &syst : systs) {
      std::cout << "Looping through syst....." << syst->ShortName() << std::endl;


      SystShifts pm1SigmaShift;
      pm1SigmaShift.SetShift(syst, +1.);
      TH1 *hUp1 = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(predBundle.pot,
                                                                                      EExposureType::kPOT,
                                                                                      kBinDensity);
//      std::cout << "Up integral: " << hUp1->Integral() << std::endl;
      up1Shifts.push_back(hUp1);

      pm1SigmaShift.SetShift(syst, -1.);
      TH1 *hDn1 = predBundle.pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(predBundle.pot,
                                                                                      EExposureType::kPOT,
                                                                                      kBinDensity);
//      std::cout << "Down integral: " << hDn1->Integral() << std::endl;
      dn1Shifts.push_back(hDn1);

    } // systs to create error bands


		// do the plotting
      c.cd();
      c.Clear();

      auto * xAxisEHad = new TGaxis(0.001, 0.5, 0.8, 0.501, 0., 0.8, 10, "");
      xAxisEHad->SetLabelOffset(-0.015); // default is 0.005
      xAxisEHad->SetLabelFont(42);
      xAxisEHad->SetTitle("E_{had}^{vis} (GeV)");
      xAxisEHad->SetTitleOffset(1.2);
      xAxisEHad->CenterTitle();
      xAxisEHad->SetTitleFont(42);

      //pavetext to print out the events
      TPaveText ptEnuEvents(0.7, 0.60, 0.85, 0.67, "ARC NDC");
      ptEnuEvents.SetFillColor(0);
      ptEnuEvents.SetFillStyle(0);
      ptEnuEvents.SetBorderSize(0);
      ptEnuEvents.SetTextSize(0.032);
      ptEnuEvents.SetTextFont(102);

      // contains all systs
      PlotWithSystErrorBand((IPrediction *&) predBundle.pred, systs, calc2020BF.get(), predBundle.pot, kGray + 2, kGray);
      c.SaveAs(ndfit::FullFilename(outDirPlot, "profiled_error_bands_plot_" + predBundle.name + ".png").c_str());
      c.Clear();
      // 2D profile




      /// Plot comparison and ratio on save canvas
      SplitCanvas(0.25, p1, p2);
      // EHadVis
      //create the histograms for the PlotWithSystErrorBand() function
      std::cout << "Producing EHadVis plots for " << predBundle.name << "......" << std::endl;
      TH1 * hCVPred = predBundle.pred->PredictSyst(calc2020BF.get(), SystShifts::Nominal()).ToTH1(predBundle.pot,
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


      p1->cd();
      auto ErrorBand = PlotWithSystErrorBand(hCVPred, up1Shifts, dn1Shifts, kGray + 2, kGray);
      hCVPred->Draw("same hist e");
		if (predCounter != 0) {
			TH1D * hCVPredClone = (TH1D*) hCVPred->Clone("hCVPredClone");
			hCVPredClone->SetLineColor(kGreen + 2);
			hCVPredClone->SetFillColor(kGreen + 2);
		}
		hCVPred->GetYaxis()->SetTitle("10^{6} Events / GeV");
      hCVPred->GetYaxis()->SetTitleSize(0.036);
      hCVPred->GetYaxis()->SetTitleOffset(1.1);
      hCVPred->SetMaximum(hCVPred->GetMaximum() * 2.0);
      hCVPred->GetXaxis()->SetLabelSize(0.0);
      hCVPred->GetXaxis()->SetTitleSize(0.0);
			hCVPred->GetXaxis()->SetRangeUser(0., 0.8);

      TLegend leg(0.45, 0.65, 0.9, 0.9);
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      const std::string errorBands = "#pm1#sigma #pi^{#pm} unc.";
      const std::string cvPred = "NOvA 2024 MC";
      std::string legCVText;
      leg.AddEntry(hCVPred, cvPred.c_str(), "l");
      up1Shifts.at(0)->SetFillColor(kGray);
      up1Shifts.at(0)->SetLineColor(kGray);
      leg.AddEntry(up1Shifts.at(0), Form("%s", errorBands.c_str()), "f");
      leg.Draw("same");
      TLatex latex;
      latex.SetTextSize(0.04);
      latex.SetTextAlign(13);
      latex.DrawLatexNDC(.15, .85, (Form("%s", beamType.c_str())));
      latex.DrawLatexNDC(.15, .8, Form("%s", quantileString.c_str()));
      latex.Draw("same");
//    ptEnuEvents.Draw("same");
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
      hUnity->SetYTitle("MC Ratio");
      hUnity->GetYaxis()->CenterTitle();
      xAxisEHad->Draw("same");


      ndfit::visuals::DetectorLabel(predBundle.det);
      for (const auto &ext: {".png", ".pdf"}) // ".root"
        c.SaveAs(ndfit::FullFilename(outDirPlot, "plot_" + predBundle.name + "_EHadVis_xsec24_errorbands" + ext).c_str());

    // write captions here...
    if (saveCaptions) {
      std::ofstream ofile;
      ofile.open(ndfit::FullFilename(outDirPlot, "plot_" + predBundle.name + "_EHadVis_xsec24_errorbands" + ".txt").c_str());
      ofile << "Near Detector " << beamType << " Prod5.1 Ana2024 Monte Carlo prediction in the hadronic energy fraction:  " << quantileString
                << ". The variable in this plot is reconstructed hadronic visible energy (in dark grey)."
                   " The light grey band is the 1 sigma error from all NOvA cross-section uncertainties for the Ana2024 analysis."
                   " This includes the new RES and DIS uncertainties from the ND fitting work."
                   " The top distribution is the number of events, and the bottom is the ratio from the MC." << std::endl;
      ofile.close();
    }

 	predCounter++;
  } //predBundle in preds



}
