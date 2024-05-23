/* plot_preds_nd_ehadvis_allnumu_npi_resvpvn_dishadro_errorbands.C
 *
 * preparation for Ana2024.
 *
 * produce plots of the ND AllNumu EHadVis pred
 * of the "new" RES and DIS systematics,
 * with the true Npi^+- > 0 component
 * for Neutrino 2024 (talk) and poster.
 * Contains ONLY the RESvpn + DISHadro systs.
 * Plots on the same Canvas.
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
void plot_preds_nd_ehadvis_allnumu_npi_resvpvn_dishadro_errorbands(const std::string beam = "fhc",
				const bool saveCaptions = false)
// =====================================================================================================
{

  std::string outDirPlot = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/neutrino24-poster/plot_preds_nd_ehadvis_allnumu_npi_resvpvn_dishadro_errorbands";

  //load all systs that exist in the preds ROOT file
  NewRESDISSysts();

	std::string hc;

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
			beam + "_nd_" + pairPredName.first,
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


	// canvas for drawing BOTH preds together on same axis
  TCanvas c("c","c", 600,600); // 900, 600
  TPad * p1, * p2; //p1 upper, p2 lower

	auto * xAxisEHad = new TGaxis(0.001, 0.5, 0.8, 0.501, 0., 0.8, 10, "");
	xAxisEHad->SetLabelOffset(-0.015); // default is 0.005
	xAxisEHad->SetLabelFont(42);
	xAxisEHad->SetTitle("Hadronic Visible Energy (GeV)");
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

	TLegend leg(0.45, 0.6, 0.9, 0.85);
	leg.SetFillColor(0);
	leg.SetFillStyle(0);

	// set scale factors here.
  const double scaleFactor = 1e-6;
	double histMax = -5;

	// Assumption: we know the predictions are the both from AllNumu (Q5) sample.
	const ndfit::Quantiles q = ndfit::visuals::GetQuantileEnum(preds[0].name);
	const std::string quantileString = ndfit::visuals::GetQuantileString(q);
	const std::string beamType = ndfit::visuals::GetHornCurrent(preds[0].name);
	hc = beamType;

	TLatex latex;
	latex.SetTextSize(0.04);
	latex.SetTextAlign(13);


	assert (preds[1].name == "fhc_nd_pred_interp_Q5_chg_pi");

	// =================================== Create Histograms =================================================



	//NOTE: the histograms and error bands TH1s are outside the loop.
	std::vector<TH1*> up1Shifts_q5, dn1Shifts_q5;
	std::vector<TH1*> up1Shifts_chg_pi, dn1Shifts_chg_pi;
	for (const ana::ISyst* &syst : systs){
		SystShifts pm1SigmaShift;

		pm1SigmaShift.SetShift(syst, +1.);
		TH1 * hUp1_q5     = preds[0].pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(preds[0].pot);
		TH1 * hUp1_chg_pi = preds[1].pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(preds[1].pot);
		up1Shifts_q5.push_back(hUp1_q5);
		up1Shifts_chg_pi.push_back(hUp1_chg_pi);

		pm1SigmaShift.SetShift(syst, -1.);
		TH1 * hDn1_q5     = preds[0].pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(preds[0].pot);
		TH1 * hDn1_chg_pi = preds[1].pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(preds[1].pot);
		dn1Shifts_q5.push_back(hDn1_q5);
		dn1Shifts_chg_pi.push_back(hDn1_chg_pi);
	} // systs


	TH1 * hPredQ5        = preds[0].pred->PredictSyst(calc2020BF.get(), SystShifts::Nominal()).ToTH1(preds[0].pot);
	TH1 * hPredQ5_chg_pi = preds[1].pred->PredictSyst(calc2020BF.get(), SystShifts::Nominal()).ToTH1(preds[1].pot);


	TH1D * hPredQ5_chg_pi_Clone = (TH1D*) hPredQ5_chg_pi->Clone("hPredQ5_chg_pi_Clone");
	hPredQ5_chg_pi_Clone->SetLineColor(kGreen + 4); // slightly darker line.
	hPredQ5_chg_pi_Clone->SetFillColor(kGreen + 2);
	hPredQ5_chg_pi_Clone->SetFillColorAlpha(kGreen + 2, 0.5);

	// create the ratios -- ONLY for the Q5 prediction.
	TH1 *hUnity = (TH1F *) hPredQ5->Clone("hEUnity");
	// do the necessary removal of info.
	for (TH1* h : {(TH1*) hPredQ5, (TH1*) hPredQ5_chg_pi, (TH1*) hPredQ5_chg_pi_Clone, (TH1*) hUnity}) {
		h->GetYaxis()->SetTitleSize(0.036);
		h->GetYaxis()->SetTitleOffset(1.25);
		h->GetXaxis()->SetLabelSize(0.0);
		h->GetXaxis()->SetTitleSize(0.0);
		h->SetTitle(";;");
		h->GetXaxis()->SetRangeUser(0., 0.8);
		h->Scale(scaleFactor); // scale first, then divide.
	}
	hUnity->Divide(hPredQ5); // scale first, then divide.

	// scale the vectors of TH1s.
	for (const auto &vecHist : {up1Shifts_q5, up1Shifts_chg_pi, dn1Shifts_q5, dn1Shifts_chg_pi}){
		for (TH1 *histShift: vecHist)
			histShift->Scale(scaleFactor);
	}

	/// only create the ratios for the error bands of the Q5 AllNumu pred.
	std::vector<TH1*> up1ShiftsRatio = up1Shifts_q5;
	std::vector<TH1*> dn1ShiftsRatio = dn1Shifts_q5;
	for (auto &hist: up1ShiftsRatio)
		hist->Divide(hPredQ5);
	for (auto &hist: dn1ShiftsRatio)
		hist->Divide(hPredQ5);

	
	/// aesthetics for the Canvas
	hPredQ5->SetLineColor(kGray + 2);
	hPredQ5->SetLineWidth(3);
	hPredQ5->SetTitle(";;");
	hPredQ5->GetXaxis()->SetLabelSize(0);

	hPredQ5_chg_pi_Clone->SetTitle("; ; ");
	hPredQ5_chg_pi_Clone->GetXaxis()->SetLabelSize(0);

	up1Shifts_q5.at(0)->SetFillColor(kGray);
	leg.AddEntry(hPredQ5, "NOvA 2024 MC", "l");
	leg.AddEntry(up1Shifts_q5.at(0), "#pm1#sigma #pi^{#pm} unc.", "f");
	leg.AddEntry(hPredQ5_chg_pi_Clone, "True N#pi^{#pm} > 0", "f");



	///  =================== Draw ! ===================
	c.cd();

	/// Plot comparison and ratio on save canvas
	SplitCanvas(0.25, p1, p2);
	p1->cd();

	hPredQ5->Draw("same hist e");
	hPredQ5_chg_pi_Clone->Draw("same hist e");

	hPredQ5->GetYaxis()->SetTitle("10^{6} Events / GeV");
	hPredQ5->SetMaximum(hPredQ5->GetMaximum() * 1.8);

	// draw the error bands and CV of the two preds.
	auto tgQ5 = PlotWithSystErrorBand(hPredQ5, up1Shifts_q5, dn1Shifts_q5, kGray + 2, kGray);
//	auto tg_chg_pi = PlotWithSystErrorBand(hPredQ5_chg_pi, up1Shifts_chg_pi, dn1Shifts_chg_pi, kGreen + 4, kGreen + 2);

	latex.DrawLatexNDC(.15, .85, hc.c_str());
	latex.DrawLatexNDC(.15, .8,  quantileString.c_str());
	leg.Draw("same");
	latex.Draw("same");

	// post-hist drawing.
	Simulation();
	NeutrinoLabel(ndfit::NeutrinoType::kNumu, beamType == "Antineutrino Beam");
	ndfit::visuals::DetectorLabel(caf::kNEARDET);


	// Draw Ratio on p2.
	p2->cd();
	p2->SetGridy(1);

	// Draw Q5 error band.
	PlotWithSystErrorBand(hUnity, up1ShiftsRatio, dn1ShiftsRatio, kGray + 2, kGray);

	hUnity->SetTitle(";;");
	hUnity->GetXaxis()->SetTitleOffset(1.);
	hUnity->GetXaxis()->SetTitleSize(0.045);
	hUnity->SetXTitle(""); // set from the TAxis object
	hUnity->GetYaxis()->CenterTitle();
	hUnity->GetYaxis()->SetRangeUser(0.5, 1.5);
	hUnity->GetXaxis()->SetRangeUser(0., 0.8);
	hUnity->GetYaxis()->SetTitleSize(0.02);
	hUnity->GetYaxis()->SetLabelSize(0.02);
	hUnity->GetYaxis()->SetTitleOffset(1.5);
	hUnity->GetXaxis()->CenterTitle();
	hUnity->GetYaxis()->CenterTitle();

	hUnity->SetYTitle("MC Ratio");
	xAxisEHad->Draw("same");






	const std::string plotname = Form("plot_nd_allnumu_npi_%s_EHadVis_resvpvn_dishadro_errorbands", beam.c_str());
	for (const auto &ext: {".png", ".pdf"}) // ".root"
		c.SaveAs(ndfit::FullFilename(outDirPlot, plotname + ext).c_str());

	// write captions here...
	if (saveCaptions) {
		std::ofstream ofile;
		ofile.open(ndfit::FullFilename(outDirPlot, plotname + ".txt").c_str());
		ofile << "Near Detector " << hc << " Prod5.1 Ana2024 Monte Carlo prediction in hadronic visible energy. "
						 " The prediction is broken down by the total ND selected sample (grey) and the component that contains"
						 " at least one true charged pion, pre-FSI (green). "
						 " The light grey band is the 1 sigma error from NOvA's RESvp/vn Nu and Nubar systematics and "
						 "DIS Nu and Nubar hadronization systematics, implemented for the first time in the Ana2024 analysis."
						 " The top distribution is the number of events; the bottom is the ratio relative to the nominal MC." << std::endl;
		ofile.close();
	}


}
