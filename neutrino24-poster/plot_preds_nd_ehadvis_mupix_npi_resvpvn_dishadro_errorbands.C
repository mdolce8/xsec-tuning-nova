/* plot_preds_nd_ehadvis_mupix_npi_resvpvn_dishadro_errorbands.C
 *
 * preparation for Ana2024.
 *
 * produce plots of the ND Mupix EHadVis pred
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
				{"pred_interp_mupix", "pred_interp_mupix_nd_mupix_%s_resvpvn_dishadro_ehadvis.root"},
				{"pred_interp_mupix_chg_pi", "pred_interp_mupix_chg_pi_nd_mupix_%s_resvpvn_dishadro_ehadvis.root"}
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
void plot_preds_nd_ehadvis_mupix_npi_resvpvn_dishadro_errorbands(const std::string beam = "fhc",
				const bool saveCaptions = false)
// =====================================================================================================
{

  std::string outDirPlot = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/neutrino24-poster/plot_preds_nd_ehadvis_mupix_npi_resvpvn_dishadro_errorbands";

  //load all systs that exist in the preds ROOT file
  NewRESDISSysts();

	std::string hc;

  // Load the ND predictions.
  const std::string& inputDir = "/exp/nova/data/users/mdolce/preds+spectra/ana2024/neutrino24-poster/generate_nd_mupix_ehadvis_charged_vs_neutral_hadrons/test";
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

	auto * xAxisEHad = new TGaxis(0.001, 0.9, 0.8, 0.901, 0., 0.8, 10, "");
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

	TLegend leg(0.45, 0.7, 0.9, 0.88);
	leg.SetFillColor(0);
  leg.SetFillStyle(0);
  TLegend legUnc(0.45, 0.65, 0.9, 0.68);
  legUnc.SetFillColor(0);
  legUnc.SetFillStyle(0);

	// set scale factors here.
  const double scaleFactor = 1e-4;
	double histMax = -5;

	// Assumption: we know the predictions are the both from Mupix sample.
	const ndfit::TopologyType t = ndfit::visuals::GetTopologyName("fhc_nd_pred_interp_MuPiEtc", "p5p1");
	const std::string topoLtx = ndfit::visuals::GetTopologyLatexName(t);
	const std::string beamType = ndfit::visuals::GetHornCurrent(preds[0].name);
	hc = beamType;

	TLatex latex;
	latex.SetTextSize(0.04);
	latex.SetTextAlign(13);

	std::cout << "First  prediction is: " << preds[0].name << std::endl;
	std::cout << "Second prediction is: " << preds[1].name << std::endl;
	assert (preds[1].name == Form("%s_nd_pred_interp_mupix_chg_pi", beam.c_str()));

	// =================================== Create Histograms =================================================

	TH1 * hPred_MuPiX        = preds[0].pred->PredictSyst(calc2020BF.get(), SystShifts::Nominal()).ToTH1(preds[0].pot);
	TH1 * hPred_MuPiX_chgpi = preds[1].pred->PredictSyst(calc2020BF.get(), SystShifts::Nominal()).ToTH1(preds[1].pot);

	//NOTE: the histograms and error bands TH1s are outside the loop.
	std::vector<TH1*> up1Shifts_MuPiX, dn1Shifts_MuPiX;
	std::vector<TH1*> up1Shifts_chg_pi, dn1Shifts_chg_pi;
	for (const ana::ISyst* &syst : systs){
		SystShifts pm1SigmaShift;

		pm1SigmaShift.SetShift(syst, +1., true);
		TH1D * hUp1_q5     = preds[0].pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(preds[0].pot);
		TH1D * hUp1_chg_pi = preds[1].pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(preds[1].pot);
		up1Shifts_MuPiX.push_back(hUp1_q5);
		up1Shifts_chg_pi.push_back(hUp1_chg_pi);

		pm1SigmaShift.SetShift(syst, -1., true);
		TH1D * hDn1_q5     = preds[0].pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(preds[0].pot);
		TH1D * hDn1_chg_pi = preds[1].pred->PredictSyst(calc2020BF.get(), pm1SigmaShift).ToTH1(preds[1].pot);
		dn1Shifts_MuPiX.push_back(hDn1_q5);
		dn1Shifts_chg_pi.push_back(hDn1_chg_pi);

		std::cout << "syst: " << syst->ShortName() << std::endl;
		std::cout << "compare Up1 vs PredMuPiX: " << hUp1_q5->Integral() << ", " << hPred_MuPiX->Integral() << std::endl;
		std::cout << "compare Dn1 vs PredMuPiX: " << hDn1_q5->Integral() << ", " << hPred_MuPiX->Integral() << std::endl;
	} // systs



	TH1D * hPred_MuPiX_chg_pi_Clone = (TH1D*) hPred_MuPiX_chgpi->Clone("hPred_MuPiX_chgpi_Clone");
	hPred_MuPiX_chg_pi_Clone->SetLineColor(kGreen + 4); // slightly darker line.
	hPred_MuPiX_chg_pi_Clone->SetFillColor(kGreen + 1);
	hPred_MuPiX_chg_pi_Clone->SetFillColorAlpha(kGreen + 1, 0.5);

	// create the ratios -- ONLY for the MuPiX prediction.
	TH1 *hUnity = (TH1F *) hPred_MuPiX->Clone("hEUnity");
	// do the necessary removal of info.

	// for aesthetics
	TH1D hBlank("hBlank", "", 1, 0., 0.8);
	hBlank.SetBinContent(1, 1.);

	for (TH1* h : {(TH1*) hPred_MuPiX, (TH1*) hPred_MuPiX_chgpi, (TH1*) hPred_MuPiX_chg_pi_Clone, (TH1*) hUnity, (TH1*) &hBlank}) {
		h->GetYaxis()->SetTitleSize(0.036);
		h->GetYaxis()->SetTitleOffset(1.25);
		h->GetXaxis()->SetLabelSize(0.0);
		h->GetXaxis()->SetTitleSize(0.0);
		h->SetTitle(";;");
		h->GetXaxis()->SetRangeUser(0., 0.8);
		h->Scale(scaleFactor); // scale first, then divide.
	}
	hUnity->Divide(hPred_MuPiX); // scale first, then divide.

	// scale the vectors of TH1s.
	for (const auto &vecHist : {up1Shifts_MuPiX, up1Shifts_chg_pi, dn1Shifts_MuPiX, dn1Shifts_chg_pi}){
		for (TH1 *histShift: vecHist)
			histShift->Scale(scaleFactor);
	}
//	for (TH1* hShift : up1Shifts_q5) {
//		hShift->Scale(scaleFactor);
//	}
//	for (TH1* hShift : dn1Shifts_q5) {
//		hShift->Scale(scaleFactor);
//	}



	/// only create the ratios for the error bands of the Q5 Mupix pred.
	std::vector<TH1*> up1ShiftsRatio = up1Shifts_MuPiX;
	std::vector<TH1*> dn1ShiftsRatio = dn1Shifts_MuPiX;
	for (auto &hist: up1ShiftsRatio)
		hist->Divide(hPred_MuPiX);
	for (auto &hist: dn1ShiftsRatio)
		hist->Divide(hPred_MuPiX);

	
	/// aesthetics for the Canvas
	hPred_MuPiX->SetLineColor(kGray + 2);
	hPred_MuPiX->SetLineWidth(3);
	hPred_MuPiX->SetTitle(";;");
	hPred_MuPiX->GetXaxis()->SetLabelSize(0);

	hPred_MuPiX_chg_pi_Clone->SetTitle("; ; ");
	hPred_MuPiX_chg_pi_Clone->GetXaxis()->SetLabelSize(0);

	up1Shifts_MuPiX.at(0)->SetFillColor(kGray);
	up1Shifts_chg_pi.at(0)->SetFillColor(kGreen + 4);
	leg.AddEntry(hPred_MuPiX, "2024 sim.", "l");
//	leg.AddEntry(up1Shifts_q5.at(0), "#pm1#sigma #pi^{#pm} unc.", "f");
	leg.AddEntry(hPred_MuPiX_chg_pi_Clone, "True N#pi^{#pm} > 0", "f");
  legUnc.SetHeader("Cross section unc.", "c");
  legUnc.AddEntry(up1Shifts_chg_pi.at(0), "#frac{#sigma_{RES}(#nup)}{#sigma_{RES}(#nun)} + DIS Hadro.", "f");



	///  =================== Draw ! ===================
	c.cd();

	/// Plot comparison and ratio on save canvas
	SplitCanvas(0.25, p1, p2);
	p1->cd();

	hBlank.SetLineColor(0);
	hBlank.Draw("same");
	hBlank.SetMaximum(hPred_MuPiX->GetMaximum() * 1.8);
	// draw the error bands and CV of the two preds.

	auto tg_chg_pi = PlotWithSystErrorBand(hPred_MuPiX_chgpi, up1Shifts_chg_pi, dn1Shifts_chg_pi, kGreen + 4, kGreen + 3, 1, false, 0.8, false);
//	auto tgQ5 = PlotWithSystErrorBand(hPredQ5, up1Shifts_q5, dn1Shifts_q5, kGray + 2, kGray, 1, false, 0.8, false);

	hPred_MuPiX->Draw("same hist");
	hPred_MuPiX_chg_pi_Clone->Draw("same hist");


	hBlank.GetYaxis()->SetTitle("10^{4} Events / GeV");
	hBlank.GetYaxis()->CenterTitle();
	hPred_MuPiX->SetMaximum(hPred_MuPiX->GetMaximum() * 1.8);

	latex.DrawLatexNDC(.15, .85, hc.c_str());
	latex.DrawLatexNDC(.15, .8,  topoLtx.c_str());
	leg.Draw("same");
  legUnc.Draw("same");
	latex.Draw("same");

	// post-hist drawing.
	Simulation();
	NeutrinoLabel(ndfit::NeutrinoType::kNumu, beamType == "Antineutrino Beam");
	ndfit::visuals::DetectorLabel(caf::kNEARDET);


	// Draw Ratio on p2.
	p2->cd();
	p2->SetGridy(1);

	// Draw Q5 error band.
	PlotWithSystErrorBand(hUnity, up1ShiftsRatio, dn1ShiftsRatio, kGreen + 4, kGreen + 3);

	hUnity->SetTitle(";;");
	hUnity->GetXaxis()->SetTitleOffset(1.);
	hUnity->GetXaxis()->SetTitleSize(0.045);
	hUnity->SetXTitle(""); // set from the TAxis object
	hUnity->GetYaxis()->CenterTitle();
	hUnity->GetYaxis()->SetRangeUser(0.9, 1.1);
	hUnity->GetXaxis()->SetRangeUser(0., 0.8);
	hUnity->GetYaxis()->SetTitleSize(0.02);
	hUnity->GetYaxis()->SetLabelSize(0.02);
	hUnity->GetYaxis()->SetTitleOffset(1.8);
	hUnity->GetXaxis()->CenterTitle();
	hUnity->GetYaxis()->CenterTitle();

	hUnity->SetYTitle("Shifted / Nominal");
	xAxisEHad->Draw("same");



 c.Update();


	const std::string plotname = Form("plot_nd_mupix_npi_%s_EHadVis_resvpvn_dishadro_errorbands", beam.c_str());
	for (const auto &ext: {".png", ".pdf"}) // ".root"
		c.SaveAs(ndfit::FullFilename(outDirPlot, plotname + ext).c_str());

	// write captions here...
	if (saveCaptions) {
		std::ofstream ofile;
		ofile.open(ndfit::FullFilename(outDirPlot, plotname + ".txt").c_str());
		ofile << "Near Detector " << hc << " Prod5.1 Ana2024 Monte Carlo prediction in hadronic visible energy. "
						 " The prediction is broken down by the charge-pion and prong topological sample (grey) and the component that contains"
						 " at least one true charged pion, pre-FSI (green). "
						 " The light grey band is the 1 sigma error from NOvA's RESvp/vn Nu and Nubar systematics and "
						 "DIS Nu and Nubar hadronization systematics, implemented for the first time in the Ana2024 analysis."
						 " The top distribution is the number of events; the bottom is the ratio relative to the nominal MC." << std::endl;
		ofile.close();
	}


}
