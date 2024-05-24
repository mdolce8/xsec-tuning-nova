/*
 *  apply_crpa_lfg_lowE_nue_ratio.C
 *
 * 	Apply the ratio from the cRPA paper to our
 *  NOvA FD FHC LowE Nue (appeared + background) spectra.
 *
 *  This time we will use:
 *  RR_CRPA_C_LFG_O_e_mu
 *  (this assumes the cRPA effect is applied
 *  to numu as well).
 *
 *
 *
 *  May. 2024
 *  M. Dolce
 */

#include <iostream>


#include "3FlavorAna/NDFit/NDFitHelper.h"
#include "3FlavorAna/NDFit/InitializeFit.h"

#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Analysis/Exposures.h"

#include "OscLib/OscCalcPMNSOpt.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "Utilities/rootlogon.C"


using namespace ana;

// =====================================================================================================
void apply_crpa_lfg_lowE_nue_ratio(const std::string& beam = "fhc")
// =====================================================================================================
{

	const std::string outDir = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/cRPA/apply_crpa_lfg_lowE_nue_ratio/";

	// ROOT file Summed TH2
	TFile * fnue = TFile::Open(Form("/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/cRPA/plot_fd_enu_theta_lowE_nue/th2_fd_%s_enu_theta_lowE_nue_summed.root", beam.c_str()), "read");

	// TH2 sum of the nue_app and nuebar_app portions -- "Signal"
	TH2D * h2Signal = (TH2D*) fnue->Get("Signal");
	std::cout << "h2Sum events: " << h2Signal->Integral() << std::endl;

	// make a clone for the reweighted histogram
	TH2D * h2SumRwgt = (TH2D*) h2Signal->Clone("nue_sig_cRPA");
	h2SumRwgt->Reset("ICESM"); // maintain the binning, just clear the content
	assert (h2SumRwgt->Integral() == 0);

	// I think this is the right ratio to read in at the moment?
	TFile * fcRPA = TFile::Open("/exp/nova/data/users/mdolce/xsec-tuning-nova/RR_CRPA_C_LFG_O_e_mu.root", "read");

	// TH2 of the cRPA / LFG ratio for nue
	TH2D * h2RPARatio = (TH2D*) fcRPA->Get("RR_CRPA_C_LFG_O_e_mu");
	std::cout << "h2Ratio events: " << h2RPARatio->Integral() << std::endl;


	// LowE Sample: 0.5 GeV - 1.5 GeV. ----> 10 bins of 100 MeV each.
	// But the cRPA paper ratio starts at 0 GeV, of 100 MeV bins,
	// So we only want to apply the cRPA ratio of the middle 10 bins!


	// NOTE: _this_ is the right way: start at [1, GetNbins()].
	// NOTE: the bin widths are identical, so this should be easy...?
	// NOTE: must do Y loop first, because we want the last row (i.e. same y bin).
	for (unsigned int binIdxY = 1; binIdxY <= h2Signal->GetNbinsY(); binIdxY++){
		for (unsigned int binIdxX = 1; binIdxX <= h2Signal->GetNbinsX(); binIdxX++){

			// standard case, apply ratio as normal.
			// We want to skip the first 5 ROWS (Y bins) of the distribution. Bc LowE starts at E = 0.5 GeV
			if (binIdxX <= h2RPARatio->GetNbinsX() && binIdxY <= h2RPARatio->GetNbinsY()) {
				const double contentRatio = h2RPARatio->GetBinContent(binIdxX, binIdxY + 5); // from the cRPA ratio paper.
				const double contentSumNue = h2Signal->GetBinContent(binIdxX, binIdxY);  // NOvA's FD LowE sample

				const double contentRwgt = contentRatio * contentSumNue; // this is the value after re-weighting to cRPA.
				h2SumRwgt->SetBinContent(binIdxX, binIdxY, contentRwgt);

				if (binIdxX % 5 == 0) std::cout << "Applying weight as expected. To bin: " << binIdxX << ", " << binIdxY <<
																				". To values: (" << h2Signal->GetXaxis()->GetBinCenter(binIdxX) << ", " << h2Signal->GetYaxis()->GetBinCenter(binIdxY) <<
																				"  Ratio to be applied to Nue App. = " << contentRatio << std::endl;
			}

			// special case: address bin values above 1.2 GeV....
			else if (binIdxY > h2RPARatio->GetNbinsY()){
				if (binIdxX % 15 == 0) { // one for each row.
					std::cout << "Neutrino Energy: " << h2RPARatio->GetYaxis()->GetBinCenter(binIdxY) << " GeV" << std::endl;
					std::cout << "Applying weight to Angle (deg) " << h2Signal->GetXaxis()->GetBinCenter(binIdxX) << std::endl;
					std::cout << "binIdxY == " << binIdxY << std::endl;
				}

				// want the bin content from the top row. That would be this value.
				// NOTE: we need to scan the x values still. Only the Y bin is constant
				const double contentRatioTopRow = h2RPARatio->GetBinContent(binIdxX, h2RPARatio->GetNbinsY());
				const double contentSumNue = h2Signal->GetBinContent(binIdxX, binIdxY);

				const double contentRwgt = contentSumNue * contentRatioTopRow;
				h2SumRwgt->SetBinContent(binIdxX, binIdxY, contentRwgt);
			}

			else {
				std::cout << "binIdxX, binIDxY" << binIdxX << ", " << binIdxY << std::endl;
				std::cerr << "I don't know what to do...exit." << std::endl;
				exit(0);
			}

		} // binIdxY
	} // binIdxX


	// Draw the Re-weighted histogram now. No need to scale, already scaled to 2024 POT.
	TCanvas c;

	h2SumRwgt->SetTitle("Approx. cRPA");
	h2SumRwgt->Draw("same hist colz");
//	h2SumRwgt->GetZaxis()->SetRangeUser(0.5, 1.5); // this works. but leave out for now. Want this for a new TH2, the ratio of the ratio...?
	TLatex latex;
	latex.SetTextColor(kGray);
	latex.DrawLatexNDC(0.62, 0.6, Form(beam == "fhc" ? "Neutrino Beam" : "Antineutrino Beam"));
	latex.SetTextSize(0.85);
	TLatex ltx2;
	ltx2.SetTextColor(kGray);
	ltx2.DrawLatexNDC(0.62, 0.5, "Asimov A");
	ltx2.SetTextSize(0.85);
	Simulation();

	std::cout << " events integral: " << h2SumRwgt->Integral() << std::endl;

	c.SaveAs(Form("%s/plot_crpa_lfg_fd_%s_prod5.1_enu_theta_lowE_nue.png", outDir.c_str(),  beam.c_str()));


	// Request to make projection onto Y axis.
	c.Clear();
	auto hProjLowE_Rwgt = h2SumRwgt->ProjectionY();
	auto hProjLowE_Signal = h2Signal->ProjectionY();
	hProjLowE_Rwgt->SetLineWidth(3);
	hProjLowE_Rwgt->SetLineColor(kBlue);
	hProjLowE_Rwgt->Draw("same hist");
	hProjLowE_Signal->SetLineColor(kRed);
	hProjLowE_Signal->SetLineWidth(3);
	hProjLowE_Signal->Draw("same hist");
	latex.DrawLatexNDC(0.62, 0.6, Form(beam == "fhc" ? "Neutrino Beam" : "Antineutrino Beam"));
	ltx2.DrawLatexNDC(0.62, 0.5, "Asimov A");
	TLatex ltx3;
	ltx3.SetTextColor(kGray);
	ltx3.SetTextSize(0.85);
	ltx3.DrawLatexNDC(0.62, 0.7, "(#theta_e, E_#nu) Projection");

	TLegend leg(0.15, 0.6, 0.5, 0.85);
	leg.SetFillColor(0);
	leg.SetFillStyle(0);
	leg.AddEntry(hProjLowE_Signal, "NOvA Signal", "l");
	leg.AddEntry(hProjLowE_Rwgt, "~cRPA Effect", "l");
	leg.Draw("same");

	c.SaveAs(Form("%s/plot_crpa_lfg_fd_%s_prod5.1_enu_theta_lowE_nue_projection.png", outDir.c_str(),  beam.c_str()));


}