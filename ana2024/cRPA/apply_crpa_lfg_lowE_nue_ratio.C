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
void apply_crpa_lfglowE_nue_ratio(const std::string& beam = "fhc")
// =====================================================================================================
{

	const std::string outDir = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/cRPA/apply_crpa_lfg_lowE_nue_ratio/";

	// ROOT file Summed TH2
	TFile * fnue = TFile::Open(Form("/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/cRPA/plot_fd_enu_theta_lowE_nue/th2_fd_%s_enu_theta_lowE_nue_summed.root", beam.c_str()), "read");

	// TH2 sum of the nue_app and nuebar_app portions.
	TH2D * h2Sum = (TH2D*) fnue->Get("nue_app");
	std::cout << "h2Sum events: " << h2Sum->Integral() << std::endl;

	// make a clone for the reweighted histogram
	TH2D * h2SumRwgt = (TH2D*) h2Sum->Clone("nue_app_cRPA");
	h2SumRwgt->Reset("ICESM"); // maintain the binning, just clear the content
	assert (h2SumRwgt->Integral() == 0);

	// I think this is the right ratio to read in at the moment?
	TFile * fcRPA = TFile::Open("/exp/nova/data/users/mdolce/xsec-tuning-nova/RR_CRPA_C_LFG_O_e_mu.root", "read");

	// TH2 of the cRPA / LFG ratio for nue
	TH2D * h2Ratio = (TH2D*) fcRPA->Get("RR_CRPA_C_LFG_O_e_ae");
	std::cout << "h2Ratio events: " << h2Ratio->Integral() << std::endl;


	// NOTE: _this_ is the right way: start at [1, GetNbins()].
	// NOTE: the bin widths are identical, so this should be easy...?
	// NOTE: must do Y loop first, because we want the last row (i.e. same y bin).
	for (unsigned int binIdxY = 1; binIdxY <= h2Sum->GetNbinsY(); binIdxY++){
		for (unsigned int binIdxX = 1; binIdxX <= h2Sum->GetNbinsX(); binIdxX++){

			// standard case, apply ratio as normal.
			if (binIdxX <= h2Ratio->GetNbinsX() && binIdxY <= h2Ratio->GetNbinsY()) {
				const double contentRatio = h2Ratio->GetBinContent(binIdxX, binIdxY);
				const double contentSumNue = h2Sum->GetBinContent(binIdxX, binIdxY);

				const double contentRwgt = contentRatio * contentSumNue; // this is the value after re-weighting to cRPA.
				h2SumRwgt->SetBinContent(binIdxX, binIdxY, contentRwgt);

				if (binIdxX % 5 == 0) std::cout << "Applying weight as expected. To bin: " << binIdxX << ", " << binIdxY <<
				". To values: (" << h2Sum->GetXaxis()->GetBinCenter(binIdxX) << ", " << h2Sum->GetYaxis()->GetBinCenter(binIdxY) <<
				"  Ratio to be applied to Nue App. = " << contentRatio << std::endl;
			}

			// special case: address bin values above 1.2 GeV....
			else if (binIdxY > h2Ratio->GetNbinsY()){
				if (binIdxX % 15 == 0) { // one for each row.
					std::cout << "Neutrino Energy: " << h2Ratio->GetYaxis()->GetBinCenter(binIdxY) << " GeV" << std::endl;
					std::cout << "Applying weight to Angle (deg) " << h2Sum->GetXaxis()->GetBinCenter(binIdxX) << std::endl;
					std::cout << "binIdxY == " << binIdxY << std::endl;
				}

				// want the bin content from the top row. That would be this value.
				// NOTE: we need to scan the x values still. Only the Y bin is constant
				const double contentRatioTopRow = h2Ratio->GetBinContent(binIdxX, h2Ratio->GetNbinsY());
				const double contentSumNue = h2Sum->GetBinContent(binIdxX, binIdxY);

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
	latex.DrawLatexNDC(0.62, 0.6, Form(beam == "fhc" ? "Neutrino Beam" : "AntiNeutrino Beam"));
	latex.SetTextSize(0.85);
	TLatex ltx2;
	ltx2.SetTextColor(kGray);
	ltx2.DrawLatexNDC(0.62, 0.5, "Asimov A");
	ltx2.SetTextSize(0.85);
	Simulation();

	std::cout << " events integral: " << h2SumRwgt->Integral() << std::endl;

	c.SaveAs(Form("%s/plot_crpa_lfg_fd_%s_prod5.1_enu_theta_nue.png", outDir.c_str(),  beam.c_str()));

}