/*
 *  apply_crpa_lfg_ratio.C
 *
 * 	Apply the ratio from the cRPA paper to our
 *  NOvA FD Nue (appeared + background) spectra.
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

// TODO: read in ROOT file of the Summed TH2 and the RR ratio from paper.
// TODO: figure out how to apply the ratio to the prediction.

using namespace ana;

// =====================================================================================================
void apply_crpa_lfg_ratio(const std::string& beam)
// =====================================================================================================
{

	const std::string outDir = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/cRPA/apply_crpa_lfg_nue_ratio/";

	// ROOT file Summed TH2
	TFile fnue("/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/cRPA/plot_fd_enu_theta_nue/th2_fd_fhc_enu_theta_nue_summed.root", "read");

	// I think this is the right ratio to read in at the moment?
	TFile fcRPA("/exp/nova/data/users/mdolce/xsec-tuning-nova/RR_CRPA_C_LFG_O_e_ae.root", "read");


	// TH2 sum of the nue_app and nuebar_app portions.
	TH2D * h2Sum = (TH2D*) fnue.Get("nue_app");

	// make a clone for the reweighted histogram
	TH2D * h2SumRwgt = (TH2D*) h2Sum->Clone("nue_app_cRPA");
	h2SumRwgt->Reset("ICESM"); // maintain the binning, just clear the content
	assert (h2SumRwgt->Integral() == 0);

	// TH2 of the cRPA / LFG ratio for nue
	TH2D * h2Ratio = (TH2D*) fcRPA.Get("CC_RPA_LFG_O_e_ae.root");

	// TODO: is this right? start at 0 ?
	for (unsigned int binIdxX = 0; binIdxX <= h2Sum->GetNbinsX(); binIdxX++){
		for (unsigned int binIdxY = 0; binIdxY <= h2Sum->GetNbinsY(); binIdxY++){

			const double contentRatio = h2Ratio->GetBinContent(binIdxX, binIdxY);
			const double contentSumNue = h2Sum->GetBinContent(binIdxX, binIdxY);

			const double contentRwgt = contentRatio * contentSumNue; // this is the value after re-weighting to cRPA.

			h2SumRwgt->SetBinContent(binIdxX, binIdxY, contentRwgt);

			// TODO: still leaves us with the bin values above 1.2 GeV....

		} // binIdxY
	} // binIdxX


	// Draw the Reweighted histogram now.
	TCanvas c;

	h2SumRwgt->Scale(beam == "fhc" ? kAna2024FHCPOT : kAna2024RHCPOT);

	h2SumRwgt->Draw("same hist colz");
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

	c.SaveAs(Form("%s/plot_fd_%s_prod5.1_enu_theta_nue.png", outDir.c_str(),  beam.c_str()));


}