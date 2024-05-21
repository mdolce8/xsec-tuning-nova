/*
 *  plot_fd_enu_theta_lowE_nue.C
 *
 *  macro to plot inclusive FD Nue LowE spectra
 *  after each round of cuts.
 *
 *
 *  May 2024
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

// NOTE: LowE sample is ONLY for FHC !

// =====================================================================================================
void plot_fd_enu_theta_lowE_nue(const std::string& beam = "fhc")
// =====================================================================================================
{

  // Asimov A. The 2020 best fit.
  auto calc = new osc::OscCalcPMNSOpt();
  calc->SetL(810);
  calc->SetRho(2.84);
  calc->SetDmsq21(7.53e-5);
  calc->SetTh12(asin(sqrt(0.307)));
  calc->SetDmsq32(2.41e-3);
  calc->SetTh23(asin(sqrt(0.57)));
  calc->SetdCP(0.82*M_PI);
  calc->SetTh13(asin(sqrt(2.18e-2)));

  struct Component
  {
      Flavors::Flavors_t flav;
      Current::Current_t curr;
      Sign::Sign_t sign;
  };


	// nuebarcc_app is small for FHC, but let's include anyway.
  std::map<std::string, Component> flavors = {
					{"nuecc", {Flavors::kAllNuE, Current::kCC, Sign::kBoth}},
          {"beam_nuecc", {Flavors::kNuEToNuE, Current::kCC, Sign::kBoth}},
          {"app_nuecc", {Flavors::kNuMuToNuE, Current::kCC, Sign::kNu}},
          {"app_nuebarcc", {Flavors::kNuMuToNuE, Current::kCC, Sign::kAntiNu}},
  };


	// create Signal, Signal + Beam, ALL
	TH2D h2Signal = TH2D("nue_sig", "FD #nu_{e} Signal", 18, 0., 180, 10, 0.5,1.5);
	TH2D h2SigBeam = TH2D("nue_sig_beam", "FD #nu_{e} Sig + Beam", 18, 0., 180, 10, 0.5, 1.5);
	TH2D h2All = TH2D("nue_all", "FD #nu_{e} All", 18, 0., 180, 10, 0.5, 1.5);

  // dir of the FD Numu Data ROOT files
  const std::string inputDir = "/exp/nova/data/users/mdolce/preds+spectra/ana2024/crpa/";
  const std::string outDir = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/cRPA/plot_fd_enu_theta_lowE_nue/";

	// lowE file. All Spectra are in this file.
	const std::string fname = Form("pred_nxp_fd_%s_prod5.1_enu_theta_lowE_nue.root", beam.c_str());
	const std::string fPath = inputDir + "/" + fname;

  // loop over the flavors.
  for (const auto& pairFlavor : flavors){
		std::cout << "looping through flavor......." << pairFlavor.first << std::endl;

    // these are the enum Flavors category name
    const std::string sName = Form("pred_nxp_enu_theta_nue_%s_all", pairFlavor.first.c_str()); // this is a dir.
    TFile * f = TFile::Open(fPath.c_str());
    Spectrum s = *ana::Spectrum::LoadFrom(f, sName);
    TH2 * h2Flavor = s.ToTH2(kAna2024FHCPOT); // there should be NO RHC.

		// perform the correct operation for appropriate case
		if (pairFlavor.first == "nuecc"){
			std::cout << "Found 'flavor':" << pairFlavor.first << std::endl;
			h2All.Add(h2Flavor);
		}

		else if (pairFlavor.first == "app_nuecc" || pairFlavor.first == "app_nuebarcc"){
			std::cout << "Found 'flavor':" << pairFlavor.first << std::endl;
			h2Signal.Add(h2Flavor);
			h2SigBeam.Add(h2Flavor);
		}

		else if (pairFlavor.first == "beam_nuecc"){
			std::cout << "Found 'flavor':" << pairFlavor.first << std::endl;
			h2SigBeam.Add(h2Flavor);
		}

		else {
			std::cerr << "Error. Unknown flavor: " << pairFlavor.first << std::endl;
			exit(0);
		}

		// draw the flavor
		TCanvas c;
    h2Flavor->Draw("same hist colz");

    TLatex latex;
		latex.SetTextColor(kGray);
		latex.DrawLatexNDC(0.62, 0.6, "Neutrino Beam");
    latex.SetTextSize(0.85);
    TLatex ltx2;
		ltx2.SetTextColor(kGray);
		ltx2.DrawLatexNDC(0.62, 0.5, "Asimov A");
    ltx2.SetTextSize(0.85);
    std::cout << " events integral: " << h2Flavor->Integral() << std::endl;
    Simulation();

    c.SaveAs(Form("%s/plot_fd_%s_prod5.1_enu_theta_lowE_nue_flavor_%s.png", outDir.c_str(),  beam.c_str(), pairFlavor.first.c_str()));

  } // flavors


	// Now plot the histograms that were actually requested:
	std::unordered_map<std::string, TH2*> map_FlavorTH2
	{
		{"Signal", (TH2D*) &h2Signal},
		{"Sig_Beam", (TH2D*) &h2SigBeam},
		{"All", (TH2D*) &h2All},
	};

	// save to ROOT file.
	TFile ofile(Form("%s/th2_fd_%s_enu_theta_lowE_nue_summed.root", outDir.c_str(), beam.c_str()), "recreate");

	for (const auto pairHist : map_FlavorTH2) {
		TCanvas c;
		pairHist.second->Draw("colz");
		c.SetRightMargin(0.1);
		Simulation();
		pairHist.second->SetYTitle("E_{#nu} (GeV)");
		pairHist.second->SetXTitle("#theta (deg)");

		TLatex latex;
		latex.SetTextColor(kGray);
		latex.DrawLatexNDC(0.62, 0.6, Form(beam == "fhc" ? "Neutrino Beam" : "AntiNeutrino Beam"));
		latex.SetTextSize(0.85);
		TLatex ltx2;
		ltx2.SetTextColor(kGray);
		ltx2.DrawLatexNDC(0.62, 0.5, "Asimov A");
		ltx2.SetTextSize(0.85);

		std::cout << " Summed events integral: " << pairHist.second->Integral() << std::endl;

		c.SaveAs(Form("%s/plot_fd_%s_enu_theta_lowE_nue_summed_%s.png", outDir.c_str(), beam.c_str(), pairHist.first.c_str()));

		pairHist.second->SetDirectory(&ofile);
		pairHist.second->Write(pairHist.second->GetName());
	}

	ofile.Write();
	ofile.Close();
}