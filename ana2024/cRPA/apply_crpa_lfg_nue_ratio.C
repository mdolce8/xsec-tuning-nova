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


	// we NEED nue_app signal for RHC.
	// nuebar_app is small for FHC, but let's include anyway.
  std::map<std::string, Component> flavors = {
          //{"nuecc", {Flavors::kAllNuE, Current::kCC, Sign::kBoth}},
//          {"beam_nuecc", {Flavors::kNuEToNuE, Current::kCC, Sign::kBoth}},
          {"app_nuecc", {Flavors::kNuMuToNuE, Current::kCC, Sign::kNu}},
          {"app_nuebarcc", {Flavors::kNuMuToNuE, Current::kCC, Sign::kAntiNu}},
//          {"numucc", {Flavors::kAllNuMu, Current::kCC, Sign::kBoth}},
//          {"nc", {Flavors::kAll, Current::kNC, Sign::kBoth}}
  };


	// create sum of the nue_app and nuebar_app portions.
	TH2D h2Sum = TH2D(Form("nue_app_%s", beam.c_str()), "FD #nu_{e}", 18, 0., 180, 30, 0., 3.6);

  // dir of the FD Numu Data ROOT files
  const std::string inputDir = "/exp/nova/data/users/mdolce/preds+spectra/ana2024/crpa/";
  const std::string outDir = "/exp/nova/data/users/mdolce/xsec-tuning-nova/plots/ana2024/cRPA/plot_fd_enu_theta_nue/";


  // loop over the flavors.
  for (const auto& pairFlavor : flavors){
		std::cout << "looping through flavor......." << pairFlavor.first << std::endl;
    const std::string fname = Form("pred_nxp_fd_%s_prod5.1_enu_theta_nue.root", beam.c_str());
    const std::string fPath = inputDir + "/" + fname;


    // these are the "app_nuecc" category
    const std::string sName = Form("pred_nxp_enu_theta_nue_%s_all", pairFlavor.first.c_str()); // this is a dir.

    TFile * f = TFile::Open(fPath.c_str());


    Spectrum s = *ana::Spectrum::LoadFrom(f, sName);

    // do plotting
    TCanvas c;

    TH2 * h2Flavor = s.ToTH2(beam == "fhc" ? kAna2024FHCPOT : kAna2024RHCPOT);

		// add this flavor histogram to the h2Sum
		h2Sum.Add(h2Flavor);

    h2Flavor->Draw("same hist colz");


//    h2->SetMaximum(h2->GetMaximum() * 1.5);

    TLatex latex;
		latex.SetTextColor(kGray);
		latex.DrawLatexNDC(0.62, 0.6, Form(beam == "fhc" ? "Neutrino Beam" : "AntiNeutrino Beam"));
    latex.SetTextSize(0.85);
    TLatex ltx2;
		ltx2.SetTextColor(kGray);
		ltx2.DrawLatexNDC(0.62, 0.5, "Asimov A");
    ltx2.SetTextSize(0.85);

    std::cout << " events integral: " << h2Flavor->Integral() << std::endl;

//    TLegend leg(0.1, 0.62, 0.9, 0.9);
//    leg.SetFillStyle(0);
//    leg.AddEntry(h2, Form("%s", cutPair.first.c_str()), "l");
//    leg.Draw("same");

    Simulation();

    c.SaveAs(Form("%s/plot_fd_%s_prod5.1_enu_theta_nue.png", outDir.c_str(),  beam.c_str()));

  } // flavors


	TCanvas c;
	h2Sum.Draw("colz");
	c.SetRightMargin(0.1);
	Simulation();
	h2Sum.SetYTitle("E_{#nu} (GeV)");
	h2Sum.SetXTitle("#theta (deg)");

	TLatex latex;
	latex.SetTextColor(kGray);
	latex.DrawLatexNDC(0.62, 0.6, Form(beam == "fhc" ? "Neutrino Beam" : "AntiNeutrino Beam"));
	latex.SetTextSize(0.85);
	TLatex ltx2;
	ltx2.SetTextColor(kGray);
	ltx2.DrawLatexNDC(0.62, 0.5, "Asimov A");
	ltx2.SetTextSize(0.85);

	std::cout << " Summed events integral: " << h2Sum.Integral() << std::endl;


	c.SaveAs(Form("%s/plot_fd_%s_enu_theta_nue_summed.png", outDir.c_str(), beam.c_str()));



}