/*
 * compare_numu_event_cuts_2020_vs_2024.C:
 *    Try to create some printouts of the events and determine which
 *    events are passing which cuts.
 *
 *    Author: M. Dolce
 *    Date:  April 2024
 *
 */

#include <3FlavorAna/Cuts/NumuCuts2024.h>
#include <3FlavorAna/Cuts/NumuCuts2020.h>
#include <CAFAna/Core/Loaders.h>
#include <CAFAna/Prediction/PredictionNoExtrap.h>
#include <TCanvas.h>
#include <OscLib/OscCalcPMNSOpt.h>
#include "3FlavorAna/Cuts/QuantileCuts2020.h"
#include "3FlavorAna/NDFit/Samples/UsefulCutsVars.h"
#include "3FlavorAna/Vars/HistAxes.h"
#include "3FlavorAna/Vars/NumuVars.h"

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"


//  ------ cafe -bq -l <file-number> --stride too?
// NOTE: this uses Prod5.1 !

using namespace ana;

// =====================================================================================================
void compare_numu_event_cuts_2020_vs_2024(const std::string& beam,        // fhc or rhc
                                          const std::string& outDir      // $data/xsec-tuning-nova/ana2024/box-opening ("." for grid: -o $scratch/data )
)
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

  std::cout << "Ana2024 Box Opening........" << std::endl;
  std::cout << "Plotting Prod5.1 FD Numu Quantile(s) Data with 2020 cuts in Reco Enu........" << std::endl;


  // we are only looking at p1-10 data right now -- no new data.

  std::string defMC_Non, defMC_Flux, defMC_Tau;
  if (beam == "fhc") {
    defMC_Non = "prod_sumdecaf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_numu2020"; // 3F concat -- 150 files
    defMC_Flux = "prod_sumdecaf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_fluxswap_fhc_nova_v08_full_v1_numu2020";
    defMC_Tau = "prod_sumdecaf_R20-11-25-prod5.1reco.j_fd_genie_N1810j0211a_tauswap_fhc_nova_v08_full_v1_numu2020";
  }
  else if (beam == "rhc") {
    defMC_Non = "prod_sumdecaf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_v1_numu2020"; // 3F concat -- 150 files
    defMC_Flux = "prod_sumdecaf_R20-11-25-prod5.1reco.j.l_fd_genie_N1810j0211a_fluxswap_rhc_nova_v08_full_v1_numu2020";
    defMC_Tau = "prod_sumdecaf_R20-11-25-prod5.1reco.j_fd_genie_N1810j0211a_tauswap_rhc_nova_v08_full_v1_numu2020";
  }
  else {std::cerr << "Unknown 'beam'. exit..." << std::endl; exit(1);}

  Loaders loader;
  loader.SetLoaderPath(defMC_Non, caf::kFARDET, ana::Loaders::kMC, DataSource::kBeam, Loaders::kNonSwap);
  loader.SetLoaderPath(defMC_Flux, caf::kFARDET, ana::Loaders::kMC, DataSource::kBeam, ana::Loaders::kFluxSwap);
  loader.SetSpillCut(kStandardSpillCuts);

//  std::vector < Cut > cutQuantiles = GetNumuEhadFracQuantCuts2020(beam != "fhc");

  std::unordered_map<std::string, ana::Var> vars;


  int cut_FD24, cut_FD20, good_events = 0;

  // Create Vars of the weights that include print statements
  vars.try_emplace("FD Cut",
                   ([&cut_FD20, &cut_FD24, &good_events](const caf::SRProxy *sr) {
                       if (kNumu2024FD(sr) && !kNumu2020FD(sr))
                       {
                         std::cout << "2024 passed = " << kNumu2024FD(sr) << ", 2020 failed = " << kNumu2020FD(sr) << ": " << sr->hdr.run << "/" << sr->hdr.subrun << "/" << sr->hdr.evt << std::endl;
                         cut_FD24++;
                       }
                       else if (!kNumu2024FD(sr) && kNumu2020FD(sr))
                       {
                         std::cout << "2024 failed = " << kNumu2024FD(sr) << ", 2020 passed = " << kNumu2020FD(sr) << ": "  << sr->hdr.run << "/" << sr->hdr.subrun << "/" << sr->hdr.evt << std::endl;
                         cut_FD20++;
                       }
                       else {
                         std::cerr << "IDK what happened. 2020 == " << kNumu2024FD(sr) << ". 2024 == " << kNumu2020FD(sr) << std::endl;
                       }
                       return -5.;
                   }) // Var lambda

  ); // map try_emplace

  HistAxis haxis("Cut outcome combinations", ana::Binning::Simple(4, 0, 4), vars.at("FD Cut"));


//  Spectrum s(loader.GetLoader(caf::kFARDET, ana::Loaders::kMC, ana::DataSource::kBeam, ana::Loaders::kNonSwap), haxis, kNoCut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2024);

//  Spectrum sflux(loader.GetLoader(caf::kFARDET, ana::Loaders::kMC, ana::DataSource::kBeam, ana::Loaders::kFluxSwap), haxis, kNoCut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2024);

  auto * pnxp = new PredictionNoExtrap(loader, haxis, kNoCut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2024);

  loader.Go();

  TFile ofile(Form("%s/compare_numu_event_cuts_2020_vs_2024.root", outDir.c_str()), "recreate");
  Spectrum s = pnxp->PredictComponent(calc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth);
  s.SaveTo(&ofile, "numucc_all");



  std::cout << "good_events: " << good_events << std::endl;
  std::cout << "cut_FD24: " << cut_FD24 << std::endl;
  std::cout << "cut_FD20: " << cut_FD20 << std::endl;

  TH1D h = TH1D("h", beam.c_str(), 4, 0, 4);
  h.SetBinContent(0, good_events);
  h.SetBinContent(1, cut_FD20);
  h.SetBinContent(2, cut_FD24);

  h.GetXaxis()->SetBinLabel(0, "pass both");
  h.GetXaxis()->SetBinLabel(0, "pass 2020");
  h.GetXaxis()->SetBinLabel(0, "pass 2024");

  TCanvas c;
  c.SetBottomMargin(0.15);
  h.Draw("hist e");

  h.Write("h_cuts_outcome");
  c.SaveAs(Form("%s/%s_cuts.png", outDir.c_str(), beam.c_str()));

  ofile.Close();

}
