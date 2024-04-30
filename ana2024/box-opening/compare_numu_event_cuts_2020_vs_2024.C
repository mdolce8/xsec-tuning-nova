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
#include <CAFAna/Analysis/Exposures.h>
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

  std::ofstream fboth;
  std::ofstream f24;
  std::ofstream f20;

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


  std::unordered_map<std::string, ana::Var> vars;

  fboth.open(Form("%s/pass_both_%s.txt", outDir.c_str(), beam.c_str()));
  f24.open(Form("%s/pass_%s_kNumu2024FD.txt", outDir.c_str(), beam.c_str()));
  f20.open(Form("%s/pass_%s_kNumu2020FD.txt", outDir.c_str(), beam.c_str()));

  // these are for checks later...
  int cut_FD24 = 0, cut_FD20 = 0, pass_both = 0, fail_both = 0, unclear = 0;

  vars.try_emplace("kNumu2020FD_Pass",
                   ([&f20, &cut_FD20](const caf::SRProxy *sr) {
                     int pass = 0;
                     if (!kNumu2024FD(sr) && kNumu2020FD(sr)) {
                       cut_FD20++;
                       pass = 1;
                       f20 << sr->hdr.run << "/" << sr->hdr.subrun << "/" << sr->hdr.evt << std::endl;
                     }
                     return pass;
                    }
  ));

  vars.try_emplace("kNumu2024FD_Pass",
                   ([&f24, &cut_FD24](const caf::SRProxy *sr) {
                     int pass =0;
                     if (kNumu2024FD(sr) && !kNumu2020FD(sr)){
                       cut_FD24++;
                       pass = 1;
                       f24 << sr->hdr.run << "/" << sr->hdr.subrun << "/" << sr->hdr.evt << std::endl;
                     }
                     return pass;
                   }
  ));

  vars.try_emplace("Pass_Both",
                   ([&fboth, &pass_both](const caf::SRProxy *sr) {
                       int pass = 0;
                       if (kNumu2024FD(sr) && kNumu2020FD(sr)){
                          pass_both++;
                          pass = 1;
                          fboth << sr->hdr.run << "/" << sr->hdr.subrun << "/" << sr->hdr.evt << std::endl;
                       }
                       return pass;
                   }
                   ));

  std::unordered_map<std::string, HistAxis> map_haxis;
  std::unordered_map<std::string, PredictionNoExtrap*> map_pnxp;

  // two bins: pass or not.
  for (const auto& varPair : vars)
    map_haxis.try_emplace(varPair.first, varPair.first, ana::Binning::Simple(2, 0, 2), varPair.second);

  for (const auto& haxisPair : map_haxis)
    map_pnxp.try_emplace(haxisPair.first, new PredictionNoExtrap(loader, haxisPair.second, kNoCut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2024));

  loader.Go();


  std::unordered_map<std::string, double> map_integrals;

  TFile ofile(Form("%s/compare_%s_numu_event_cuts_2020_vs_2024.root", outDir.c_str(), beam.c_str()), "recreate");
  for (const auto& pnxpPair : map_pnxp) {
    Spectrum s = pnxpPair.second->PredictComponent(calc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth);
    s.SaveTo(&ofile, Form("%s_numucc_all", pnxpPair.first.c_str()));

    TH1D * h = s.ToTH1(beam == "fhc" ? kAna2020FHCPOT : kAna2020RHCPOT);
    h->Write(Form("h_%s_numucc_all", pnxpPair.first.c_str()));

    // write the integrals in here.
    auto integral = h->Integral(1,2); // only the passing events
    std::cout << pnxpPair.first << " integral: " << integral << std::endl;
    map_integrals.try_emplace(pnxpPair.first, integral);
  }

  fboth.close();
  f20.close();
  f24.close();

  std::cout << "good_events: " << pass_both << std::endl;
  std::cout << "cut_FD24: " << cut_FD24 << std::endl;
  std::cout << "cut_FD20: " << cut_FD20 << std::endl;
  std::cout << "fail_both: " << fail_both << std::endl;
  std::cout << "unclear: " << unclear << std::endl;


  // ---- save information into a simple TH1 and ROOT file. ----
  TH1D h = TH1D("h", Form("%s kNumu2020FD vs kNumu2024FD", beam.c_str()), 5, 0, 5);
  h.SetBinContent(1, pass_both);
  h.SetBinContent(2, cut_FD20);
  h.SetBinContent(3, cut_FD24);
  h.SetBinContent(4, unclear);

  h.GetXaxis()->SetBinLabel(1, "pass both");
  h.GetXaxis()->SetBinLabel(2, "pass 2020");
  h.GetXaxis()->SetBinLabel(3, "pass 2024");
  h.GetXaxis()->SetBinLabel(4, "unclear");

  TCanvas c;
  c.SetBottomMargin(0.15);
  h.Draw("hist e");

  h.Write("h_cuts_outcome");
  c.SaveAs(Form("%s/%s_cuts.png", outDir.c_str(), beam.c_str()));

  ofile.Close();

}
