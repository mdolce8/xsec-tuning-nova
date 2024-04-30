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
  loader.SetLoaderPath(defMC_Non, caf::kFARDET, ana::Loaders::kMC);
  loader.SetLoaderPath(defMC_Flux, caf::kFARDET, ana::Loaders::kMC, DataSource::kBeam, ana::Loaders::kFluxSwap);
  loader.SetSpillCut(kStandardSpillCuts);

  std::vector < Cut > cutQuantiles = GetNumuEhadFracQuantCuts2020(beam != "fhc");

  // todo: do we need spill cuts too?
  std::unordered_map<std::string, ana::Var> vars;
  std::unordered_map<std::string, ana::Cut> map_cut_names
          {
//                  {"kNumu", kStandardSpillCuts},
                  {"kNumu2020DecafCut", kNumu2020FDDecafCut},
                  {"kNumu2024DecafCut", kNumu2024FDDecafCut},
          };
  std::unordered_map<std::string, int> map_cut_status{};

  int cut_FD24, cut_FD20, good_events = 0;

  const Cut kMyNumu2020CosRejLoose(
          [](const caf::SRProxy* sr)
          {
              return kNumuContPID(sr) > 0.4 && kCVNm_looseptp(sr) > 0.;
          });

  // Create Vars of the weights that include print statements
  vars.try_emplace("FD Cut",
                   ([&cut_FD20, &cut_FD24, &good_events](const caf::SRProxy *sr) {
                       if (kNumu2024FD(sr)) {
                         good_events++;
                       }
                       else if (kNumu2024FD(sr) && !kNumu2020FD(sr))
                       {
                         std::cout << "2024 passed, 2020 failed: " << sr->hdr.run << "/" << sr->hdr.subrun << "/" << sr->hdr.evt << std::endl;
                         cut_FD24++;
                       }
                       else if (!kNumu2024FD(sr) && kNumu2020FD(sr))
                       {
                         std::cout << "2024 failed, 2020 passed: " << sr->hdr.run << "/" << sr->hdr.subrun << "/" << sr->hdr.evt << std::endl;
                         cut_FD20++;
                       }
                       else {
                         std::cerr << "IDK what happened." << std::endl;
                       }
                       return -5.;
                   }) // Var lambda

  ); // map try_emplace

  HistAxis haxis("label", ana::Binning::Simple(10, 1,10), vars.at("FD Cut"));

//  Spectrum s(loader.GetLoader(caf::kFARDET, ana::Loaders::kMC, ana::DataSource::kBeam, ana::Loaders::kNonSwap), haxis, kNoCut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2024);

//  Spectrum sflux(loader.GetLoader(caf::kFARDET, ana::Loaders::kMC, ana::DataSource::kBeam, ana::Loaders::kFluxSwap), haxis, kNoCut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2024);

  auto * pnxp = new PredictionNoExtrap(loader, haxis, kNoCut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2024);

  TH1D h = TH1D("h", "h", 10, 1, 10);
  h.SetBinContent(1, good_events);
  h.SetBinContent(2, cut_FD20);
  h.SetBinContent(3, cut_FD24);

  h.Draw("hist e");



}
