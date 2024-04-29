/*
 * compare_numu_event_cuts_2020_vs_2024.C:
 *    Try to create some printouts of the events and determine which
 *    events are passing which cuts.
 *
 *    Author: M. Dolce
 *    Date:  April 2024
 *
 */

#include "3FlavorAna/Cuts/QuantileCuts2020.h"
#include "3FlavorAna/Cuts/NumuCuts2020.h"
#include "3FlavorAna/NDFit/Samples/UsefulCutsVars.h"
#include "3FlavorAna/Vars/HistAxes.h"
#include "3FlavorAna/Vars/NumuVars.h"

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Cuts/SpillCuts.h"

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

  std::string defData;
  if (beam == "fhc") defData = "tbezerra_prod_sumrestricteddecaf_R20-11-25-prod5.1reco_fd_numi_fhc_p1-10_v1_goodruns_numu2024"; // 3F concat -- 600 files
  else if (beam == "rhc") defData = ""; // 3F concat
  else {std::cerr << "Unknown 'beam'. exit..." << std::endl; exit(1);}

  SpectrumLoader dataLoader(defData);
  dataLoader.SetSpillCut(kStandardSpillCuts);

  std::vector<Cut> cutQuantiles = GetNumuEhadFracQuantCuts2020(beam != "fhc");

  std::unordered_map<std::string, ana::Var> vars;
  std::unordered_map<std::string, int> map_cut_status{};


  // Create Vars of the weights that include print statements
  vars.try_emplace(syst->ShortName() + "_" + shift.second,
                   ([threshold, syst, shift](const caf::SRProxy *sr) {
                       double weight = 1.0;
                       auto mutableSR = const_cast<caf::SRProxy *>(sr);
                       syst->Shift(shift.first, mutableSR, weight);
                       if (weight > threshold || weight < 0.0) {
                         std::cout << shift.second << "-sigma shift is bananas from systematic shift: " << syst->ShortName() << std::endl;
                         std::cout << "The Weight is == " << weight << std::endl;
                         std::cout << "Event = " << sr->hdr.run << "/" << sr->hdr.subrun << "/" << sr->hdr.evt
                                   << std::endl;
                         // if you want you can dump out truth information from sr->mc.nu[0] as well
                       }
                       return weight;
                   }) // Var lambda

  ); // map emplace

}
