/*
 *  Plot the values of the weights for the new parameters
 *  ( DISHadronQ{0,1}{Nu,NuBar}Syst )
 *  and look more carefully.
 *  Just looking at the spectra they seem identical.
 *
 * Feb 2024
 * M. Dolce
 */

#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TLatex.h"

#include "3FlavorAna/Cuts/NumuCuts2020.h"
#include "3FlavorAna/Systs/3FlavorAna2020Systs.h"
#include "3FlavorAna/Vars/NumuVars.h"
#include "3FlavorAna/NDFit/Samples/CutsPngCVNOpt.h"
#include "3FlavorAna/NDFit/Systs/LoadSysts.h"

#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Weight.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionXSecTuning.h"
#include "CAFAna/Systs/RESSysts.h"
#include "CAFAna/Systs/DISSysts.h"
#include "CAFAna/Weights/GenieWeights.h"
#include "CAFAna/Weights/PPFXWeights.h"


#include "OscLib/OscCalcPMNSOpt.h"



using namespace ana;

// put title on the right.
void TitleSide(const std::string& title)
{
  TLatex* prelim = new TLatex(.93, .95, Form("%s", title.c_str()));
//  prelim->SetTextColor(kBlue);
  prelim->SetNDC();
  prelim->SetTextSize(1.7/40.);
  prelim->SetTextAngle(270);
  prelim->SetTextAlign(12);
  prelim->Draw();
}


 // for less files do: cafe -bq -l <file-number>

using namespace ana;

// =====================================================================================================
void validate_DISHadronSyst_weights(
                       const double threshold = 5.,              // value of weight you want to dump in the printouts
											 const bool gridSubmission=false
)
// =====================================================================================================
{  

  std::cout << " -- Using prod5.1 data and MC definitions." << std::endl;
  std::cout << " -- GENIE Skew Fix Weight Applied. Using Prod5.1 MC defs." << std::endl;

  const std::string  outDir = "/nova/ana/users/mdolce/preds+spectra/ana2024/validate_DISHadronSyst_weights/test/";

  std::cout << "**********************************************************" << std::endl;

// -------------------------------------------------------------------------------------------
// 		Definitions:

  // mc
	const std::string& fhcMCDef = "prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_numu2020"; // 3F concat
          //"def_snapshot prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_ndphysics_contain_v1"; // from ND Physics concat. stopped using March 8 2022.



	// mc
	const std::string& rhcMCDef = "prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_v1_numu2020prod5.1"; // 3F concat
//          "def_snapshot prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_ndphysics_contain_v1"; // from ND Physics concat. stopped using March 8.


  if ( fhcMCDef.empty() or  rhcMCDef.empty() ) throw std::runtime_error( "MC SAM Definition is empty" );


  Loaders fhcLoaders, rhcLoaders;
  fhcLoaders.SetLoaderPath(fhcMCDef  , caf::kNEARDET, Loaders::kMC   );
  fhcLoaders.SetSpillCut(kStandardSpillCuts );

	rhcLoaders.SetLoaderPath(rhcMCDef  , caf::kNEARDET, Loaders::kMC   );
	rhcLoaders.SetSpillCut(kStandardSpillCuts );

  //Define the space for two variables;
	double q0_min=0; double q0_max=0.8; int binsq0=40;
	double q3_min=0; double q3_max=2.0; int binsq3=40;

  HistAxis haRecoq( "Reco. |#vec{q}| (GeV)"   , Binning::Simple( binsq3, q3_min, q3_max ), kRecoQmag);

  HistAxis haEHad ("Reco. E_{had, vis} (GeV)", Binning::Simple( binsq0, q0_min, q0_max ), kNumuHadVisE );



	const Cut kEvElseFHC = !onlyMuOpt(0.5) && !kMuPrOpt(0.5,0.5)  && !kMuPrEtcOpt(0.5,0.5) && !kMu1PiOpt(0.5,0.7);
	const Cut kEvElseRHC = !onlyMuOpt(0.5) && !kMuPrOpt(0.5,0.5)  && !kMuEtcOpt(0.5,0.5) && !kMu1PiOpt(0.5,0.7);


  // Store here the cuts you will actually use to fill your spectra
  std::vector< std::pair< Cut, std::string > > fhcTopologicalCuts, rhcTopologicalCuts;

  // create an AllNumu as well.
//  SelCuts.emplace_back(std::make_pair(kNumu2020ND, "allNumu"));
  fhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && onlyMuOpt(0.5), "Muon"));
  rhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && onlyMuOpt(0.5), "Muon"));
  fhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && kMuPrOpt(0.5,0.5), "MuPr"));
  rhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && kMuPrOpt(0.5,0.5), "MuPr"));
  fhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && kMu1PiOpt(0.5,0.7), "MuPiEtc"));
  rhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && kMu1PiOpt(0.5,0.7), "MuPiEtc"));
  // some cuts are different for FHC or RHC
  fhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && kMuPrEtcOpt(0.5,0.5), "MuPrEtc"));
  rhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && kMuEtcOpt(0.5,0.5), "MuEtc"));
  fhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && kEvElseFHC, "EvElse"));
  rhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND && kEvElseRHC, "EvElse"));
  // AllNumu too...
  fhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND, "AllNumu"));
  rhcTopologicalCuts.emplace_back(std::make_pair(kNumu2020ND, "AllNumu"));


  // there are the problematic systs....
  std::vector<const ana::ISyst*> systs;
  systs.push_back(&kDISNuHadronQ1Syst);
  systs.push_back(&kDISNuBarHadronQ0Syst);
  systs.push_back(&kRESvpvnRatioNuXSecSyst);
  systs.push_back(&kRESvpvnRatioNubarXSecSyst);


  std::unordered_map<std::string, ana::Var> vars;
  std::unordered_map<double, std::string> shifts
  {
          {+3., "+3"},
          {-3., "-3"},
          {+2.5, "+2.5"},
  };




  for (const auto & syst : systs) {
    for (const auto & shift : shifts) {

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
    } // shift
  } // systs

  std::unordered_map<std::string, Spectrum*> weightSpectraMap;

  const auto wgtBins = ana::Binning::Simple(30, -0.1, 3., {"Weight"});

  // These topologies are the known problematic ones. Start here....(RHC MuPiEtc, FHC EvElse, FHC MuPiEtc)
  // loops through each var, which is a var for each syst's weight!
  for (auto & var : vars){

    // RHC
//    weightSpectraMap.try_emplace(Form("RHC_Muon_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, rhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && onlyMuOpt(0.5)));
//
//    weightSpectraMap.try_emplace(Form("RHC_MuPr_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, rhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && kMuPrOpt(0.5, 0.5)));
//
//    weightSpectraMap.try_emplace(Form("RHC_MuPiEtc_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, rhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && kMu1PiOpt(0.5, 0.7)));
//
//    weightSpectraMap.try_emplace(Form("RHC_MuEtc_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, rhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && kMuEtcOpt(0.5, 0.5)));
//
//    weightSpectraMap.try_emplace(Form("RHC_EvElse_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, rhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && kEvElseRHC));

    weightSpectraMap.try_emplace(Form("RHC_AllNumu_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, rhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND));


    // FHC
//    weightSpectraMap.try_emplace(Form("FHC_Muon_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, fhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && onlyMuOpt(0.5)));
//
//    weightSpectraMap.try_emplace(Form("FHC_MuPr_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, fhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && kMuPrOpt(0.5, 0.5)));
//
//    weightSpectraMap.try_emplace(Form("FHC_MuPiEtc_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, fhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && kMu1PiOpt(0.5, 0.7)));
//
//    weightSpectraMap.try_emplace(Form("FHC_MuPrEtc_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, fhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && kMuPrEtcOpt(0.5, 0.5)));
//
//    weightSpectraMap.try_emplace(Form("FHC_EvElse_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, fhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND && kEvElseFHC));

    weightSpectraMap.try_emplace(Form("FHC_AllNumu_%s", var.first.c_str()), new Spectrum("Weight", wgtBins, fhcLoaders.GetLoader(caf::kNEARDET, Loaders::kMC), var.second, kNumu2020ND));


  } // vars


  fhcLoaders.Go();
	rhcLoaders.Go();



  std::string out_dir;
	if (gridSubmission) out_dir = "."; // for grid: " -o /pnfs/nova/scratch/users/mdolce/ ? "

	else { // for local, provide
		out_dir = outDir + "/" + out_dir;
    std::cout << "local submission...." << std::endl;
	}

	if ( gSystem->AccessPathName( out_dir.c_str() ) ) gSystem->mkdir( out_dir.c_str(), true );

  const std::string fileName = Form("validate_DISHadronSyst_pm3sigma_weights.root");
  const std::string fileNamePath= out_dir + "/" + fileName;
  TFile ofile(fileNamePath.c_str(), "recreate");

  // the string of the map should also contain the systematic name...
  for (const auto & wgtSpecPair : weightSpectraMap){
    TString hc = wgtSpecPair.first.c_str();
    double pot = -5;
    if (hc.Contains("FHC")) pot = kAna2020FHCPOT;
    if (hc.Contains("RHC")) pot = kAna2020RHCPOT;
//    wgtSpecPair.second->SaveTo(&ofile, wgtSpecPair.first.c_str());
    TH1 * h = wgtSpecPair.second->ToTH1(pot);
    h->SetName(Form("h_%s", wgtSpecPair.first.c_str()));
    h->SetDirectory(&ofile);
    h->Write(h->GetName());

    TCanvas c;
    c.Clear(); c.SetCanvasSize(400,400);
    h->Draw("hist");
    TitleSide(wgtSpecPair.first);
    c.SaveAs((out_dir + "/" + "plot_weights_" + wgtSpecPair.first + ".png").c_str());

  } // weightSpectraMap

  ofile.Close();
  std::cout << fileName << " created." << std::endl;


	std::cout << "done." << std::endl;

}
