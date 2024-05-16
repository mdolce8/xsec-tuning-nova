/*
 *  generate_fd_lowE_nue_spectra_crpa_study.C:
 *    Create Spectrum objects from Prod5.1 MC
 *    in Ev & theta for the lowE Nue sample.
 *
 *
 *    Author: M. Dolce
 *    Date:  May 2024
 *
 */

#include <CAFAna/Prediction/PredictionNoExtrap.h>
#include <CAFAna/Core/Loaders.h>
#include <OscLib/OscCalcPMNSOpt.h>
#include "3FlavorAna/Cuts/QuantileCuts2024.h"
#include "3FlavorAna/Cuts/NueCuts2024.h"
#include "3FlavorAna/NDFit/Samples/UsefulCutsVars.h"
#include "3FlavorAna/Vars/HistAxes.h"
#include "3FlavorAna/Vars/NumuVars.h"

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Weights/XsecTunes.h"
#include "CAFAna/Weights/PPFXWeights.h"

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"


//  ------ cafe -bq -l <file-number> --stride too?
// NOTE: this uses Prod5.1 !

using namespace ana;

// =====================================================================================================
void generate_fd_lowE_nue_spectra_crpa_study(const std::string& beam,        // fhc or rhc
																				const std::string& outDir,      // $data/preds+spectra/ana2024/crpa/ ("." for grid: -o $scratch/data )
																			 	const bool gridSubmission = false
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

  struct Component
  {
      Flavors::Flavors_t flav;
      Current::Current_t curr;
      Sign::Sign_t sign;
  };

  std::map<std::string, Component> flavors = {
          //{"nuecc", {Flavors::kAllNuE, Current::kCC, Sign::kBoth}},
          {"beam_nuecc", {Flavors::kNuEToNuE, Current::kCC, Sign::kBoth}},
          {"app_nuecc", {Flavors::kNuMuToNuE, Current::kCC, Sign::kNu}},
          {"app_nuebarcc", {Flavors::kNuMuToNuE, Current::kCC, Sign::kAntiNu}},
//          {"numucc", {Flavors::kAllNuMu, Current::kCC, Sign::kBoth}},
//          {"nc", {Flavors::kAll, Current::kNC, Sign::kBoth}}
  };


  std::cout << "Ana2024 CRPA study........" << std::endl;
  std::cout << "Making spectra for FD Nue sample in Ev vs. theta space........" << std::endl;


  // we are only looking at p1-10. no new periods

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
//  loader.SetLoaderPath(defMC_Tau, caf::kFARDET, ana::Loaders::kMC, DataSource::kBeam, ana::Loaders::kTauSwap);

  loader.SetSpillCut(kStandardSpillCuts);


	const ana::NuTruthVar kTrueElectronThetaST([](const caf::SRNeutrinoProxy * nu) {
			if (abs(nu->pdg) != 12 || !nu->iscc)
				return -5.0;
			int nprims = nu->prim.size();
			for (int iprim = 0; iprim < nprims; iprim++)
			{
				if (abs(nu->prim[iprim].pdg) == 11)
				{

					TVector3 edir = nu->prim[iprim].p.Vect();
					TVector3 beamdir = ana::NuMIBeamDirection(caf::kNEARDET);

					float theta = edir.Unit().Dot(beamdir.Unit());
					return (float) TMath::ACos(theta) * 180 / TMath::Pi();
				}
			}
			return -5.0;
	});

	const ana::Var kTrueElectronTheta = ana::VarFromNuTruthVar(kTrueElectronThetaST);



  std::map<std::string, const PredictionNoExtrap*> predNxp;

	const ana::Binning bins_theta = ana::Binning::Simple(18, 0., 180.);
	const ana::Binning bins_ENu = ana::Binning::Simple(40, 0., 4.0);
	HistAxis ha_ENu_Theta("#theta (deg)", bins_theta, kTrueElectronTheta,
												"E_{#nu} (GeV)", bins_ENu, kCCE);

	// any relevant appeared Nue.
	const Cut kNueAll = kNue2024FDAllSamples || kNue2024FDLE;

	// Cut is both Core and Peripheral OR LowE sample.
  predNxp.try_emplace("pred_nxp_enu_theta_nue", new PredictionNoExtrap(loader, ha_ENu_Theta, , kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2024));

  loader.Go();



  std::string out_dir;
  if (gridSubmission) out_dir = "."; // for grid: " -o /pnfs/nova/scratch/users/mdolce/<path>/<outDir> "
  else {out_dir = outDir;}   // for local
  if ( gSystem->AccessPathName( out_dir.c_str() ) ) gSystem->mkdir( out_dir.c_str(), true );


  // save the spectra to each ROOT file
  int quantileCount = 1;
  for (const auto &prednxp : predNxp){
    std::string fileName = Form("pred_nxp_fd_%s_prod5.1_enu_theta_nue.root",  beam.c_str());
    const std::string& finalOutDir = out_dir + "/" + fileName;
    TFile ofile(Form("%s", finalOutDir.c_str()), "recreate");

    // save the PredNoExtrap info
    std::cout << "saving PredNxp: " << std::endl;
    prednxp.second->SaveTo(&ofile, prednxp.first);

    // save the Spectrum info.
    std::cout << "saving Spectrum: " << std::endl;
    std::map<std::string, Spectrum> specs_preds;
    for (auto flavor : flavors){
      auto tmp_spec = prednxp.second->PredictComponent(calc, flavor.second.flav, flavor.second.curr, flavor.second.sign);
      std::string key = prednxp.first+"_"+flavor.first+"_all";
      specs_preds.insert({key, tmp_spec});
    } // flavors

    for(auto spec: specs_preds)
      spec.second.SaveTo(&ofile, spec.first.c_str());

    std::cout << "saving PredNxp: " << std::endl;

    ofile.Close();
    std::cout << "Wrote file: " << finalOutDir << std::endl;

    quantileCount++;
  } // preds

}
