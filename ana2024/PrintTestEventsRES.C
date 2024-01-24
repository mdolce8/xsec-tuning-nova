/*
 * PrintTestEventsRES.C:
 *
 *   Print out test events for the NOvARwgt test suite
 *   of the Prod5.1 RES systematics from the ND fitting.
 *
 *      Author: M. Dolce <mdolce@fnal.gov>
 *
 *
 *  To run this script, many libraries are needed:
 *  $ cafe --libdir /cvmfs/nova.opensciencegrid.org/externals/genie/v3_04_00/Linux64bit+3.10-2.17-e20-inclxx-prof/lib  --lib GRwFwk --lib GRwIO --lib GRwClc  --lib GFwMsg --lib GFwReg --lib GFwAlg --lib GFwInt --lib GFwGHEP --lib GFwNum --lib GFwUtl --lib GFwParDat --lib GFwEG --lib GFwNtp --lib GPhXSIg --lib GPhPDF --lib GPhNuclSt --lib GPhCmn --lib GPhDcy --lib GPhHadTransp --lib GPhHadnz --lib GPhDeEx --lib GPhAMNGXS --lib GPhAMNGEG --lib GPhChmXS --lib GPhCohXS --lib GPhCohEG --lib GPhDISXS --lib GPhDISEG --lib GPhDfrcXS --lib GPhDfrcEG --lib GPhHELptnXS --lib GPhHELptnEG --lib GPhIBDXS --lib GPhIBDEG --lib GPhHadTens --lib GPhMNucXS --lib GPhMNucEG --lib GPhMEL --lib GPhNuElXS --lib GPhNuElEG --lib GPhQELXS --lib GPhQELEG --lib GPhResXS --lib GPhResEG --lib GPhStrXS --lib GPhStrEG --lib GPhHEDISXS --lib GPhHEDISEG --lib GPhNDcy --lib GPhNNBarOsc --lib GPhBDMXS --lib GPhBDMEG --lib GPhHNL --lib GTlGeo --lib GTlFlx   --libdir $CLHEP_LIB_DIR --lib CLHEP  PrintTestEventsRES.C
 */

#include <iostream>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"

#include "nugen/NuReweight/ReweightLabels.h"

#include "3FlavorAna/Vars/NumuVars.h"
#include "3FlavorAna/NDFit/Systs/LoadSysts.h"

#include "CAFAna/Analysis/Prod51Loaders.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/ISyst.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Systs/XSecSystLists.h"
#include "CAFAna/Systs/RESSysts.h"
#include "CAFAna/Systs/DISSysts.h"
#include "CAFAna/Vars/TruthVars.h"
// disable the following in really old CAFAna as needed
#include "CAFAna/Weights/XsecTunes.h"
#include "CAFAna/Weights/GenieWeights.h"

#include "StandardRecord/Proxy/SRProxy.h"


namespace rwgt
{
	std::string ModeToEnumName(int mode)
	{
		switch(mode)
		{
			case caf::kQE:
				return "kScQuasiElastic";

			case caf::kMEC:
				return "kScMEC";

			case caf::kDIS:
				return "kScDeepInelastic";

			case caf::kRes:
				return "kScResonant";

			case caf::kCoh:
				return "kScCoherent";

			default:
				std::cout << "Unknown mode: " << mode << std::endl;
				abort();

		}
	}

// TODO: these particles are post-FSI. They should be pre-FSI for the "new" systs.
// Need the Pre-FSI energies, use prefsi using prim.
    namespace fsi_unc
    {
        double VarPostFSIMultiplicity(const std::vector<int>& pdgs,
                               const caf::SRNeutrinoProxy *nu)
        {
          unsigned int nPart = 0;
          for (const auto &part : nu->prim)
          {
            if (std::find(pdgs.begin(), pdgs.end(), part.pdg) != pdgs.end())
              nPart++;
          }
          return nPart;
        };

        double VarPreFSIMultiplicity(const std::vector<int>& pdgs,
                                     const caf::SRNeutrinoProxy  *nu)
        {
          unsigned nPart = 0;
          for (const auto &part : nu->prefsi)
          {
            if (std::find(pdgs.begin(), pdgs.end(), part.pdg) != pdgs.end())
              nPart++;
          }
          return nPart;
        }



        double VarPostFSIKE(const std::vector<int>& pdgs,
                     const caf::SRNeutrinoProxy *nu){

          double KE = 0;
          for (const auto &part : nu->prim) //prefsi , prim is postFSI (this gets fed into GEANT).
          {
            if (std::find(pdgs.begin(), pdgs.end(), part.pdg) != pdgs.end())
              KE += (1. - 1./part.p.Gamma()) * part.p.E;
          }
          // no particles altogether.  not interesting
          if (KE == 0)
            return -1.;
          return KE;
        }

        double VarPreFSIKE(const std::vector<int>& pdgs,
                           const caf::SRNeutrinoProxy *nu){

          double KE = 0;
          for (const auto &part :nu->prefsi)
          {
            if (std::find(pdgs.begin(), pdgs.end(), part.pdg) != pdgs.end())
              KE +=(1. -1./part.p.Gamma()) * part.p.E;
          }

          // no particle, nothing here.
          if (KE == 0)
            return -1;
          return KE;
        }
    }


	std::string KnobToKnobName(rwgt::EReweightLabel knob)
	{
		return genie::rew::GSyst::AsString( static_cast<genie::rew::GSyst_t>(knob) );
	}

	const ana::Var loudXsecWgt([](const caf::SRProxy * sr)
  {
    if (sr->mc.nnu != 1)
      return 1.;

    // won't have any syst weights
    if (!sr->mc.nu[0].isvtxcont)
      return 1.;

    // don't go too crazy
    static std::size_t n = 0;
    n++;
    if (n>50)
    {
      std::cout << "that's plenty of events. exiting...." << std::endl;
                      exit(0);
    }

    const caf::SRNeutrinoProxy & nu = sr->mc.nu[0];
    double El = ana::kTrueMuonE(sr);
    double q0 = ana::kTrueQ0(sr);
    double q3 = ana::kTrueQ3(sr);
    double W = sqrt(nu.W2);
    double maqeWgt = 1;
    auto& srMutable = const_cast<caf::SRProxy&>(*sr);
    auto& nuMutable = const_cast<caf::SRNeutrinoProxy&>(nu);

    TLorentzVector k1(nu.p.px, nu.p.py, nu.p.pz, nu.p.E);
    TLorentzVector k2(nu.prim[0].p.px, nu.prim[0].p.py,  nu.prim[0].p.pz, nu.prim[0].p.E);

    // we want to rotate everything so that we can assume the neutrino was incident in the +z-direction
    const TVector3 zhat(0, 0, 1);
    auto r1 = k1.Vect();
    auto rotationAxis = r1.Cross(zhat);
    auto rotationAngle = r1.Angle(zhat);
    k1.Rotate(rotationAngle, rotationAxis);
    k2.Rotate(rotationAngle, rotationAxis);

    auto q = k1 - k2;

    // Mike & Maria's new RES systematics
    std::vector<const ana::ISyst*> systsRES{};
    systsRES.push_back(&ana::kRESvpvnNuRatioXSecSyst);
    systsRES.push_back(&ana::kRESvpvnNuBarRatioXSecSyst);
    systsRES.push_back(&ana::kRESDeltaScaleSyst);
    systsRES.push_back(&ana::kRESOtherScaleSyst);




    if (nu.mode == caf::kRes) {
      std::cout << "ev.nupdg = " << nu.pdg << ";" << std::endl;
      std::cout << "ev.isCC = " << (nu.iscc ? "true" : "false") << ";" << std::endl;
      std::cout << "ev.reaction = novarwgt::" << ModeToEnumName(nu.mode) << ";" << std::endl;
      std::cout << "ev.struckNucl = " << nu.hitnuc << ";" << std::endl;
      std::cout << "ev.A = " << nu.tgtA << ";" << std::endl;
      std::cout << "ev.Z = " << nu.tgtZ << ";" << std::endl;
      std::cout << "ev.Enu = " << nu.E << ";" << std::endl;
      std::cout << "ev.q = {" << q.Px() << ", " << q.Py() << ", " << q.Pz() << ", " << q0 << "};"
                << std::endl;  // note: TLorentzVector initializer is (x, y, z, t)
      std::cout << "ev.y = " << nu.y << ";" << std::endl;
      std::cout << "ev.W = " << W << ";" << std::endl;
      std::cout << "ev.resnum = " << nu.resnum << ";" << std::endl;
      auto wgtCV = ana::kXSecCVWgt2020GSFProd51(sr);
      std::cout << "---------------−" << std::endl;
      std::cout << "CV kXSecCVWgt2020GSFProd51 wgt = " << wgtCV << std::endl;
      std::cout << "ev.npiplus = " << nu.npiplus << ";" << std::endl;
      std::cout << "ev.npizero = " << nu.npizero << ";" << std::endl;
      std::cout << "ev.npiminus = " << nu.npiminus << ";" << std::endl;
      std::cout << "ev.nproton = " << nu.nproton << ";" << std::endl;
      std::cout << "ev.nneutron = " << nu.nneutron << ";" << std::endl;



      // printout the particle mult. and energy.
      std::map<std::string, double> nutruthVars
              {
//                    {"chgpi_mult", fsi_unc::MultiplicityVar({211, -211})},
                      {"pi+_KE (Pre-FSI)", fsi_unc::VarPreFSIKE({211}, &nu)},
                      {"pi0_KE (Pre-FSI)", fsi_unc::VarPreFSIKE({111}, &nu)},
                      {"pi-_KE (Pre-FSI)", fsi_unc::VarPreFSIKE({-211}, &nu)},
                      {"n_KE (Pre-FSI)", fsi_unc::VarPreFSIKE({2112}, &nu)},
                      {"p_KE (Pre-FSI)", fsi_unc::VarPreFSIKE({2212}, &nu)},

                      {"pi+_mult (Pre-FSI)", fsi_unc::VarPreFSIMultiplicity({211}, &nu)},
                      {"pi0_mult (Pre-FSI)", fsi_unc::VarPreFSIMultiplicity({111}, &nu)},
                      {"pi-_mult (Pre-FSI)", fsi_unc::VarPreFSIMultiplicity({-211}, &nu)},
                      {"n_mult (Pre-FSI)", fsi_unc::VarPreFSIMultiplicity({2112}, &nu)},
                      {"p_mult (Pre-FSI)", fsi_unc::VarPreFSIMultiplicity({2212}, &nu)},

                      {"pi+_KE (Post-FSI)", fsi_unc::VarPostFSIKE({211}, &nu)},
                      {"pi0_KE (Post-FSI)", fsi_unc::VarPostFSIKE({111}, &nu)},
                      {"pi-_KE (Post-FSI)", fsi_unc::VarPostFSIKE({-211}, &nu)},
                      {"n_KE (Post-FSI)", fsi_unc::VarPostFSIKE({2112}, &nu)},
                      {"p_KE (Post-FSI)", fsi_unc::VarPostFSIKE({2212}, &nu)},

                      {"pi+_mult (Post-FSI)", fsi_unc::VarPostFSIMultiplicity({211}, &nu)},
                      {"pi0_mult (Post-FSI)", fsi_unc::VarPostFSIMultiplicity({111}, &nu)},
                      {"pi-_mult (Post-FSI)", fsi_unc::VarPostFSIMultiplicity({-211}, &nu)},
                      {"n_mult (Post-FSI)", fsi_unc::VarPostFSIMultiplicity({2112}, &nu)},
                      {"p_mult (Post-FSI)", fsi_unc::VarPostFSIMultiplicity({2212}, &nu)},
//                      {"Enu", ana::kTrueE_NT},
//                      {"Q2", ana::kTrueQ2_NT},
//                      {"W", ana::kTrueW_NT},
//                      {"z", ana::kTruePartonZ_NT},
//                      {"wgt", ana::NuTruthVarFromNuTruthWeight(ana::kXSecCVWgt2020GSF_NT)},
              }; // nutruthVars

      for (const auto & partPair : nutruthVars){
        std::cout << "{" << partPair.first << ", " << partPair.second << "};" << std::endl;
      }



      if (nu.iscc) {
        std::cout << "Expected syst weights:" << std::endl;
        std::cout << "{" << std::endl;
        for (const auto &wgtr: systsRES) {
          std::cout << "starting loop over...." << wgtr->ShortName() << std::endl;
          for (const auto &sigma: {-1., -0.5, 0.5, 1.}) {
            double wgt = 1.0;
            std::cout << "Nominal weight = " << wgt << std::endl;
            //							wgtr->Shift(sigma, &nuMutable, wgt);  // newer CAFAna works with this one
            wgtr->Shift(sigma, &srMutable, wgt);  // this is the version needed for old CAFAna
            std::string name = wgtr->ShortName();

            // 	// let's not litter the output with weights that aren't doing anything.
            // 	// we'll just test the ones that have real effects
            // 	// (even though in principle it would be better to test that the ones
            // 	//  that aren't _supposed_ to be doing anything really aren't)
            if (wgt == 1.0) {
              std::cout << "wgt=1.0. continue" << std::endl;
              continue;
            }
            std::cout << "  {novarwgt::GetSystKnobByName(\"" << name << "\"), "
                      << "{" << sigma << ", " << wgt << "}},"
                      << std::endl;
          } //sigma
        } //systs
        std::cout << "}" << std::endl;
      } // CC
    } // RES

    std::cout << std::endl;
    auto wgtVar = 1.0;
    return wgtVar;
});
}

void PrintTestEventsRES()
{
	ana::Prod51NomLoaders loaders(ana::kNumuConcat, ana::Loaders::kFHC);
	auto & loader = loaders.GetLoader(caf::kNEARDET, ana::Loaders::kMC);

	// dummy spectrum just to trigger the printout
	ana::Spectrum dummy("dummy", ana::Binning::Simple(1, 0, 10), loader,  rwgt::loudXsecWgt, ana::kNoCut);

	loaders.Go();
}
