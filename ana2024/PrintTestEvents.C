/*
 * PrintTestEvents.C:
 *
 *   Print out test events for the NOvARwgt test suite
 *   using the NOvA 2018 files and CAFAna version of the tests
 *   (the authoritative pre-NOvARwgt incarnation of the weights).
 *
 *   For GENIE 2 testing events, run it in an old pre-NOvARwgt tag (tested with R19-04-17-2019ana.a).
 *
 *      Author: J. Wolcott <jwolcott@fnal.gov>
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

    const auto & nu = sr->mc.nu[0];
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
    std::vector<const ana::ISyst*> systsDIS;
    systsDIS.push_back(&ana::kDISNuHadroQ1Syst);
    systsDIS.push_back(&ana::kDISNuBarHadroQ0Syst);

    std::vector<const ana::ISyst*> systsRES{};
    systsRES.push_back(&ana::kRESvpvnNuRatioXSecSyst);
    systsRES.push_back(&ana::kRESvpvnNuBarRatioXSecSyst);
    systsRES.push_back(&ana::kRESDeltaScaleSyst);
    systsRES.push_back(&ana::kRESOtherScaleSyst);

      if (nu.mode == caf::kRes || nu.mode == caf::kDIS) {
        std::cout << "ev.nupdg = " << nu.pdg << ";" << std::endl;
        std::cout << "ev.isCC = " << (nu.iscc ? "true" : "false") << ";" << std::endl;
        std::cout << "ev.reaction = novarwgt::" << ModeToEnumName(nu.mode) << ";" << std::endl;
        std::cout << "ev.struckNucl = " << nu.hitnuc << ";" << std::endl;
        std::cout << "ev.A = " << nu.tgtA << ";" << std::endl;
        std::cout << "ev.Enu = " << nu.E << ";" << std::endl;
        std::cout << "ev.q = {" << q.Px() << ", " << q.Py() << ", " << q.Pz() << ", " << q0 << "};"
                  << std::endl;  // note: TLorentzVector initializer is (x, y, z, t)
        std::cout << "ev.y = " << nu.y << ";" << std::endl;
        std::cout << "ev.W = " << W << ";" << std::endl;

        if (nu.mode == caf::kDIS && nu.prefsi.size() == 2) {
          for (const auto &wgtr: systsDIS) {
            std::cout << "...DIS Interaction with 2 Final State Hadrons..." << std::endl;
            for (const auto sigma: {-1, 0, 1, 2, 3}) {
              double wgt = 1.0;
              std::cout << "Nominal weight = " << wgt << std::endl;
              if (nu.hitnuc == 2112 || nu.hitnuc == 2212) {
                std::cout << "neutrino PDG: " << nu.pdg << std::endl;
                std::cout << "ev.npiplus = " << nu.npiplus << ";" << std::endl;
                std::cout << "ev.npizero = " << nu.npizero << ";" << std::endl;
                std::cout << "ev.npiminus = " << nu.npiminus << ";" << std::endl;
                std::cout << "ev.nproton = " << nu.nproton << ";" << std::endl;
                std::cout << "ev.nneutron = " << nu.nneutron << ";" << std::endl;

                wgtr->Shift(sigma, &srMutable, wgt);  // this is the version needed for old CAFAna
                std::cout << "  {novarwgt::GetSystKnobByName(\"" << wgtr->ShortName() << "\"), "
                          << "{" << sigma << ", " << wgt << "}},"
                          << std::endl;
              } // hit nuc
            } // sigma
          } // systs
        } // DIS & Mult. 2.

        // move the RES chunk below into here, it is much cleaner...

      } // RES or DIS


    if (nu.mode == caf::kRes && nu.iscc) {
      std::cout << "Expected syst weights:" << std::endl;
      std::cout << "{" << std::endl;
      for (const auto &wgtr: systsRES) {
        std::cout << "starting loop over...." << wgtr->ShortName() << std::endl;
        for (const auto &sigma: {-0.5, 0.5}) {
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
        }
      }
      std::cout << "}" << std::endl;
    }
  std::cout << std::endl;
     auto wgtVar = 1.0;
     return wgtVar;
});
}

void PrintTestEvents()
{
	ana::Prod51NomLoaders loaders(ana::kNumuConcat, ana::Loaders::kFHC);
	auto & loader = loaders.GetLoader(caf::kNEARDET, ana::Loaders::kMC);

	// dummy spectrum just to trigger the printout
	ana::Spectrum dummy("dummy", ana::Binning::Simple(1, 0, 10), loader,  rwgt::loudXsecWgt, ana::kNoCut);

	loaders.Go();
}
