/*
 * plot_specific_systematic_errorbands_with_spectra.C
 *
 * Macro to plot error bands of specific systematics.
 * Goal is to identify what is wrong with PlotWithSystErrorBand()
 * and, if desired, plot data, MC, and samples on top.
 * in the different topologies.
 *
 * July, 2021
 * M Dolce
 */

#include <iostream>

#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Fit/Bayesian1DMarginal.h"
#include "CAFAna/Fit/BayesianSurface.h"
#include "CAFAna/Fit/MCMCSamples.h"

#include "3FlavorAna/Systs/DummyRockScaleSyst.h"
#include "3FlavorAna/MCMC/VarsWithPriors.h"
#include "3FlavorAna/MCMC/MCMC3FShared.h"
#include "3FlavorAna/NDFit/InitializeFit.h"
#include "3FlavorAna/NDFit/LoadTopoPreds.h"
#include "3FlavorAna/NDFit/Systs/LoadSysts.h"
#include "3FlavorAna/NDFit/NDFitHelper.h"


#include "OscLib/OscCalcAnalytic.h"

#include "Utilities/func/MathUtil.h"
#include "Utilities/rootlogon.C"

#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THashList.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TProfile.h"



namespace test
{
    void YouAreHere(const int & location)
    {
      std::cout << "You are at step number: " << location << "." << std::endl;
    }
}



//helpful info: ShortName() -- "MaCCRES" , "ppfx_hadp_beam_nd_pc00". (found in chi2-fhc-shifts-summary.txt)
 using namespace ana;

// =====================================================================================================
void plot_specific_systematic_errorbands_with_spectra(
																				              const std::string& samplesFilename,     // $ana/mcmc/samples/data-samples OR fakedata-samples
																				              const std::string& dataFilename,        // $ana/mcmc/data OR fake-data
																				              const std::string& outDirPrefix,        // $ana/mcmc/plot/evaluate-systematics-plots/specific-syst-errorbands/
																				              const std::string& predType,            // "prod5" or "pionrej" at the moment
																				              const std::string& systShortName,       // the ShortName() for the desired syst to plot error band
																				              bool plotRepSample=false,               // plot the 'Rep Sample' on canvas
																				              //todo: add data flag too...
																				              bool isDataFake=false)                  // using REAL data or Fake Data?
// =====================================================================================================
{

  std::string dataStatusString;
  isDataFake ? dataStatusString="Fake Data" : dataStatusString="Data";
  bool loadFakeDataFlag;
  isDataFake ? loadFakeDataFlag= true : loadFakeDataFlag= false;
	std::cout << "Plotting the NOvA CV with error band for systematic: " << systShortName << std::endl;
	const std::string plotName = Form("plot_NOvA_systematic_%s_errorband", systShortName.c_str());

	std::cout << "Also plotting the NOvA......" << dataStatusString << std::endl;
	if (plotRepSample) std::cout << "Also plotting the MCMC 'Representative Sample' (Max-Likelihood) with: " << dataStatusString << std::endl;
	std::cout << "Categorization: " << predType << std::endl;



  //load all systs that exist in the preds ROOT file
  getAllXsecSysts_2020();
  getNeutronSyst_2020();
  getMECtuneSysts();
  getAllFluxSysts_2020();

  mcmc::FileContents fileContents = mcmc::LoadFromFile(samplesFilename, isDataFake);
  std::unique_ptr<ana::MCMCSamples> warmup = std::move(fileContents.warmup);
  std::unique_ptr<ana::MCMCSamples> samples = std::move(fileContents.samples);

  std::string dataDir;
  isDataFake ? dataDir="spectra" : dataDir="data";
  std::map<std::string, ana::Spectrum> dataSpectra = ndfit::LoadData(dataFilename, dataDir.c_str());

  if (isDataFake)
  {
    fileContents = mcmc::LoadFromFile(dataFilename, isDataFake);
    std::unique_ptr<ana::SystShifts> systTruePulls = std::move(fileContents.trueSystPulls);
    std::unique_ptr<osc::IOscCalcAdjustable> calcTruth = std::move(fileContents.trueOscParams);
  }


  TString systDir = static_cast<TString>(systShortName);
  TString outDirPath = outDirPrefix+"/"+systDir;
	gSystem->MakeDirectory(outDirPath);
	const std::string& outDirPathString = static_cast<std::string>(outDirPath);
	std::string repSampleName;
	if (plotRepSample) repSampleName = "_with_rep_sample";


	// load the ND preds
	std::vector<ana::FitPredictions> preds = ndfit::LoadNDTopologicalPreds(predType, false);



	/// Create calc and syst objects
	auto bestFitIdx = samples->BestFitSampleIdx(); // sample with the largest LogLikelihood
	auto calc = std::make_unique<osc::OscCalcAnalytic>();
	ndfit::Calculator2020BestFit(*calc);
	for (const auto & var : samples->Vars())
		var->SetValue(calc.get(), samples->SampleValue(var, bestFitIdx));

	/// Initialize SystShifts objects for all of our needs...
	auto bestShifts = std::make_unique<ana::SystShifts>(); // the Min LL pull
	auto cvShifts = std::make_unique<ana::SystShifts>(); // cv -- 0 sigma
	auto  plusOneSigmaShifts = std::make_unique<ana::SystShifts>(); // +1sigma (for error bands)
	auto minusOneSigmaShifts = std::make_unique<ana::SystShifts>(); // -1sigma (for error bands)
	auto desiredSystOnePlus  = std::make_unique<ana::SystShifts>(); // for the systShortName
	auto desiredSystOneMinus = std::make_unique<ana::SystShifts>(); // for the systShortName
	std::cout << "systs in MCMC samples: " << samples->Systs().size() << std::endl;
	for (const auto & syst : samples->Systs()) {
		bestShifts->SetShift(syst, samples->SampleValue(syst, bestFitIdx));
		minusOneSigmaShifts->SetShift(syst, -1.0);
		plusOneSigmaShifts->SetShift(syst, 1.0);
		cvShifts->SetShift(syst, 0);
		if ( syst->ShortName() == systShortName )
        {
		  desiredSystOneMinus->SetShift(syst, -1.0);
		  desiredSystOnePlus->SetShift(syst, +1.0);
        }
//		else {
//			continue;
//		}
		//todo: not sure how to resolve this so it only kicks in after looping through every syst first.
//		else {
//			std::cerr << "Systematic ShortName() " << systShortName << " not found. Abort." << std::endl;
//			exit(1);
//		}
	}




  //load the expts to get the LL for each topology
  auto expts = mcmc::BuildExperiments(dataSpectra, preds);
  auto expt = mcmc::BuildMultiExperiment(mcmc::ExptPtrs(expts));

  // total ChiSq. Value printed at end of macro -- easier to read.
  double chiSq = expt->ChiSq(calc.get(),*bestShifts);

  //save the ChiSq for each topology into here
  std::unordered_map<std::string, double> chiSqMCMCMap {};
  std::unordered_map<std::string, double> chiSqNominalMap {};



  TCanvas c("c","c", 900,600); //("c", "c", 1000, 1000); // try 4x3 ?
  TPad * p1, * p2; //p1 upper, p2 lower

  /// Ehad axis stuff
	TGaxis * xAxisEHad = new TGaxis(0.001, 0.5, 0.8, 0.501, 0., 0.8, 10, "");
	xAxisEHad->SetLabelOffset(-0.015); // //	  std::cout << xAxisEHad->GetLabelOffset() --> 0.005
	xAxisEHad->SetLabelFont(42);
	TLatex * mcmcTxtEHad = new TLatex(0.82, 1.5, "#frac{MCMC Rep. Sample}{NOvA CV}");
	mcmcTxtEHad->SetTextAlign(11); mcmcTxtEHad->SetTextSize(0.02);
	mcmcTxtEHad->SetTextAngle(270); mcmcTxtEHad->SetTextColor(kAzure - 4);

	/// q3 reco axis stuff
	TGaxis * xAxisq3 = new TGaxis(0.001, 0.5, 2.0, 0.501, 0., 2.0, 10, "");
	xAxisq3->SetLabelOffset(-0.015); // //	  std::cout << xAxisEHad->GetLabelOffset() --> 0.005
	xAxisq3->SetLabelFont(42);
	TLatex * mcmcTxtq3 = new TLatex(2.05, 1.5, "#frac{MCMC Rep. Sample}{NOvA CV}");
	mcmcTxtq3->SetTextAlign(11); mcmcTxtq3->SetTextSize(0.02);
	mcmcTxtq3->SetTextAngle(270); mcmcTxtq3->SetTextColor(kAzure - 4);


  const double scaleFactor = 1./10000;
  int sampleType = 0;
  //I don't think this is needed...
  const std::vector<const ISyst *> systs = bestShifts->ActiveSysts();

  std::cout << "Plotting the ratio and fitted predictions with:  "<< dataStatusString << std::endl;
  std::cout << "Split into the categories: " << predType << std::endl;
  for (const auto &predBundle : preds) {
    std::vector<TH1*> q3PlusSigmaHistVector, q3MinusSigmaHistVector, EHadVisPlusSigmaVector, EHadVisMinusSigmaVector;

    // retrieve the ChiSq()s for each topological sample.
    double tmpChiSqMCMC = expts.at(sampleType)->ChiSq(calc.get(), *bestShifts);
    double tmpChiSqNom  = expts.at(sampleType)->ChiSq(calc.get(), *cvShifts);
    chiSqMCMCMap.try_emplace(predBundle.name, tmpChiSqMCMC);
    chiSqNominalMap.try_emplace(predBundle.name, tmpChiSqNom);

    c.cd();
    c.Clear();

    TPaveText ptTopology(0.7, 0.68, 0.85, 0.75, "ARC NDC");
    ptTopology.SetFillStyle(0);
    ptTopology.SetFillColor(0);
    ptTopology.SetBorderSize(0);

	  const std::string beamType = visuals::GetHornCurrent(predBundle.name);
	  const std::string topologyType = visuals::GetTopologyName(predBundle.name, predType);
	  ptTopology.AddText(Form("%s", beamType.c_str()));
	  ptTopology.AddText(Form("%s", topologyType.c_str()));


    //pavetext to print out the events for each topology
    TPaveText ptq3Events(0.7, 0.60, 0.85, 0.67, "ARC NDC"), ptEHadEvents(0.7, 0.60, 0.85, 0.67, "ARC NDC");
    ptq3Events.SetFillColor(0); ptq3Events.SetFillStyle(0); ptq3Events.SetBorderSize(0); ptq3Events.SetTextSize(0.032); ptq3Events.SetTextFont(102);
    ptEHadEvents.SetFillColor(0); ptEHadEvents.SetFillStyle(0); ptEHadEvents.SetBorderSize(0); ptEHadEvents.SetTextSize(0.032); ptEHadEvents.SetTextFont(102);

	  double weed = dataSpectra.at(predBundle.name).POT();
	  std::cout << "POT for predBundle: " << predBundle.name << ": " << weed << std::endl;

    // 2D profiled view the two variables: q3_Reco & EHadVis.
    PlotWithSystErrorBand((IPrediction *&) predBundle.pred, systs,calc.get(),weed,kGray+2, kGray);
    c.SaveAs(ndfit::FullFilename(outDirPathString, "profiled_error_bands_plot_" + predBundle.name + ".png").c_str());
    c.Clear();
    // 2D profile




    /// Plot comparison and ratio on save canvas
    SplitCanvas(0.25,p1,p2);
    // q3-Reco
    //create the histograms for the PlotWithSystErrorBand() function
    std::cout << "Producing q3-Reco plots for " << predBundle.name << "......" << std::endl;
    TH1 * hq3Data = dataSpectra.at(predBundle.name).ToTH2(weed)->ProjectionX();
    TH1 * hq3BestFit = predBundle.pred->PredictSyst(calc.get(), *bestShifts).ToTH2(weed)->ProjectionX();
//    TH1 * hDesiredSyst = predBundle.pred->PredictSyst(calc.get(), *d)

    TH1 * hq3CVPred = predBundle.pred->PredictSyst(calc.get(), *cvShifts).ToTH2(weed)->ProjectionX();
    TH1 * hq3PlusSigmaDesiredSyst = predBundle.pred->PredictSyst(calc.get(), *desiredSystOnePlus).ToTH2(weed)->ProjectionX();
    q3PlusSigmaHistVector.push_back(hq3PlusSigmaDesiredSyst);
    TH1 * hq3MinusSigmaDesiredSyst = predBundle.pred->PredictSyst(calc.get(), *desiredSystOneMinus).ToTH2(weed)->ProjectionX();
    q3MinusSigmaHistVector.push_back(hq3MinusSigmaDesiredSyst);

    //get the events BEFORE the re-scaling
    int dataq3Events = 0.0, mcmcq3Events = 0.0;
    dataq3Events = hq3Data->Integral();
    mcmcq3Events = hq3BestFit->Integral();
    ptq3Events.AddText("Data/MCMC:");
    ptq3Events.AddText(Form("%d/%d", dataq3Events, mcmcq3Events));
    std:: cout << "Events: " << dataq3Events << "/" << mcmcq3Events << std::endl;

    visuals::PredPreDrawAesthetics(hq3CVPred, scaleFactor, false);
    visuals::PredPreDrawAesthetics(hq3MinusSigmaDesiredSyst, scaleFactor, false);
    visuals::PredPreDrawAesthetics(hq3PlusSigmaDesiredSyst, scaleFactor, false);
    visuals::PredPreDrawAesthetics(hq3BestFit, scaleFactor, true);
    visuals::DataPreDrawAesthetics(hq3Data, scaleFactor);

    p1->cd();
    /// should plot a band of only the systematic I requested
    auto q3ErrorBand = PlotWithSystErrorBand(hq3CVPred, q3PlusSigmaHistVector, q3MinusSigmaHistVector, kGray + 2,kGray);
    hq3CVPred->Draw("same hist e");
    if (plotRepSample) hq3BestFit->Draw("same hist e");
    hq3Data->Draw("same hist p"); // draw as points (to distinguish with data)
    visuals::ComparisonPlotPostDrawAesthetics(hq3CVPred, visuals::maxRecoq3BinContent.at(predBundle.name));

    TLegend leg(0.6, 0.75, 0.9, 0.9);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    const std::string systLegName = "CV Prediction & 1#sigma: " + systShortName;
    leg.AddEntry(hq3CVPred, Form("%s", systLegName.c_str()), "l");
    leg.AddEntry(hq3Data, Form( "%s" , dataStatusString.c_str()), "p");
    if (plotRepSample) {
    	leg.AddEntry(hq3BestFit, "MCMC 'Rep. Sample'", "l");
	    ptq3Events.Draw("same");
    }
    leg.Draw("same");
    ptTopology.Draw("same");
    if (isDataFake) FakeData();
    else Preliminary();

    //q3-Reco ratio
    p2->cd();
    p2->SetGridy(1);
    TH1 * hq3Unity = (TH1F *) hq3CVPred->Clone("hq3Unity");
    hq3Unity->Divide(hq3CVPred);
    TH1 * hq3DataRatio = (TH1F *) hq3Data->Clone("hq3DataRatio");
    hq3DataRatio->Divide(hq3CVPred);
    TH1 * hq3MCMCRatio = (TH1F *) hq3BestFit->Clone("hq3MCMCRatio");
    hq3MCMCRatio->Divide(hq3CVPred);

    q3PlusSigmaHistVector[0]->Divide(hq3CVPred);
    q3MinusSigmaHistVector[0]->Divide(hq3CVPred);
    PlotWithSystErrorBand(hq3Unity, q3PlusSigmaHistVector, q3MinusSigmaHistVector, kGray + 2, kGray);

    hq3Unity->GetXaxis()->CenterTitle();
    hq3Unity->GetXaxis()->SetTitleOffset(1.);
    hq3Unity->GetXaxis()->SetTitleSize(0.045);
    hq3Unity->SetXTitle("#vec{q}_{3} Reco (GeV)");
    hq3Unity->GetYaxis()->CenterTitle();
    hq3Unity->GetYaxis()->SetRangeUser(0.5, 1.5);
    hq3Unity->GetYaxis()->SetTitleSize(0.02);
    hq3Unity->GetYaxis()->SetLabelSize(0.02);
    hq3Unity->GetYaxis()->SetTitleOffset(1.5);
    hq3Unity->SetYTitle(Form("#frac{%s}{NOvA CV}", dataStatusString.c_str()));
    hq3Unity->GetYaxis()->CenterTitle();
    xAxisq3->Draw("same");


    hq3DataRatio->SetMarkerColor(kBlack);
    hq3DataRatio->SetMarkerStyle(kFullCircle);
    hq3DataRatio->Draw("hist same pe");

    visuals::RatioPlotDrawAesthetics(hq3DataRatio,dataStatusString,"#vec{q}_{3} Reco (GeV)");
    hq3Unity->Draw("same hist pe");

	  if (plotRepSample) {
		  hq3MCMCRatio->Draw("same e");
		  mcmcTxtq3->Draw();
	  }
	  c.SaveAs(ndfit::FullFilename(outDirPathString, plotName + "_q3Reco_" + predBundle.name + repSampleName + ".png").c_str());
//	  for ( auto & ext : {".png", ".pdf", ".root"} )
//      c.SaveAs(ndfit::FullFilename(bestfitDir, plotName + "_q3Reco_" + predBundle.name + ext).c_str());




    c.Clear();
    c.cd();
    SplitCanvas(0.25,p1,p2);




      // EHadVis comparison
    std::cout << "Producing EHadVis plots for " << predBundle.name << "......" << std::endl;
    TH1 * hEHadVisData = dataSpectra.at(predBundle.name).ToTH2(weed)->ProjectionY();
    TH1 * hEHadVisBestFit = predBundle.pred->PredictSyst(calc.get(), *bestShifts).ToTH2(
            weed)->ProjectionY();


    TH1 * hEHadVisCVPred = predBundle.pred->PredictSyst(calc.get(), *cvShifts).ToTH2(weed)->ProjectionY();
    TH1 * hEHadVisPlusSigmaDesiredSyst = predBundle.pred->PredictSyst(calc.get(), *desiredSystOnePlus).ToTH2(weed)->ProjectionY();
    EHadVisPlusSigmaVector.push_back(hEHadVisPlusSigmaDesiredSyst);
    TH1 * hEHadVisMinusSigmaDesiredSyst = predBundle.pred->PredictSyst(calc.get(), *desiredSystOneMinus).ToTH2(weed)->ProjectionY();
    EHadVisMinusSigmaVector.push_back(hEHadVisMinusSigmaDesiredSyst);

    //get the events BEFORE the re-scaling
    int dataEHadEvents = 0.0, mcmcEHadEvents = 0.0;
    dataEHadEvents = hEHadVisData->Integral();
    mcmcEHadEvents = hEHadVisBestFit->Integral();
    ptEHadEvents.AddText("Data/MCMC Events:");
    ptEHadEvents.AddText(Form("%d/%d", dataEHadEvents, mcmcEHadEvents));
    std:: cout << "Events: " << dataEHadEvents << "/" << mcmcEHadEvents << std::endl;

    visuals::DataPreDrawAesthetics(hEHadVisData, scaleFactor);
    visuals::PredPreDrawAesthetics(hEHadVisBestFit, scaleFactor, true);
    visuals::PredPreDrawAesthetics(hEHadVisCVPred, scaleFactor, false);
    visuals::PredPreDrawAesthetics(hEHadVisPlusSigmaDesiredSyst, scaleFactor, false);
    visuals::PredPreDrawAesthetics(hEHadVisMinusSigmaDesiredSyst, scaleFactor, false);

    p1->cd();
    auto eHadErrorBand = PlotWithSystErrorBand(hEHadVisCVPred, EHadVisPlusSigmaVector, EHadVisMinusSigmaVector, kGray + 2,kGray);
    hEHadVisCVPred->Draw("same hist e");
    if (plotRepSample) {
    	hEHadVisBestFit->Draw("same hist e");
	    ptEHadEvents.Draw("same");
    }
    hEHadVisData->Draw("same hist p"); // draw as points (to distinguish with data)
    visuals::ComparisonPlotPostDrawAesthetics(hEHadVisCVPred, visuals::maxEHadBinContent.at(predBundle.name));

    leg.Draw("same");
    ptTopology.Draw("same");
    if (isDataFake) FakeData();
    else Preliminary();

    //EHadVis ratio
    p2->cd();
    p2->SetGridy(1);
    TH1 * hEHadUnity = (TH1F *) hEHadVisCVPred->Clone("hEHadUnity");
    hEHadUnity->Divide(hEHadVisCVPred);
    TH1 * hEHadVisDataRatio = (TH1F *) hEHadVisData->Clone("hEHadVisDataRatio");
	  hEHadVisDataRatio->Divide(hEHadVisCVPred);
    TH1 * hEHadMCMCRatio = (TH1F *) hEHadVisBestFit->Clone("hEHadMCMCRatio");
	  hEHadMCMCRatio->Divide(hEHadVisCVPred);

    EHadVisPlusSigmaVector[0]->Divide(hEHadVisCVPred);
    EHadVisMinusSigmaVector[0]->Divide(hEHadVisCVPred);
    PlotWithSystErrorBand(hEHadUnity, EHadVisPlusSigmaVector, EHadVisMinusSigmaVector, kGray + 2, kGray);

    hEHadUnity->GetXaxis()->CenterTitle();
    hEHadUnity->GetXaxis()->SetTitleOffset(1.);
    hEHadUnity->GetXaxis()->SetTitleSize(0.045);
    hEHadUnity->SetXTitle("E_{had}^{vis} (GeV)");
    hEHadUnity->GetYaxis()->CenterTitle();
    hEHadUnity->GetYaxis()->SetRangeUser(0.5, 1.5);
    hEHadUnity->GetYaxis()->SetTitleSize(0.02);
    hEHadUnity->GetYaxis()->SetLabelSize(0.02);
    hEHadUnity->GetYaxis()->SetTitleOffset(1.5);
    hEHadUnity->SetYTitle(Form("#frac{%s}{NOvA CV}", dataStatusString.c_str()));
    hEHadUnity->GetYaxis()->CenterTitle();
    xAxisEHad->Draw("same");


	  hEHadVisDataRatio->SetMarkerColor(kBlack);
	  hEHadVisDataRatio->SetMarkerStyle(kFullCircle);
	  hEHadVisDataRatio->Draw("hist same pe");

      visuals::RatioPlotDrawAesthetics(hEHadVisDataRatio, dataStatusString, "E_{had}^{vis} (GeV)");
      hEHadUnity->Draw("same hist");


	  if (plotRepSample) {
		  hEHadMCMCRatio->Draw("same e");
		  mcmcTxtEHad->Draw();
	  }


	  c.SaveAs(ndfit::FullFilename(outDirPathString, plotName + "_EHadVis_" + predBundle.name + repSampleName + ".png").c_str());
//	  for ( auto & ext : {".png", ".pdf", ".root"} )
//      c.SaveAs(ndfit::FullFilename(bestfitDir, plotName + "_EHadVis_" + predBundle.name + ext).c_str());

    sampleType++;
  } //predBundle in preds

//  std::cout << "X^2 is: " << chiSq << std::endl;
//  std::cout << "========== Nominal X^2 ===========" << std::endl;
//  for (auto const & pair : chiSqNominalMap) std::cout << pair.first << " --> " << pair.second << std::endl;
//  std::cout << "========== MCMC X^2 ===========" << std::endl;
//  for (auto const & pair : chiSqMCMCMap) std::cout << pair.first << " --> " << pair.second << std::endl;

}
