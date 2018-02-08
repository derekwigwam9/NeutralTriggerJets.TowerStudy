// 'SimplePlot.C'
// Derek Anderson
// 01.17.2018
//
// Use this to plot a summed histogram from
// embedding (the output of 'SimpleSum.C')
// against the equivalent histogram from data.
//
// NOTE: 'E' for embedding, 'D' for data


#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;


// constants
static const Bool_t   doNorm(true);
static const Bool_t   doIntNorm(false);
static const Bool_t   doRebin(true);
static const UInt_t   nRebin(3);
static const Double_t normer(20700.);
static const Double_t range[2] = {-3.15, 3.15};



void SimplePlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning (simple) plot script..." << endl;

  // io parameters
  const TString sOut("fTrgPhys.dataXembed.et9vz55hadXpi0.d8m2y2018.root");
  const TString sInD("pp200r9.fTrgCheck.et9vz55pi0.root");
  const TString sInE("pp200r9embed.fTrgPhys.et9vz55had.d7m2y2018.root");
  const TString sHistD("hClustPhi");
  const TString sHistE("hNorm");

  // histogram parameters
  const TString sTitle("Trigger (physics) #varphi");
  const TString sTitleX("#varphi^{trg}");
  const TString sTitleY("(1/N^{trg}) dN^{trg}/d#varphi^{trg}");
  const TString sTitleR("embedding / data");
  const TString sCanvas("cTrgPhi");


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fInD = new TFile(sInD.Data(), "read");
  TFile *fInE = new TFile(sInE.Data(), "read");
  if (!fInD || !fInE) {
    cerr << "PANIC: couldn't open an input file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hData  = (TH1D*) fInD -> Get(sHistD.Data()) -> Clone();
  TH1D *hEmbed = (TH1D*) fInE -> Get(sHistE.Data()) -> Clone();
  fOut -> cd();
  cout << "    Grabbed histograms." << endl;


  // rebin histograms (if needed)
  if (doRebin) {
    hData  -> Rebin(nRebin);
    hEmbed -> Rebin(nRebin);
    cout << "    Histograms rebinned." << endl;
  }

  // normalize data (if needed)
  if (doNorm) {
    if (doIntNorm) {
      const Double_t integral = hData -> Integral();
      hData  -> Scale(1. / integral);
      cout << "    Normalized data.\n"
           << "      normalization = " << integral
           << endl;
    }
    else {
      hData  -> Scale(1. / normer);
      cout << "    Normalized data.\n"
           << "      normalization - " << normer
           << endl;
    }
  }

  // calculate ratio
  const UInt_t   nBins = hData -> GetNbinsX();
  const Double_t xBin1 = hData -> GetBinLowEdge(1);
  const Double_t xBin2 = hData -> GetBinLowEdge(nBins + 1);

  TH1D *hRatio = new TH1D("hRatio", "", nBins, xBin1, xBin2);
  hRatio -> Sumw2();
  hRatio -> Divide(hEmbed, hData , 1., 1.);
  cout << "    Divided histograms." << endl;


  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.03);
  const Float_t fLabR(0.055);
  const Float_t fSizR(0.055);
  const Float_t fOffR(0.75);
  const UInt_t  fCol[3] = {899, 859, 819};
  const UInt_t  fMar[3] = {8, 29, 8};

  // data
  hData -> SetLineColor(fCol[0]);
  hData -> SetMarkerColor(fCol[0]);
  hData -> SetMarkerStyle(fMar[0]);
  hData -> SetTitle(sTitle.Data());
  hData -> SetTitleFont(fTxt);
  hData -> GetXaxis() -> SetTitle(sTitleX.Data());
  hData -> GetXaxis() -> SetTitleFont(fTxt);
  hData -> GetXaxis() -> CenterTitle(fCnt);
  hData -> GetXaxis() -> SetRangeUser(range[0], range[1]);
  hData -> GetXaxis() -> SetLabelSize(fLab);
  hData -> GetYaxis() -> SetTitle(sTitleY.Data());
  hData -> GetYaxis() -> SetTitleFont(fTxt);
  hData -> GetYaxis() -> CenterTitle(fCnt);
  hData -> GetYaxis() -> SetLabelSize(fLab);

  // embedding
  hEmbed -> SetLineColor(fCol[1]);
  hEmbed -> SetMarkerColor(fCol[1]);
  hEmbed -> SetMarkerStyle(fMar[1]);
  hEmbed -> SetTitle(sTitle.Data());
  hEmbed -> SetTitleFont(fTxt);
  hEmbed -> GetXaxis() -> SetTitle(sTitleX.Data());
  hEmbed -> GetXaxis() -> SetTitleFont(fTxt);
  hEmbed -> GetXaxis() -> CenterTitle(fCnt);
  hEmbed -> GetXaxis() -> SetRangeUser(range[0], range[1]);
  hEmbed -> GetXaxis() -> SetLabelSize(fLab);
  hEmbed -> GetYaxis() -> SetTitle(sTitleY.Data());
  hEmbed -> GetYaxis() -> SetTitleFont(fTxt);
  hEmbed -> GetYaxis() -> CenterTitle(fCnt);
  hEmbed -> GetYaxis() -> SetLabelSize(fLab);

  // ratio
  hRatio -> SetLineColor(fCol[2]);
  hRatio -> SetMarkerColor(fCol[2]);
  hRatio -> SetMarkerStyle(fMar[2]);
  hRatio -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatio -> GetXaxis() -> SetTitleFont(fTxt);
  hRatio -> GetXaxis() -> CenterTitle(fCnt);
  hRatio -> GetXaxis() -> SetRangeUser(range[0], range[1]);
  hRatio -> GetXaxis() -> SetTitleSize(fSizR);
  hRatio -> GetXaxis() -> SetTitleOffset(fOffR);
  hRatio -> GetXaxis() -> SetLabelSize(fLabR);
  hRatio -> GetYaxis() -> SetTitle(sTitleR.Data());
  hRatio -> GetYaxis() -> SetTitleFont(fTxt);
  hRatio -> GetYaxis() -> CenterTitle(fCnt);
  hRatio -> GetYaxis() -> SetLabelSize(fLabR);
  hRatio -> GetYaxis() -> SetTitleSize(fSizR);
  hRatio -> GetYaxis() -> SetTitleOffset(fOffR);
  cout << "    Set styles." << endl;


  // make legend
  const UInt_t  fColL(0);
  const UInt_t  fStyL(0);
  const Float_t xyLeg[4] = {0.1, 0.1, 0.3, 0.3};

  TLegend *legend = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legend -> SetFillColor(fColL);
  legend -> SetFillStyle(fStyL);
  legend -> SetLineColor(fColL);
  legend -> SetLineStyle(fStyL);
  legend -> SetTextFont(fTxt);
  legend -> AddEntry(hData, "data");
  legend -> AddEntry(hEmbed, "embedding");
  cout << "    Made legend." << endl;


  // make line
  const UInt_t fColI(1);
  const UInt_t fStyI(2);

  TLine *one = new TLine(range[0], 1., range[1], 1.);
  one -> SetLineColor(fColI);
  one -> SetLineStyle(fStyI);
  cout << "    Made line." << endl;


  // make plots
  const UInt_t  fLog(1);
  const UInt_t  fGrid(0);
  const UInt_t  fTicks(1);
  const UInt_t  width(750);
  const UInt_t  height(850);
  const Float_t margin(0.);
  const Float_t xyPadR[4] = {0., 0., 1., 0.35};
  const Float_t xyPadD[4] = {0., 0.35, 1., 1.};

  TCanvas *canvas = new TCanvas(sCanvas.Data(), "", width, height);
  TPad    *pRatio = new TPad("pRatio", "", xyPadR[0], xyPadR[1], xyPadR[2], xyPadR[3]);
  TPad    *pDistr = new TPad("pDistribution", "", xyPadD[0], xyPadD[1], xyPadD[2], xyPadD[3]);
  pRatio -> SetGrid(fGrid, fGrid);
  pRatio -> SetTicks(fTicks, fTicks);
  pRatio -> SetTopMargin(margin);
  pDistr -> SetLogy(fLog);
  pDistr -> SetGrid(fGrid, fGrid);
  pDistr -> SetTicks(fTicks, fTicks);
  pDistr -> SetBottomMargin(margin);
  canvas -> cd();
  pRatio -> Draw();
  pDistr -> Draw();
  pRatio -> cd();
  hRatio -> Draw();
  one    -> Draw();
  pDistr -> cd();
  hData  -> Draw();
  hEmbed -> Draw("same");
  legend -> Draw();
  canvas -> Write();
  canvas -> Close();
  cout << "    Drawn plots." << endl;


  // save and close
  fOut   -> cd();
  hData  -> Write();
  hEmbed -> Write();
  hRatio -> Write();
  fOut   -> Close();
  fInD   -> cd();
  fInD   -> Close();
  fInE   -> cd();
  fInE   -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
