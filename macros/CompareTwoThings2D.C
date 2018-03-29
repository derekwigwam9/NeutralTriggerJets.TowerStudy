// 'CompareTwoThings2D.C'
// Derek Anderson
// 09.27.2017
//
// Use this to compare a couple of histograms.


#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// i/o parameters
static const TString sOutput("hotTowerMap.mergedListsDataVsEmbedding.d24m10y2017.root");
static const TString sInputA("hotTwrMap.mergedListsData.d24m10y2017.root");
static const TString sInputB("hotTwrMap.mergedListsEmbedding.d24m10y2017.root");
static const TString sHistA("hTwrPhiVsEta");
static const TString sHistB("hTwrPhiVsEta");
// histogram parameters
static const Float_t pRangeX[2] = {-1., 1.};
static const Float_t pRangeY[2] = {-1. * TMath::Pi(), 1. * TMath::Pi()};
static const TString sNameA("hTwrPhiVsEta_data");
static const TString sNameB("hTwrPhiVsEta_embedding");
static const TString sTitleA("Combined hot tower lists, data");
static const TString sTitleB("Combined hot tower lists, embedding");
static const TString sTitleX("#eta^{twr}");
static const TString sTitleY("#varphi^{twr}");
// label parameters
static const TString sLabel1("pp collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("E_{T}^{trg} #in (9,20) GeV, |#eta^{trg}| < 0.9");
static const TString sLabel3("p_{T}^{trk} > 0.2 GeV, |#eta^{trk}| < 1");
static const TString sLabel4("Anti-k_{T}, R = 0.3, full jets");
// canvas parameters
static const TString sCanvas("cTwrPhiVsEta");



void CompareTwoThings2D() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning comparison script..." << endl;


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInputA = new TFile(sInputA.Data(), "read");
  TFile *fInputB = new TFile(sInputB.Data(), "read");
  if (!fOutput || !fInputA || !fInputB) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fInputA);
    assert(fInputB);
  }
  cout << "    Files opened." << endl;

  TH1D *hInputA = (TH1D*) fInputA -> Get(sHistA.Data());
  TH1D *hInputB = (TH1D*) fInputB -> Get(sHistB.Data());
  if (!hInputA || !hInputB) {
    cerr << "PANIC: couldn't grab histogram!" << endl;
    assert(hInputA);
    assert(hInputB);
  }
  cout << "    Histograms grabbed." << endl;


  const Int_t    txt = 42;
  const Int_t    cnt = 1;
  const Double_t lab = 0.02;
  hInputA -> SetTitleFont(txt);
  hInputA -> SetTitle(sTitleA.Data());
  hInputA -> SetName(sNameA.Data());
  hInputA -> GetXaxis() -> SetRangeUser(pRangeX[0], pRangeX[1]);
  hInputA -> GetXaxis() -> SetLabelSize(lab);
  hInputA -> GetXaxis() -> CenterTitle(cnt);
  hInputA -> GetXaxis() -> SetTitleFont(txt);
  hInputA -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputA -> GetYaxis() -> SetRangeUser(pRangeY[0], pRangeY[1]);
  hInputA -> GetYaxis() -> SetLabelSize(lab);
  hInputA -> GetYaxis() -> CenterTitle(cnt);
  hInputA -> GetYaxis() -> SetTitleFont(txt);
  hInputA -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputA -> GetZaxis() -> SetLabelSize(lab);
  hInputB -> SetTitleFont(txt);
  hInputB -> SetTitle(sTitleB.Data());
  hInputB -> SetName(sNameB.Data());
  hInputB -> GetXaxis() -> SetRangeUser(pRangeX[0], pRangeX[1]);
  hInputB -> GetXaxis() -> SetLabelSize(lab);
  hInputB -> GetXaxis() -> CenterTitle(cnt);
  hInputB -> GetXaxis() -> SetTitleFont(txt);
  hInputB -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputB -> GetYaxis() -> SetRangeUser(pRangeY[0], pRangeY[1]);
  hInputB -> GetYaxis() -> SetLabelSize(lab);
  hInputB -> GetYaxis() -> CenterTitle(cnt);
  hInputB -> GetYaxis() -> SetTitleFont(txt);
  hInputB -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputB -> GetZaxis() -> SetLabelSize(lab);
  cout << "    Styles set." << endl;


  // rebin histograms
  const Bool_t rebinA  = false;
  const Bool_t rebinB  = false;
  const UInt_t nRebinA = 5;
  const UInt_t nRebinB = 5;
  if (rebinA) {
    hInputA -> Rebin(nRebinA);
    cout << "    Input A rebinned." << endl;
  }
  if (rebinB) {
    hInputB -> Rebin(nRebinB);
    cout << "    Input B rebinned." << endl;
  }


  // scale histograms
  const Bool_t   scaleA = false;
  const Bool_t   scaleB = false;
  const Double_t aScale = 1429435.;
  const Double_t bScale = 1336827.;
  if (scaleA) {
    hInputA -> Scale(1. / aScale);
    cout << "    Input A scaled." << endl;
  }
  if (scaleB) {
    hInputB -> Scale(1. / bScale);
    cout << "    Input B scaled." << endl;
  }


  const UInt_t   cL  = 0;
  const UInt_t   sL  = 0;
  const Double_t x1P = 0.3;
  const Double_t x2P = 0.5;
  const Double_t y1P = 0.1;
  const Double_t y2P = 0.3;
  TPaveText *pLabel = new TPaveText(x1P, y1P, x2P, y2P, "NDC NB");
  pLabel -> SetFillColor(cL);
  pLabel -> SetFillStyle(sL);
  pLabel -> SetLineColor(cL);
  pLabel -> SetLineStyle(sL);
  pLabel -> SetTextFont(txt);
  pLabel -> AddText(sLabel1.Data());
  pLabel -> AddText(sLabel2.Data());
  pLabel -> AddText(sLabel3.Data());
  //pLabel -> AddText(sLabel4.Data());
  cout << "    Label created." << endl;


  // make plots
  fOutput -> cd();

  const Int_t   wC       = 1600;
  const Int_t   hC       = 800;
  const Int_t   grd      = 0;
  const Int_t   log      = 1;
  const Float_t xPadA[2] = {0., 0.5};
  const Float_t xPadB[2] = {0.5, 1.};
  const Float_t yPad[2]  = {0., 1.};
  TCanvas *cPlot  = new TCanvas(sCanvas.Data(), "", wC, hC);
  TPad    *pPlotA = new TPad("pPlotA", "", xPadA[0], yPad[0], xPadA[1], yPad[1]);
  TPad    *pPlotB = new TPad("pPlotB", "", xPadB[0], yPad[0], xPadB[1], yPad[1]);
  pPlotA  -> SetGrid(grd, grd);
  pPlotA  -> SetLogz(log);
  pPlotB  -> SetGrid(grd, grd);
  pPlotB  -> SetLogz(log);
  pPlotA  -> Draw();
  pPlotB  -> Draw();
  pPlotA  -> cd();
  hInputA -> Draw("colz");
  pPlotB  -> cd();
  hInputB -> Draw("colz");
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Plot drawn." << endl;


  fOutput -> cd();
  hInputA -> Write();
  hInputB -> Write();
  fOutput -> Close();
  fInputB -> cd();
  fInputB -> Close();
  fInputA -> cd();
  fInputA -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
