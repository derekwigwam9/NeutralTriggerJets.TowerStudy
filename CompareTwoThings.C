// 'CompareTwoThings.C'
// Derek Anderson
// 06.28.2017
//
// Use this to compare a couple of histograms.


#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// i/o parameters
static const TString sOutput("dataXembedding.eTtrg.comparison.allFiles.d2m11y2017.root");
static const TString sInputA("scaleFactors.eTtrg.thirdJetMaker.allFiles.d2m11y2017.root");
static const TString sInputB("scaleFactors.eTtrg.thirdJetMaker.allFiles.d2m11y2017.root");
static const TString sHistA("hData");
static const TString sHistB("hSumNormalized");
// histogram parameters
static const Float_t pRange[2] = {4., 25.};
static const TString sNameA("hTrgEt_data");
static const TString sNameB("hTrgEt_embed");
static const TString sTitle("Trigger E_{T}");
static const TString sTitleX("E_{T}^{trg}");
static const TString sTitleY("dN^{trg}_{eff}/dE_{T}^{trg}");
static const TString sTitleYR("embedding / data");
// legend parameters
static const TString sLegendA("data");
static const TString sLegendB("embedding");
// label parameters
static const TString sLabel1("pp collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("E_{T}^{trg} #in (9,20) GeV, |#eta^{trg}| < 0.9");
static const TString sLabel3("p_{T}^{trk} > 0.2 GeV, |#eta^{trk}| < 1");
static const TString sLabel4("Anti-k_{T}, R = 0.3, full jets");
// canvas parameters
static const TString sCanvas("cTrgEt");
static const TString sCanvasR("cTrgEtRatio");



void CompareTwoThings() {

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


  const Int_t    cA  = 810;
  const Int_t    cB  = 890;
  const Int_t    mA  = 7;
  const Int_t    mB  = 4;
  const Int_t    txt = 42;
  const Int_t    cnt = 1;
  const Double_t lab = 0.02;
  hInputA -> SetLineColor(cA);
  hInputA -> SetMarkerColor(cA);
  hInputA -> SetMarkerStyle(mA);
  hInputA -> SetTitleFont(txt);
  hInputA -> SetTitle(sTitle.Data());
  hInputA -> SetName(sNameA.Data());
  hInputA -> GetXaxis() -> SetLabelSize(lab);
  hInputA -> GetXaxis() -> CenterTitle(cnt);
  hInputA -> GetXaxis() -> SetTitleFont(txt);
  hInputA -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputA -> GetYaxis() -> SetLabelSize(lab);
  hInputA -> GetYaxis() -> CenterTitle(cnt);
  hInputA -> GetYaxis() -> SetTitleFont(txt);
  hInputA -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputB -> SetLineColor(cB);
  hInputB -> SetMarkerColor(cB);
  hInputB -> SetMarkerStyle(mB);
  hInputB -> SetTitleFont(txt);
  hInputB -> SetTitle(sTitle.Data());
  hInputB -> SetName(sNameB.Data());
  hInputB -> GetXaxis() -> SetLabelSize(lab);
  hInputB -> GetXaxis() -> CenterTitle(cnt);
  hInputB -> GetXaxis() -> SetTitleFont(txt);
  hInputB -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputB -> GetYaxis() -> SetLabelSize(lab);
  hInputB -> GetYaxis() -> CenterTitle(cnt);
  hInputB -> GetYaxis() -> SetTitleFont(txt);
  hInputB -> GetYaxis() -> SetTitle(sTitleY.Data());
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


  // make ratio
  const UInt_t   nBins = hInputA -> GetNbinsX();
  const Double_t xBin1 = hInputA -> GetBinLowEdge(1);
  const Double_t xBin2 = hInputA -> GetBinLowEdge(nBins + 1);
  TH1D *hRatio = new TH1D("hRatio", "", nBins, xBin1, xBin2);
  hRatio -> Divide(hInputB, hInputA, 1., 1.);

  // set ratio's style
  const UInt_t  cR   = 850;
  const UInt_t  mR   = 4;
  const UInt_t  tR   = 42;
  const UInt_t  cntR = 1;
  const Float_t labR = 0.035;
  const Float_t ttlR = 0.065;
  const Float_t offR = 0.65;
  hRatio -> SetLineColor(cR);
  hRatio -> SetMarkerColor(cR);
  hRatio -> SetMarkerStyle(mR);
  hRatio -> GetXaxis() -> SetLabelSize(labR);
  hRatio -> GetXaxis() -> SetTitleSize(ttlR);
  hRatio -> GetXaxis() -> SetTitleOffset(offR);
  hRatio -> GetXaxis() -> SetTitleFont(tR);
  hRatio -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatio -> GetXaxis() -> CenterTitle(cntR);
  hRatio -> GetYaxis() -> SetLabelSize(labR);
  hRatio -> GetYaxis() -> SetTitleSize(ttlR);
  hRatio -> GetYaxis() -> SetTitleOffset(offR);
  hRatio -> GetYaxis() -> SetTitleFont(tR);
  hRatio -> GetYaxis() -> SetTitle(sTitleYR.Data());
  hRatio -> GetYaxis() -> CenterTitle(cntR);


  const Int_t    cL  = 0;
  const Int_t    fL  = 0;
  const Int_t    sL  = 0;
  const Double_t x1L = 0.1;
  const Double_t x2L = 0.3;
  const Double_t y1L = 0.1;
  const Double_t y2L = 0.3;
  TLegend *lLegend = new TLegend(x1L, y1L, x2L, y2L);
  lLegend -> SetFillColor(cL);
  lLegend -> SetFillStyle(sL);
  lLegend -> SetLineColor(cL);
  lLegend -> SetLineStyle(sL);
  lLegend -> SetTextFont(txt);
  lLegend -> AddEntry(hInputA, sLegendA.Data());
  lLegend -> AddEntry(hInputB, sLegendB.Data());
  cout << "    Legend created." << endl;

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


  // make line and set range
  const UInt_t cO = 1;
  const UInt_t sO = 2;
  TLine *lOne = new TLine(pRange[0], 1., pRange[1], 1.);
  lOne    -> SetLineColor(cO);
  lOne    -> SetLineStyle(sO);
  hInputA -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
  hInputB -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
  hRatio  -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);


  // make plots
  fOutput -> cd();

  const Int_t wC  = 800;
  const Int_t hC  = 800;
  const Int_t grd = 0;
  const Int_t log = 0;
  TCanvas *cPlot = new TCanvas(sCanvas.Data(), "", wC, hC);
  cPlot   -> SetGrid(grd, grd);
  cPlot   -> SetLogy(log);
  hInputA -> Draw();
  hInputB -> Draw("same");
  lLegend -> Draw();
  //pLabel  -> Draw();
  cPlot   -> Write();
  cPlot   -> Close();

  const Float_t border   = 0.;
  const Float_t xPad[2]  = {0., 1.};
  const Float_t yPadR[2] = {0., 0.3};
  const Float_t yPadD[2] = {0.3, 1.};
  TCanvas *cRatio = new TCanvas(sCanvasR.Data(), "", wC, hC);
  TPad    *pRatio = new TPad("pRatio", "", xPad[0], yPadR[0], xPad[1], yPadR[1]);
  TPad    *pDist  = new TPad("pDist", "", xPad[0], yPadD[0], xPad[1], yPadD[1]);
  pRatio  -> SetGrid(grd, grd);
  pRatio  -> SetLogy(log);
  pRatio  -> SetTopMargin(border);
  pDist   -> SetGrid(grd, grd);
  pDist   -> SetLogy(log);
  pDist   -> SetBottomMargin(border);
  pRatio  -> Draw();
  pDist   -> Draw();
  pRatio  -> cd();
  hRatio  -> Draw();
  lOne    -> Draw();
  pDist   -> cd();
  hInputA -> Draw();
  hInputB -> Draw("same");
  lLegend -> Draw();
  cRatio  -> Write();
  cRatio  -> Close();
  cout << "    Plot drawn." << endl;


  fOutput -> cd();
  hInputA -> Write();
  hInputB -> Write();
  hRatio  -> Write();
  fOutput -> Close();
  fInputB -> cd();
  fInputB -> Close();
  fInputA -> cd();
  fInputA -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
