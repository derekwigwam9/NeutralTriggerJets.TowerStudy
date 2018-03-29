// 'CompareFourThringsWithRatio.C'
// Derek Anderson
// 07.25.2017
//
// Use this to compare four histograms and make some ratios.


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
static const TString sOutput("dataXembedding.eTwrRec.withVsWithoutConsistentHotTowers.d9m10y2017.root");
static const TString sInputA("dataXembedding.eTwrRec.smallerRangeBiggerBins.d25m9y2017.root");
static const TString sInputB("dataXembedding.eTwrRec.smallerRangeBiggerBins.d25m9y2017.root");
static const TString sInputC("dataXembedding.eTwrRec.pi0xHighTwr.d27m9y2017.root");
static const TString sInputD("dataXembedding.eTwrRec.pi0xHighTwr.d27m9y2017.root");
static const TString sHistA("hDataPerTrigger");
static const TString sHistB("hSumPerTrigger");
static const TString sHistC("hDataPerTrigger");
static const TString sHistD("hSumPerTrigger");
// histogram parameters
static const TString sNameA("hTwrEneWo_data");
static const TString sNameB("hTwrEneWo_embed");
static const TString sNameC("hTwrEneW_data");
static const TString sNameD("hTwrEneW_embed");
static const TString sTitle("Tower energy (raw)");
static const TString sTitleX("E^{twr}_{raw}");
static const TString sTitleY("(1/N^{trg}_{eff}) dN^{twr}/dE^{twr}_{raw}");
static const TString sTitleYR("embedding / data");
// legend parameters
static const TString sLegendA("data (inconsistent)");
static const TString sLegendB("embedding (inconsistent)");
static const TString sLegendC("data (consistent)");
static const TString sLegendD("embedding (consistent)");
static const TString sLegendAB("inconsistent");
static const TString sLegendCD("consistent");
// label parameters
static const TString sLabel1("pp collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("E_{T}^{trg} #in (9, 20) GeV, |#eta^{trg}| < 0.9");
static const TString sLabel3("p_{T}^{trk} > 0.2 GeV/c, |#eta^{trk}| < 1");
static const TString sLabel4("Identified #pi^{0} vs. high tower");
// canvas parameters
static const TString sCanvas("cTwrEne");
static const TString sCanvasR("cTwrEneRatio");



void CompareFourThingsWithRatio() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning comparison script..." << endl;


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInputA = new TFile(sInputA.Data(), "read");
  TFile *fInputB = new TFile(sInputB.Data(), "read");
  TFile *fInputC = new TFile(sInputC.Data(), "read");
  TFile *fInputD = new TFile(sInputD.Data(), "read");
  if (!fOutput || !fInputA || !fInputB || !fInputC || !fInputD) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fInputA);
    assert(fInputB);
    assert(fInputC);
    assert(fInputD);
  }
  cout << "    Files opened." << endl;

  TH1D *hInputA = (TH1D*) fInputA -> Get(sHistA.Data());
  TH1D *hInputB = (TH1D*) fInputB -> Get(sHistB.Data());
  TH1D *hInputC = (TH1D*) fInputC -> Get(sHistC.Data());
  TH1D *hInputD = (TH1D*) fInputD -> Get(sHistD.Data());
  if (!hInputA || !hInputB || !hInputC || !hInputD) {
    cerr << "PANIC: couldn't grab histogram!" << endl;
    assert(hInputA);
    assert(hInputB);
    assert(hInputC);
    assert(hInputD);
  }
  cout << "    Histograms grabbed." << endl;


  const Int_t    cA     = 810;
  const Int_t    cB     = 800;
  const Int_t    cC     = 890;
  const Int_t    cD     = 880;
  const Int_t    mA     = 7;
  const Int_t    mB     = 4;
  const Int_t    mC     = 1;
  const Int_t    mD     = 25;
  const Int_t    txt    = 42;
  const Int_t    cnt    = 1;
  const Double_t lab    = 0.02;
  const Double_t rng[2] = {0., 30.};
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
  hInputA -> GetXaxis() -> SetRangeUser(rng[0], rng[1]);
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
  hInputB -> GetXaxis() -> SetRangeUser(rng[0], rng[1]);
  hInputB -> GetYaxis() -> SetLabelSize(lab);
  hInputB -> GetYaxis() -> CenterTitle(cnt);
  hInputB -> GetYaxis() -> SetTitleFont(txt);
  hInputB -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputC -> SetLineColor(cC);
  hInputC -> SetMarkerColor(cC);
  hInputC -> SetMarkerStyle(mC);
  hInputC -> SetTitleFont(txt);
  hInputC -> SetTitle(sTitle.Data());
  hInputC -> SetName(sNameC.Data());
  hInputC -> GetXaxis() -> SetLabelSize(lab);
  hInputC -> GetXaxis() -> CenterTitle(cnt);
  hInputC -> GetXaxis() -> SetTitleFont(txt);
  hInputC -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputC -> GetXaxis() -> SetRangeUser(rng[0], rng[1]);
  hInputC -> GetYaxis() -> SetLabelSize(lab);
  hInputC -> GetYaxis() -> CenterTitle(cnt);
  hInputC -> GetYaxis() -> SetTitleFont(txt);
  hInputC -> GetYaxis() -> SetTitle(sTitleY.Data());
  hInputD -> SetLineColor(cD);
  hInputD -> SetMarkerColor(cD);
  hInputD -> SetMarkerStyle(mD);
  hInputD -> SetTitleFont(txt);
  hInputD -> SetTitle(sTitle.Data());
  hInputD -> SetName(sNameD.Data());
  hInputD -> GetXaxis() -> SetLabelSize(lab);
  hInputD -> GetXaxis() -> CenterTitle(cnt);
  hInputD -> GetXaxis() -> SetTitleFont(txt);
  hInputD -> GetXaxis() -> SetTitle(sTitleX.Data());
  hInputD -> GetXaxis() -> SetRangeUser(rng[0], rng[1]);
  hInputD -> GetYaxis() -> SetLabelSize(lab);
  hInputD -> GetYaxis() -> CenterTitle(cnt);
  hInputD -> GetYaxis() -> SetTitleFont(txt);
  hInputD -> GetYaxis() -> SetTitle(sTitleY.Data());
  cout << "    Styles set." << endl;


  // rebin input
  const UInt_t nRebinA = 5;
  const UInt_t nRebinB = 5;
  const UInt_t nRebinC = 5;
  const UInt_t nRebinD = 5;
  const Bool_t rebinA  = false;
  const Bool_t rebinB  = false;
  const Bool_t rebinC  = false;
  const Bool_t rebinD  = false;
  if (rebinA) {
    hInputA -> Rebin(nRebinA);
    hInputA -> Scale(1. / (Double_t) nRebinA);
    cout << "    Rebinned input A." << endl;
  }
  if (rebinB) {
    hInputB -> Rebin(nRebinB);
    hInputB -> Scale(1. / (Double_t) nRebinB);
    cout << "    Rebinned input B." << endl;
  }
  if (rebinC) {
    hInputC -> Rebin(nRebinC);
    hInputC -> Scale(1. / (Double_t) nRebinC);
    cout << "    Rebinned input C." << endl;
  }
  if (rebinD) {
    hInputD -> Rebin(nRebinD);
    hInputD -> Scale(1. / (Double_t) nRebinD);
    cout << "    Rebinned input D." << endl;
  }

  // scale input
  const Bool_t scaleByEntries = false;
  if (scaleByEntries) {
    const UInt_t nA = hInputA -> GetEntries();
    const UInt_t nB = hInputB -> GetEntries();
    const UInt_t nC = hInputC -> GetEntries();
    const UInt_t nD = hInputD -> GetEntries();
    hInputA -> Scale(1. / (Double_t) nA);
    hInputB -> Scale(1. / (Double_t) nB);
    hInputC -> Scale(1. / (Double_t) nC);
    hInputD -> Scale(1. / (Double_t) nD);
    cout << "    Scaled input." << endl;
  }

  // scale by bin-width
  const Bool_t scaleBinA = false;
  const Bool_t scaleBinB = false;
  const Bool_t scaleBinC = false;
  const Bool_t scaleBinD = false;
  if (scaleBinA) {
    const Float_t binA = hInputA -> GetBinWidth(17);
    hInputA -> Scale(1. / binA);
  }
  if (scaleBinB) {
    const Float_t binB = hInputB -> GetBinWidth(17);
    hInputB -> Scale(1. / binB);
  }
  if (scaleBinC) {
    const Float_t binC = hInputC -> GetBinWidth(17);
    hInputC -> Scale(1. / binC);
  }
  if (scaleBinD) {
    const Float_t binD = hInputD -> GetBinWidth(17);
    hInputD -> Scale(1. / binD);
  }


  // calculate ratios
  const UInt_t nAB   = hInputA -> GetNbinsX();
  const UInt_t nCD   = hInputC -> GetNbinsX();
  const Double_t ab1 = hInputA -> GetBinLowEdge(1);
  const Double_t cd1 = hInputC -> GetBinLowEdge(1);
  const Double_t ab2 = hInputA -> GetBinLowEdge(nAB + 1);
  const Double_t cd2 = hInputC -> GetBinLowEdge(nCD + 1);
  TH1D *hRatioAB = new TH1D("hRatioAB", "", nAB, ab1, ab2);
  TH1D *hRatioCD = new TH1D("hRatioCD", "", nCD, cd1, cd2);
  hRatioAB -> Sumw2();
  hRatioCD -> Sumw2();
  hRatioAB -> Divide(hInputB, hInputA, 1., 1.);
  hRatioCD -> Divide(hInputD, hInputC, 1., 1.);

  const Int_t    cAB = 810;
  const Int_t    cCD = 890;
  const Int_t    mAB = 4;
  const Int_t    mCD = 25;
  const Int_t    txt = 42;
  const Int_t    cnt = 1;
  const Double_t lab = 0.02;
  hRatioAB -> SetLineColor(cAB);
  hRatioAB -> SetMarkerColor(cAB);
  hRatioAB -> SetMarkerStyle(mAB);
  hRatioAB -> SetTitleFont(txt);
  hRatioAB -> SetTitle(sTitle.Data());
  hRatioAB -> SetName(sNameA.Data());
  hRatioAB -> GetXaxis() -> SetLabelSize(lab);
  hRatioAB -> GetXaxis() -> CenterTitle(cnt);
  hRatioAB -> GetXaxis() -> SetTitleFont(txt);
  hRatioAB -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatioAB -> GetXaxis() -> SetRangeUser(rng[0], rng[1]);
  hRatioAB -> GetYaxis() -> SetLabelSize(lab);
  hRatioAB -> GetYaxis() -> CenterTitle(cnt);
  hRatioAB -> GetYaxis() -> SetTitleFont(txt);
  hRatioAB -> GetYaxis() -> SetTitle(sTitleYR.Data());
  hRatioCD -> SetLineColor(cCD);
  hRatioCD -> SetMarkerColor(cCD);
  hRatioCD -> SetMarkerStyle(mCD);
  hRatioCD -> SetTitleFont(txt);
  hRatioCD -> SetTitle(sTitle.Data());
  hRatioCD -> SetName(sNameA.Data());
  hRatioCD -> GetXaxis() -> SetLabelSize(lab);
  hRatioCD -> GetXaxis() -> CenterTitle(cnt);
  hRatioCD -> GetXaxis() -> SetTitleFont(txt);
  hRatioCD -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatioCD -> GetXaxis() -> SetRangeUser(rng[0], rng[1]);
  hRatioCD -> GetYaxis() -> SetLabelSize(lab);
  hRatioCD -> GetYaxis() -> CenterTitle(cnt);
  hRatioCD -> GetYaxis() -> SetTitleFont(txt);
  hRatioCD -> GetYaxis() -> SetTitle(sTitleYR.Data());
  cout << "    Ratios calculated." << endl;


  // make line
  const Int_t   cI = 1;
  const Int_t   sI = 2;
  TLine *lOne = new TLine(rng[0], 1., rng[1], 1.);
  lOne -> SetLineColor(cI);
  lOne -> SetLineStyle(sI);
  cout << "     Line made." << endl;


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
  lLegend -> AddEntry(hInputC, sLegendC.Data());
  lLegend -> AddEntry(hInputD, sLegendD.Data());
  cout << "    Legend created." << endl;

  TLegend *lLegendR = new TLegend(x1L, y1L, x2L, y2L);
  lLegendR -> SetFillColor(cL);
  lLegendR -> SetFillStyle(sL);
  lLegendR -> SetLineColor(cL);
  lLegendR -> SetLineStyle(sL);
  lLegendR -> SetTextFont(txt);
  lLegendR -> AddEntry(hRatioAB, sLegendAB.Data());
  lLegendR -> AddEntry(hRatioCD, sLegendCD.Data());
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
  pLabel -> AddText(sLabel4.Data());
  cout << "    Label created." << endl;


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
  hInputC -> Draw("same");
  hInputD -> Draw("same");
  lLegend -> Draw();
  pLabel  -> Draw();
  cPlot   -> Write();
  cPlot   -> Close();

  const UInt_t  mar   = 0;
  const Float_t x[2]  = {0., 1.};
  const Float_t yR[2] = {0., 0.25};
  const Float_t yP[2] = {0.25, 1.};
  TCanvas *cRatio = new TCanvas(sCanvasR.Data(), "", wC, hC);
  TPad    *pRatio = new TPad("pRatio", "", x[0], yR[0], x[1], yR[1]);
  TPad    *pPlot  = new TPad("pPlot", "", x[0], yP[0], x[1], yP[1]);
  pRatio   -> SetGrid(grd, grd);
  pRatio   -> SetTopMargin(mar);
  pPlot    -> SetGrid(grd, grd);
  pPlot    -> SetLogy(log);
  pPlot    -> SetBottomMargin(mar);
  pRatio   -> Draw();
  pPlot    -> Draw();
  pRatio   -> cd();
  hRatioAB -> Draw();
  hRatioCD -> Draw("same");
  lOne     -> Draw();
  lLegendR -> Draw();
  pPlot    -> cd();
  hInputA  -> Draw();
  hInputB  -> Draw("same");
  hInputC  -> Draw("same");
  hInputD  -> Draw("same");
  lLegend  -> Draw();
  pLabel   -> Draw();
  cRatio   -> Write();
  cRatio   -> Close();
  cout << "    Plots drawn." << endl;


  fOutput -> cd();
  hInputA -> Write();
  hInputB -> Write();
  hInputC -> Write();
  hInputD -> Write();
  fOutput -> Close();
  fInputA -> cd();
  fInputA -> Close();
  fInputB -> cd();
  fInputB -> Close();
  fInputC -> cd();
  fInputC -> Close();
  fInputD -> cd();
  fInputD -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
