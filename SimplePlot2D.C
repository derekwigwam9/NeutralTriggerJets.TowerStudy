// 'SimplePlot2D.C'
// Derek Anderson
// 02.22.2018
//
// Use this to plot a summed histogram from
// embedding (the output of 'SimpleSum.C')
// against the equivalent histogram from data.
//
// NOTE: 'E' for embedding, 'D' for data


#include <iostream>
#include "TH2.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;


// constants
static const Bool_t   doNorm(true);
static const Bool_t   doIntNorm(false);
static const Bool_t   doRebinX(false);
static const Bool_t   doRebinY(true);
static const UInt_t   nRebinX(1);
static const UInt_t   nRebinY(3);
static const Double_t normer(20700.);
static const Double_t rangeX[2] = {-1., 1.};
static const Double_t rangeY[2] = {-3.15, 3.15};



void SimplePlot2D() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning (simple) plot script..." << endl;

  // io parameters
  const TString sOut("phiVsEtaTrg.pi0xHad.d22m2y2018.root");
  const TString sInD("pp200r9.fTrgCheck.et9vz55pi0.root");
  const TString sInE("pp200r9embed.phiVsEtaTrg.et9vz55had.d22m2y2018.root");
  const TString sHistD("hTrgPhiVsEtaTwr");
  const TString sHistE("hNorm");

  // histogram parameters
  const TString sTitleD("Trigger #varphi vs. #eta, data (#pi^{0})");
  const TString sTitleE("Trigger #varphi vs. #eta, embedding (h^{#pm})");
  const TString sTitleX("#eta^{trg}");
  const TString sTitleY("#varphi^{trg}");
  const TString sCanvas("cTrgPhiVsEta");


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
  TH2D *hData  = (TH2D*) fInD -> Get(sHistD.Data()) -> Clone();
  TH2D *hEmbed = (TH2D*) fInE -> Get(sHistE.Data()) -> Clone();
  fOut -> cd();
  cout << "    Grabbed histograms." << endl;


  // rebin histograms (if needed)
  if (doRebinX) {
    hData  -> RebinX(nRebinX);
    hEmbed -> RebinX(nRebinX);
    cout << "    X-axis rebinned." << endl;
  }
  if (doRebinY) {
    hData  -> RebinY(nRebinY);
    hEmbed -> RebinY(nRebinY);
    cout << "    Y-axis rebinned." << endl;
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


  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.025);
  const Float_t fLabR(0.055);
  const Float_t fSizR(0.055);
  const Float_t fOffR(0.75);

  // data
  hData -> SetTitle(sTitleD.Data());
  hData -> SetTitleFont(fTxt);
  hData -> GetXaxis() -> SetTitle(sTitleX.Data());
  hData -> GetXaxis() -> SetTitleFont(fTxt);
  hData -> GetXaxis() -> CenterTitle(fCnt);
  hData -> GetXaxis() -> SetRangeUser(rangeX[0], rangeX[1]);
  hData -> GetXaxis() -> SetLabelSize(fLab);
  hData -> GetYaxis() -> SetTitle(sTitleY.Data());
  hData -> GetYaxis() -> SetTitleFont(fTxt);
  hData -> GetYaxis() -> CenterTitle(fCnt);
  hData -> GetYaxis() -> SetRangeUser(rangeY[0], rangeY[1]);
  hData -> GetYaxis() -> SetLabelSize(fLab);
  hData -> GetZaxis() -> SetLabelSize(fLab);

  // embedding
  hEmbed -> SetTitle(sTitleE.Data());
  hEmbed -> SetTitleFont(fTxt);
  hEmbed -> GetXaxis() -> SetTitle(sTitleX.Data());
  hEmbed -> GetXaxis() -> SetTitleFont(fTxt);
  hEmbed -> GetXaxis() -> CenterTitle(fCnt);
  hEmbed -> GetXaxis() -> SetRangeUser(rangeX[0], rangeX[1]);
  hEmbed -> GetXaxis() -> SetLabelSize(fLab);
  hEmbed -> GetYaxis() -> SetTitle(sTitleY.Data());
  hEmbed -> GetYaxis() -> SetTitleFont(fTxt);
  hEmbed -> GetYaxis() -> CenterTitle(fCnt);
  hEmbed -> GetYaxis() -> SetRangeUser(rangeY[0], rangeY[1]);
  hEmbed -> GetYaxis() -> SetLabelSize(fLab);
  hEmbed -> GetZaxis() -> SetLabelSize(fLab);
  cout << "    Set styles." << endl;


  // make plots
  const UInt_t  fLog(1);
  const UInt_t  fGrid(0);
  const UInt_t  fTicks(1);
  const UInt_t  width(1500);
  const UInt_t  height(750);
  const Float_t margin(0.);
  const Float_t xyPadD[4] = {0., 0., 0.5, 1.};
  const Float_t xyPadE[4] = {0.5, 0., 1., 1.};

  TCanvas *canvas = new TCanvas(sCanvas.Data(), "", width, height);
  TPad    *pData  = new TPad("pData", "", xyPadD[0], xyPadD[1], xyPadD[2], xyPadD[3]);
  TPad    *pEmbed = new TPad("pEmbed", "", xyPadE[0], xyPadE[1], xyPadE[2], xyPadE[3]);
  pData  -> SetLogz(fLog);
  pData  -> SetGrid(fGrid, fGrid);
  pData  -> SetTicks(fTicks, fTicks);
  pEmbed -> SetLogz(fLog);
  pEmbed -> SetGrid(fGrid, fGrid);
  pEmbed -> SetTicks(fTicks, fTicks);
  canvas -> cd();
  pData  -> Draw();
  pEmbed -> Draw();
  pData  -> cd();
  hData  -> Draw("colz");
  pEmbed -> cd();
  hEmbed -> Draw("colz");
  canvas -> Write();
  canvas -> Close();
  cout << "    Drawn plots." << endl;


  // save and close
  fOut   -> cd();
  hData  -> Write();
  hEmbed -> Write();
  fOut   -> Close();
  fInD   -> cd();
  fInD   -> Close();
  fInE   -> cd();
  fInE   -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
