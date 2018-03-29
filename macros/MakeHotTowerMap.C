// 'MakeHotTowerMap.C'
// Derek Anderson
// 10.11.2017
//
// Use this to produce an eta-phi map of hot/bad towers.
// 'sInputList' should be the name of a text file with
// the hot/bad tower list, and 'sInputHist' should be
// the name of a root file with the histograms specified
// by 'sTwrEtaVsID' and 'sTwrPhiVsID'.

#include <vector>
#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"


// i/o parameters
static const TString sInputList("pp200r12.mineAndKoljasListCombine.d24m10y2017.list");
static const TString sInputHist("/star/scratch/DancingWithDerek/pp200r9.hotTowerMapWithDeDxForReal.root");
static const TString sOutput("hotTwrMap.mergedListsEmbedding.d24m10y2017.root");


void MakeHotTowerMap() {

  gErrorIgnoreLevel = kError;
  cout << "\n Beginning map script..." << endl;

  // constants
  const UInt_t   nTow(4802);
  const TString  sTwrEtaVsID("hTwrEtaVsId_e0");
  const TString  sTwrPhiVsID("hTwrPhiVsId_e0");
  const Double_t tow1(-1.);
  const Double_t tow2(4801.);


  // open files and grab histograms
  TFile *fOutput     = new TFile(sOutput.Data(), "recreate");
  TFile *fInputHist  = new TFile(sInputHist.Data(), "read");
  if (!fOutput || !fInputHist) {
    cerr << "PANIC: couldn't open file!\n"
         << "       fInput  -- " << fInput << "\n"
         << "       fOutput -- " << fOutput
         << endl;
    return;
  }
  cout << "    Files opened." << endl;

  TH2D *hTwrEtaVsID = (TH2D*) fInputHist -> Get(sTwrEtaVsID.Data());
  TH2D *hTwrPhiVsID = (TH2D*) fInputHist -> Get(sTwrPhiVsID.Data());
  cout << "    Histograms grabbed." << endl;


  // load hot towers
  ifstream hotTwrList(sInputList.Data());
  if (!hotTwrList) {
    cerr << "PANIC: couldn't open hot tower list!" << endl;
    return;
  }
  vector<UInt_t> hotTwrs;
  hotTwrs.clear();

  UInt_t iHotTwr(0);
  UInt_t hotTwrID(0);
  while (hotTwrList) {
    hotTwrList >> hotTwrID;
    hotTwrs.push_back(hotTwrID);
    iHotTwr++;
  }
  hotTwrList.close();
  cout << "    Hot towers loaded:\n"
       << "      " << iHotTwr << " hot towers."
       << endl;


  // create eta-phi map
  const UInt_t   nBinsH = hTwrEtaVsID -> GetNbinsY();
  const UInt_t   nBinsF = hTwrPhiVsID -> GetNbinsY();
  const Double_t hBin1  = hTwrEtaVsID -> GetYaxis() -> GetBinLowEdge(1);
  const Double_t hBin2  = hTwrEtaVsID -> GetYaxis() -> GetBinLowEdge(nBinsH + 1);
  const Double_t fBin1  = hTwrPhiVsID -> GetYaxis() -> GetBinLowEdge(1);
  const Double_t fBin2  = hTwrPhiVsID -> GetYaxis() -> GetBinLowEdge(nBinsF + 1);
  TH1D *hHotTwrID    = new TH1D("hHotTwrID", "", nTow, tow1, tow2);
  TH2D *hTwrPhiVsEta = new TH2D("hTwrPhiVsEta", "", nBinsH, hBin1, hBin2, nBinsF, fBin1, fBin2);

  const UInt_t nTwrs = (UInt_t) hotTwrs.size();
  for (UInt_t iTwr = 0; iTwr < nTwrs; iTwr++) {

    const UInt_t twrId  = hotTwrs.at(iTwr);
    const UInt_t iBinId = hHotTwrID   -> FindBin(twrId);
    const UInt_t iBinIH = hTwrEtaVsID -> GetXaxis() -> FindBin(twrId);
    const UInt_t iBinIF = hTwrPhiVsID -> GetXaxis() -> FindBin(twrId);
    hHotTwrID -> SetBinContent(iBinId, 1.);

    Bool_t   goodEta = false;
    Bool_t   goodPhi = false;
    Double_t twrEta  = -2.;
    Double_t twrPhi  = -6.;
    for (UInt_t iBinH = 1; iBinH < nBinsH + 1; iBinH++) {
      const Double_t binCnt = hTwrEtaVsID -> GetBinContent(iBinIH, iBinH);
      if (binCnt != 0.) {
        twrEta = hTwrEtaVsID -> GetYaxis() -> GetBinCenter(iBinH);
        goodEta = true;
        break;
      }
    }
    for (UInt_t iBinF = 1; iBinF < nBinsF + 1; iBinF++) {
      const Double_t binCnt = hTwrPhiVsID -> GetBinContent(iBinIF, iBinF);
      if (binCnt != 0.) {
        twrPhi = hTwrPhiVsID -> GetYaxis() -> GetBinCenter(iBinF);
        goodPhi = true;
        break;
      }
    }

    const UInt_t iBinX = hTwrPhiVsEta -> GetXaxis() -> FindBin(twrEta);
    const UInt_t iBinY = hTwrPhiVsEta -> GetYaxis() -> FindBin(twrPhi);
    const UInt_t iBinG = hTwrPhiVsEta -> GetBin(iBinX, iBinY);
    if (goodEta && goodPhi) hTwrPhiVsEta -> SetBinContent(iBinG, 1.);

  }  // end tower loop
  cout << "    Map made." << endl;


  // change names and close files
  hTwrEtaVsID -> SetName("hTwrEtaVsID");
  hTwrPhiVsID -> SetName("hTwrPhiVsID");

  fOutput      -> cd();
  hHotTwrID    -> Write();
  hTwrEtaVsID  -> Write();
  hTwrPhiVsID  -> Write();
  hTwrPhiVsEta -> Write();
  fOutput      -> Close();
  fInputHist   -> cd();
  fInputHist   -> Close();
  cout << "  Map script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
