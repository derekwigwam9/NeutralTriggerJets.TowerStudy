// 'HotTowerStudy.C'
// Derek Anderson
// 07.26.2017
//
// This takes the output of 'DataTowerStudy.C' or
// 'EmbeddingTowerStudy.C' and uses it to create
// a list of hot towers.

#include <vector>
#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;

// global constants
static const UInt_t  nEneBins(6);
// i/o parameters
static const TString sList("hotTwrTest.list");
static const TString sInput("pp200r9.hotTowerMapWithDeDxForReal.root");
static const TString sOutput("hotTwrTest.root");
static const TString sBins[nEneBins] = {"_e0", "_e2", "_e5", "_e10", "_e20", "_e35"};



void HotTowerStudy() {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning hot tower script..." << endl;


  // constants
  const TString  sTwrID("hTwrId");
  const TString  sTwrEne("hTwrEne");
  const TString  sTwrIDvsEne("hTwrEneVsId");
  const TString  sTwrEtaVsID("hTwrEtaVsId_e0");
  const TString  sTwrPhiVsID("hTwrPhiVsId_e0");
  const Double_t nSigMax(3.);


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput  = new TFile(sInput.Data(), "read");
  if (!fOutput || !fInput) {
    cerr << "PANIC: couldn't open file!\n"
         << "       fInput  -- " << fInput << "\n"
         << "       fOutput -- " << fOutput
         << endl;
    return;
  }
  cout << "    Files opened." << endl;

  TH1D *hTwrID[nEneBins];
  TH1D *hTwrEne[nEneBins];
  TH2D *hTwrIDvsEne[nEneBins];
  for (UInt_t iEneBin = 0; iEneBin < nEneBins; iEneBin++) {

    TString sId(sTwrID.Data());
    TString sEne(sTwrEne.Data());
    TString sIdVsEne(sTwrIDvsEne.Data());
    sId      += sBins[iEneBin].Data();
    sEne     += sBins[iEneBin].Data();
    sIdVsEne += sBins[iEneBin].Data();

    hTwrID[iEneBin]      = (TH1D*) fInput -> Get(sId.Data());
    hTwrEne[iEneBin]     = (TH1D*) fInput -> Get(sEne.Data());
    hTwrIDvsEne[iEneBin] = (TH2D*) fInput -> Get(sIdVsEne.Data());
    if (!hTwrID[iEneBin] || !hTwrEne[iEneBin] || !hTwrIDvsEne[iEneBin]) {
      cerr << "PANIC: couldn't grab histogram " << iEneBin << "!\n"
           << "       hTwrID      -- " << hTwrID[iEneBin] << "\n"
           << "       hTwrEne     -- " << hTwrEne[iEneBin] << "\n"
           << "       hTwrIdVsEne -- " << hTwrIDvsEne[iEneBin]
           << endl;
      return;
    }
  }  // end bin loop
  TH1D *hTwrEtaVsID = (TH1D*) fInput -> Get(sTwrEtaVsID.Data());
  TH1D *hTwrPhiVsID = (TH1D*) fInput -> Get(sTwrPhiVsID.Data());
  cout << "    Histograms grabbed." << endl;


  // create hot tower list per bin
  UInt_t         nHotTwrs(0);
  Double_t       freqMean[nEneBins];
  Double_t       freqSigma[nEneBins];
  vector<UInt_t> twrBin[nEneBins];
  for (UInt_t iEneBin = 0; iEneBin < nEneBins; iEneBin++) {

    const UInt_t nTwr = hTwrID[iEneBin] -> GetNbinsX();
    nHotTwrs = 0;

    // calculate mean
    Double_t fSum  = 0.;
    Double_t fMean = 0.;
    for (UInt_t iTwr = 1; iTwr < nTwr + 1; iTwr++) {
      const Double_t freq = hTwrID[iEneBin] -> GetBinContent(iTwr);
      fSum += freq;
    }
    fMean = fSum / (Double_t) nTwr;

    // calculate std. dev.
    Double_t fSum2  = 0.;
    Double_t fSigma = 0.;
    for (UInt_t iTwr = 1; iTwr < nTwr + 1; iTwr++) {
      const Double_t freq  = hTwrID[iEneBin] -> GetBinContent(iTwr);
      const Double_t fDif  = freq - fMean;
      const Double_t fDif2 = fDif * fDif;
      fSum2 += fDif2;
    }
    fSigma = TMath::Sqrt(fSum2 / (Double_t) (nTwr - 1));

    // determine hot towers
    for (UInt_t iTwr = 1; iTwr < nTwr + 1; iTwr++) {
      const Double_t freq = hTwrID[iEneBin] -> GetBinContent(iTwr);
      const Double_t fDif = TMath::Abs(freq - fMean);
      const Double_t nSig = fDif / fSigma;
      if (nSig > nSigMax) {
        twrBin[iEneBin].push_back(iTwr);
        nHotTwrs++;
      }
    }

    freqMean[iEneBin]  = fMean;
    freqSigma[iEneBin] = fSigma;
    cout << "      Bin " << iEneBin << ": mean = " << fMean
         << ", sigma = " << fSigma << ".\n"
         << "             Found " << nHotTwrs << " hot towers."
         << endl;

  }  // end bin loop
  cout << "    Lists created." << endl;


  // consolidate lists
  vector<UInt_t> twrTot;
  for (UInt_t iEneBin = 0; iEneBin < nEneBins; iEneBin++) {
    for (UInt_t iTwrBin = 0; iTwrBin < twrBin[iEneBin].size(); iTwrBin++) {

      // determine if tower should be added
      Bool_t addTwr = false;
      Bool_t isNew  = true;
      if (iEneBin == 0)
        addTwr = true;
      else {
        for (UInt_t iTwrTot = 0; iTwrTot < twrTot.size(); iTwrTot++) {
          if (twrBin[iEneBin][iTwrBin] ==  twrTot[iTwrTot]) {
            isNew = false;
            break;
          }
        }  // end total list loop
        if (isNew) addTwr = true;
      }

      // add to list if necessary
      if (addTwr) twrTot.push_back(twrBin[iEneBin][iTwrBin]);

    }  // end per bin list loop
  }  // end bin loop
  cout << "    Consolidated list: \n"
       << "      " << twrTot.size() << " towers total." << endl;


  // sort list
  vector<UInt_t> finalList;

  TH1D *hTowerList = new TH1D("hTowerList", "Hot Towers", 4800, 0., 4800.);
  for (UInt_t iTwr = 0; iTwr < twrTot.size(); iTwr++) {
    hTowerList -> Fill(twrTot[iTwr]);
  }

  const UInt_t nBins = hTowerList -> GetNbinsX();
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {
    const UInt_t bin = (UInt_t) hTowerList -> GetBinContent(iBin);
    if (bin != 0) finalList.push_back((hTowerList -> GetBinLowEdge(iBin)) - 1);
  }
  cout << "    Sorted list: \n"
       << "      " << finalList.size() << " towers total." << endl;


  // create eta-phi map
  const UInt_t   nBinsH = hTwrEtaVsID -> GetNbinsY();
  const UInt_t   nBinsF = hTwrPhiVsID -> GetNbinsY();
  const Double_t hBin1  = hTwrEtaVsID -> GetYaxis() -> GetBinLowEdge(1);
  const Double_t hBin2  = hTwrEtaVsID -> GetYaxis() -> GetBinLowEdge(nBinsH + 1);
  const Double_t fBin1  = hTwrPhiVsID -> GetYaxis() -> GetBinLowEdge(1);
  const Double_t fBin2  = hTwrPhiVsID -> GetYaxis() -> GetBinLowEdge(nBinsF + 1);
  TH2D *hHotTwrPhiVsEta = new TH2D("hHotTwrPhiVsEta", "", nBinsH, hBin1, hBin2, nBinsF, fBin1, fBin2);

  const UInt_t nTwrs = (UInt_t) finalList.size();
  for (UInt_t iTwr = 0; iTwr < nTwrs; iTwr++) {

    const UInt_t twrId  = finalList.at(iTwr);
    const UInt_t iBinIH = hTwrEtaVsID -> GetXaxis() -> FindBin(twrId);
    const UInt_t iBinIF = hTwrPhiVsID -> GetXaxis() -> FindBin(twrId);

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

    const UInt_t iBinX = hHotTwrPhiVsEta -> GetXaxis() -> FindBin(twrEta);
    const UInt_t iBinY = hHotTwrPhiVsEta -> GetYaxis() -> FindBin(twrPhi);
    const UInt_t iBinG = hHotTwrPhiVsEta -> GetBin(iBinX, iBinY);
    if (goodEta && goodPhi) hHotTwrPhiVsEta -> SetBinContent(iBinG, 1.);

  }
  cout << "    Eta-phi map created." << endl;


  // write list to file
  ofstream towers(sList.Data());
  for (UInt_t iTwr = 0; iTwr < finalList.size(); iTwr++) {
    towers << finalList[iTwr];
    if ((iTwr + 1) != finalList.size())
      towers << endl;
  }
  towers.close();

  fOutput         -> cd();
  hTowerList      -> Write();
  hHotTwrPhiVsEta -> Write();
  fOutput         -> Close();
  fInput          -> cd();
  fInput          -> Close();
  cout << "  Hot tower script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
