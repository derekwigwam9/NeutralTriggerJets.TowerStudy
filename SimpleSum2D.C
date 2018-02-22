// 'SimpleSum2D.C'
// Derek Anderson
// 02.22.2018
//
// Use this to sum a set of histograms with
// a given set of weights from embedding.


#include <iostream>
#include "TH2.h"
#include "TFile.h"
#include "TString.h"

using namespace std;


// constants
static const UInt_t   nHist(7);
static const UInt_t   nTotal(10);
static const Bool_t   doIntNorm(true);
static const Double_t weights[nTotal] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};




void SimpleSum2D() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning sum script..." << endl;

  // io parameters
  const TString  sOut("pp200r9embed.phiVsEtaTrg.et9vz55had.d22m2y2018.root");
  const TString  sIn[nHist]   = {"pp200r9pt5u.fTrgCheck.et9vz55had.root", "pp200r9pt7u.fTrgCheck.et9vz55had.root", "pp200r9pt9u.fTrgCheck.et9vz55had.root", "pp200r9pt11u.fTrgCheck.et9vz55had.root", "pp200r9pt15u.fTrgCheck.et9vz55had.root", "pp200r9pt25u.fTrgCheck.et9vz55had.root", "pp200r9pt35u.fTrgCheck.et9vz55had.root"};
  const TString  sHist[nHist] = {"hTrgPhiVsEtaTwr", "hTrgPhiVsEtaTwr", "hTrgPhiVsEtaTwr", "hTrgPhiVsEtaTwr", "hTrgPhiVsEtaTwr", "hTrgPhiVsEtaTwr", "hTrgPhiVsEtaTwr"};
  const Double_t norms[nHist] = {10., 89., 468., 1495., 6800., 19287., 9099.};


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fIn[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    fIn[iHist] = new TFile(sIn[iHist].Data(), "read");
    if (!fIn[iHist]) {
      cerr << "PANIC: couldn't open input file no. " << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // get histograms
  TH2D *hHist[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    hHist[iHist] = (TH2D*) fIn[iHist] -> Get(sHist[iHist].Data());
    if (!hHist[iHist]) {
      cerr << "PANIC: couldn't grab histogram no. " << iHist << "!" << endl;
      return;
    }
  }
  fOut -> cd();
  cout << "    Grabbed histograms." << endl;


  // scale histograms
  Double_t normer(0.);
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    hHist[iHist] -> Scale(weights[(nTotal - nHist) + iHist]);
    normer += norms[iHist] * weights[(nTotal - nHist) + iHist];
  }
  cout << "    Scaled histograms." << endl;

  // sum histograms
  const UInt_t  nBinsX = hHist[0] -> GetNbinsX();
  const UInt_t  nBinsY = hHist[0] -> GetNbinsY();
  const Float_t xBin1  = hHist[0] -> GetXaxis() -> GetBinLowEdge(1);
  const Float_t xBin2  = hHist[0] -> GetXaxis() -> GetBinLowEdge(nBinsX + 1);
  const Float_t yBin1  = hHist[0] -> GetYaxis() -> GetBinLowEdge(1);
  const Float_t yBin2  = hHist[0] -> GetYaxis() -> GetBinLowEdge(nBinsY + 1);

  TH2D *hSum  = new TH2D("hSum", "", nBinsX, xBin1, xBin2, nBinsY, yBin1, yBin2);
  TH2D *hNorm = new TH2D("hNorm", "", nBinsX, xBin1, xBin2, nBinsY, yBin1, yBin2);
  hSum  -> Sumw2();
  hNorm -> Sumw2();
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    hSum  -> Add(hHist[iHist]);
    hNorm -> Add(hHist[iHist]);
  }

  // normalize sum
  if (doIntNorm) {
    const Double_t integral = hNorm -> Integral();
    hNorm -> Scale(1. / integral);
    cout << "    Normalized histogram.\n"
         << "      normalization = " << integral
         << endl;
  }
  else {
    hNorm -> Scale(1. / normer);
    cout << "    Normalized histogram.\n"
         << "      normalization = " << normer
         << endl;
  }


  // save and close
  fOut  -> cd();
  hSum  -> Write();
  hNorm -> Write();
  fOut  -> Close();
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    fIn[iHist] -> cd();
    fIn[iHist] -> Close();
  }
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
