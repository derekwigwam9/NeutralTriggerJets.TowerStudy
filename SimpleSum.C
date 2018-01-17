// 'SimpleSum.C'
// Derek Anderson
// 01.17.2018
//
// Use this to sum a set of histograms with
// a given set of weights from embedding.


#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TString.h"

using namespace std;


// constants
static const UInt_t   nHist(7);
static const UInt_t   nTotal(10);
static const Bool_t   doIntNorm(true);
static const Double_t weights[nTotal] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};




void SimpleSum() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning sum script..." << endl;

  // io parameters
  const TString  sOut("ePhadr.other.embedding.summed.d17m1y2018.root");
  const TString  sIn[nHist]   = {"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt5.match.root", "/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt7.match.root", "/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt9.match.root", "/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt11.match.root", "/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt15.match.root", "/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt25.match.root", "/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt35.match.root"};
  const TString  sHist[nHist] = {"hEbp_hadr", "hEbp_hadr", "hEbp_hadr", "hEbp_hadr", "hEbp_hadr", "hEbp_hadr", "hEbp_hadr"};
  const Double_t norms[nHist] = {72937., 138140., 158352., 119083., 115052., 88558., 27347.};


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
  TH1D *hHist[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    hHist[iHist] = (TH1D*) fIn[iHist] -> Get(sHist[iHist].Data());
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
  const UInt_t  nBins = hHist[0] -> GetNbinsX();
  const Float_t xBin1 = hHist[0] -> GetBinLowEdge(1);
  const Float_t xBin2 = hHist[0] -> GetBinLowEdge(nBins + 1);

  TH1D *hSum  = new TH1D("hSum", "", nBins, xBin1, xBin2);
  TH1D *hNorm = new TH1D("hNorm", "", nBins, xBin1, xBin2);
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
