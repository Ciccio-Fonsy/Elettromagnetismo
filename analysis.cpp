
#include "analysis.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iostream>

int main() {
  TApplication theApp("App", 0, 0);
  // Apri il file ROOT
  TCanvas*     canvas = new TCanvas("c1", "Canvas per Istogramma", 800, 600);

  TFile* file = TFile::Open("IstogrammiParticelle.root");

  if (!file || file->IsZombie()) {
    std::cerr << "Errore nell'apertura del file ROOT." << std::endl;
    return 1;
  }

  // Caricamento degli istogrammi
  TH1D* hMassInvariant        = (TH1D*)file->Get("hMassInvariant");
  TH1D* hMassSameSign         = (TH1D*)file->Get("hMassSameSign");
  TH1D* hMassOppositeSign     = (TH1D*)file->Get("hMassOppositeSign");
  TH1D* hMassPionKaonOpposite = (TH1D*)file->Get("hMassPionKaonOpposite");
  TH1D* hMassPionKaonSame     = (TH1D*)file->Get("hMassPionKaonSame");
  TH1D* hMassKStarDecay       = (TH1D*)file->Get("hMassKStarDecay");
  TH1D* hType                 = (TH1D*)file->Get("hType");
  TH1D* hEnergy               = (TH1D*)file->Get("hEnergy");
  TH1D* hTheta                = (TH1D*)file->Get("hTheta");
  TH1D* hPhi                  = (TH1D*)file->Get("hPhi");
  TH1D* hPout                 = (TH1D*)file->Get("hPout");
  TH1D* hPtrasv               = (TH1D*)file->Get("hPtrasv");

  if (!hMassInvariant || !hMassOppositeSign || !hMassSameSign
      || !hMassPionKaonOpposite || !hMassPionKaonSame || !hMassKStarDecay
      || !hType || !hEnergy || !hTheta || !hPhi || !hPout || !hPtrasv) {
    std::cerr << "Errore nel caricamento degli istogrammi." << std::endl;
    file->Close();
    return 2;
  }

  // costruisco un vettore di istogrammi
  std::array<TH1D*, 12> h_array = {hMassInvariant,
                                   hMassOppositeSign,
                                   hMassSameSign,
                                   hMassPionKaonOpposite,
                                   hMassPionKaonSame,
                                   hMassKStarDecay,
                                   hType,
                                   hEnergy,
                                   hTheta,
                                   hPhi,
                                   hPout,
                                   hPtrasv};

  // Numero di ingressi
  for (auto h : h_array) {
    std::cout << "Numero di ingressi in " << h->GetName() << ": "
              << h->GetEntries() << std::endl;
  }
  std::cout << '\n';

  // verifico distibuzione particelle
  for (int i = 1; i <= hType->GetNbinsX(); i++) {
    std::cout << "Contenuto del bin " << i << ": " << hType->GetBinContent(i)
              << " +/- " << hType->GetBinError(i) << "\tPrecntuale: "
              << hType->GetBinContent(i) / hType->GetEntries() * 100 << " +/- "
              << hType->GetBinError(i) / hType->GetEntries() * 100 << "%"
              << std::endl;
  }
  std::cout << '\n';

  // Fit della distribuzione angolare con una funzione uniforme
  hPhi->Fit("pol0", "Q", "", 0, 2 * M_PI);
  TF1* f_uniform_phi = hPhi->GetFunction("pol0");

  if (f_uniform_phi) {
    std::cout << "Parametri del fit uniforme di phi:" << std::endl;
    std::cout << "  Costante:\t" << f_uniform_phi->GetParameter(0) << " +/- "
              << f_uniform_phi->GetParError(0) << std::endl;
    std::cout << "  Chi2/NDF:\t"
              << f_uniform_phi->GetChisquare() / f_uniform_phi->GetNDF()
              << std::endl;
    std::cout << "  Probabilità:\t" << f_uniform_phi->GetProb() << std::endl;
  } else {
    std::cerr << "Errore nel recupero della funzione di fit." << std::endl;
  }
  std::cout << '\n';

  // Fit della distribuzione angolare con una funzione uniforme
  hTheta->Fit("pol0", "Q", "", 0, M_PI);
  TF1* f_uniform_theta = hTheta->GetFunction("pol0");

  if (f_uniform_theta) {
    std::cout << "Parametri del fit uniforme di theta:" << std::endl;
    std::cout << "  Costante:\t" << f_uniform_theta->GetParameter(0) << " +/- "
              << f_uniform_theta->GetParError(0) << std::endl;
    std::cout << "  Chi2/NDF:\t"
              << f_uniform_theta->GetChisquare() / f_uniform_theta->GetNDF()
              << std::endl;
    std::cout << "  Probabilità:\t" << f_uniform_theta->GetProb() << std::endl;
  } else {
    std::cerr << "Errore nel recupero della funzione di fit." << std::endl;
  }
  std::cout << '\n';

  // Fit della distribuzione del modulo dell'impulso con una funzione
  // esponenziale
  hPout->Fit("expo", "Q");
  TF1* f_exponential = hPout->GetFunction("expo");

  if (f_exponential) {
    std::cout << "Parametri del fit esponenziale:" << std::endl;
    std::cout << "  Costante:\t" << f_exponential->GetParameter(0) << " +/- "
              << f_exponential->GetParError(0) << std::endl;
    std::cout << "  Slope:\t" << f_exponential->GetParameter(1) << " +/- "
              << f_exponential->GetParError(1) << std::endl;
    std::cout << "  Chi2/NDF:\t"
              << f_exponential->GetChisquare() / f_exponential->GetNDF()
              << std::endl;
    std::cout << "  Probabilità:\t" << f_exponential->GetProb() << std::endl;
  } else {
    std::cerr << "Errore nel recupero della funzione di fit." << std::endl;
  }
  std::cout << '\n';
  std::cout << "Media pout:\n";
  std::cout << hPout->GetMean() << '\n';

  // sottraggo gli istogrammi delle masse invarianti con carica opposta e stessa
  // carica
  TH1F* hMassInvariantSub1 = new TH1F(
      "hMassInvariantSub1",
      "Massa invariante con carica opposta - stessa carica", 600, 0, 6);
  hMassInvariantSub1->Add(hMassSameSign, 1);
  hMassInvariantSub1->Add(hMassOppositeSign, -1);
  for (int i = 0; i < hMassInvariantSub1->GetNbinsX(); ++i) {
    hMassInvariantSub1->SetBinContent(
        i, std::abs(hMassInvariantSub1->GetBinContent(i)));
  }
  hMassInvariantSub1->GetXaxis()->SetTitle("Massa invariante [GeV/c^{2}]");
  hMassInvariantSub1->GetYaxis()->SetTitle("Eventi");
  hMassInvariantSub1->SetLineColor(kGreen);
  hMassInvariantSub1->Draw("HISTO");
  TH1F* hMassInvariantSub2 = new TH1F(
      "hMassInvariantSub2",
      "Massa invariante con carica opposta - stessa carica per pioni e kaoni",
      600, 0, 6);
  hMassInvariantSub2->Add(hMassPionKaonOpposite, 1);
  hMassInvariantSub2->Add(hMassPionKaonSame, -1);
  for (int i = 0; i < hMassInvariantSub2->GetNbinsX(); ++i) {
    hMassInvariantSub2->SetBinContent(
        i, std::abs(hMassInvariantSub2->GetBinContent(i)));
  }
  hMassInvariantSub2->GetXaxis()->SetTitle("Massa invariante [GeV/c^{2}]");
  hMassInvariantSub2->GetYaxis()->SetTitle("Eventi");
  hMassInvariantSub2->SetLineColor(kRed);
  hMassInvariantSub2->Draw("SAME,HISTO");
  hMassKStarDecay->Draw("SAME,HISTO");

  canvas->Update();

  while (gROOT->GetListOfCanvases()->FindObject("c1")) {
    gSystem->ProcessEvents(); // Process any events (including canvas events)
    gSystem->Sleep(100);      // Add a small delay to prevent 100% CPU usage
  }

  std::cout << "Canvas closed. Exiting program." << std::endl;

  return 0;
}
