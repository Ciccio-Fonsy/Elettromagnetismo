
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
  TH1F* hMassInvariant        = (TH1F*)file->Get("hMassInvariant");
  TH1F* hMassOppositeSign     = (TH1F*)file->Get("hMassOppositeSign");
  TH1F* hMassSameSign         = (TH1F*)file->Get("hMassSameSign");
  TH1F* hMassPionKaonOpposite = (TH1F*)file->Get("hMassPionKaonOpposite");
  TH1F* hMassPionKaonSame     = (TH1F*)file->Get("hMassPionKaonSame");
  TH1F* hMassKStarDecay       = (TH1F*)file->Get("hMassKStarDecay");
  TH1F* hType                 = (TH1F*)file->Get("hType");
  TH1F* hEnergy               = (TH1F*)file->Get("hEnergy");
  TH1F* hTheta                = (TH1F*)file->Get("hTheta");
  TH1F* hPhi                  = (TH1F*)file->Get("hPhi");
  TH1F* hPout                 = (TH1F*)file->Get("hPout");
  TH1F* hPtrasv               = (TH1F*)file->Get("hPtrasv");

  if (!hMassInvariant || !hMassOppositeSign || !hMassSameSign
      || !hMassPionKaonOpposite || !hMassPionKaonSame || !hMassKStarDecay
      || !hType || !hEnergy || !hTheta || !hPhi || !hPout || !hPtrasv) {
    std::cerr << "Errore nel caricamento degli istogrammi." << std::endl;
    file->Close();
    return 2;
  }

  // costruisco un vettore di istogrammi
  std::vector<TH1F*> h_vector;
  h_vector.push_back(hMassInvariant);
  h_vector.push_back(hMassOppositeSign);
  h_vector.push_back(hMassSameSign);
  h_vector.push_back(hMassPionKaonOpposite);
  h_vector.push_back(hMassPionKaonSame);
  h_vector.push_back(hMassKStarDecay);
  h_vector.push_back(hType);
  h_vector.push_back(hEnergy);
  h_vector.push_back(hTheta);
  h_vector.push_back(hPhi);
  h_vector.push_back(hPout);
  h_vector.push_back(hPtrasv);

  // Numero di ingressi
  for (auto h : h_vector) {
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

  // sottraggo gli istogrammi delle masse invarianti con carica opposta e stessa
  // carica
  TH1F* hMassInvariantSub1 = new TH1F(
      "hMassInvariantSub1",
      "Massa invariante con carica opposta - stessa carica", 300, 0, 6);
  hMassInvariantSub1->Add(hMassSameSign, 1);
  hMassInvariantSub1->Add(hMassOppositeSign, -1);
  hMassInvariantSub1->GetXaxis()->SetTitle("Massa invariante [GeV/c^{2}]");
  hMassInvariantSub1->GetYaxis()->SetTitle("Eventi");
  hMassInvariantSub1->SetLineColor(kGreen);
  hMassInvariantSub1->Draw("APE");
  TH1F* hMassInvariantSub2 = new TH1F(
      "hMassInvariantSub2",
      "Massa invariante con carica opposta - stessa carica per pioni e kaoni",
      300, 0, 6);
  hMassInvariantSub2->Add(hMassPionKaonOpposite, 1);
  hMassInvariantSub2->Add(hMassPionKaonSame, -1);
  hMassInvariantSub2->GetXaxis()->SetTitle("Massa invariante [GeV/c^{2}]");
  hMassInvariantSub2->GetYaxis()->SetTitle("Eventi");
  hMassInvariantSub2->SetLineColor(kRed);
  hMassInvariantSub2->Draw("SAME,APE");
  hMassKStarDecay->Draw("SAME,APE");

  canvas->Update();

  while (gROOT->GetListOfCanvases()->FindObject("c1")) {
    gSystem->ProcessEvents(); // Process any events (including canvas events)
    gSystem->Sleep(100);      // Add a small delay to prevent 100% CPU usage
  }

  std::cout << "Canvas closed. Exiting program." << std::endl;

  return 0;
}
