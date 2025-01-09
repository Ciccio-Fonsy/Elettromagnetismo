
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
  TCanvas*     canvas = new TCanvas("c1", "Analysis", 800, 600);

  TFile* file = TFile::Open("IstogrammiParticelle.root");

  if (!file || file->IsZombie()) {
    std::cerr << "Error during ROOT file opening." << std::endl;
    return 1;
  }

  // Caricamento degli istogrammi
  TH1D* hInvariantMass        = (TH1D*)file->Get("hInvariantMass");
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

  if (!hInvariantMass || !hMassOppositeSign || !hMassSameSign
      || !hMassPionKaonOpposite || !hMassPionKaonSame || !hMassKStarDecay
      || !hType || !hEnergy || !hTheta || !hPhi || !hPout || !hPtrasv) {
    std::cerr << "Error during histograms loading." << std::endl;
    file->Close();
    return 2;
  }

  // costruisco un vettore di istogrammi
  std::array<TH1D*, 12> h_array = {hInvariantMass,
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
    std::cout << "Entries in " << h->GetName() << ": " << h->GetEntries()
              << std::endl;
  }
  std::cout << '\n';

  // verifico distibuzione particelle
  for (int i = 1; i <= hType->GetNbinsX(); i++) {
    std::cout << "Number of " << hType->GetXaxis()->GetBinLabel(i) << ": "
              << std::setw(11) << hType->GetBinContent(i)
              << "  Percentage: " << std::setw(7)
              << hType->GetBinContent(i) / hType->GetEntries() * 100 << "%"
              << std::endl;
  }
  std::cout << '\n';

  // Fit della distribuzione angolare con una funzione uniforme
  TF1* f_uniform_phi = new TF1("f_uniform_phi", "[0]", 0, 2 * M_PI);

  f_uniform_phi->Print();
  std::cout << "\u03C6 uniform fit results:" << std::endl;
  hPhi->Fit("f_uniform_phi", "Q");
  std::cout << "  Constant:\t" << f_uniform_phi->GetParameter(0) << " +/- "
            << f_uniform_phi->GetParError(0) << std::endl;
  std::cout << "  Chi2/NDF:\t"
            << f_uniform_phi->GetChisquare() / f_uniform_phi->GetNDF()
            << std::endl;
  std::cout << "  Probability:\t" << f_uniform_phi->GetProb() << std::endl;
  std::cout << '\n';

  // Fit della distribuzione angolare con una funzione uniforme
  TF1* f_uniform_theta = new TF1("f_uniform_theta", "[0]", 0, M_PI);

  f_uniform_theta->Print();
  std::cout << "\u03B8 uniform fit results:" << std::endl;
  hTheta->Fit("f_uniform_theta", "Q");
  std::cout << "  Constant:\t" << f_uniform_theta->GetParameter(0) << " +/- "
            << f_uniform_theta->GetParError(0) << std::endl;
  std::cout << "  Chi2/NDF:\t"
            << f_uniform_theta->GetChisquare() / f_uniform_theta->GetNDF()
            << std::endl;
  std::cout << "  Probability:\t" << f_uniform_theta->GetProb() << std::endl;
  std::cout << '\n';

  // Fit della distribuzione del modulo dell'impulso con una funzione
  // esponenziale
  TF1* f_exponential = new TF1("f_exponential", "[0]*exp(-[1]*x)", 0, 6);

  f_exponential->Print();
  std::cout << "Momentum exponential fit results:" << std::endl;
  hPout->Fit("f_exponential", "Q");
  std::cout << "  Intercept:\t" << f_exponential->GetParameter(0) << " +/- "
            << f_exponential->GetParError(0) << std::endl;
  std::cout << "  Slope:\t" << f_exponential->GetParameter(1) << " +/- "
            << f_exponential->GetParError(1) << std::endl;
  std::cout << "  Chi2/NDF:\t"
            << f_exponential->GetChisquare() / f_exponential->GetNDF()
            << std::endl;
  std::cout << "  Probability:\t" << f_exponential->GetProb() << std::endl;
  std::cout << "Average momentum: " << hPout->GetMean() << std::endl;
  std::cout << '\n';

  // sottraggo gli istogrammi delle masse invarianti con carica opposta e stessa
  // carica
  TH1F* hInvariantMassSub1 = new TH1F(
      "hInvariantMassSub1",
      "Invariant mass of K* decays (blue), difference between concordant and discordant charge particles (green) and concordant and discordant #pi/K couples (red)",
      600, 0, 6);
  hInvariantMassSub1->Add(hMassSameSign, 1);
  hInvariantMassSub1->Add(hMassOppositeSign, -1);
  for (int i = 0; i < hInvariantMassSub1->GetNbinsX(); ++i) {
    hInvariantMassSub1->SetBinContent(
        i, std::abs(hInvariantMassSub1->GetBinContent(i)));
  }
  hInvariantMassSub1->GetXaxis()->SetTitle("Invariant mass [GeV/c^{2}]");
  hInvariantMassSub1->GetYaxis()->SetTitle("Events");
  hInvariantMassSub1->SetLineColor(kGreen);
  hInvariantMassSub1->SetStats(0);
  hInvariantMassSub1->Draw("HISTO");

  TH1F* hInvariantMassSub2 = new TH1F(
      "hInvariantMassSub2",
      "Invariant mass of K* decays (blue), difference between concordant and discordant charge particles (green) and concordant and discordant #pi/K couples (red)",
      600, 0, 6);
  hInvariantMassSub2->Add(hMassPionKaonOpposite, 1);
  hInvariantMassSub2->Add(hMassPionKaonSame, -1);
  for (int i = 0; i < hInvariantMassSub2->GetNbinsX(); ++i) {
    hInvariantMassSub2->SetBinContent(
        i, std::abs(hInvariantMassSub2->GetBinContent(i)));
  }
  hInvariantMassSub2->GetXaxis()->SetTitle("Invariant mass [GeV/c^{2}]");
  hInvariantMassSub2->GetYaxis()->SetTitle("Events");
  hInvariantMassSub2->SetLineColor(kRed);
  hInvariantMassSub2->SetStats(0);
  hInvariantMassSub2->Draw("SAME,HISTO");

  hMassKStarDecay->SetStats(0);
  hMassKStarDecay->Draw("SAME,HISTO");

  canvas->Update();

  TF1* f_gaussian = new TF1("f_gaussian", "gaus(0)", 0.7, 1.1);

  // Eseguiamo il fit sull'istogramma risultato dalla sottrazione
  hInvariantMassSub1->Fit("f_gaussian", "Q", nullptr, 0.7, 1.1);

  // Stampiamo i parametri del fit
  std::cout << "Concordant - discordant particles gaussian fit results:"
            << std::endl;
  std::cout << "K* mass: " << f_gaussian->GetParameter(1) << " +/- "
            << f_gaussian->GetParError(1) << std::endl;
  std::cout << "K* width: " << f_gaussian->GetParameter(2) << " +/- "
            << f_gaussian->GetParError(2) << std::endl;
  std::cout << "Chi2/NDF: " << f_gaussian->GetChisquare() / f_gaussian->GetNDF()
            << std::endl;
  std::cout << "Probability: " << f_gaussian->GetProb() << std::endl;
  std::cout << '\n';

  // Eseguiamo il fit sull'istogramma risultato dalla sottrazione
  hInvariantMassSub2->Fit("f_gaussian", "Q", nullptr, 0.7, 1.1);

  // Stampiamo i parametri del fit
  std::cout << "Concordant - discordant \u03C0/K couples gaussian fit results:"
            << std::endl;
  std::cout << "K* mass: " << f_gaussian->GetParameter(1) << " +/- "
            << f_gaussian->GetParError(1) << std::endl;
  std::cout << "K* width: " << f_gaussian->GetParameter(2) << " +/- "
            << f_gaussian->GetParError(2) << std::endl;
  std::cout << "Chi2/NDF: " << f_gaussian->GetChisquare() / f_gaussian->GetNDF()
            << std::endl;
  std::cout << "Probability: " << f_gaussian->GetProb() << std::endl;
  std::cout << '\n';

  while (gROOT->GetListOfCanvases()->FindObject("c1")) {
    gSystem->ProcessEvents(); // Process any events (including canvas events)
    gSystem->Sleep(100);      // Add a small delay to prevent 100% CPU usage
  }

  std::cout << "Canvas closed. Exiting program." << std::endl;

  return 0;
}
