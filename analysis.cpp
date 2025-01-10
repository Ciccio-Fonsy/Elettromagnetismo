
#include "analysis.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iomanip>
#include <iostream>

int main() {
  TApplication theApp("App", 0, 0);
  // Apri il file ROOT
  TCanvas*     canvas  = new TCanvas("c1", "Analysis", 800, 600);
  TCanvas*     canvas2 = new TCanvas("c2", "Analysis", 800, 600);
  canvas->Divide(3, 2);
  canvas->cd(1);

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

  hType->GetXaxis()->SetBinLabel(1, "\u03C0+");
  hType->GetXaxis()->SetBinLabel(2, "\u03C0-");

  // costruisco un vettore di istogrammi
  std::array<TH1D*, 12> h_array = {hType,
                                   hInvariantMass,
                                   hMassOppositeSign,
                                   hMassSameSign,
                                   hMassPionKaonOpposite,
                                   hMassPionKaonSame,
                                   hMassKStarDecay,
                                   hEnergy,
                                   hTheta,
                                   hPhi,
                                   hPout,
                                   hPtrasv};

  std::array<double, 12> expected_array = {
    1.01e+07, 5.05e+08, 2.55e+08, 2.50e+08, 4.40e+07, 4.40e+07,
    1.00e+05, 1.01e+07, 1.01e+07, 1.01e+07, 1.01e+07, 1.01e+07,
  };

  // Numero di ingressi
  for (int i = 0; i < 12; ++i) {
    std::cout << "Entries in " << std::setw(21) << h_array[i]->GetName() << ": "
              << std::setw(11) << h_array[i]->GetEntries()
              << ", expected: " << expected_array[i] << std::endl;
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
  TF1* f_uniform = new TF1("f_uniform", "[0] * 100000");
  f_uniform->Print();
  std::cout << '\n';

  std::cout << "\u03C6 uniform fit results:" << std::endl;
  hPhi->Fit("f_uniform", "Q");
  std::cout << "  Constant:\t" << f_uniform->GetParameter(0) << " ± "
            << f_uniform->GetParError(0) << std::endl;
  std::cout << "  \u03C7²:\t\t" << f_uniform->GetChisquare() << std::endl;
  std::cout << "  NDF:\t\t" << f_uniform->GetNDF() << std::endl;
  std::cout << "  \u03C7²/NDF:\t"
            << f_uniform->GetChisquare() / f_uniform->GetNDF() << std::endl;
  std::cout << "  Probability:\t" << f_uniform->GetProb() << std::endl;
  std::cout << '\n';
  hPhi->Draw("P SAME");

  canvas->cd(2);
  std::cout << "\u03B8 uniform fit results:" << std::endl;
  hTheta->Fit("f_uniform", "Q");
  std::cout << "  Constant:\t" << f_uniform->GetParameter(0) << " ± "
            << f_uniform->GetParError(0) << std::endl;
  std::cout << "  \u03C7²:\t\t" << f_uniform->GetChisquare() << std::endl;
  std::cout << "  NDF:\t\t" << f_uniform->GetNDF() << std::endl;
  std::cout << "  \u03C7²/NDF:\t"
            << f_uniform->GetChisquare() / f_uniform->GetNDF() << std::endl;
  std::cout << "  Probability:\t" << f_uniform->GetProb() << std::endl;
  std::cout << '\n';
  hTheta->Draw("P SAME");

  canvas->cd(3);
  // Fit della distribuzione del modulo dell'impulso con una funzione
  // esponenziale
  TF1* f_exponential = new TF1("f_exponential", "[0]*exp(-[1]*x)");

  f_exponential->Print();
  std::cout << "Momentum exponential fit results:" << std::endl;
  hPout->Fit("f_exponential", "Q");
  std::cout << "  Intercept:\t" << f_exponential->GetParameter(0) << " ± "
            << f_exponential->GetParError(0) << std::endl;
  std::cout << "  Slope:\t" << f_exponential->GetParameter(1) << " ± "
            << f_exponential->GetParError(1) << std::endl;
  std::cout << "  Average:\t" << f_exponential->Mean(0, 6) << " ± "
            << f_exponential->GetParError(1) * f_exponential->Mean(0, 6)
                   * f_exponential->Mean(0, 6)
            << std::endl;
  std::cout << "  \u03C7²:\t\t" << f_exponential->GetChisquare() << std::endl;
  std::cout << "  NDF:\t\t" << f_exponential->GetNDF() << std::endl;
  std::cout << "  \u03C7²/NDF:\t"
            << f_exponential->GetChisquare() / f_exponential->GetNDF()
            << std::endl;
  std::cout << "  Probability:\t" << f_exponential->GetProb() << std::endl;
  std::cout << '\n';
  hPout->Draw("P SAME");

  canvas->cd(4);
  TF1* f_gaussian = new TF1("f_gaussian", "gaus(0)", 0.7, 1.1);
  f_gaussian->SetLineColor(kBlack);

  // sottraggo gli istogrammi delle masse invarianti con carica opposta e stessa
  // carica
  TH1F* hInvariantMassSub1 = new TH1F(
      "hInvariantMassSub1",
      "Difference in invariant mass between concordant and discordant charge particles",
      hMassSameSign->GetNbinsX(), 0, 6);
  hInvariantMassSub1->Add(hMassOppositeSign, 1);
  hInvariantMassSub1->Add(hMassSameSign, -1);
  hInvariantMassSub1->GetXaxis()->SetTitle("Invariant mass [GeV/c^{2}]");
  hInvariantMassSub1->GetYaxis()->SetTitle("Events");
  hInvariantMassSub1->SetLineColor(kGreen);

  // Eseguiamo il fit sull'istogramma risultato dalla sottrazione
  hInvariantMassSub1->Fit("f_gaussian", "Q", nullptr, 0.7, 1.1);

  hInvariantMassSub1->Draw("HIST SAME");

  // Stampiamo i parametri del fit
  std::cout << "Concordant - discordant particles gaussian fit results:"
            << std::endl;
  std::cout << "  Amplitude:\t" << f_gaussian->GetParameter(0) << " ± "
            << f_gaussian->GetParError(0) << std::endl;
  std::cout << "  K* mass:\t" << f_gaussian->GetParameter(1) << " ± "
            << f_gaussian->GetParError(1) << std::endl;
  std::cout << "  K* width:\t" << f_gaussian->GetParameter(2) << " ± "
            << f_gaussian->GetParError(2) << std::endl;
  std::cout << "  \u03C7²:\t\t" << f_gaussian->GetChisquare() << std::endl;
  std::cout << "  NDF:\t\t" << f_gaussian->GetNDF() << std::endl;
  std::cout << "  \u03C7²/NDF:\t"
            << f_gaussian->GetChisquare() / f_gaussian->GetNDF() << std::endl;
  std::cout << "  Probability:\t" << f_gaussian->GetProb() << std::endl;
  std::cout << '\n';

  canvas->cd(5);

  TH1F* hInvariantMassSub2 = new TH1F(
      "hInvariantMassSub2",
      "Difference in invariant mass between concordant and discordant #pi/K couples",
      hMassPionKaonSame->GetNbinsX(), 0, 6);
  hInvariantMassSub2->Add(hMassPionKaonOpposite, 1);
  hInvariantMassSub2->Add(hMassPionKaonSame, -1);
  hInvariantMassSub2->GetXaxis()->SetTitle("Invariant mass [GeV/c^{2}]");
  hInvariantMassSub2->GetYaxis()->SetTitle("Events");
  hInvariantMassSub2->SetLineColor(kRed);

  // Eseguiamo il fit sull'istogramma risultato dalla sottrazione
  hInvariantMassSub2->Fit("f_gaussian", "Q", nullptr, 0.7, 1.1);

  hInvariantMassSub2->Draw("HIST SAME");

  // Stampiamo i parametri del fit
  std::cout << "Concordant - discordant \u03C0/K couples gaussian fit results:"
            << std::endl;
  std::cout << "  Amplitude:\t" << f_gaussian->GetParameter(0) << " ± "
            << f_gaussian->GetParError(0) << std::endl;
  std::cout << "  K* mass:\t" << f_gaussian->GetParameter(1) << " ± "
            << f_gaussian->GetParError(1) << std::endl;
  std::cout << "  K* width:\t" << f_gaussian->GetParameter(2) << " ± "
            << f_gaussian->GetParError(2) << std::endl;
  std::cout << "  \u03C7²:\t\t" << f_gaussian->GetChisquare() << std::endl;
  std::cout << "  NDF:\t\t" << f_gaussian->GetNDF() << std::endl;
  std::cout << "  \u03C7²/NDF:\t"
            << f_gaussian->GetChisquare() / f_gaussian->GetNDF() << std::endl;
  std::cout << "  Probability:\t" << f_gaussian->GetProb() << std::endl;
  std::cout << '\n';

  canvas->cd(6);

  // Eseguiamo il fit sull'istogramma risultato dalla sottrazione
  hMassKStarDecay->Fit("f_gaussian", "Q", nullptr, 0.7, 1.1);

  hMassKStarDecay->Draw("HIST SAME");

  // Stampiamo i parametri del fit
  std::cout << "True K* decays gaussian fit results:" << std::endl;
  std::cout << "  Amplitude:\t" << f_gaussian->GetParameter(0) << " ± "
            << f_gaussian->GetParError(0) << std::endl;
  std::cout << "  K* mass:\t" << f_gaussian->GetParameter(1) << " ± "
            << f_gaussian->GetParError(1) << std::endl;
  std::cout << "  K* width:\t" << f_gaussian->GetParameter(2) << " ± "
            << f_gaussian->GetParError(2) << std::endl;
  std::cout << "  \u03C7²:\t\t" << f_gaussian->GetChisquare() << std::endl;
  std::cout << "  NDF:\t\t" << f_gaussian->GetNDF() << std::endl;
  std::cout << "  \u03C7²/NDF:\t"
            << f_gaussian->GetChisquare() / f_gaussian->GetNDF() << std::endl;
  std::cout << "  Probability:\t" << f_gaussian->GetProb() << std::endl;
  std::cout << '\n';

  canvas2->cd();
  hInvariantMassSub1->Draw("HIST");
  hInvariantMassSub2->Draw("HIST SAME");
  hMassKStarDecay->Draw("HIST SAME");

  while (gROOT->GetListOfCanvases()->FindObject("c1")
         && gROOT->GetListOfCanvases()->FindObject("c2")) {
    gSystem->ProcessEvents(); // Process any events (including canvas events)
    gSystem->Sleep(100);      // Add a small delay to prevent 100% CPU usage
  }

  std::cout << "Canvas closed. Exiting program." << std::endl;

  return 0;
}
