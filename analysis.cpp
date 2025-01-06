#include "analysis.hpp"

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h>

#include <iostream>

void analyze_histograms() {
  // Apri il file ROOT
  TFile* file = TFile::Open("IstogrammiParticelle.root");

  if (!file || file->IsZombie()) {
    std::cerr << "Errore nell'apertura del file ROOT." << std::endl;
    return;
  }

  // Caricamento degli istogrammi
  TH1F* hMassInvariant        = (TH1F*)file->Get("hMassInvariant");
  TH1F* hMassOppositeSign     = (TH1F*)file->Get("hMassOppositeSign");
  TH1F* hMassSameSign         = (TH1F*)file->Get("hMassSameSign");
  TH1F* hMassPionKaonOpposite = (TH1F*)file->Get("hMassPionKaonOpposite");
  TH1F* hMassPionKaonSame     = (TH1F*)file->Get("hMassPionKaonSame");
  TH1F* hMassKStarDecay       = (TH1F*)file->Get("hMassKStarDecay");
  TH1F* h_type                = (TH1F*)file->Get("h_type");
  TH1F* h_energy              = (TH1F*)file->Get("h_energy");
  TH1F* h_theta               = (TH1F*)file->Get("h_theta");
  TH1F* h_phi                 = (TH1F*)file->Get("h_phi");
  TH1F* h_pout                = (TH1F*)file->Get("h_pout");
  TH1F* h_Ptrasv              = (TH1F*)file->Get("h_Ptrasv");

  if (!hMassInvariant || !hMassOppositeSign || !hMassSameSign
      || !hMassPionKaonOpposite || !hMassPionKaonSame || !hMassKStarDecay
      || !h_type || !h_energy || !h_theta || !h_phi || !h_pout || !h_Ptrasv) {
    std::cerr << "Errore nel caricamento degli istogrammi." << std::endl;
    file->Close();
    return;
  }

  // costruisco un vettore di istogrammi
  std::vector<TH1F*> h_vector;
  h_vector.push_back(hMassInvariant);
  h_vector.push_back(hMassOppositeSign);
  h_vector.push_back(hMassSameSign);
  h_vector.push_back(hMassPionKaonOpposite);
  h_vector.push_back(hMassPionKaonSame);
  h_vector.push_back(hMassKStarDecay);
  h_vector.push_back(h_type);
  h_vector.push_back(h_energy);
  h_vector.push_back(h_theta);
  h_vector.push_back(h_phi);
  h_vector.push_back(h_pout);
  h_vector.push_back(h_Ptrasv);

  // Numero di ingressi
  for (auto h : h_vector) {
    std::cout << "Numero di ingressi in " << h->GetName() << ": "
              << h->GetEntries() << std::endl;
  }
  std::cout << '\n';

  // verifico distibuzione particelle
  for (int i = 1; i <= h_type->GetNbinsX(); i++) {
    std::cout << "Contenuto del bin " << i << ": " << h_type->GetBinContent(i)
              << "\t+/- " << h_type->GetBinError(i) << "\tPrecntuale: "
              << h_type->GetBinContent(i) / h_type->GetEntries() * 100
              << "\t+/- " << h_type->GetBinError(i) / h_type->GetEntries() * 100
              << "%" << std::endl;
  }
  std::cout << '\n';

  // Fit della distribuzione angolare con una funzione uniforme
  h_phi->Fit("pol0",
             "Q"); // "pol0" è una funzione predefinita per un fit costante
  TF1* f_uniform = h_phi->GetFunction("pol0");

  if (f_uniform) {
    std::cout << "Parametri del fit uniforme:" << std::endl;
    std::cout << "  Costante: " << f_uniform->GetParameter(0) << " +/- "
              << f_uniform->GetParError(0) << std::endl;
    std::cout << "  Chi2/NDF: "
              << f_uniform->GetChisquare() / f_uniform->GetNDF() << std::endl;
    std::cout << "  Probabilità: " << f_uniform->GetProb() << std::endl;
  } else {
    std::cerr << "Errore nel recupero della funzione di fit." << std::endl;
  }
  std::cout << '\n';

  // Fit della distribuzione del modulo dell'impulso con una funzione
  // esponenziale
  h_pout->Fit("expo", "Q");
  TF1* f_exponential = h_pout->GetFunction("expo");

  if (f_exponential) {
    std::cout << "Parametri del fit esponenziale:" << std::endl;
    std::cout << "  Costante: " << f_exponential->GetParameter(0) << " +/- "
              << f_exponential->GetParError(0) << std::endl;
    std::cout << "  Slope: " << f_exponential->GetParameter(1) << " +/- "
              << f_exponential->GetParError(1) << std::endl;
    std::cout << "  Chi2/NDF: "
              << f_exponential->GetChisquare() / f_exponential->GetNDF()
              << std::endl;
    std::cout << "  Probabilità: " << f_exponential->GetProb() << std::endl;
  } else {
    std::cerr << "Errore nel recupero della funzione di fit." << std::endl;
  }
  std::cout << '\n';

  // Chiusura del file
  file->Close();
}

int main() { analyze_histograms(); }
