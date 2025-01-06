#include "analysis.hpp"
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <iostream>

void analyze_histograms() {
    // Apri il file ROOT
    TFile *file = TFile::Open("IstogrammiParticelle.root");
    
    if (!file || file->IsZombie()) {
        std::cerr << "Errore nell'apertura del file ROOT." << std::endl;
        return;
    }

  
    // Caricamento degli istogrammi
    TH1F *hMassInvariant = (TH1F*)file->Get("hMassInvariant");
    TH1F *hMassOppositeSign = (TH1F*)file->Get("hMassOppositeSign");
    TH1F *hMassSameSign = (TH1F*)file->Get("hMassSameSign");
    TH1F *hMassPionKaonOpposite = (TH1F*)file->Get("hMassPionKaonOpposite");
    TH1F *hMassPionKaonSame = (TH1F*)file->Get("hMassPionKaonSame");
    TH1F *hMassKStarDecay = (TH1F*)file->Get("hMassKStarDecay");
    TH1F *h_type = (TH1F*)file->Get("h_type");
    TH1F *h_energy = (TH1F*)file->Get("h_energy");
    TH1F *h_theta = (TH1F*)file->Get("h_theta");
    TH1F *h_phi = (TH1F*)file->Get("h_phi");
    TH1F *h_pout = (TH1F*)file->Get("h_pout");
    TH1F *h_Ptrasv = (TH1F*)file->Get("h_Ptrasv");



    if (!h_phi || !h_pout || !hMassInvariant) {
        std::cerr << "Errore nel caricamento degli istogrammi." << std::endl;
        file->Close();
        return;
    }

    // Numero di ingressi
    std::cout << "Numero di ingressi in h_phi: " << h_phi->GetEntries() << std::endl;
    std::cout << "Numero di ingressi in h_pout: " << h_pout->GetEntries() << std::endl;
    std::cout << "Numero di ingressi in h_mass_invariant: " << hMassInvariant->GetEntries() << std::endl;

    // Fit della distribuzione angolare con una funzione uniforme
  h_phi->Fit("pol0"); // "pol0" è una funzione predefinita per un fit costante
  TF1* f_uniform = h_phi->GetFunction("pol0");

  if (f_uniform) {
    f_uniform->Print();
    std::cout << "Chi2/NDF per il fit uniforme: "
              << f_uniform->GetChisquare() / f_uniform->GetNDF() << std::endl;
    std::cout << "Probabilità del fit uniforme: " << f_uniform->GetProb()
              << std::endl;
  } else {
    std::cerr << "Errore nel recupero della funzione di fit." << std::endl;
  }

  // Fit della distribuzione del modulo dell'impulso con una funzione
  // esponenziale
  h_pout->Fit("expo");
  TF1* f_exponential = h_pout->GetFunction("expo");

  if (f_exponential) {
    f_exponential->Print();
    std::cout << "Chi2/NDF per il fit esponenziale: "
              << f_exponential->GetChisquare() / f_exponential->GetNDF() << std::endl;
    std::cout << "Probabilità del fit esponenziale: " << f_exponential->GetProb()
              << std::endl;
  } else {
    std::cerr << "Errore nel recupero della funzione di fit." << std::endl;
  }

  // Chiusura del file
  file->Close();
}
