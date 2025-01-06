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
    TF1 *f_uniform = new TF1("f_uniform", "[0]", 0, 2 * M_PI);
    h_phi->Fit("f_uniform");
    f_uniform->Print();

    // Fit della distribuzione del modulo dell'impulso con una funzione esponenziale
    TF1 *f_exponential = new TF1("f_exponential", "[0]*exp(-[1]*x)", 0, 10);
    h_pout->Fit("f_exponential");
    f_exponential->Print();

    // Chi-quadro e probabilità del fit
    std::cout << "Chi2/NDF per il fit uniforme: " << f_uniform->GetChisquare() / f_uniform->GetNDF() << std::endl;
    std::cout << "Probabilità del fit uniforme: " << f_uniform->GetProb() << std::endl;

    std::cout << "Chi2/NDF per il fit esponenziale: " << f_exponential->GetChisquare() / f_exponential->GetNDF() << std::endl;
    std::cout << "Probabilità del fit esponenziale: " << f_exponential->GetProb() << std::endl;

    // Chiusura del file
    file->Close();
}
