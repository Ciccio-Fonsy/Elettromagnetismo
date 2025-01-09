#include "generator.hpp"

#include "particle.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>

#include <iostream>

void     createInstances(int n_event, int n_particles_event, TRandom3* rand,
                         std::array<Particle, 120>& event_particles,
                         TH1D* hMassKStarDecay, TH1D* hInvariantMass,
                         TH1D* hMassOppositeSign, TH1D* hMassSameSign,
                         TH1D* hMassPionKaonOpposite, TH1D* hMassPionKaonSame,
                         TH1D* hType, TH1D* hEnergy, TH1D* hTheta, TH1D* hPhi,
                         TH1D* hPout, TH1D* hPtrasv);
void     randomParticlePosition(TRandom3* rand, Particle& particle);
Particle createRandomParticle(TRandom3* rand);
void     fillHistogram(const std::array<Particle, 120>& particle,
                       TH1D* hInvariantMass, TH1D* hMassOppositeSign,
                       TH1D* hMassSameSign, TH1D* hMassPionKaonOpposite,
                       TH1D* hMassPionKaonSame, TH1D* hType, TH1D* hEnergy,
                       TH1D* hTheta, TH1D* hPhi, TH1D* hPout, TH1D* hPtrasv,
                       int n_particles);
void     saveHistograms(const std::array<TH1D*, 12>& histograms,
                        const std::string&           filename);

int main() {
  // Crea un'applicazione ROOT
  TApplication theApp("App", nullptr, nullptr);

  std::cout << "Generating canvas..." << std::endl;
  TCanvas* canvas = new TCanvas("c1", "Generation", 800, 600);
  canvas->Divide(4, 3); // Dividi il canvas in 3 righe e 4 colonne

  // name, mass, charge, whidth
  Particle::addParticleType('Q', 0.13957, 1);
  Particle::addParticleType('q', 0.13957, -1);
  Particle::addParticleType('P', 0.93827, 1);
  Particle::addParticleType('p', 0.93827, -1);
  Particle::addParticleType('K', 0.49367, 1);
  Particle::addParticleType('k', 0.49367, -1);
  Particle::addParticleType('*', 0.89166, 0, 0.050);

  int                       n_event           = 100000;
  int                       n_particles_event = 100;
  std::array<Particle, 120> event_particles;

  // secondo Massimo "si dovrebbero cambiare i nomi dei grafici con dei
  // numeri così alla fine si può usare un ciclo for per scrivere i grafici
  // a fine codice e per il metodo sumw :)"

  // range tra 0 e 3, il numero dei bin è , non 60 se no è troopo largo,
  // 180.

  TH1D* hInvariantMass =
      new TH1D("hInvariantMass", "Invariant Mass distribution", 600, 0, 6);
  hInvariantMass->Sumw2();
  TH1D* hMassOppositeSign = new TH1D(
      "hMassOppositeSign", "Invariant Mass - Discordant Charge", 600, 0, 6);
  hMassOppositeSign->Sumw2();
  TH1D* hMassSameSign = new TH1D(
      "hMassSameSign", "Invariant Mass - Concordant Charge", 600, 0, 6);
  hMassSameSign->Sumw2();
  TH1D* hMassPionKaonOpposite = new TH1D(
      "hMassPionKaonOpposite", "Invariant Mass- Pi/K discordant", 600, 0, 6);
  hMassPionKaonOpposite->Sumw2();
  TH1D* hMassPionKaonSame = new TH1D(
      "hMassPionKaonSame", "Invariant Mass - Pi/K concordant", 600, 0, 6);
  hMassPionKaonSame->Sumw2();
  TH1D* hMassKStarDecay =
      new TH1D("hMassKStarDecay", "Invariant Mass - K* decay", 100, 0.7, 1.07);
  hMassKStarDecay->Sumw2();

  TH1D* hType = new TH1D("hType", "Particle type distribution", 7, 0, 7);
  hType->GetXaxis()->SetBinLabel(1, "\u03C0+");
  hType->GetXaxis()->SetBinLabel(2, "\u03C0-");
  hType->GetXaxis()->SetBinLabel(3, "P+");
  hType->GetXaxis()->SetBinLabel(4, "P-");
  hType->GetXaxis()->SetBinLabel(5, "K+");
  hType->GetXaxis()->SetBinLabel(6, "K-");
  hType->GetXaxis()->SetBinLabel(7, "K*");

  TH1D* hEnergy = new TH1D("hEnergy", "Energy distribution", 600, 0, 6);
  TH1D* hTheta  = new TH1D("hTheta", "Theta distribution", 314, 0, M_PI);
  TH1D* hPhi    = new TH1D("hPhi", "Phi distribution", 628, 0, M_PI * 2);
  TH1D* hPout   = new TH1D("hPout", "Momentum distribution", 600, 0, 6);
  TH1D* hPtrasv =
      new TH1D("hPtrasv", "Transverse momentum distribution", 600, 0, 6);

  TRandom3 rand;
  rand.SetSeed(0); // set seed for random number generator non so se ho
                   // incluso la classe giusta

  createInstances(n_event, n_particles_event, &rand, event_particles,
                  hMassKStarDecay, hInvariantMass, hMassOppositeSign,
                  hMassSameSign, hMassPionKaonOpposite, hMassPionKaonSame,
                  hType, hEnergy, hTheta, hPhi, hPout, hPtrasv);

  // Collect all histograms into a vector
  std::array<TH1D*, 12> histograms = {hType,
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

  // Save histograms to a file
  saveHistograms(histograms, "IstogrammiParticelle.root");

  std::cout << "Drawing histograms..." << std::endl;

  for (int i = 0; i < 12; ++i) {
    canvas->cd(i + 1);
    histograms[i]->Draw("APE");
  }

  canvas->Update(); // Aggiorna il canvas per visualizzare i grafici

  // Start the ROOT application event loop, but check periodically if the canvas
  // is closed
  while (gROOT->GetListOfCanvases()->FindObject("c1")) {
    gSystem->ProcessEvents(); // Process any events (including canvas events)
    gSystem->Sleep(100);      // Add a small delay to prevent 100% CPU usage
  }

  std::cout << "Canvas closed. Exiting program." << std::endl;

  return 0;
}
// g++ main.cpp $(root-config --cflags --libs) -o main
