#include "particle.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>

#include <iostream>

void randomParticlePosition(TRandom3* rand, Particle& particle) {
  double phi   = rand->Uniform(0, 2 * M_PI);
  double theta = rand->Uniform(0, M_PI);
  double pout  = rand->Exp(1);
  double px    = pout * std::sin(theta) * std::cos(phi);
  double py    = pout * std::sin(theta) * std::sin(phi);
  double pz    = pout * std::cos(theta);

  particle.set_p(px, py, pz);
}

Particle createRandomParticle(TRandom3* rand) {
  char name;
  int  c = rand->Integer(200);
  if (c < 80) {
    name = 'Q';
  } else if (c < 160) {
    name = 'q';
  } else if (c < 170) {
    name = 'K';
  } else if (c < 180) {
    name = 'k';
  } else if (c < 189) {
    name = 'P';
  } else if (c < 198) {
    name = 'p';
  } else {
    name = '*';
  }

  Particle particle(name);
  randomParticlePosition(rand, particle);
  return particle;
}

void fillHistogram(const std::array<Particle, 120>& particle,
                   TH1D* hInvariantMass, TH1D* hMassOppositeSign,
                   TH1D* hMassSameSign, TH1D* hMassPionKaonOpposite,
                   TH1D* hMassPionKaonSame, TH1D* hType, TH1D* hEnergy,
                   TH1D* hTheta, TH1D* hPhi, TH1D* hPout, TH1D* hPtrasv,
                   int n_particles) {
  for (int i = 0; i < n_particles; ++i) {
    // altri istogtammi credo
    double x     = particle[i].get_px();
    double y     = particle[i].get_py();
    double z     = particle[i].get_pz();
    double r     = std::sqrt(x * x + y * y + z * z);
    double theta = std::acos(z / r);
    double phi   = std::atan2(y, x);
    if (phi < 0) phi = phi + 2 * M_PI;

    hType->Fill(particle[i].get_index());
    hEnergy->Fill(particle[i].get_energy());

    hTheta->Fill(theta);
    hPhi->Fill(phi);
    hPout->Fill(r);
    hPtrasv->Fill(std::sqrt(x * x + y * y));

    for (int j = i + 1; j < n_particles; ++j) {
      double mass_inv = particle[i].invMass(particle[j]);

      hInvariantMass->Fill(mass_inv);

      // calcolo la massa invariante tra le particelle i e j di base e la metto
      // negli istogrammi solo se le parti soddisfano le cond
      //  Istogrammi basati su carica discordante (cariche opposte)
      if (particle[i].get_charge() * particle[j].get_charge() < 0) {
        hMassOppositeSign->Fill(mass_inv);
      }

      // Istogrammi basati su carica concordante (cariche uguali)
      else if (particle[i].get_charge() * particle[j].get_charge() > 0) {
        hMassSameSign->Fill(mass_inv);
      }

      // FATTI DA MASSIMO, DUNQUE POCO AFFIDABILI(però dai su insomma si spera
      // che siano fatti decentemente non doveva fare molto)

      // Massa invariante tra particelle di tipo Pion+/Kaon- e Pion-/Kaon+
      if ((particle[i].get_name() == 'Q' && particle[j].get_name() == 'k')
          || (particle[i].get_name() == 'q' && particle[j].get_name() == 'K')
          || (particle[i].get_name() == 'K' && particle[j].get_name() == 'q')
          || (particle[i].get_name() == 'k' && particle[j].get_name() == 'Q')) {
        hMassPionKaonOpposite->Fill(mass_inv);
      }

      // Massa invariante tra particelle di tipo Pion+/Kaon+ e Pion-/Kaon-
      else if ((particle[i].get_name() == 'Q' && particle[j].get_name() == 'K')
               || (particle[i].get_name() == 'q'
                   && particle[j].get_name() == 'k')
               || (particle[i].get_name() == 'K'
                   && particle[j].get_name() == 'Q')
               || (particle[i].get_name() == 'k'
                   && particle[j].get_name() == 'q')) {
        hMassPionKaonSame->Fill(mass_inv);
      }
    }
  }
}

void createInstances(int n_event, int n_particles_event, TRandom3* rand,
                     std::array<Particle, 120>& event_particles,
                     TH1D* hMassKStarDecay, TH1D* hInvariantMass,
                     TH1D* hMassOppositeSign, TH1D* hMassSameSign,
                     TH1D* hMassPionKaonOpposite, TH1D* hMassPionKaonSame,
                     TH1D* hType, TH1D* hEnergy, TH1D* hTheta, TH1D* hPhi,
                     TH1D* hPout, TH1D* hPtrasv) {
  std::cout << "Generating particles..." << std::endl;
  for (int i = 0; i < n_event; ++i) {
    int particle_count = 0;

    for (int j = 0; j < n_particles_event; ++j) {
      Particle new_particle = createRandomParticle(rand);

      if (new_particle.get_name() == '*') {
        char dau1_name = 0;
        char dau2_name = 0;

        switch (rand->Integer(2)) {
        case 0:
          dau1_name = 'Q';
          dau2_name = 'k';
          break;
        case 1:
          dau1_name = 'q';
          dau2_name = 'K';
          break;
        }

        event_particles[particle_count] = Particle(dau1_name);

        event_particles[particle_count + 1] = Particle(dau2_name);

        new_particle.decay2body(event_particles[particle_count],
                                event_particles[particle_count + 1]);

        hMassKStarDecay->Fill(event_particles[particle_count].invMass(
            event_particles[particle_count + 1]));

        particle_count += 2;
      } else {
        event_particles[particle_count] = new_particle;

        ++particle_count;
      }
    }

    // riempimento istogrammi
    fillHistogram(event_particles, hInvariantMass, hMassOppositeSign,
                  hMassSameSign, hMassPionKaonOpposite, hMassPionKaonSame,
                  hType, hEnergy, hTheta, hPhi, hPout, hPtrasv, particle_count);
  }
}

void saveHistograms(const std::array<TH1D*, 12>& histograms,
                    const std::string&           filename) {
  std::cout << "Saving data in " << filename << "..." << std::endl;
  // Create a new TFile object to write the histograms
  TFile file(filename.c_str(), "RECREATE");
  if (!file.IsOpen()) {
    std::cerr << "Error: Could not open file " << filename << " for writing."
              << std::endl;
    return;
  }

  // Loop over the histograms and write each to the file
  for (auto hist : histograms) {
    if (hist) { hist->Write(); }
  }

  // Close the file
  file.Close();
}

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
      new TH1D("hMassKStarDecay", "Invariant Mass - K* decay", 40, 0.7, 1.1);
  hMassKStarDecay->Sumw2();

  TH1D* hType = new TH1D("hType", "Particle type distribution", 7, 0, 7);
  hType->SetStats(0);
  hType->GetXaxis()->SetBinLabel(1, "#pi+");
  hType->GetXaxis()->SetBinLabel(2, "#pi-");
  hType->GetXaxis()->SetBinLabel(3, "P+");
  hType->GetXaxis()->SetBinLabel(4, "P-");
  hType->GetXaxis()->SetBinLabel(5, "K+");
  hType->GetXaxis()->SetBinLabel(6, "K-");
  hType->GetXaxis()->SetBinLabel(7, "K*");

  TH1D* hEnergy = new TH1D("hEnergy", "Energy distribution", 600, 0, 6);
  TH1D* hTheta  = new TH1D("hTheta", "Theta distribution", 314, 0, M_PI);
  TH1D* hPhi    = new TH1D("hPhi", "Phi distribution", 628, 0, 2 * M_PI);
  TH1D* hPout   = new TH1D("hPout", "Momentum distribution", 600, 0, 6);
  TH1D* hPtrasv =
      new TH1D("hPtrasv", "Transverse momentum distribution", 600, 0, 6);

  TRandom3 rand;
  rand.SetSeed(0);

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

  std::cout << "Drawing histograms..." << std::endl;

  canvas->cd(1);
  hType->Draw("E");
  for (int i = 1; i < 12; ++i) {
    canvas->cd(i + 1);
    histograms[i]->Draw("P");
  }

  saveHistograms(histograms, "IstogrammiParticelle.root");

  canvas->Update(); // Aggiorna il canvas per visualizzare i grafici

  while (gROOT->GetListOfCanvases()->FindObject("c1")) {
    gSystem->ProcessEvents(); // Process any events (including canvas events)
    gSystem->Sleep(100);      // Add a small delay to prevent 100% CPU usage
  }

  std::cout << "Canvas closed. Exiting program." << std::endl;

  return 0;
}
