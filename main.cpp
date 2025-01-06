
#include "TH1D.h"
#include "TRandom3.h"

#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"

#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TStyle.h>
// #include <TSystem.h>

#include <cmath>
#include <iostream>
#include <random>

void     randomParticlePosition(TRandom3* rand, Particle& particle);
Particle createRandomParticle(TRandom3* rand);
void     fillInstogram(const Particle (&particle)[120], TH1D* hMassOppositeSign,
                       TH1D* hMassSameSign, TH1D* hMassPionKaonOpposite,
                       TH1D* hMassPionKaonSame, TH1D* hMassKStarDecay, TH1D* h_type,
                       TH1D* h_energy, TH1D* h_theta, TH1D* h_phi, TH1D* h_pout,
                       TH1D* h_Ptrasv, TFile* file, int n_particles);
void     saveHistograms(const std::vector<TH1D*>& histograms,
                        const std::string&        filename);

int main() {
  // Crea un'applicazione ROOT
  TApplication theApp("App", 0, 0);

  std::cout << "Creazione del canvas..." << std::endl;
  TCanvas* canvas = new TCanvas("c1", "Canvas per Istogramma", 800, 600);
  canvas->Divide(4, 3); // Dividi il canvas in 2 righe e 2 colonne
  if (canvas == nullptr) {
    std::cerr << "Errore nella creazione del canvas!" << std::endl;
    return -1;
  }

  TRandom3 rand;
  rand.SetSeed(0); // set seed for random number generator non so se ho
                   // incluso la classe giusta
                   // name, mass, charge, whidth
  Particle::addParticleType('Q', 0.13957, 1);
  Particle::addParticleType('q', 0.13957, -1);
  Particle::addParticleType('P', 0.93827, 1);
  Particle::addParticleType('p', 0.93827, -1);
  Particle::addParticleType('K', 0.49367, 1);
  Particle::addParticleType('k', 0.49367, -1);
  Particle::addParticleType('*', 0.89166, 0, 0.050);

  int      n_event           = 1000;
  int      n_particles_event = 100;
  Particle event_particles[120];

  // MASSIMO HA AGGIUNTO SUMW E TFILE  "si dovrebbe mettere così, non so cosa
  // serva, però nel caso aggiustate"

  TFile* file = new TFile("IstogrammiParticelle.root", "RECREATE");

  // secondo Massimo "si dovrebbero cambiare i nomi dei grafici con dei
  // numeri così alla fine si può usare un ciclo for per scrivere i grafici
  // a fine codice e per il metodo sumw :)"

  // range tra 0 e 3, il numero dei bin è , non 60 se no è troopo largo,
  // 180.

  TH1D* hMassInvariant =
      new TH1D("hMassInvariant", "Mass Invariant distribution", 100, 0, 100001);
  hMassInvariant->Sumw2();
  TH1D* hMassOppositeSign = new TH1D(
      "hMassOppositeSign", "Massa Invariante - Carica Discorde", 300, 0, 6);
  hMassOppositeSign->Sumw2();
  TH1D* hMassSameSign = new TH1D(
      "hMassSameSign", "Massa Invariante - Carica Concorde", 300, 0, 6);
  hMassSameSign->Sumw2();
  TH1D* hMassPionKaonOpposite =
      new TH1D("hMassPionKaonOpposite",
               "Massa Invariante - Pione+/Kaone- e Pione-/Kaone+", 300, 0, 6);
  hMassPionKaonOpposite->Sumw2();
  TH1D* hMassPionKaonSame =
      new TH1D("hMassPionKaonSame",
               "Massa Invariante - Pione+/Kaone+ e Pione-/Kaone-", 300, 0, 6);
  hMassPionKaonSame->Sumw2();
  TH1D* hMassKStarDecay =
      new TH1D("hMassKStarDecay", "Massa Invariante - Decadimento Risonanza K*",
               100, 0.7, 1.05);
  hMassKStarDecay->Sumw2();

  TH1D* h_type   = new TH1D("h_type", "Particle type distribution", 7, 0, 7);
  TH1D* h_energy = new TH1D("energy", "Energy distribution", 300, 0, 6);
  TH1D* h_theta  = new TH1D("theta", "Theta distribution", 100, 0, M_PI);
  TH1D* h_phi    = new TH1D("phi", "Phi distribution", 100, 0, M_PI * 2);
  TH1D* h_pout   = new TH1D("pout", "Pout distribution", 200, 0, 6);
  TH1D* h_Ptrasv = new TH1D("Ptrasve", "Ptrasv distribution", 200, 0, 6);

  double total_energy = 0;
  double total_px     = 0;
  double total_py     = 0;
  double total_pz     = 0;

  std::cout << "Generazione delle istanze..." << std::endl;
  for (int i = 0; i < n_event; ++i) {
    int    particle_count = 0;
    double mass_invariant = 0;

    for (int j = 0; j < n_particles_event; ++j) {
      Particle new_particle = createRandomParticle(&rand);

      if (new_particle.get_name() == '*') {
        char dau1_name = 0;
        char dau2_name = 0;

        switch (rand.Integer(2)) {
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

        total_energy += (event_particles[particle_count].get_energy()
                         + event_particles[particle_count + 1].get_energy());
        total_px     += (event_particles[particle_count].get_px()
                     + event_particles[particle_count + 1].get_px());
        total_py     += (event_particles[particle_count].get_py()
                     + event_particles[particle_count + 1].get_py());
        total_pz     += (event_particles[particle_count].get_pz()
                     + event_particles[particle_count + 1].get_pz());

        hMassKStarDecay->Fill(event_particles[particle_count].invMass(
            event_particles[particle_count + 1]));

        particle_count += 2;
      } else {
        event_particles[particle_count] = new_particle;

        total_energy += new_particle.get_energy();
        total_px     += new_particle.get_px();
        total_py     += new_particle.get_py();
        total_pz     += new_particle.get_pz();

        ++particle_count;
      }
    }

    mass_invariant = std::sqrt(total_energy * total_energy - total_px * total_px
                               - total_py * total_py - total_pz * total_pz);
    hMassInvariant->Fill(mass_invariant);
    // riempimento istogrammi
    fillInstogram(event_particles, hMassOppositeSign, hMassSameSign,
                  hMassPionKaonOpposite, hMassPionKaonSame, hMassKStarDecay,
                  h_type, h_energy, h_theta, h_phi, h_pout, h_Ptrasv, file,
                  n_particles_event);
  }

  // Collect all histograms into a vector
  std::vector<TH1D*> histograms = {hMassInvariant,
                                   hMassOppositeSign,
                                   hMassSameSign,
                                   hMassPionKaonOpposite,
                                   hMassPionKaonSame,
                                   hMassKStarDecay,
                                   h_type,
                                   h_energy,
                                   h_theta,
                                   h_phi,
                                   h_pout,
                                   h_Ptrasv};

  // Save histograms to a file
  saveHistograms(histograms, "IstogrammiParticelle.root");

  double mass_invariant =
      std::sqrt(total_energy * total_energy - total_px * total_px
                - total_py * total_py - total_pz * total_pz);

  std::cout << "Disegno dei grafici..." << std::endl;

  // Disegna su pannello 1: hMassInvariant
  canvas->cd(1); // Seleziona il pannello 1
  hMassInvariant->Draw("APE");

  // Disegna su pannello 2: hMassOppositeSign
  canvas->cd(2); // Seleziona il pannello 2
  hMassOppositeSign->Draw("APE");

  // Disegna su pannello 3: hMassSameSign
  canvas->cd(3); // Seleziona il pannello 3
  hMassSameSign->Draw("APE");

  // Disegna su pannello 4: hMassPionKaonOpposite
  canvas->cd(4); // Seleziona il pannello 4
  hMassPionKaonOpposite->Draw("APE");

  // Disegna su pannello 5: hMassPionKaonSame
  canvas->cd(5); // Seleziona il pannello 5
  hMassPionKaonSame->Draw("APE");

  // Disegna su pannello 6: hMassKStarDecay
  canvas->cd(6); // Seleziona il pannello 6
  hMassKStarDecay->Draw("APE");

  // Disegna su pannello 7: h_type
  canvas->cd(7); // Seleziona il pannello 7
  h_type->Draw("APE");

  // Disegna su pannello 8: h_energy
  canvas->cd(8); // Seleziona il pannello 8
  h_energy->Draw("APE");

  // Disegna su pannello 9: h_theta
  canvas->cd(9); // Seleziona il pannello 9
  h_theta->Draw("APE");

  // Disegna su pannello 10: h_phi
  canvas->cd(10); // Seleziona il pannello 10
  h_phi->Draw("APE");

  // Disegna su pannello 11: h_pout
  canvas->cd(11); // Seleziona il pannello 11
  h_pout->Draw("APE");

  // Disegna su pannello 12: h_Ptrasv
  canvas->cd(12); // Seleziona il pannello 12
  h_Ptrasv->Draw("APE");

  canvas->Update(); // Aggiorna il canvas per visualizzare i grafici

  // Mantieni la finestra ROOT aperta
  theApp.Run();

  gApplication->Run();
  return 0;
}

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

void fillInstogram(const Particle (&particle)[120], TH1D* hMassOppositeSign,
                   TH1D* hMassSameSign, TH1D* hMassPionKaonOpposite,
                   TH1D* hMassPionKaonSame, TH1D* hMassKStarDecay, TH1D* h_type,
                   TH1D* h_energy, TH1D* h_theta, TH1D* h_phi, TH1D* h_pout,
                   TH1D* h_Ptrasv, TFile* file, int n_particles) {
  for (int i = 0; i < n_particles; ++i) {
    // altri istogtammi credo
    double x = particle[i].get_px();
    double y = particle[i].get_py();
    double z = particle[i].get_pz();

    h_type->Fill(particle[i].get_index());
    h_energy->Fill(particle[i].get_energy());
    h_theta->Fill(std::acos(z / std::sqrt(x * x + y * y + z * z)));

    double phi = std::atan2(y, x);
    if (phi < 0) {
      phi += 2 * M_PI; // Convert the range from [-pi, pi] to [0, 2*pi]
    }

    h_phi->Fill(phi);
    h_pout->Fill(std::sqrt(particle[i].get_px() * particle[i].get_px()
                           + particle[i].get_py() * particle[i].get_py()
                           + particle[i].get_pz() * particle[i].get_pz()));
    h_Ptrasv->Fill(std::sqrt(particle[i].get_px() * particle[i].get_px()
                             + particle[i].get_py() * particle[i].get_py()));

    for (int j = i + 1; j < n_particles; ++j) {
      // Calcola la massa invariante tra particella i e particella
      // j(assoulutamente inutile)
      double px_sum   = particle[i].get_px() + particle[j].get_px();
      double py_sum   = particle[i].get_py() + particle[j].get_py();
      double pz_sum   = particle[i].get_pz() + particle[j].get_pz();
      double E_sum    = particle[i].get_energy() + particle[j].get_energy();
      double mass_inv = std::sqrt(E_sum * E_sum - px_sum * px_sum
                                  - py_sum * py_sum - pz_sum * pz_sum);

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
          || (particle[i].get_name() == 'q' && particle[j].get_name() == 'K')) {
        hMassPionKaonOpposite->Fill(mass_inv);
      }

      // Massa invariante tra particelle di tipo Pion+/Kaon+ e Pion-/Kaon-
      else if ((particle[i].get_name() == 'Q' && particle[j].get_name() == 'K')
               || (particle[i].get_name() == 'q'
                   && particle[j].get_name() == 'k')) {
        hMassPionKaonSame->Fill(mass_inv);
      }

      /* massa invariante fra le particelle generate in ogni evento che derivano
      dal decadimento della risonanza K* (quelle in “coda”, per intenderci. N.B.
      considerate esclusivamente coppie di particelle figlie che provengono
      dalla stessa “madre”, i.e. non mischiare “figlie” di “madri” diverse,
      poiché solo nel primo caso si osserverà in questo istogramma il
      caratteristico picco di massa invariante). Questo istogramma fa da
      istogramma di benchmark affinché si possa essere certi che la gestione del
      decadimento, le formule di massa invariante/energia e la virtualizzazione
      delle classi ParticleType/ResonanceType siano state correttamente
      implementate. Se tutto è ok in questa fase, dovreste osservare una
      gaussiana con media la massa della K*, e RMS la larghezza della K* come
      impostata all’inizio del programma.    */
    }
  }
  hMassOppositeSign->Write();
  hMassSameSign->Write();
  hMassPionKaonOpposite->Write();
  hMassPionKaonSame->Write();
  hMassKStarDecay->Write();
  // file->Close();
}

void saveHistograms(const std::vector<TH1D*>& histograms,
                    const std::string&        filename) {
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

// g++ main.cpp $(root-config --cflags --libs) -o main
