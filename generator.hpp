#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include "particle.hpp"

#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>

void createInstances(int n_event, int n_particles_event, TRandom3* rand,
                     std::array<Particle, 120>& event_particles,
                     double& total_energy, double& total_px, double& total_py,
                     double& total_pz, TH1D* hMassKStarDecay,
                     TH1D* hMassInvariant, TH1D* hMassOppositeSign,
                     TH1D* hMassSameSign, TH1D* hMassPionKaonOpposite,
                     TH1D* hMassPionKaonSame, TH1D* hType, TH1D* hEnergy,
                     TH1D* hTheta, TH1D* hPhi, TH1D* hPout, TH1D* hPtrasv) {
  std::cout << "Generazione delle istanze..." << std::endl;
  for (int i = 0; i < n_event; ++i) {
    int    particle_count = 0;
    double mass_invariant = 0;

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
    fillHistogram(event_particles, hMassOppositeSign, hMassSameSign,
                  hMassPionKaonOpposite, hMassPionKaonSame, hType, hEnergy,
                  hTheta, hPhi, hPout, hPtrasv, particle_count);
  }
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

void fillHistogram(const std::array<Particle, 120>& particle,
                   TH1D* hMassOppositeSign, TH1D* hMassSameSign,
                   TH1D* hMassPionKaonOpposite, TH1D* hMassPionKaonSame,
                   TH1D* hType, TH1D* hEnergy, TH1D* hTheta, TH1D* hPhi,
                   TH1D* hPout, TH1D* hPtrasv, int n_particles) {
  for (int i = 0; i < n_particles; ++i) {
    // altri istogtammi credo
    double x = particle[i].get_px();
    double y = particle[i].get_py();
    double z = particle[i].get_pz();

    hType->Fill(particle[i].get_index());
    hEnergy->Fill(particle[i].get_energy());
    hTheta->Fill(std::acos(z / std::sqrt(x * x + y * y + z * z)));

    double phi = std::atan2(y, x);
    if (phi < 0) {
      phi += 2 * M_PI; // Convert the range from [-pi, pi] to [0, 2*pi]
    }

    hPhi->Fill(phi);
    hPout->Fill(std::sqrt(particle[i].get_px() * particle[i].get_px()
                          + particle[i].get_py() * particle[i].get_py()
                          + particle[i].get_pz() * particle[i].get_pz()));
    hPtrasv->Fill(std::sqrt(particle[i].get_px() * particle[i].get_px()
                            + particle[i].get_py() * particle[i].get_py()));

    for (int j = i + 1; j < n_particles; ++j) {
      double mass_inv = particle[i].invMass(particle[j]);

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
}

void saveHistograms(const std::array<TH1D*, 12>& histograms,
                    const std::string&           filename) {
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

#endif