#include <iostream>
#include <random>
#include <cmath>
#include <TRandom3.h>

#include "TH1D.h"
#include "TRandom3.h"

#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"



void randomParticlePosition(Particle& particle);
Particle createRandomParticle();

int main(){
    gRandom->SetSeed(); // set seed for random number generator non so se ho incluso la classe giusta
    
    Particle::addParticleType('P', 0.93827, 1, 0);
    Particle::addParticleType('p', 0.93827, -1, 0);
    Particle::addParticleType('K', 0.49367, 1, 0);
    Particle::addParticleType('k', 0.49367, -1, 0);
    Particle::addParticleType('Q', 0.13957, 1, 0);
    Particle::addParticleType('q', 0.13957, -1, 0);
    Particle::addParticleType('*', 0.89166, 0, 0.050);

    int n_event = 1E5;
    int n_particles_event = 100;
    Particle event_particles[n_particles_event+20];
    

    //MASSIMO HA AGGIUNTO SUMW E TFILE  "si dovrebbe mettere così, non so cosa serva, però nel caso aggiustate" 

    TFile*file = Tfile::open("IstogrammiParticelle.root","RECREATE") 

    //secondo Massimo "si dovrebbero cambiare i nomi dei grafici con dei numeri così alla fine si può usare un ciclo for per scrivere i grafici a fine codice e per il metodo sumw :)"

    //range tra 0 e 3, il numero dei bin è , non 60 se no è troopo largo, 180. 

    TH1D* h_mass_invariant = new TH1D("h_mass_invariant", "Mass Invariant distribution", 100, 0, 15);
    h_mass_invariant->Sumw2();
    TH1D* hMassOppositeSign = new TH1D("hMassOppositeSign", "Massa Invariante - Carica Discorre", 100, 0, 2);
    hMassOppositeSign->Sumw2();
    TH1D* hMassSameSign = new TH1D("hMassSameSign", "Massa Invariante - Carica Concorde", 100, 0, 2);
    hMassSameSign->Sumw2();
    TH1D* hMassPionKaonOpposite = new TH1D("hMassPionKaonOpposite", "Massa Invariante - Pione+/Kaone- e Pione-/Kaone+", 100, 0, 2);
    hMassPionKaonOpposite->Sumw2();
    TH1D* hMassPionKaonSame = new TH1D("hMassPionKaonSame", "Massa Invariante - Pione+/Kaone+ e Pione-/Kaone-", 100, 0, 2);
    hMassPionKaonSame->Sumw2();
    TH1D* hMassKStarDecay = new TH1D("hMassKStarDecay", "Massa Invariante - Decadimento Risonanza K*", 100, 0.7, 1.0);
    hMassKStarDecay->Sumw2();

    TH1D* h_type = new TH1D("h_type", "Particle type distribution", 7, 0, 7);
    TH1D* h_energy = new TH1D("energy", "Energy distribution", 100, 0, 2);
    TH1D* h_theta = new TH1D("theta", "Theta distribution", 100, -M_PI/2, M_PI/2);
    TH1D* h_phi = new TH1D("phi", "Phi distribution", 100, 0, 2*M_PI);
    TH1D* h_pout = new TH1D("pout", "Pout distribution", 100, 0, 1);
    TH1D* h_Ptrasv = new TH1D("Ptrasve", "Ptrasv distribution", 100, 0, 1);

    double total_energy = 0;
    double total_px = 0;
    double total_py = 0;
    double total_pz = 0;

    for (int i=0; i <= n_event; ++i)
    {
           int particle_count=0;
           double mass_invariant = 0;
        
        for(int j=0; j <= n_particles_event; ++j)
        {
            Particle new_particle = createRandomParticle();
           
            if (new_particle.get_name() == '*')
            {
                char dau1_name = 0;
                char dau2_name = 0;
                TRandom3 rand;
                switch (    rand.Integer(2))
                {
                case 0:
                dau1_name = 'Q';
                dau2_name = 'k';
                break;
                case 1:
                dau1_name = 'q';
                dau2_name = 'K';
                break;}
                
                event_particles[particle_count] = Particle(dau1_name);
                
                event_particles[particle_count + 1] = Particle(dau2_name);                
                
                new_particle.decay2body(event_particles[particle_count + 1], event_particles[particle_count + 1]);
               
                total_energy += (event_particles[particle_count].get_energy() + event_particles[particle_count + 1].get_energy());
                total_px += (event_particles[particle_count].get_px() + event_particles[particle_count + 1].get_px());
                total_py += (event_particles[particle_count].get_py() + event_particles[particle_count + 1].get_py());
                total_pz += (event_particles[particle_count].get_pz() + event_particles[particle_count + 1].get_pz());
                
                hMassKStarDecay->Fill(event_particles[particle_count].invMass(event_particles[particle_count + 1]));

                particle_count += 2;
            } else {

                event_particles[particle_count] = new_particle;

                total_energy += new_particle.get_energy();
                total_px += new_particle.get_px();  
                total_py += new_particle.get_py();
                total_pz += new_particle.get_pz();

                ++particle_count;
            }
        }

        mass_invariant = std::sqrt(total_energy*total_energy - total_px*total_px - total_py*total_py - total_pz*total_pz);
        h_mass_invariant->Fill(mass_invariant);
        //riempimento istogrammi



    }

    double mass_invariant = std::sqrt(total_energy*total_energy - total_px*total_px - total_py*total_py - total_pz*total_pz);

    return 0;
}

void randomParticlePosition(Particle& particle)
{
    TRandom3 rand;
    double phi =     rand.Uniform(0, 2*M_PI);
    double theta =     rand.Uniform(-M_PI/2, M_PI/2);
    double pout =     rand.Exp(1);
    double px = pout * std::sin(theta) * std::cos(phi);
    double py = pout * std::sin(theta) * std::sin(phi);
    double pz = pout * std::cos(theta);

   particle.set_p(px, py, pz);
}

Particle createRandomParticle()
{
    char name;
    TRandom3 rand;
    int c =     rand.Integer(200);
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
    randomParticlePosition(particle);
    return particle;
}

void fillInstogram(Particle particle[120], TH1D* hMassOppositeSign, TH1D* hMassSameSign, TH1D* hMassPionKaonOpposite, TH1D* hMassPionKaonSame, TH1D* hMassKStarDecay, TH1D* h_type, TH1D* h_energy, TH1D* h_theta, TH1D* h_phi, TH1D* h_pout, TH1D* h_Ptrasv)
{




    for (int i = 0; i < n_particles; ++i) {
        //altri istogtammi credo
        h_type->Fill(particle[i].get_name());
        h_energy->Fill(particle[i].get_energy());
        h_theta->Fill(std::acos(particle[i].get_pz()/std::sqrt(particle[i].get_px()*particle[i].get_px() + particle[i].get_py()*particle[i].get_py() + particle[i].get_pz()*particle[i].get_pz())));
        h_phi->Fill(std::atan(particle[i].get_py()/particle[i].get_px()));
        h_pout->Fill(std::sqrt(particle[i].get_px()*particle[i].get_px() + particle[i].get_py()*particle[i].get_py() + particle[i].get_pz()*particle[i].get_pz()));
        h_Ptrasv->Fill(std::sqrt(particle[i].get_px()*particle[i].get_px() + particle[i].get_py()*particle[i].get_py()));
    
            for (int j = i + 1; j < n_particles; ++j) {
                // Calcola la massa invariante tra particella i e particella j(assoulutamente inutile)
                double px_sum = particles[i].get_px() + particles[j].get_px();
                double py_sum = particles[i].get_py() + particles[j].get_py();
                double pz_sum = particles[i].get_pz() + particles[j].get_pz();
                double E_sum = particles[i].get_energy() + particles[j].get_energy();
                double mass_inv = std::sqrt(E_sum * E_sum - px_sum * px_sum - py_sum * py_sum - pz_sum * pz_sum);

//calcolo la massa invariante tra le particelle i e j di base e la metto negli istogrammi solo se le parti soddisfano le cond
                // Istogrammi basati su carica discordante (cariche opposte)
                if (particles[i].get_charge() * particles[j].get_charge() < 0) {
                    hMassOppositeSign->Fill(mass_inv);  
                }

                // Istogrammi basati su carica concordante (cariche uguali)
                else if (particles[i].get_charge() * particles[j].get_charge() > 0) {
                    hMassSameSign->Fill(mass_inv);
                }
                
                //FATTI DA MASSIMO, DUNQUE POCO AFFIDABILI(però dai su insomma si spera che siano fatti decentemente non doveva fare molto)
                
                // Massa invariante tra particelle di tipo Pion+/Kaon- e Pion-/Kaon+
                else if ((particles[i].type == 0 && particles[j].type == 3) || (particles[i].type == 1 && particles[j].type == 2)) {
                    hMassPionKaonOpposite->Fill(mass_inv);
                }

                // Massa invariante tra particelle di tipo Pion+/Kaon+ e Pion-/Kaon-
                else if ((particles[i].type == 0 && particles[j].type == 2) || (particles[i].type == 1 && particles[j].type == 3)) {
                   hMassPionKaonSame->Fill(mass_inv); 
                }
                
                


                /* massa invariante fra le particelle generate in ogni evento che derivano dal decadimento della
                risonanza K* (quelle in “coda”, per intenderci. N.B. considerate esclusivamente coppie di particelle
                figlie che provengono dalla stessa “madre”, i.e. non mischiare “figlie” di “madri” diverse, poiché
                solo nel primo caso si osserverà in questo istogramma il caratteristico picco di massa invariante).
                Questo istogramma fa da istogramma di benchmark affinché si possa essere certi che la gestione del
                decadimento, le formule di massa invariante/energia e la virtualizzazione delle classi
                ParticleType/ResonanceType siano state correttamente implementate. Se tutto è ok in questa fase,
                dovreste osservare una gaussiana con media la massa della K*, e RMS la larghezza della K*
                come impostata all’inizio del programma.    */
            }
        } 
    h_mass_invariant->write();
    hMassOppositeSign->write();
    hMassSameSign->write();
    hMassPionKaonOpposite->write();
    hMassPionKaonSame->write();
    hMassKStarDecay->write();
    file->close();

}
