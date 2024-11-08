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
    
    
    TH1D* h_mass_invariant = new TH1D("h_mass_invariant", "Mass Invariant distribution", 100, 0, 15);
    TH1D* h_type = new TH1D("h_type", "Particle type distribution", 7, 0, 7);
    TH1D* h_charge = new TH1D("h_charge", "Charge distribution", 3, -1, 2);
    TH1D* h_width = new TH1D("h_width", "Width distribution", 100, 0, 0.1);
    TH1D* h_energy = new TH1D("energy", "Energy distribution", 100, 0, 2);
    TH1D* h_theta = new TH1D("theta", "Theta distribution", 100, -M_PI/2, M_PI/2);
    TH1D* h_phi = new TH1D("phi", "Phi distribution", 100, 0, 2*M_PI);
    TH1D* h_pout = new TH1D("pout", "Pout distribution", 100, 0, 1);
    TH1D* h_Ptrasv = new TH1D("Ptrasve", "Ptrasv distribution", 100, 0, 1);




    for (int i=0; i <= n_event; ++i)
    {
           int particle_count=0;
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
                new_particle.decay2body(event_particles[particle_count], event_particles[particle_count + 1]);
                particle_count += 2;
            } else {
                event_particles[particle_count] = new_particle;
                ++particle_count;
            }
        }
        double mass_invariant = 0;
        for (int j=0; j <= particle_count; ++j)
        {
            mass_invariant += event_particles[j].get_mass();
        }
    }

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