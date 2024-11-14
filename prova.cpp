#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"
#include "particle.hpp"

int main() {
  ParticleType  particle('P', 5.68, -2);
  ResonanceType resonance('K', 17.3, 0, 6.7);

  std::array<ParticleType*, 2> particle_types;

  particle_types[0] = &particle;
  particle_types[1] = &resonance;

  for (int i = 0; i < 2; ++i) particle_types[i]->print();

  Particle::addParticleType(particle.get_name(), particle.get_mass(), particle.get_charge());

  Particle new_particle(particle.get_name());

  new_particle.printParticle();
}