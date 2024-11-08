#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP

#include <iomanip>
#include <iostream>

class ParticleType {
 public:
  char   get_name() const;
  double get_mass() const;
  int    get_charge() const;
  virtual double get_width() const { return 0; }

  virtual void print() const;

  bool operator==(const ParticleType& other) const;
  bool operator!=(const ParticleType& other) const;

  ParticleType(char name, double mass, int charge);

 private:
  const char   name_;
  const double mass_;
  const int    charge_;
};

inline char   ParticleType::get_name() const { return name_; }
inline double ParticleType::get_mass() const { return mass_; }
inline int    ParticleType::get_charge() const { return charge_; }

inline void ParticleType::print() const {
  std::cout << "Particle " << get_name() << ": mass = " << std::setw(6)
            << get_mass() << ", charge = " << std::setw(6) << get_charge()
            << '\n';
}

inline ParticleType::ParticleType(char name, double mass, int charge)
    : name_(name)
    , mass_(mass)
    , charge_(charge) {}

inline bool ParticleType::operator==(const ParticleType& other) const {
  return get_charge() == other.get_charge() && get_mass() == other.get_mass();
}

inline bool ParticleType::operator!=(const ParticleType& other) const {
  return !(*this == other);
}

#endif