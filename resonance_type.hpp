#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP

#include "particle_type.hpp"

#include <iomanip>
#include <iostream>

class ResonanceType : public ParticleType {
 public:
  double get_width() const override;

  void print() const override;

  bool operator==(const ResonanceType& other) const;
  bool operator!=(const ResonanceType& other) const;

  ResonanceType(char name, double mass, int charge, double width);

 private:
  const double width_;
};

inline double ResonanceType::get_width() const { return width_; }

inline void ResonanceType::print() const {
  std::cout << "Particle " << get_name() << ": mass = " << std::setw(6)
            << get_mass() << ", charge = " << std::setw(6) << get_charge()
            << ", width = " << std::setw(6) << get_width() << '\n';
}

inline ResonanceType::ResonanceType(char name, double mass, int charge,
                                    double width)
    : ParticleType(name, mass, charge)
    , width_(width) {}

inline bool ResonanceType::operator==(const ResonanceType& other) const {
  return get_charge() == other.get_charge() && get_mass() == other.get_mass()
      && get_width() == other.get_width();
}

inline bool ResonanceType::operator!=(const ResonanceType& other) const {
  return !(*this == other);
}
#endif