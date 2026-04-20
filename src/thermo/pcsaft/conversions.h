/**
 * @file conversions.h
 * @brief Conversion utilities between legacy pcsaft:: and DMThermo::Core types.
 */

#ifndef THERMO_LEGACY_PCSAFT_CONVERSIONS_H
#define THERMO_LEGACY_PCSAFT_CONVERSIONS_H

#include "thermo/core/component.h"
#include "thermo/core/mixture.h"

namespace pcsaft {
class Component;
class Mixture;
struct AssociationParams;
enum class AssociationScheme;
} // namespace pcsaft

namespace DMThermo {
namespace Core {

AssociationScheme convertScheme(pcsaft::AssociationScheme legacy);
pcsaft::AssociationScheme convertScheme(AssociationScheme scheme);

Component fromLegacy(const pcsaft::Component& legacy);
pcsaft::Component toLegacy(const Component& component);

Mixture fromLegacy(const pcsaft::Mixture& legacy);
pcsaft::Mixture toLegacy(const Mixture& mixture, const std::vector<double>& composition);

Component getComponentByName(const std::string& name);
Mixture createMixture(const std::vector<std::string>& names);
Mixture createMixture(const std::vector<std::string>& names, const BinaryParameters& kij);

} // namespace Core
} // namespace DMThermo

#endif // THERMO_LEGACY_PCSAFT_CONVERSIONS_H
