/**
 * @file eos_factory.h
 * @brief Factory for creating equation of state objects
 */

#ifndef THERMO_FACTORY_EOS_FACTORY_H
#define THERMO_FACTORY_EOS_FACTORY_H

#include <memory>
#include <string>
#include <functional>
#include <unordered_map>
#include "thermo/eos.h"
#include "thermo/core/mixture.h"
#include "thermo/data/databanks.h"
#include "thermo/data/diagnostics.h"
#include "mixture_factory.h"

namespace DMThermo {
namespace Factory {

/**
 * @brief Factory for creating EOS implementations
 *
 * Supports registration of custom EOS implementations and provides
 * convenience methods for built-in EOS types.
 *
 * Example usage:
 * @code
 * // Create PC-SAFT EOS for a mixture
 * auto eos = EOSFactory::createPCSaft(mixture);
 *
 * // Or use generic creation by name
 * auto eos = EOSFactory::create("PC-SAFT", mixture);
 *
 * // Register custom EOS
 * EOSFactory::registerEOS("CustomEOS", [](const Core::Mixture& mix) {
 *     return std::make_shared<CustomEOS>(mix);
 * });
 * @endcode
 */
class EOSFactory {
public:
    /// Creator function type
    using EOSCreator = std::function<EOSPtr(const Core::Mixture&)>;

    // =========================================================================
    // Generic Creation
    // =========================================================================

    /**
     * @brief Create EOS by name
     *
     * @param eos_type EOS type name (e.g., "PC-SAFT", "Peng-Robinson")
     * @param mixture Mixture definition
     * @return Shared pointer to EOS
     * @throws std::invalid_argument if EOS type not registered
     */
    static EOSPtr create(
        const std::string& eos_type,
        const Core::Mixture& mixture
    );

    /**
     * @brief Create EOS by name with mixture pointer
     */
    static EOSPtr create(
        const std::string& eos_type,
        Core::MixturePtr mixture
    );

    // =========================================================================
    // Specific EOS Creation (Convenience Methods)
    // =========================================================================

    /**
     * @brief Create PC-SAFT EOS from new Mixture
     *
     * @param mixture Mixture definition
     * @return PC-SAFT EOS instance
     */
    static EOSPtr createPCSaft(const Core::Mixture& mixture);
    static EOSPtr createPCSaft(Core::MixturePtr mixture);

    /**
     * @brief Create Peng-Robinson EOS
     *
     * @param mixture Mixture definition
     * @return Peng-Robinson EOS instance
     */
    static EOSPtr createPengRobinson(const Core::Mixture& mixture);
    static EOSPtr createPengRobinson(Core::MixturePtr mixture);

    /**
     * @brief Create Soave-Redlich-Kwong EOS
     *
     * @param mixture Mixture definition
     * @return SRK EOS instance
     */
    static EOSPtr createSRK(const Core::Mixture& mixture);
    static EOSPtr createSRK(Core::MixturePtr mixture);

    // =========================================================================
    // Registration
    // =========================================================================

    /**
     * @brief Register custom EOS creator
     *
     * @param name EOS type name
     * @param creator Function to create EOS from mixture
     */
    static void registerEOS(const std::string& name, EOSCreator creator);

    /**
     * @brief Check if EOS type is registered
     *
     * @param name EOS type name
     * @return true if registered
     */
    static bool isRegistered(const std::string& name);

    /**
     * @brief Get list of registered EOS types
     *
     * @return Vector of registered type names
     */
    static std::vector<std::string> registeredTypes();

    // =========================================================================
    // Databank-driven Creation
    // =========================================================================

    /**
     * @brief Build a mixture from CSV databanks and create an EOS by name.
     *
     * @param eos_type EOS type name (e.g., "PC-SAFT", "Peng-Robinson", "SRK")
     * @param databanks CSV-backed databanks
     * @param identifiers List of CAS/name identifiers
     * @param options Mixture build options
     * @param diag Optional diagnostics sink
     */
    static EOSPtr createFromDatabanks(
        const std::string& eos_type,
        const Data::Databanks& databanks,
        const std::vector<std::string>& identifiers,
        MixtureBuildOptions options = {},
        Data::Diagnostics* diag = nullptr
    );

private:
    static std::unordered_map<std::string, EOSCreator>& registry();
    static void initializeBuiltinEOS();
    static bool initialized_;
};

} // namespace Factory
} // namespace DMThermo

#endif // THERMO_FACTORY_EOS_FACTORY_H
