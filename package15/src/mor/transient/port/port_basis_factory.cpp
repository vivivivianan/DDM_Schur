#include "port_basis_factory.hpp"

namespace mor::transient::port {

std::unique_ptr<IPortBasisBuilder> makePortBasisBuilder(
    const RandomizedTransferPortOptions& options)
{
    return std::make_unique<RandomizedTransferBasisBuilder>(options);
}

} // namespace mor::transient::port
