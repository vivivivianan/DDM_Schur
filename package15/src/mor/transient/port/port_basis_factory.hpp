#pragma once

#include "IPortBasisBuilder.hpp"
#include "randomized_transfer_port.hpp"

#include <memory>

namespace mor::transient::port {

std::unique_ptr<IPortBasisBuilder> makePortBasisBuilder(
    const RandomizedTransferPortOptions& options);

} // namespace mor::transient::port
