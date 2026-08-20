#include "local_subdomain_model.hpp"

namespace mor::local {

namespace {

template <typename T>
std::size_t vectorBytes(const std::vector<T>& values)
{
    return values.capacity() * sizeof(T);
}

} // namespace

std::size_t estimateModelBytes(const Model& model)
{
    std::size_t bytes = sizeof(Model)
        + vectorBytes(model.interfaceGlobalDofs)
        + vectorBytes(model.interfaceEntries)
        + model.subdomains.capacity() * sizeof(SubdomainModel);
    for (const SubdomainModel& local : model.subdomains) {
        bytes += vectorBytes(local.interiorGlobalDofs)
            + vectorBytes(local.interfaceGlobalDofs)
            + vectorBytes(local.interfaceIndices)
            + vectorBytes(local.referenceInterior)
            + vectorBytes(local.basis)
            + vectorBytes(local.singularValues)
            + vectorBytes(local.retainedEnergy)
            + vectorBytes(local.reducedInterior)
            + vectorBytes(local.reducedInteriorInterface)
            + vectorBytes(local.reducedInterfaceInterior)
            + vectorBytes(local.interiorReferenceImage)
            + vectorBytes(local.interfaceReferenceImage);
    }
    return bytes;
}

} // namespace mor::local
