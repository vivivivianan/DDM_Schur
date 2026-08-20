#include "parametric_reduced_schur.hpp"

#include "mor/pod_basis.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#include <mkl_lapacke.h>
#endif

namespace mor::parametric {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed(const Clock::time_point& start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

std::uint64_t hashIntegers(const std::vector<int>& values)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    for (int value : values) {
        const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (std::size_t i = 0; i < sizeof(value); ++i) {
            hash ^= bytes[i];
            hash *= UINT64_C(1099511628211);
        }
    }
    return hash;
}

std::vector<double> projectEntries(
    const std::vector<ddm_schur::Entry>& entries,
    const std::vector<double>& leftBasis,
    int leftRows,
    int leftRank,
    const std::vector<double>& rightBasis,
    int rightRows,
    int rightRank)
{
    std::vector<double> product(
        static_cast<std::size_t>(leftRows * rightRank), 0.0);
    for (const ddm_schur::Entry& entry : entries) {
        for (int column = 0; column < rightRank; ++column) {
            product[static_cast<std::size_t>(entry.row + column * leftRows)] +=
                entry.value * rightBasis[static_cast<std::size_t>(
                    entry.col + column * rightRows)];
        }
    }
    std::vector<double> reduced(
        static_cast<std::size_t>(leftRank * rightRank), 0.0);
#ifdef USE_MKL_PARDISO
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
        leftRank, rightRank, leftRows, 1.0,
        leftBasis.data(), leftRows,
        product.data(), leftRows,
        0.0, reduced.data(), leftRank);
#else
    for (int column = 0; column < rightRank; ++column) {
        for (int rowMode = 0; rowMode < leftRank; ++rowMode) {
            double value = 0.0;
            for (int row = 0; row < leftRows; ++row) {
                value += leftBasis[static_cast<std::size_t>(row + rowMode * leftRows)]
                    * product[static_cast<std::size_t>(row + column * leftRows)];
            }
            reduced[static_cast<std::size_t>(rowMode + column * leftRank)] = value;
        }
    }
#endif
    return reduced;
}

std::vector<double> projectVector(const std::vector<double>& basis,
                                  int rows,
                                  int rank,
                                  const std::vector<int>& globalDofs,
                                  const std::vector<double>& vector)
{
    std::vector<double> result(static_cast<std::size_t>(rank), 0.0);
    for (int mode = 0; mode < rank; ++mode) {
        double value = 0.0;
        for (int row = 0; row < rows; ++row) {
            value += basis[static_cast<std::size_t>(row + mode * rows)]
                * vector[static_cast<std::size_t>(
                    globalDofs[static_cast<std::size_t>(row)])];
        }
        result[static_cast<std::size_t>(mode)] = value;
    }
    return result;
}

std::vector<double> projectSources(const std::vector<double>& basis,
                                   int rows,
                                   int rank,
                                   const std::vector<int>& globalDofs,
                                   const SourceParameterization& sources)
{
    const int channels = static_cast<int>(sources.channels.size());
    std::vector<double> result(
        static_cast<std::size_t>(rank * channels), 0.0);
    for (int channel = 0; channel < channels; ++channel) {
        const std::vector<double>& load =
            sources.channels[static_cast<std::size_t>(channel)].rhsPerWatt;
        for (int mode = 0; mode < rank; ++mode) {
            double value = 0.0;
            for (int row = 0; row < rows; ++row) {
                value += basis[static_cast<std::size_t>(row + mode * rows)]
                    * load[static_cast<std::size_t>(
                        globalDofs[static_cast<std::size_t>(row)])];
            }
            result[static_cast<std::size_t>(mode + channel * rank)] = value;
        }
    }
    return result;
}

std::vector<double> effectiveOffset(const SparseMatrix& matrix,
                                    const std::vector<double>& rhs,
                                    const std::vector<double>& reference)
{
    const std::vector<double> image = matrix.multiply(reference);
    std::vector<double> result(rhs.size(), 0.0);
    for (std::size_t row = 0; row < rhs.size(); ++row) {
        result[row] = rhs[row] - image[row];
    }
    return result;
}

struct ProjectedEntries {
    std::vector<ddm_schur::Entry> gammaGamma;
    struct Local {
        std::vector<ddm_schur::Entry> ii;
        std::vector<ddm_schur::Entry> iGamma;
        std::vector<ddm_schur::Entry> gammaI;
    };
    std::vector<Local> locals;
};

ProjectedEntries splitEntries(const SparseMatrix& matrix,
                              const ddm_schur::InterfacePartition& partition)
{
    const int size = matrix.size();
    std::vector<int> domainSlot(static_cast<std::size_t>(size), -1);
    std::vector<int> interiorIndex(static_cast<std::size_t>(size), -1);
    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        const auto& domain = partition.domains[slot];
        for (std::size_t local = 0; local < domain.interiorGlobalDofs.size(); ++local) {
            const int global = domain.interiorGlobalDofs[local];
            domainSlot[static_cast<std::size_t>(global)] = static_cast<int>(slot);
            interiorIndex[static_cast<std::size_t>(global)] = static_cast<int>(local);
        }
        for (int global : domain.interfaceGlobalDofs) {
            domainSlot[static_cast<std::size_t>(global)] = static_cast<int>(slot);
        }
    }
    ProjectedEntries result;
    result.locals.resize(partition.domains.size());
    matrix.forEachEntry([&](int row, int column, double value) {
        const int gammaRow = partition.globalToInterface[static_cast<std::size_t>(row)];
        const int gammaColumn = partition.globalToInterface[static_cast<std::size_t>(column)];
        if (gammaRow >= 0 && gammaColumn >= 0) {
            result.gammaGamma.push_back({gammaRow, gammaColumn, value});
        } else if (gammaRow < 0 && gammaColumn < 0) {
            const int slot = domainSlot[static_cast<std::size_t>(row)];
            if (slot < 0 || slot != domainSlot[static_cast<std::size_t>(column)]) {
                throw std::runtime_error("Parametric interior block partition is inconsistent.");
            }
            result.locals[static_cast<std::size_t>(slot)].ii.push_back({
                interiorIndex[static_cast<std::size_t>(row)],
                interiorIndex[static_cast<std::size_t>(column)], value});
        } else if (gammaRow < 0) {
            const int slot = domainSlot[static_cast<std::size_t>(row)];
            result.locals[static_cast<std::size_t>(slot)].iGamma.push_back({
                interiorIndex[static_cast<std::size_t>(row)], gammaColumn, value});
        } else {
            const int slot = domainSlot[static_cast<std::size_t>(column)];
            result.locals[static_cast<std::size_t>(slot)].gammaI.push_back({
                gammaRow, interiorIndex[static_cast<std::size_t>(column)], value});
        }
    });
    return result;
}

using DenseSymmetricFactor = ParametricDenseFactorWorkspace;

void symmetrizeDense(std::vector<double>& matrix, int size)
{
    for (int column = 0; column < size; ++column) {
        for (int row = column + 1; row < size; ++row) {
            const std::size_t lower = static_cast<std::size_t>(row + column * size);
            const std::size_t upper = static_cast<std::size_t>(column + row * size);
            const double average = 0.5 * (matrix[lower] + matrix[upper]);
            matrix[lower] = average;
            matrix[upper] = average;
        }
    }
}

void factorDenseSymmetricInPlace(
    DenseSymmetricFactor& factor,
    int size,
    const std::string& context,
    const std::vector<double>* originalMatrix = nullptr)
{
    symmetrizeDense(factor.values, size);
#ifdef USE_MKL_PARDISO
    lapack_int info = LAPACKE_dpotrf(
        LAPACK_COL_MAJOR, 'L', size, factor.values.data(), size);
    if (info == 0) {
        factor.cholesky = true;
        return;
    }
    factor.cholesky = false;
    if (originalMatrix == nullptr) {
        throw std::runtime_error(
            "Parametric reduced Cholesky factorization failed for " + context
            + " and no unfactored fallback matrix was supplied.");
    }
    std::copy(originalMatrix->begin(), originalMatrix->end(), factor.values.begin());
    symmetrizeDense(factor.values, size);
    factor.pivots.resize(static_cast<std::size_t>(size));
    static_assert(sizeof(lapack_int) == sizeof(int),
                  "Parametric workspace requires LP64 LAPACK integers.");
    info = LAPACKE_dsytrf(
        LAPACK_COL_MAJOR, 'L', size, factor.values.data(), size,
        reinterpret_cast<lapack_int*>(factor.pivots.data()));
    if (info != 0) {
        throw std::runtime_error(
            "Parametric reduced symmetric factorization failed for " + context
            + " with info=" + std::to_string(info));
    }
#else
    factor.cholesky = true;
    for (int column = 0; column < size; ++column) {
        double diagonal = factor.values[static_cast<std::size_t>(column + column * size)];
        for (int k = 0; k < column; ++k) {
            const double value = factor.values[static_cast<std::size_t>(column + k * size)];
            diagonal -= value * value;
        }
        if (!(diagonal > 0.0)) {
            throw std::runtime_error(
                "Parametric reduced matrix is not SPD for " + context + '.');
        }
        factor.values[static_cast<std::size_t>(column + column * size)] =
            std::sqrt(diagonal);
        for (int row = column + 1; row < size; ++row) {
            double value = factor.values[static_cast<std::size_t>(row + column * size)];
            for (int k = 0; k < column; ++k) {
                value -= factor.values[static_cast<std::size_t>(row + k * size)]
                    * factor.values[static_cast<std::size_t>(column + k * size)];
            }
            factor.values[static_cast<std::size_t>(row + column * size)] =
                value / factor.values[static_cast<std::size_t>(column + column * size)];
        }
    }
#endif
}

void solveDenseSymmetric(const DenseSymmetricFactor& factor,
                         int size,
                         std::vector<double>& rhs,
                         int rightHandSides)
{
#ifdef USE_MKL_PARDISO
    const lapack_int info = factor.cholesky
        ? LAPACKE_dpotrs(
            LAPACK_COL_MAJOR, 'L', size, rightHandSides,
            factor.values.data(), size, rhs.data(), size)
        : LAPACKE_dsytrs(
            LAPACK_COL_MAJOR, 'L', size, rightHandSides,
            factor.values.data(), size,
            reinterpret_cast<const lapack_int*>(factor.pivots.data()),
            rhs.data(), size);
    if (info != 0) {
        throw std::runtime_error(
            "Parametric reduced dense solve failed with info="
            + std::to_string(info));
    }
#else
    for (int column = 0; column < rightHandSides; ++column) {
        double* vector = rhs.data() + static_cast<std::size_t>(column * size);
        for (int row = 0; row < size; ++row) {
            for (int k = 0; k < row; ++k) {
                vector[row] -= factor.values[static_cast<std::size_t>(row + k * size)]
                    * vector[k];
            }
            vector[row] /= factor.values[static_cast<std::size_t>(row + row * size)];
        }
        for (int row = size - 1; row >= 0; --row) {
            for (int k = row + 1; k < size; ++k) {
                vector[row] -= factor.values[static_cast<std::size_t>(k + row * size)]
                    * vector[k];
            }
            vector[row] /= factor.values[static_cast<std::size_t>(row + row * size)];
        }
    }
#endif
}

void gemmSubtract(const std::vector<double>& left,
                  int rows,
                  int inner,
                  const std::vector<double>& right,
                  int columns,
                  std::vector<double>& target)
{
#ifdef USE_MKL_PARDISO
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        rows, columns, inner, -1.0,
        left.data(), rows, right.data(), inner,
        1.0, target.data(), rows);
#else
    for (int column = 0; column < columns; ++column) {
        for (int k = 0; k < inner; ++k) {
            const double coefficient = right[static_cast<std::size_t>(k + column * inner)];
            for (int row = 0; row < rows; ++row) {
                target[static_cast<std::size_t>(row + column * rows)] -=
                    left[static_cast<std::size_t>(row + k * rows)] * coefficient;
            }
        }
    }
#endif
}

void affineDenseInto(const std::vector<double>& constant,
                     const std::vector<double>& linear,
                     const std::vector<double>& harmonic,
                     const std::vector<double>& coefficients,
                     std::vector<double>& result)
{
    if (result.size() != constant.size()) {
        result.resize(constant.size());
    }
    std::copy(constant.begin(), constant.end(), result.begin());
    for (std::size_t index = 0; index < result.size(); ++index) {
        result[index] += coefficients[0] * linear[index];
        for (std::size_t group = 1; group < coefficients.size(); ++group) {
            result[index] += coefficients[group]
                * harmonic[(group - 1) * result.size() + index];
        }
    }
}

void reducedRhsInto(const std::vector<double>& constant,
                    const std::vector<double>& linear,
                    const std::vector<double>& harmonic,
                    const std::vector<double>& coefficients,
                    const std::vector<double>& sources,
                    int rank,
                    int channels,
                    const std::vector<double>& powers,
                    std::vector<double>& result)
{
    if (result.size() != constant.size()) {
        result.resize(constant.size());
    }
    std::copy(constant.begin(), constant.end(), result.begin());
    for (int row = 0; row < rank; ++row) {
        result[static_cast<std::size_t>(row)] +=
            coefficients[0] * linear[static_cast<std::size_t>(row)];
        for (std::size_t group = 1; group < coefficients.size(); ++group) {
            result[static_cast<std::size_t>(row)] += coefficients[group]
                * harmonic[(group - 1) * static_cast<std::size_t>(rank)
                    + static_cast<std::size_t>(row)];
        }
    }
    for (int channel = 0; channel < channels; ++channel) {
        for (int row = 0; row < rank; ++row) {
            result[static_cast<std::size_t>(row)] +=
                powers[static_cast<std::size_t>(channel)]
                * sources[static_cast<std::size_t>(row + channel * rank)];
        }
    }
}

class BinaryWriter {
public:
    explicit BinaryWriter(const std::filesystem::path& path) : out_(path, std::ios::binary)
    {
        if (!out_) {
            throw std::runtime_error("Cannot create parametric ROM model: " + path.string());
        }
    }

    template <typename T> void scalar(const T& value)
    {
        out_.write(reinterpret_cast<const char*>(&value), sizeof(value));
    }

    void string(const std::string& value)
    {
        const std::uint64_t size = static_cast<std::uint64_t>(value.size());
        scalar(size);
        out_.write(value.data(), static_cast<std::streamsize>(value.size()));
    }

    template <typename T> void vector(const std::vector<T>& values)
    {
        const std::uint64_t size = static_cast<std::uint64_t>(values.size());
        scalar(size);
        if (!values.empty()) {
            out_.write(reinterpret_cast<const char*>(values.data()),
                       static_cast<std::streamsize>(values.size() * sizeof(T)));
        }
    }

private:
    std::ofstream out_;
};

class BinaryReader {
public:
    explicit BinaryReader(const std::filesystem::path& path) : in_(path, std::ios::binary)
    {
        if (!in_) {
            throw std::runtime_error("Cannot open parametric ROM model: " + path.string());
        }
    }

    template <typename T> T scalar()
    {
        T value{};
        in_.read(reinterpret_cast<char*>(&value), sizeof(value));
        require();
        return value;
    }

    std::string string()
    {
        const std::uint64_t size = scalar<std::uint64_t>();
        if (size > UINT64_C(1) << 30) {
            throw std::runtime_error("Invalid string in parametric ROM model.");
        }
        std::string value(static_cast<std::size_t>(size), '\0');
        in_.read(value.data(), static_cast<std::streamsize>(size));
        require();
        return value;
    }

    template <typename T> std::vector<T> vector()
    {
        const std::uint64_t size = scalar<std::uint64_t>();
        if (size > (UINT64_C(1) << 34)) {
            throw std::runtime_error("Invalid array in parametric ROM model.");
        }
        std::vector<T> values(static_cast<std::size_t>(size));
        if (!values.empty()) {
            in_.read(reinterpret_cast<char*>(values.data()),
                     static_cast<std::streamsize>(values.size() * sizeof(T)));
            require();
        }
        return values;
    }

private:
    void require()
    {
        if (!in_) {
            throw std::runtime_error("Parametric ROM model is truncated.");
        }
    }
    std::ifstream in_;
};

void writeLocal(BinaryWriter& writer, const ParametricLocalBlock& block)
{
    writer.scalar(block.domainId);
    writer.scalar(block.rank);
    writer.scalar(block.dofOrderingHash);
    writer.vector(block.globalDofs);
    writer.vector(block.basis);
    writer.vector(block.singularValues);
    writer.vector(block.aiiConstant);
    writer.vector(block.aiiLinear);
    writer.vector(block.aiiHarmonic);
    writer.vector(block.aiGammaConstant);
    writer.vector(block.aiGammaLinear);
    writer.vector(block.aiGammaHarmonic);
    writer.vector(block.aGammaIConstant);
    writer.vector(block.aGammaILinear);
    writer.vector(block.aGammaIHarmonic);
    writer.vector(block.rhsConstant);
    writer.vector(block.rhsLinear);
    writer.vector(block.rhsHarmonic);
    writer.vector(block.sourceChannels);
}

ParametricLocalBlock readLocal(BinaryReader& reader)
{
    ParametricLocalBlock block;
    block.domainId = reader.scalar<int>();
    block.rank = reader.scalar<int>();
    block.dofOrderingHash = reader.scalar<std::uint64_t>();
    block.globalDofs = reader.vector<int>();
    block.basis = reader.vector<double>();
    block.singularValues = reader.vector<double>();
    block.aiiConstant = reader.vector<double>();
    block.aiiLinear = reader.vector<double>();
    block.aiiHarmonic = reader.vector<double>();
    block.aiGammaConstant = reader.vector<double>();
    block.aiGammaLinear = reader.vector<double>();
    block.aiGammaHarmonic = reader.vector<double>();
    block.aGammaIConstant = reader.vector<double>();
    block.aGammaILinear = reader.vector<double>();
    block.aGammaIHarmonic = reader.vector<double>();
    block.rhsConstant = reader.vector<double>();
    block.rhsLinear = reader.vector<double>();
    block.rhsHarmonic = reader.vector<double>();
    block.sourceChannels = reader.vector<double>();
    return block;
}

std::vector<double> leadingDenseBlocks(
    const std::vector<double>& source,
    int oldRows,
    int oldColumns,
    int newRows,
    int newColumns,
    std::size_t blocks)
{
    if (source.size() != blocks * static_cast<std::size_t>(oldRows * oldColumns)) {
        throw std::runtime_error("Parametric model dense block dimensions are inconsistent.");
    }
    std::vector<double> result;
    result.reserve(blocks * static_cast<std::size_t>(newRows * newColumns));
    for (std::size_t block = 0; block < blocks; ++block) {
        const std::size_t offset = block * static_cast<std::size_t>(oldRows * oldColumns);
        for (int column = 0; column < newColumns; ++column) {
            for (int row = 0; row < newRows; ++row) {
                result.push_back(source[offset + static_cast<std::size_t>(
                    row + column * oldRows)]);
            }
        }
    }
    return result;
}

} // namespace

double modelOrthogonalityError(const std::vector<double>& basis,
                               int rows,
                               int columns)
{
    double maximum = 0.0;
    for (int left = 0; left < columns; ++left) {
        for (int right = 0; right < columns; ++right) {
#ifdef USE_MKL_PARDISO
            const double product = cblas_ddot(
                rows,
                basis.data() + static_cast<std::size_t>(left * rows), 1,
                basis.data() + static_cast<std::size_t>(right * rows), 1);
#else
            double product = 0.0;
            for (int row = 0; row < rows; ++row) {
                product += basis[static_cast<std::size_t>(row + left * rows)]
                    * basis[static_cast<std::size_t>(row + right * rows)];
            }
#endif
            maximum = std::max(maximum,
                std::abs(product - (left == right ? 1.0 : 0.0)));
        }
    }
    return maximum;
}

ParametricReducedModel buildParametricReducedModel(
    const Mesh& mesh,
    const AffineFemComponents& components,
    const ddm_schur::InterfacePartition& partition,
    const ParametricTrainingSnapshots& snapshots,
    const mor::Options& options)
{
    if (snapshots.interiors.size() != partition.domains.size()) {
        throw std::runtime_error("Parametric local snapshot partition mismatch.");
    }
    ParametricReducedModel model;
    model.globalDofs = static_cast<int>(mesh.nodes.size());
    model.interfaceDofs = static_cast<int>(partition.interfaceGlobalDofs.size());
    model.sourceChannels = static_cast<int>(components.sources.channels.size());
    model.parameter = components.parameter;
    model.affineConstantHash = components.constantMatrixHash;
    model.affineLinearHash = components.linearMatrixHash;
    model.affineHarmonicHashes = components.harmonicMatrixHashes;
    model.referenceTemperature = snapshots.referenceTemperature;
    model.interfaceGlobalDofs = partition.interfaceGlobalDofs;
    model.dofs.resize(mesh.nodes.size());
    for (std::size_t row = 0; row < mesh.nodes.size(); ++row) {
        model.dofs[row] = {mesh.nodes[row].p.x, mesh.nodes[row].p.y,
                           mesh.nodes[row].p.z, mesh.nodes[row].subdomain,
                           mesh.nodes[row].sourceVertex};
    }

    const auto basisStart = Clock::now();
    const PodResult interfacePod = buildGramPod(
        snapshots.interfaceSnapshots, options.interfaceRank,
        options.energyTolerance, options.singularValueTolerance);
    model.interfaceRank = interfacePod.selectedRank;
    model.interfaceBasis = interfacePod.basis;
    model.interfaceSingularValues = interfacePod.singularValues;
    model.locals.resize(partition.domains.size());
    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        const PodResult localPod = buildGramPod(
            snapshots.interiors[slot], options.localRank,
            options.interiorEnergyTolerance,
            options.interiorSingularValueTolerance);
        ParametricLocalBlock& block = model.locals[slot];
        block.domainId = partition.domains[slot].domainId;
        block.rank = localPod.selectedRank;
        block.globalDofs = partition.domains[slot].interiorGlobalDofs;
        block.dofOrderingHash = hashIntegers(block.globalDofs);
        block.basis = localPod.basis;
        block.singularValues = localPod.singularValues;
    }
    model.basisSeconds = elapsed(basisStart);

    const auto projectionStart = Clock::now();
    const ProjectedEntries constantEntries = splitEntries(
        components.matrixConstant, partition);
    const ProjectedEntries linearEntries = splitEntries(
        components.matrixLinear, partition);
    std::vector<ProjectedEntries> harmonicEntries;
    harmonicEntries.reserve(components.matrixHarmonic.size());
    for (const SparseMatrix& matrix : components.matrixHarmonic) {
        harmonicEntries.push_back(splitEntries(matrix, partition));
    }
    model.projectionPreparationSeconds = elapsed(projectionStart);
    const auto interfaceProjectionStart = Clock::now();
    model.aGammaGammaConstant = projectEntries(
        constantEntries.gammaGamma,
        model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
        model.interfaceBasis, model.interfaceDofs, model.interfaceRank);
    model.aGammaGammaLinear = projectEntries(
        linearEntries.gammaGamma,
        model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
        model.interfaceBasis, model.interfaceDofs, model.interfaceRank);
    for (const ProjectedEntries& entries : harmonicEntries) {
        const std::vector<double> projected = projectEntries(
            entries.gammaGamma,
            model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
            model.interfaceBasis, model.interfaceDofs, model.interfaceRank);
        model.aGammaGammaHarmonic.insert(
            model.aGammaGammaHarmonic.end(), projected.begin(), projected.end());
    }

    const std::vector<double> offsetConstant = effectiveOffset(
        components.matrixConstant, components.rhsConstant,
        model.referenceTemperature);
    const std::vector<double> offsetLinear = effectiveOffset(
        components.matrixLinear, components.rhsLinear,
        model.referenceTemperature);
    std::vector<std::vector<double>> offsetHarmonic;
    for (std::size_t group = 0; group < components.matrixHarmonic.size(); ++group) {
        offsetHarmonic.push_back(effectiveOffset(
            components.matrixHarmonic[group], components.rhsHarmonic[group],
            model.referenceTemperature));
    }
    model.rhsGammaConstant = projectVector(
        model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
        model.interfaceGlobalDofs, offsetConstant);
    model.rhsGammaLinear = projectVector(
        model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
        model.interfaceGlobalDofs, offsetLinear);
    for (const std::vector<double>& offset : offsetHarmonic) {
        const std::vector<double> projected = projectVector(
            model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
            model.interfaceGlobalDofs, offset);
        model.rhsGammaHarmonic.insert(
            model.rhsGammaHarmonic.end(), projected.begin(), projected.end());
    }
    model.sourceGamma = projectSources(
        model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
        model.interfaceGlobalDofs, components.sources);
    model.interfaceProjectionSeconds = elapsed(interfaceProjectionStart);

    const auto localProjectionStart = Clock::now();
    for (std::size_t slot = 0; slot < model.locals.size(); ++slot) {
        ParametricLocalBlock& block = model.locals[slot];
        const int rows = static_cast<int>(block.globalDofs.size());
        const int rank = block.rank;
        block.aiiConstant = projectEntries(
            constantEntries.locals[slot].ii,
            block.basis, rows, rank, block.basis, rows, rank);
        block.aiiLinear = projectEntries(
            linearEntries.locals[slot].ii,
            block.basis, rows, rank, block.basis, rows, rank);
        for (const ProjectedEntries& entries : harmonicEntries) {
            const std::vector<double> projected = projectEntries(
                entries.locals[slot].ii,
                block.basis, rows, rank, block.basis, rows, rank);
            block.aiiHarmonic.insert(
                block.aiiHarmonic.end(), projected.begin(), projected.end());
        }
        block.aiGammaConstant = projectEntries(
            constantEntries.locals[slot].iGamma,
            block.basis, rows, rank,
            model.interfaceBasis, model.interfaceDofs, model.interfaceRank);
        block.aiGammaLinear = projectEntries(
            linearEntries.locals[slot].iGamma,
            block.basis, rows, rank,
            model.interfaceBasis, model.interfaceDofs, model.interfaceRank);
        for (const ProjectedEntries& entries : harmonicEntries) {
            const std::vector<double> projected = projectEntries(
                entries.locals[slot].iGamma,
                block.basis, rows, rank,
                model.interfaceBasis, model.interfaceDofs, model.interfaceRank);
            block.aiGammaHarmonic.insert(
                block.aiGammaHarmonic.end(), projected.begin(), projected.end());
        }
        block.aGammaIConstant = projectEntries(
            constantEntries.locals[slot].gammaI,
            model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
            block.basis, rows, rank);
        block.aGammaILinear = projectEntries(
            linearEntries.locals[slot].gammaI,
            model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
            block.basis, rows, rank);
        for (const ProjectedEntries& entries : harmonicEntries) {
            const std::vector<double> projected = projectEntries(
                entries.locals[slot].gammaI,
                model.interfaceBasis, model.interfaceDofs, model.interfaceRank,
                block.basis, rows, rank);
            block.aGammaIHarmonic.insert(
                block.aGammaIHarmonic.end(), projected.begin(), projected.end());
        }
        block.rhsConstant = projectVector(
            block.basis, rows, rank, block.globalDofs, offsetConstant);
        block.rhsLinear = projectVector(
            block.basis, rows, rank, block.globalDofs, offsetLinear);
        for (const std::vector<double>& offset : offsetHarmonic) {
            const std::vector<double> projected = projectVector(
                block.basis, rows, rank, block.globalDofs, offset);
            block.rhsHarmonic.insert(
                block.rhsHarmonic.end(), projected.begin(), projected.end());
        }
        block.sourceChannels = projectSources(
            block.basis, rows, rank, block.globalDofs, components.sources);
    }
    model.localProjectionSeconds = elapsed(localProjectionStart);
    model.projectionSeconds = elapsed(projectionStart);

    for (const SourceChannel& channel : components.sources.channels) {
        model.nominalPowersW.push_back(channel.nominalPowerW);
        model.minimumPowersW.push_back(channel.minimumPowerW);
        model.maximumPowersW.push_back(channel.maximumPowerW);
        model.sourceSubdomains.push_back(channel.subdomain);
        model.sourceDomainEntities.push_back(channel.domainEntity);
    }
    return model;
}

void ParametricRomWorkspace::initialize(const ParametricReducedModel& model)
{
    globalDofs = model.globalDofs;
    interfaceDofs = model.interfaceDofs;
    interfaceRank = model.interfaceRank;
    sourceChannels = model.sourceChannels;
    parameterCoefficients.assign(
        model.parameter.harmonicNeighborConductivities.size() + 1, 0.0);
    powersW.assign(static_cast<std::size_t>(sourceChannels), 0.0);
    schurFactor.values.assign(
        static_cast<std::size_t>(interfaceRank * interfaceRank), 0.0);
    schurFactor.pivots.assign(static_cast<std::size_t>(interfaceRank), 0);
    interfaceCoordinates.assign(static_cast<std::size_t>(interfaceRank), 0.0);
    locals.resize(model.locals.size());
    std::size_t largestFactor = schurFactor.values.size();
    for (std::size_t slot = 0; slot < model.locals.size(); ++slot) {
        const int rank = model.locals[slot].rank;
        ParametricLocalWorkspace& local = locals[slot];
        local.rank = rank;
        local.factor.values.assign(static_cast<std::size_t>(rank * rank), 0.0);
        local.factor.pivots.assign(static_cast<std::size_t>(rank), 0);
        local.eliminatedCoupling.assign(
            static_cast<std::size_t>(rank * interfaceRank), 0.0);
        local.rhsSolution.assign(static_cast<std::size_t>(rank), 0.0);
        local.aGammaI.assign(
            static_cast<std::size_t>(interfaceRank * rank), 0.0);
        largestFactor = std::max(largestFactor, local.factor.values.size());
    }
    factorBackup.clear();
    factorBackup.reserve(largestFactor);
    fullTemperature.assign(static_cast<std::size_t>(globalDofs), 0.0);
    interfaceTemperature.assign(static_cast<std::size_t>(interfaceDofs), 0.0);
    errorScratch.clear();
    outputScratch.assign(64U * 1024U, 0);
}

void ParametricRomWorkspace::resetForNextCase()
{
    std::fill(parameterCoefficients.begin(), parameterCoefficients.end(), 0.0);
    std::fill(powersW.begin(), powersW.end(), 0.0);
    std::fill(interfaceCoordinates.begin(), interfaceCoordinates.end(), 0.0);
    // Large matrices and reconstructed fields are overwritten in place by the
    // next solve, deliberately preserving their capacity and virtual address.
}

namespace {

template <typename T>
std::size_t vectorCapacityBytes(const std::vector<T>& vector)
{
    return vector.capacity() * sizeof(T);
}

bool workspaceMatches(const ParametricRomWorkspace& workspace,
                      const ParametricReducedModel& model)
{
    if (workspace.globalDofs != model.globalDofs
        || workspace.interfaceDofs != model.interfaceDofs
        || workspace.interfaceRank != model.interfaceRank
        || workspace.sourceChannels != model.sourceChannels
        || workspace.locals.size() != model.locals.size()
        || workspace.parameterCoefficients.size()
            != model.parameter.harmonicNeighborConductivities.size() + 1) {
        return false;
    }
    for (std::size_t slot = 0; slot < model.locals.size(); ++slot) {
        if (workspace.locals[slot].rank != model.locals[slot].rank) {
            return false;
        }
    }
    return true;
}

} // namespace

std::size_t ParametricRomWorkspace::workspaceBytes() const
{
    std::size_t bytes = vectorCapacityBytes(parameterCoefficients)
        + vectorCapacityBytes(powersW) + vectorCapacityBytes(factorBackup)
        + vectorCapacityBytes(schurFactor.values)
        + vectorCapacityBytes(schurFactor.pivots)
        + vectorCapacityBytes(interfaceCoordinates);
    for (const ParametricLocalWorkspace& local : locals) {
        bytes += vectorCapacityBytes(local.factor.values)
            + vectorCapacityBytes(local.factor.pivots)
            + vectorCapacityBytes(local.eliminatedCoupling)
            + vectorCapacityBytes(local.rhsSolution)
            + vectorCapacityBytes(local.aGammaI);
    }
    return bytes;
}

std::size_t ParametricRomWorkspace::reconstructionBytes() const
{
    return vectorCapacityBytes(fullTemperature)
        + vectorCapacityBytes(interfaceTemperature)
        + vectorCapacityBytes(errorScratch);
}

std::size_t ParametricModelMemoryBreakdown::totalBytes() const
{
    return interfaceBasisBytes + localBasisBytes + reducedBlockBytes
        + sourceBytes + metadataBytes;
}

ParametricModelMemoryBreakdown parametricModelMemoryBreakdown(
    const ParametricReducedModel& model)
{
    ParametricModelMemoryBreakdown bytes;
    bytes.interfaceBasisBytes = vectorCapacityBytes(model.interfaceBasis)
        + vectorCapacityBytes(model.interfaceSingularValues);
    bytes.reducedBlockBytes = vectorCapacityBytes(model.aGammaGammaConstant)
        + vectorCapacityBytes(model.aGammaGammaLinear)
        + vectorCapacityBytes(model.aGammaGammaHarmonic)
        + vectorCapacityBytes(model.rhsGammaConstant)
        + vectorCapacityBytes(model.rhsGammaLinear)
        + vectorCapacityBytes(model.rhsGammaHarmonic);
    bytes.sourceBytes = vectorCapacityBytes(model.sourceGamma)
        + vectorCapacityBytes(model.nominalPowersW)
        + vectorCapacityBytes(model.minimumPowersW)
        + vectorCapacityBytes(model.maximumPowersW);
    bytes.metadataBytes = vectorCapacityBytes(model.dofs)
        + vectorCapacityBytes(model.interfaceGlobalDofs)
        + vectorCapacityBytes(model.referenceTemperature)
        + vectorCapacityBytes(model.sourceSubdomains)
        + vectorCapacityBytes(model.sourceDomainEntities)
        + vectorCapacityBytes(model.parameter.harmonicNeighborConductivities)
        + vectorCapacityBytes(model.affineHarmonicHashes);
    for (const ParametricLocalBlock& local : model.locals) {
        bytes.localBasisBytes += vectorCapacityBytes(local.basis)
            + vectorCapacityBytes(local.singularValues);
        bytes.reducedBlockBytes += vectorCapacityBytes(local.aiiConstant)
            + vectorCapacityBytes(local.aiiLinear)
            + vectorCapacityBytes(local.aiiHarmonic)
            + vectorCapacityBytes(local.aiGammaConstant)
            + vectorCapacityBytes(local.aiGammaLinear)
            + vectorCapacityBytes(local.aiGammaHarmonic)
            + vectorCapacityBytes(local.aGammaIConstant)
            + vectorCapacityBytes(local.aGammaILinear)
            + vectorCapacityBytes(local.aGammaIHarmonic)
            + vectorCapacityBytes(local.rhsConstant)
            + vectorCapacityBytes(local.rhsLinear)
            + vectorCapacityBytes(local.rhsHarmonic);
        bytes.sourceBytes += vectorCapacityBytes(local.sourceChannels);
        bytes.metadataBytes += vectorCapacityBytes(local.globalDofs);
    }
    return bytes;
}

ParametricOnlineResult solveParametricRomInWorkspace(
    const ParametricReducedModel& model,
    double parameterValue,
    const std::vector<double>& powersW,
    bool allowExtrapolation,
    ParametricRomWorkspace& workspace)
{
    const auto totalStart = Clock::now();
    if (!std::isfinite(parameterValue) || !(parameterValue > 0.0)) {
        throw std::runtime_error("Parametric ROM input is not a positive finite value.");
    }
    if (powersW.size() != static_cast<std::size_t>(model.sourceChannels)) {
        throw std::runtime_error("Parametric ROM power vector has the wrong dimension.");
    }
    if (!workspaceMatches(workspace, model)) {
        workspace.initialize(model);
    } else {
        workspace.resetForNextCase();
    }
    ParametricOnlineResult result;
    result.parameterValue = parameterValue;
    result.extrapolated = parameterValue < model.parameter.minimum
        || parameterValue > model.parameter.maximum;
    if (result.extrapolated && !allowExtrapolation) {
        throw std::runtime_error(
            "Pure parametric ROM rejected an out-of-training-range parameter. Use --mor-allow-extrapolation to override explicitly.");
    }

    const auto coefficientStart = Clock::now();
    workspace.parameterCoefficients[0] = parameterValue;
    for (std::size_t group = 0;
         group < model.parameter.harmonicNeighborConductivities.size(); ++group) {
        workspace.parameterCoefficients[group + 1] =
            harmonicTheta(model.parameter, parameterValue, group);
    }
    std::copy(powersW.begin(), powersW.end(), workspace.powersW.begin());
    const int gammaRank = model.interfaceRank;
    affineDenseInto(model.aGammaGammaConstant, model.aGammaGammaLinear,
                    model.aGammaGammaHarmonic, workspace.parameterCoefficients,
                    workspace.schurFactor.values);
    reducedRhsInto(model.rhsGammaConstant, model.rhsGammaLinear,
                   model.rhsGammaHarmonic, workspace.parameterCoefficients,
                   model.sourceGamma, gammaRank, model.sourceChannels,
                   workspace.powersW, workspace.interfaceCoordinates);
    result.timing.coefficientSeconds = elapsed(coefficientStart);

    const auto localAssemblyStart = Clock::now();
    for (std::size_t slot = 0; slot < model.locals.size(); ++slot) {
        const ParametricLocalBlock& block = model.locals[slot];
        ParametricLocalWorkspace& local = workspace.locals[slot];
        affineDenseInto(block.aiiConstant, block.aiiLinear, block.aiiHarmonic,
                        workspace.parameterCoefficients, local.factor.values);
        affineDenseInto(block.aiGammaConstant, block.aiGammaLinear,
                        block.aiGammaHarmonic, workspace.parameterCoefficients,
                        local.eliminatedCoupling);
        affineDenseInto(block.aGammaIConstant, block.aGammaILinear,
                        block.aGammaIHarmonic, workspace.parameterCoefficients,
                        local.aGammaI);
        reducedRhsInto(block.rhsConstant, block.rhsLinear, block.rhsHarmonic,
                       workspace.parameterCoefficients, block.sourceChannels,
                       block.rank, model.sourceChannels, workspace.powersW,
                       local.rhsSolution);
    }
    result.timing.localAssemblySeconds = elapsed(localAssemblyStart);

    const auto localFactorStart = Clock::now();
    for (std::size_t slot = 0; slot < workspace.locals.size(); ++slot) {
        ParametricLocalWorkspace& local = workspace.locals[slot];
        workspace.factorBackup.assign(
            local.factor.values.begin(), local.factor.values.end());
        factorDenseSymmetricInPlace(
            local.factor, local.rank,
            "subdomain " + std::to_string(model.locals[slot].domainId),
            &workspace.factorBackup);
        solveDenseSymmetric(local.factor, local.rank,
                            local.eliminatedCoupling, gammaRank);
        solveDenseSymmetric(local.factor, local.rank, local.rhsSolution, 1);
    }
    result.timing.localFactorSeconds = elapsed(localFactorStart);

    const auto schurAssemblyStart = Clock::now();
    for (const ParametricLocalWorkspace& local : workspace.locals) {
        gemmSubtract(local.aGammaI, gammaRank, local.rank,
                     local.eliminatedCoupling, gammaRank,
                     workspace.schurFactor.values);
        for (int row = 0; row < gammaRank; ++row) {
            double correction = 0.0;
            for (int k = 0; k < local.rank; ++k) {
                correction += local.aGammaI[static_cast<std::size_t>(row + k * gammaRank)]
                    * local.rhsSolution[static_cast<std::size_t>(k)];
            }
            workspace.interfaceCoordinates[static_cast<std::size_t>(row)] -= correction;
        }
    }
    result.timing.reducedSchurAssemblySeconds = elapsed(schurAssemblyStart);

    const auto reducedSolveStart = Clock::now();
    workspace.factorBackup.assign(workspace.schurFactor.values.begin(),
                                  workspace.schurFactor.values.end());
    factorDenseSymmetricInPlace(workspace.schurFactor, gammaRank,
                                "global reduced Schur system",
                                &workspace.factorBackup);
    solveDenseSymmetric(workspace.schurFactor, gammaRank,
                        workspace.interfaceCoordinates, 1);
    for (ParametricLocalWorkspace& local : workspace.locals) {
        for (int row = 0; row < local.rank; ++row) {
            double correction = 0.0;
            for (int column = 0; column < gammaRank; ++column) {
                correction += local.eliminatedCoupling[
                    static_cast<std::size_t>(row + column * local.rank)]
                    * workspace.interfaceCoordinates[static_cast<std::size_t>(column)];
            }
            local.rhsSolution[static_cast<std::size_t>(row)] -= correction;
        }
    }
    result.timing.reducedSolveSeconds = elapsed(reducedSolveStart);

    const auto reconstructionStart = Clock::now();
    std::copy(model.referenceTemperature.begin(), model.referenceTemperature.end(),
              workspace.fullTemperature.begin());
    for (int row = 0; row < model.interfaceDofs; ++row) {
        const int global = model.interfaceGlobalDofs[static_cast<std::size_t>(row)];
        double value = model.referenceTemperature[static_cast<std::size_t>(global)];
        for (int mode = 0; mode < gammaRank; ++mode) {
            value += model.interfaceBasis[static_cast<std::size_t>(
                row + mode * model.interfaceDofs)]
                * workspace.interfaceCoordinates[static_cast<std::size_t>(mode)];
        }
        workspace.fullTemperature[static_cast<std::size_t>(global)] = value;
        workspace.interfaceTemperature[static_cast<std::size_t>(row)] = value;
    }
    for (std::size_t slot = 0; slot < model.locals.size(); ++slot) {
        const ParametricLocalBlock& block = model.locals[slot];
        const ParametricLocalWorkspace& local = workspace.locals[slot];
        const int rows = static_cast<int>(block.globalDofs.size());
        for (int row = 0; row < rows; ++row) {
            const int global = block.globalDofs[static_cast<std::size_t>(row)];
            double value = model.referenceTemperature[static_cast<std::size_t>(global)];
            for (int mode = 0; mode < block.rank; ++mode) {
                value += block.basis[static_cast<std::size_t>(row + mode * rows)]
                    * local.rhsSolution[static_cast<std::size_t>(mode)];
            }
            workspace.fullTemperature[static_cast<std::size_t>(global)] = value;
        }
    }
    result.timing.reconstructionSeconds = elapsed(reconstructionStart);
    result.timing.totalSeconds = elapsed(totalStart);
    result.status = "success";
    return result;
}

ParametricOnlineResult solveParametricRom(
    const ParametricReducedModel& model,
    double parameterValue,
    const std::vector<double>& powersW,
    bool allowExtrapolation)
{
    ParametricRomWorkspace workspace;
    ParametricOnlineResult result = solveParametricRomInWorkspace(
        model, parameterValue, powersW, allowExtrapolation, workspace);
    result.powersW = powersW;
    result.temperature = std::move(workspace.fullTemperature);
    result.interfaceTemperature = std::move(workspace.interfaceTemperature);
    return result;
}

void saveParametricModel(const ParametricReducedModel& model,
                         const std::filesystem::path& directory)
{
    std::filesystem::create_directories(directory);
    BinaryWriter writer(directory / "model.bin");
    writer.string("DDM_SCHUR_STAGE2B1");
    writer.scalar(model.formatVersion);
    writer.scalar(model.globalDofs);
    writer.scalar(model.interfaceDofs);
    writer.scalar(model.interfaceRank);
    writer.scalar(model.sourceChannels);
    writer.string(model.parameter.name);
    writer.string(model.parameter.units);
    writer.scalar(model.parameter.subdomain);
    writer.scalar(model.parameter.regionId);
    writer.scalar(model.parameter.minimum);
    writer.scalar(model.parameter.maximum);
    writer.scalar(model.parameter.reference);
    writer.scalar(model.parameter.ambientTemperature);
    writer.vector(model.parameter.harmonicNeighborConductivities);
    writer.scalar(model.parameter.touchesInterface);
    writer.scalar(model.parameter.selectedTetCount);
    writer.scalar(model.parameter.selectedBoundaryFaceCount);
    writer.scalar(model.parameter.definitionHash);
    writer.scalar(model.fingerprints.mesh);
    writer.scalar(model.fingerprints.system);
    writer.scalar(model.fingerprints.interfaceOrdering);
    writer.scalar(model.fingerprints.sources);
    writer.scalar(model.affineConstantHash);
    writer.scalar(model.affineLinearHash);
    writer.vector(model.affineHarmonicHashes);
    writer.scalar(model.sourceDefinitionHash);
    writer.vector(model.dofs);
    writer.vector(model.interfaceGlobalDofs);
    writer.vector(model.referenceTemperature);
    writer.vector(model.interfaceBasis);
    writer.vector(model.interfaceSingularValues);
    writer.vector(model.aGammaGammaConstant);
    writer.vector(model.aGammaGammaLinear);
    writer.vector(model.aGammaGammaHarmonic);
    writer.vector(model.rhsGammaConstant);
    writer.vector(model.rhsGammaLinear);
    writer.vector(model.rhsGammaHarmonic);
    writer.vector(model.sourceGamma);
    writer.vector(model.nominalPowersW);
    writer.vector(model.minimumPowersW);
    writer.vector(model.maximumPowersW);
    writer.vector(model.sourceSubdomains);
    writer.vector(model.sourceDomainEntities);
    const std::uint64_t localCount = static_cast<std::uint64_t>(model.locals.size());
    writer.scalar(localCount);
    for (const ParametricLocalBlock& block : model.locals) {
        writeLocal(writer, block);
    }

    std::ofstream metadata(directory / "metadata.json");
    metadata << std::setprecision(17)
        << "{\n"
        << "  \"format_version\": " << model.formatVersion << ",\n"
        << "  \"model_type\": \"parametric-local-rb-reduced-schur\",\n"
        << "  \"matrix_parameter\": \"" << model.parameter.name << "\",\n"
        << "  \"parameter_units\": \"" << model.parameter.units << "\",\n"
        << "  \"parameter_subdomain\": " << model.parameter.subdomain << ",\n"
        << "  \"parameter_region_id\": " << model.parameter.regionId << ",\n"
        << "  \"training_min\": " << model.parameter.minimum << ",\n"
        << "  \"training_max\": " << model.parameter.maximum << ",\n"
        << "  \"reference\": " << model.parameter.reference << ",\n"
        << "  \"harmonic_penalty_groups\": "
        << model.parameter.harmonicNeighborConductivities.size() << ",\n"
        << "  \"global_dofs\": " << model.globalDofs << ",\n"
        << "  \"interface_dofs\": " << model.interfaceDofs << ",\n"
        << "  \"interface_rank\": " << model.interfaceRank << ",\n"
        << "  \"source_channels\": " << model.sourceChannels << ",\n"
        << "  \"subdomains\": " << model.locals.size() << ",\n"
        << "  \"pure_online_sparse_matrix\": false,\n"
        << "  \"pure_online_pardiso\": false\n"
        << "}\n";
    std::ofstream definition(directory / "parameter_definition.json");
    definition << std::setprecision(17)
        << "{\"name\":\"" << model.parameter.name
        << "\",\"units\":\"" << model.parameter.units
        << "\",\"subdomain\":" << model.parameter.subdomain
        << ",\"region_id\":" << model.parameter.regionId
        << ",\"minimum\":" << model.parameter.minimum
        << ",\"maximum\":" << model.parameter.maximum
        << ",\"reference\":" << model.parameter.reference << "}\n";
    std::ofstream hashes(directory / "affine_operator_hashes.json");
    hashes << "{\"constant\":\"" << std::hex << model.affineConstantHash
        << "\",\"linear\":\"" << model.affineLinearHash
        << "\",\"harmonic_groups\":[";
    for (std::size_t group = 0; group < model.affineHarmonicHashes.size(); ++group) {
        if (group > 0) hashes << ',';
        hashes << "\"" << model.affineHarmonicHashes[group] << "\"";
    }
    hashes << "],\"parameter_definition\":\"" << model.parameter.definitionHash
        << "\"}\n";

    std::ofstream singular(directory / "singular_values.csv");
    singular << "space,subdomain,mode,singular_value,selected\n"
        << std::setprecision(17);
    for (std::size_t mode = 0; mode < model.interfaceSingularValues.size(); ++mode) {
        singular << "interface,-1," << mode << ','
            << model.interfaceSingularValues[mode] << ','
            << (mode < static_cast<std::size_t>(model.interfaceRank) ? 1 : 0) << '\n';
    }
    for (const ParametricLocalBlock& block : model.locals) {
        for (std::size_t mode = 0; mode < block.singularValues.size(); ++mode) {
            singular << "interior," << block.domainId << ',' << mode << ','
                << block.singularValues[mode] << ','
                << (mode < static_cast<std::size_t>(block.rank) ? 1 : 0) << '\n';
        }
    }
    std::ofstream ranks(directory / "local_ranks.csv");
    ranks << "subdomain,interior_dofs,selected_rank\n";
    for (const ParametricLocalBlock& block : model.locals) {
        ranks << block.domainId << ',' << block.globalDofs.size() << ','
            << block.rank << '\n';
    }
}

ParametricReducedModel loadParametricModel(
    const std::filesystem::path& directory,
    const Fingerprints* expectedFingerprints,
    const AffineFemComponents* expectedAffine)
{
    BinaryReader reader(directory / "model.bin");
    if (reader.string() != "DDM_SCHUR_STAGE2B1") {
        throw std::runtime_error("Not a Stage 2B.1 parametric ROM model.");
    }
    ParametricReducedModel model;
    model.formatVersion = reader.scalar<int>();
    if (model.formatVersion != 3) {
        throw std::runtime_error("Unsupported Stage 2B.1 model format.");
    }
    model.globalDofs = reader.scalar<int>();
    model.interfaceDofs = reader.scalar<int>();
    model.interfaceRank = reader.scalar<int>();
    model.sourceChannels = reader.scalar<int>();
    model.parameter.name = reader.string();
    model.parameter.units = reader.string();
    model.parameter.subdomain = reader.scalar<int>();
    model.parameter.regionId = reader.scalar<int>();
    model.parameter.minimum = reader.scalar<double>();
    model.parameter.maximum = reader.scalar<double>();
    model.parameter.reference = reader.scalar<double>();
    model.parameter.ambientTemperature = reader.scalar<double>();
    model.parameter.harmonicNeighborConductivities = reader.vector<double>();
    model.parameter.touchesInterface = reader.scalar<bool>();
    model.parameter.selectedTetCount = reader.scalar<int>();
    model.parameter.selectedBoundaryFaceCount = reader.scalar<int>();
    model.parameter.definitionHash = reader.scalar<std::uint64_t>();
    model.fingerprints.mesh = reader.scalar<std::uint64_t>();
    model.fingerprints.system = reader.scalar<std::uint64_t>();
    model.fingerprints.interfaceOrdering = reader.scalar<std::uint64_t>();
    model.fingerprints.sources = reader.scalar<std::uint64_t>();
    model.affineConstantHash = reader.scalar<std::uint64_t>();
    model.affineLinearHash = reader.scalar<std::uint64_t>();
    model.affineHarmonicHashes = reader.vector<std::uint64_t>();
    model.sourceDefinitionHash = reader.scalar<std::uint64_t>();
    model.dofs = reader.vector<DeploymentDof>();
    model.interfaceGlobalDofs = reader.vector<int>();
    model.referenceTemperature = reader.vector<double>();
    model.interfaceBasis = reader.vector<double>();
    model.interfaceSingularValues = reader.vector<double>();
    model.aGammaGammaConstant = reader.vector<double>();
    model.aGammaGammaLinear = reader.vector<double>();
    model.aGammaGammaHarmonic = reader.vector<double>();
    model.rhsGammaConstant = reader.vector<double>();
    model.rhsGammaLinear = reader.vector<double>();
    model.rhsGammaHarmonic = reader.vector<double>();
    model.sourceGamma = reader.vector<double>();
    model.nominalPowersW = reader.vector<double>();
    model.minimumPowersW = reader.vector<double>();
    model.maximumPowersW = reader.vector<double>();
    model.sourceSubdomains = reader.vector<int>();
    model.sourceDomainEntities = reader.vector<int>();
    const std::uint64_t localCount = reader.scalar<std::uint64_t>();
    model.locals.reserve(static_cast<std::size_t>(localCount));
    for (std::uint64_t index = 0; index < localCount; ++index) {
        model.locals.push_back(readLocal(reader));
    }
    model.fileBytes = std::filesystem::file_size(directory / "model.bin");

    if (expectedFingerprints != nullptr
        && (model.fingerprints.mesh != expectedFingerprints->mesh
            || model.fingerprints.system != expectedFingerprints->system
            || model.fingerprints.interfaceOrdering != expectedFingerprints->interfaceOrdering
            || model.fingerprints.sources != expectedFingerprints->sources)) {
        throw std::runtime_error(
            "Stage 2B.1 model fingerprint rejection: mesh/system/interface/source definition differs.");
    }
    if (expectedAffine != nullptr
        && (model.parameter.definitionHash != expectedAffine->parameter.definitionHash
            || model.affineConstantHash != expectedAffine->constantMatrixHash
            || model.affineLinearHash != expectedAffine->linearMatrixHash
            || model.affineHarmonicHashes != expectedAffine->harmonicMatrixHashes)) {
        throw std::runtime_error(
            "Stage 2B.1 model fingerprint rejection: affine parameter/operator differs.");
    }
    return model;
}

void truncateParametricModel(ParametricReducedModel& model,
                             int interfaceRank,
                             int localRank)
{
    const int oldGammaRank = model.interfaceRank;
    const int newGammaRank = interfaceRank > 0
        ? std::min(interfaceRank, oldGammaRank) : oldGammaRank;
    if (newGammaRank <= 0 || localRank < 0) {
        throw std::runtime_error("Parametric ROM truncation ranks must be positive or zero for unchanged.");
    }
    const std::size_t harmonicGroups =
        model.parameter.harmonicNeighborConductivities.size();
    model.interfaceBasis.resize(static_cast<std::size_t>(
        model.interfaceDofs * newGammaRank));
    model.aGammaGammaConstant = leadingDenseBlocks(
        model.aGammaGammaConstant, oldGammaRank, oldGammaRank,
        newGammaRank, newGammaRank, 1);
    model.aGammaGammaLinear = leadingDenseBlocks(
        model.aGammaGammaLinear, oldGammaRank, oldGammaRank,
        newGammaRank, newGammaRank, 1);
    model.aGammaGammaHarmonic = leadingDenseBlocks(
        model.aGammaGammaHarmonic, oldGammaRank, oldGammaRank,
        newGammaRank, newGammaRank, harmonicGroups);
    model.rhsGammaConstant.resize(static_cast<std::size_t>(newGammaRank));
    model.rhsGammaLinear.resize(static_cast<std::size_t>(newGammaRank));
    model.rhsGammaHarmonic = leadingDenseBlocks(
        model.rhsGammaHarmonic, oldGammaRank, 1,
        newGammaRank, 1, harmonicGroups);
    model.sourceGamma = leadingDenseBlocks(
        model.sourceGamma, oldGammaRank, 1,
        newGammaRank, 1, static_cast<std::size_t>(model.sourceChannels));

    for (ParametricLocalBlock& block : model.locals) {
        const int oldLocalRank = block.rank;
        const int newLocalRank = localRank > 0
            ? std::min(localRank, oldLocalRank) : oldLocalRank;
        block.basis.resize(block.globalDofs.size()
            * static_cast<std::size_t>(newLocalRank));
        block.aiiConstant = leadingDenseBlocks(
            block.aiiConstant, oldLocalRank, oldLocalRank,
            newLocalRank, newLocalRank, 1);
        block.aiiLinear = leadingDenseBlocks(
            block.aiiLinear, oldLocalRank, oldLocalRank,
            newLocalRank, newLocalRank, 1);
        block.aiiHarmonic = leadingDenseBlocks(
            block.aiiHarmonic, oldLocalRank, oldLocalRank,
            newLocalRank, newLocalRank, harmonicGroups);
        block.aiGammaConstant = leadingDenseBlocks(
            block.aiGammaConstant, oldLocalRank, oldGammaRank,
            newLocalRank, newGammaRank, 1);
        block.aiGammaLinear = leadingDenseBlocks(
            block.aiGammaLinear, oldLocalRank, oldGammaRank,
            newLocalRank, newGammaRank, 1);
        block.aiGammaHarmonic = leadingDenseBlocks(
            block.aiGammaHarmonic, oldLocalRank, oldGammaRank,
            newLocalRank, newGammaRank, harmonicGroups);
        block.aGammaIConstant = leadingDenseBlocks(
            block.aGammaIConstant, oldGammaRank, oldLocalRank,
            newGammaRank, newLocalRank, 1);
        block.aGammaILinear = leadingDenseBlocks(
            block.aGammaILinear, oldGammaRank, oldLocalRank,
            newGammaRank, newLocalRank, 1);
        block.aGammaIHarmonic = leadingDenseBlocks(
            block.aGammaIHarmonic, oldGammaRank, oldLocalRank,
            newGammaRank, newLocalRank, harmonicGroups);
        block.rhsConstant.resize(static_cast<std::size_t>(newLocalRank));
        block.rhsLinear.resize(static_cast<std::size_t>(newLocalRank));
        block.rhsHarmonic = leadingDenseBlocks(
            block.rhsHarmonic, oldLocalRank, 1,
            newLocalRank, 1, harmonicGroups);
        block.sourceChannels = leadingDenseBlocks(
            block.sourceChannels, oldLocalRank, 1,
            newLocalRank, 1, static_cast<std::size_t>(model.sourceChannels));
        block.rank = newLocalRank;
    }
    model.interfaceRank = newGammaRank;
}

} // namespace mor::parametric
