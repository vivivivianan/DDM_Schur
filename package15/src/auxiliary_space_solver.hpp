#pragma once

// Independent P2/P1 auxiliary-space prototype.  This header is included by
// main.cpp after sipg_core.hpp and linear_solvers.hpp so it deliberately uses
// the existing mesh, CSR, PARDISO, and PCG services without changing them.

struct P2P1AuxiliaryReport {
    int p2Dofs = 0;
    int p1Dofs = 0;
    std::size_t prolongationNnz = 0;
    int prolongationRowsWithOneEntry = 0;
    int prolongationRowsWithTwoEntries = 0;
    double prolongationBuildSeconds = 0.0;
    double constantError = 0.0;
    double linearXError = 0.0;
    double linearYError = 0.0;
    double linearZError = 0.0;
    std::size_t coarseNnz = 0;
    double coarseAssemblySeconds = 0.0;
    double coarseEstimatedFactorBytes = 0.0;
    double coarseSymbolicSeconds = 0.0;
    double coarseNumericalSeconds = 0.0;
    std::size_t coarseFactorBytes = 0;
    double smootherSetupSeconds = 0.0;
    std::size_t smootherBytes = 0;
    double a2SymmetryError = 0.0;
    double acSymmetryError = 0.0;
    double acPositiveEnergy = 0.0;
    double preconditionerSymmetryError = 0.0;
    double preconditionerPositiveEnergy = 0.0;
    double setupSeconds = 0.0;
    int pcgIterations = 0;
    int preconditionerCalls = 0;
    double preconditionerSeconds = 0.0;
    double solveSeconds = 0.0;
    double totalSeconds = 0.0;
    double trueResidual = std::numeric_limits<double>::quiet_NaN();
    std::size_t peakMemoryBytes = 0;
    std::string status = "not_run";
    std::string stopReason;
};

static double auxiliaryRelativeBilinearMismatch(double left, double right)
{
    return std::abs(left - right)
        / std::max({1.0, std::abs(left), std::abs(right)});
}

class P2P1AuxiliaryPreconditioner {
public:
    P2P1AuxiliaryPreconditioner(const Mesh& mesh,
                                const SparseMatrix& fine,
                                P2P1AuxiliaryReport& report,
                                bool enforceLargeCaseGates)
        : fine_(&fine), report_(&report), largeCaseGates_(enforceLargeCaseGates)
    {
        const auto setupStart = std::chrono::steady_clock::now();
        buildProlongation(mesh);
        checkSetupTime(setupStart, "prolongation_setup_gate_exceeded");
        buildInheritedCoarse(fine);
        checkCoarseMemoryEstimate();
        checkSetupTime(setupStart, "inherited_coarse_setup_gate_exceeded");
        buildElementSmoother(mesh, fine);
        checkSetupTime(setupStart, "element_smoother_setup_gate_exceeded");
        factorCoarse();
        checkSetupTime(setupStart, "coarse_factorization_setup_gate_exceeded");
        runAlgebraicSafetyTests(fine);
        checkSetupTime(setupStart, "algebraic_safety_test_setup_gate_exceeded");
        report_->status = "setup_passed";
    }

    void apply(const std::vector<double>& residual, std::vector<double>& result)
    {
        if (residual.size() != static_cast<std::size_t>(fineRows_)) {
            throw std::runtime_error("P2/P1 auxiliary preconditioner dimension mismatch.");
        }
        const auto start = std::chrono::steady_clock::now();
        applyElementSmoother(residual, result);

        std::vector<double> coarseRhs(static_cast<std::size_t>(coarseRows_), 0.0);
        for (int row = 0; row < fineRows_; ++row) {
            for (int k = pRowPtr_[static_cast<std::size_t>(row)];
                 k < pRowPtr_[static_cast<std::size_t>(row + 1)]; ++k) {
                coarseRhs[static_cast<std::size_t>(pColInd_[static_cast<std::size_t>(k)])]
                    += pValues_[static_cast<std::size_t>(k)]
                     * residual[static_cast<std::size_t>(row)];
            }
        }
        std::vector<double> coarseSolution;
        coarseSolver_->solve(coarseRhs, coarseSolution);
        parallelFor(static_cast<std::size_t>(fineRows_), [&](std::size_t row) {
            double value = 0.0;
            for (int k = pRowPtr_[row]; k < pRowPtr_[row + 1]; ++k) {
                value += pValues_[static_cast<std::size_t>(k)]
                       * coarseSolution[static_cast<std::size_t>(
                           pColInd_[static_cast<std::size_t>(k)])];
            }
            result[row] += value;
        });
        report_->preconditionerSeconds += elapsed(start);
        ++report_->preconditionerCalls;
    }

    std::size_t memoryBytes() const
    {
        return pRowPtr_.capacity() * sizeof(int)
            + pColInd_.capacity() * sizeof(int)
            + pValues_.capacity() * sizeof(double)
            + coarse_.rowPtr.capacity() * sizeof(int)
            + coarse_.colInd.capacity() * sizeof(int)
            + coarse_.values.capacity() * sizeof(double)
            + report_->smootherBytes
            + (coarseSolver_ ? coarseSolver_->memoryBytes() : 0);
    }

private:
    struct ElementFactor {
        std::array<int, 10> dofs{};
        std::array<double, 55> lower{};
    };

    const SparseMatrix* fine_ = nullptr;
    P2P1AuxiliaryReport* report_ = nullptr;
    bool largeCaseGates_ = false;
    int fineRows_ = 0;
    int coarseRows_ = 0;
    std::vector<int> pRowPtr_;
    std::vector<int> pColInd_;
    std::vector<double> pValues_;
    SparseMatrix coarse_;
    std::unique_ptr<SubdomainDirectSolver> coarseSolver_;
    std::vector<ElementFactor> elementFactors_;
    std::vector<std::vector<double>> smootherThreadBuffers_;

    static double elapsed(const std::chrono::steady_clock::time_point& start)
    {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
    }

    void fail(const std::string& reason)
    {
        report_->status = "stopped";
        report_->stopReason = reason;
        throw std::runtime_error("P2/P1 auxiliary-space stopped: " + reason);
    }

    void checkSetupTime(const std::chrono::steady_clock::time_point& setupStart,
                        const std::string& reason)
    {
        report_->setupSeconds = elapsed(setupStart);
        report_->peakMemoryBytes = std::max(report_->peakMemoryBytes, peakWorkingSetBytes());
        if (largeCaseGates_ && report_->setupSeconds > 35.0) {
            fail(reason);
        }
    }

    void buildProlongation(const Mesh& mesh)
    {
        const auto start = std::chrono::steady_clock::now();
        fineRows_ = static_cast<int>(mesh.nodes.size());
        report_->p2Dofs = fineRows_;
        std::vector<unsigned char> isVertex(static_cast<std::size_t>(fineRows_), 0);
        for (const Tet& tet : mesh.tets) {
            for (int local = 0; local < 4; ++local) {
                isVertex[static_cast<std::size_t>(tet.dof[static_cast<std::size_t>(local)])] = 1;
            }
        }
        std::vector<int> vertexToCoarse(static_cast<std::size_t>(fineRows_), -1);
        for (int global = 0; global < fineRows_; ++global) {
            if (isVertex[static_cast<std::size_t>(global)] != 0) {
                vertexToCoarse[static_cast<std::size_t>(global)] = coarseRows_++;
            }
        }
        std::vector<std::array<int, 2>> midpointEndpoints(
            static_cast<std::size_t>(fineRows_), {{-1, -1}});
        const std::array<std::array<int, 2>, 6> edges{{
            {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}
        }};
        for (const Tet& tet : mesh.tets) {
            for (int edge = 0; edge < 6; ++edge) {
                const int midpoint = tet.dof[static_cast<std::size_t>(4 + edge)];
                std::array<int, 2> endpoints{{
                    tet.dof[static_cast<std::size_t>(edges[static_cast<std::size_t>(edge)][0])],
                    tet.dof[static_cast<std::size_t>(edges[static_cast<std::size_t>(edge)][1])]}};
                std::sort(endpoints.begin(), endpoints.end());
                std::array<int, 2>& stored = midpointEndpoints[static_cast<std::size_t>(midpoint)];
                if (stored[0] < 0) {
                    stored = endpoints;
                } else if (stored != endpoints) {
                    fail("inconsistent_nodal_p2_edge_ordering");
                }
            }
        }

        pRowPtr_.assign(static_cast<std::size_t>(fineRows_ + 1), 0);
        pColInd_.reserve(static_cast<std::size_t>(fineRows_) * 2);
        pValues_.reserve(static_cast<std::size_t>(fineRows_) * 2);
        for (int row = 0; row < fineRows_; ++row) {
            if (isVertex[static_cast<std::size_t>(row)] != 0) {
                pColInd_.push_back(vertexToCoarse[static_cast<std::size_t>(row)]);
                pValues_.push_back(1.0);
                ++report_->prolongationRowsWithOneEntry;
            } else {
                const std::array<int, 2>& endpoints = midpointEndpoints[static_cast<std::size_t>(row)];
                if (endpoints[0] < 0 || endpoints[1] < 0) {
                    fail("unclassified_p2_dof");
                }
                std::array<int, 2> coarseEndpoints{{
                    vertexToCoarse[static_cast<std::size_t>(endpoints[0])],
                    vertexToCoarse[static_cast<std::size_t>(endpoints[1])]}};
                std::sort(coarseEndpoints.begin(), coarseEndpoints.end());
                if (coarseEndpoints[0] < 0 || coarseEndpoints[1] < 0) {
                    fail("invalid_p1_endpoint_mapping");
                }
                pColInd_.push_back(coarseEndpoints[0]);
                pValues_.push_back(0.5);
                pColInd_.push_back(coarseEndpoints[1]);
                pValues_.push_back(0.5);
                ++report_->prolongationRowsWithTwoEntries;
            }
            pRowPtr_[static_cast<std::size_t>(row + 1)] =
                static_cast<int>(pColInd_.size());
        }
        report_->p1Dofs = coarseRows_;
        report_->prolongationNnz = pValues_.size();
        verifyProlongation(mesh, vertexToCoarse);
        report_->prolongationBuildSeconds = elapsed(start);
    }

    void verifyProlongation(const Mesh& mesh, const std::vector<int>& vertexToCoarse)
    {
        auto errorForField = [&](int field) {
            std::vector<double> coarseValues(static_cast<std::size_t>(coarseRows_), 0.0);
            for (int row = 0; row < fineRows_; ++row) {
                const int coarse = vertexToCoarse[static_cast<std::size_t>(row)];
                if (coarse >= 0) {
                    const Vec3& p = mesh.nodes[static_cast<std::size_t>(row)].p;
                    coarseValues[static_cast<std::size_t>(coarse)] = field == 0 ? 1.0
                        : (field == 1 ? p.x : (field == 2 ? p.y : p.z));
                }
            }
            double maximum = 0.0;
            for (int row = 0; row < fineRows_; ++row) {
                double prolongated = 0.0;
                for (int k = pRowPtr_[static_cast<std::size_t>(row)];
                     k < pRowPtr_[static_cast<std::size_t>(row + 1)]; ++k) {
                    prolongated += pValues_[static_cast<std::size_t>(k)]
                        * coarseValues[static_cast<std::size_t>(
                            pColInd_[static_cast<std::size_t>(k)])];
                }
                const Vec3& p = mesh.nodes[static_cast<std::size_t>(row)].p;
                const double exact = field == 0 ? 1.0
                    : (field == 1 ? p.x : (field == 2 ? p.y : p.z));
                maximum = std::max(maximum, std::abs(prolongated - exact));
            }
            return maximum;
        };
        report_->constantError = errorForField(0);
        report_->linearXError = errorForField(1);
        report_->linearYError = errorForField(2);
        report_->linearZError = errorForField(3);
        const double coordinateScale = std::max({1.0,
            maxCoordinateMagnitude(mesh, 0), maxCoordinateMagnitude(mesh, 1),
            maxCoordinateMagnitude(mesh, 2)});
        if (report_->constantError > 1.0e-14
            || report_->linearXError > 1.0e-13 * coordinateScale
            || report_->linearYError > 1.0e-13 * coordinateScale
            || report_->linearZError > 1.0e-13 * coordinateScale) {
            fail("p1_to_p2_polynomial_reproduction_failed");
        }
    }

    static double maxCoordinateMagnitude(const Mesh& mesh, int field)
    {
        double maximum = 0.0;
        for (const Node& node : mesh.nodes) {
            const double value = field == 0 ? node.p.x : (field == 1 ? node.p.y : node.p.z);
            maximum = std::max(maximum, std::abs(value));
        }
        return maximum;
    }

    void buildInheritedCoarse(const SparseMatrix& fine)
    {
        const auto start = std::chrono::steady_clock::now();
        if (!fine.csrReady) {
            fail("fine_matrix_csr_required");
        }
#ifdef USE_MKL_PARDISO
        sparse_matrix_t fineHandle = nullptr;
        sparse_matrix_t pHandle = nullptr;
        sparse_matrix_t apHandle = nullptr;
        sparse_matrix_t coarseHandle = nullptr;
        auto destroyAll = [&]() {
            if (coarseHandle) mkl_sparse_destroy(coarseHandle);
            if (apHandle) mkl_sparse_destroy(apHandle);
            if (pHandle) mkl_sparse_destroy(pHandle);
            if (fineHandle) mkl_sparse_destroy(fineHandle);
        };
        sparse_status_t status = mkl_sparse_d_create_csr(
            &fineHandle, SPARSE_INDEX_BASE_ZERO, fineRows_, fineRows_,
            const_cast<int*>(fine.rowPtr.data()),
            const_cast<int*>(fine.rowPtr.data() + 1),
            const_cast<int*>(fine.colInd.data()),
            const_cast<double*>(fine.values.data()));
        if (status == SPARSE_STATUS_SUCCESS) {
            status = mkl_sparse_d_create_csr(
                &pHandle, SPARSE_INDEX_BASE_ZERO, fineRows_, coarseRows_,
                pRowPtr_.data(), pRowPtr_.data() + 1,
                pColInd_.data(), pValues_.data());
        }
        if (status == SPARSE_STATUS_SUCCESS) {
            status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE,
                                     fineHandle, pHandle, &apHandle);
        }
        if (status == SPARSE_STATUS_SUCCESS) {
            status = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE,
                                     pHandle, apHandle, &coarseHandle);
        }
        if (status != SPARSE_STATUS_SUCCESS) {
            destroyAll();
            fail("mkl_sparse_ptap_failed_" + std::to_string(static_cast<int>(status)));
        }
        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
        MKL_INT rows = 0, columns = 0;
        MKL_INT* starts = nullptr;
        MKL_INT* ends = nullptr;
        MKL_INT* columnsOut = nullptr;
        double* valuesOut = nullptr;
        status = mkl_sparse_d_export_csr(
            coarseHandle, &indexing, &rows, &columns,
            &starts, &ends, &columnsOut, &valuesOut);
        if (status != SPARSE_STATUS_SUCCESS || rows != coarseRows_ || columns != coarseRows_) {
            destroyAll();
            fail("mkl_sparse_ptap_export_failed");
        }
        coarse_ = SparseMatrix(coarseRows_);
        coarse_.rowPtr.resize(static_cast<std::size_t>(coarseRows_ + 1));
        coarse_.rowPtr[0] = 0;
        for (int row = 0; row < coarseRows_; ++row) {
            coarse_.rowPtr[static_cast<std::size_t>(row + 1)] = ends[row];
        }
        const std::size_t nnz = static_cast<std::size_t>(
            coarse_.rowPtr[static_cast<std::size_t>(coarseRows_)]);
        coarse_.colInd.assign(columnsOut, columnsOut + nnz);
        coarse_.values.assign(valuesOut, valuesOut + nnz);
        coarse_.csrReady = true;
        destroyAll();
#else
        coarse_ = SparseMatrix(coarseRows_);
        coarse_.triplets.reserve(fine.values.size() * 2);
        fine.forEachEntry([&](int row, int column, double value) {
            for (int kr = pRowPtr_[static_cast<std::size_t>(row)];
                 kr < pRowPtr_[static_cast<std::size_t>(row + 1)]; ++kr) {
                for (int kc = pRowPtr_[static_cast<std::size_t>(column)];
                     kc < pRowPtr_[static_cast<std::size_t>(column + 1)]; ++kc) {
                    coarse_.add(pColInd_[static_cast<std::size_t>(kr)],
                                pColInd_[static_cast<std::size_t>(kc)],
                                pValues_[static_cast<std::size_t>(kr)] * value
                              * pValues_[static_cast<std::size_t>(kc)]);
                }
            }
        });
        coarse_.finalizeCsr();
#endif
        report_->coarseNnz = coarse_.values.size();
        report_->coarseAssemblySeconds = elapsed(start);
        // Conservative sparse-Cholesky planning estimate used only for the
        // mandatory pre-factorization stop gate; actual PARDISO memory is
        // reported after phase 22.
        report_->coarseEstimatedFactorBytes = static_cast<double>(report_->coarseNnz)
            * 16.0 * 12.0;
    }

    void checkCoarseMemoryEstimate()
    {
        constexpr double limit = 7585.6 * 1024.0 * 1024.0;
        if (largeCaseGates_ && report_->coarseEstimatedFactorBytes >= limit) {
            fail("coarse_estimated_factor_memory_gate_exceeded");
        }
    }

    static double matrixValue(const SparseMatrix& matrix, int row, int column)
    {
        const int begin = matrix.rowPtr[static_cast<std::size_t>(row)];
        const int end = matrix.rowPtr[static_cast<std::size_t>(row + 1)];
        const auto first = matrix.colInd.begin() + begin;
        const auto last = matrix.colInd.begin() + end;
        const auto found = std::lower_bound(first, last, column);
        if (found == last || *found != column) {
            return 0.0;
        }
        return matrix.values[static_cast<std::size_t>(found - matrix.colInd.begin())];
    }

    void buildElementSmoother(const Mesh& mesh, const SparseMatrix& fine)
    {
        const auto start = std::chrono::steady_clock::now();
        elementFactors_.resize(mesh.tets.size());
        std::vector<std::string> errors(mesh.tets.size());
        parallelForCoarse(mesh.tets.size(), [&](std::size_t index) {
            try {
                ElementFactor& factor = elementFactors_[index];
                factor.dofs = mesh.tets[index].dof;
                for (int row = 0; row < 10; ++row) {
                    for (int column = 0; column <= row; ++column) {
                        double value = matrixValue(fine, factor.dofs[static_cast<std::size_t>(row)],
                                                   factor.dofs[static_cast<std::size_t>(column)]);
                        for (int k = 0; k < column; ++k) {
                            value -= factor.lower[packed(row, k)]
                                   * factor.lower[packed(column, k)];
                        }
                        if (row == column) {
                            if (!(value > 0.0) || !std::isfinite(value)) {
                                throw std::runtime_error("nonpositive element principal pivot");
                            }
                            factor.lower[packed(row, column)] = std::sqrt(value);
                        } else {
                            factor.lower[packed(row, column)] = value
                                / factor.lower[packed(column, column)];
                        }
                    }
                }
            } catch (const std::exception& error) {
                errors[index] = error.what();
            }
        });
        for (std::size_t index = 0; index < errors.size(); ++index) {
            if (!errors[index].empty()) {
                fail("element_cholesky_failed_at_" + std::to_string(index)
                    + ":" + errors[index]);
            }
        }
        const int workers = static_cast<int>(std::max(1u, solverParallelWorkers()));
        smootherThreadBuffers_.assign(static_cast<std::size_t>(workers),
                                      std::vector<double>(static_cast<std::size_t>(fineRows_), 0.0));
        report_->smootherBytes = elementFactors_.capacity() * sizeof(ElementFactor)
            + smootherThreadBuffers_.size() * static_cast<std::size_t>(fineRows_) * sizeof(double);
        report_->smootherSetupSeconds = elapsed(start);
    }

    static constexpr std::size_t packed(int row, int column)
    {
        return static_cast<std::size_t>(row * (row + 1) / 2 + column);
    }

    void factorCoarse()
    {
        const std::vector<MatrixEntry> entries = sparseMatrixEntries(coarse_);
        coarseSolver_ = std::make_unique<SubdomainDirectSolver>(coarseRows_, entries);
        report_->coarseSymbolicSeconds = coarseSolver_->symbolicAnalysisSeconds();
        report_->coarseNumericalSeconds = coarseSolver_->numericalFactorizationSeconds();
        report_->coarseFactorBytes = coarseSolver_->memoryBytes();
        report_->peakMemoryBytes = std::max(report_->peakMemoryBytes, peakWorkingSetBytes());
    }

    void applyElementSmoother(const std::vector<double>& residual,
                              std::vector<double>& result)
    {
        for (std::vector<double>& buffer : smootherThreadBuffers_) {
            std::fill(buffer.begin(), buffer.end(), 0.0);
        }
        const int workers = static_cast<int>(smootherThreadBuffers_.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers)
#endif
        for (long long element = 0;
             element < static_cast<long long>(elementFactors_.size()); ++element) {
#ifdef _OPENMP
            const int thread = omp_get_thread_num();
#else
            const int thread = 0;
#endif
            const ElementFactor& factor = elementFactors_[static_cast<std::size_t>(element)];
            std::array<double, 10> y{};
            std::array<double, 10> x{};
            for (int row = 0; row < 10; ++row) {
                double value = residual[static_cast<std::size_t>(
                    factor.dofs[static_cast<std::size_t>(row)])];
                for (int column = 0; column < row; ++column) {
                    value -= factor.lower[packed(row, column)]
                           * y[static_cast<std::size_t>(column)];
                }
                y[static_cast<std::size_t>(row)] = value
                    / factor.lower[packed(row, row)];
            }
            for (int row = 9; row >= 0; --row) {
                double value = y[static_cast<std::size_t>(row)];
                for (int column = row + 1; column < 10; ++column) {
                    value -= factor.lower[packed(column, row)]
                           * x[static_cast<std::size_t>(column)];
                }
                x[static_cast<std::size_t>(row)] = value
                    / factor.lower[packed(row, row)];
            }
            std::vector<double>& buffer = smootherThreadBuffers_[static_cast<std::size_t>(thread)];
            for (int local = 0; local < 10; ++local) {
                buffer[static_cast<std::size_t>(factor.dofs[static_cast<std::size_t>(local)])]
                    += x[static_cast<std::size_t>(local)];
            }
        }
        result.assign(static_cast<std::size_t>(fineRows_), 0.0);
        parallelFor(result.size(), [&](std::size_t row) {
            double value = 0.0;
            for (const std::vector<double>& buffer : smootherThreadBuffers_) {
                value += buffer[row];
            }
            result[row] = value;
        });
    }

    void runAlgebraicSafetyTests(const SparseMatrix& fine)
    {
        std::vector<double> x(static_cast<std::size_t>(fineRows_));
        std::vector<double> y(static_cast<std::size_t>(fineRows_));
        for (int i = 0; i < fineRows_; ++i) {
            x[static_cast<std::size_t>(i)] = std::sin(0.173 * static_cast<double>(i + 1));
            y[static_cast<std::size_t>(i)] = std::cos(0.097 * static_cast<double>(i + 3));
        }
        const std::vector<double> ax = fine.multiply(x);
        const std::vector<double> ay = fine.multiply(y);
        report_->a2SymmetryError = auxiliaryRelativeBilinearMismatch(
            vectorDot(x, ay), vectorDot(ax, y));

        std::vector<double> xc(static_cast<std::size_t>(coarseRows_));
        std::vector<double> yc(static_cast<std::size_t>(coarseRows_));
        for (int i = 0; i < coarseRows_; ++i) {
            xc[static_cast<std::size_t>(i)] = std::sin(0.131 * static_cast<double>(i + 1));
            yc[static_cast<std::size_t>(i)] = std::cos(0.071 * static_cast<double>(i + 2));
        }
        const std::vector<double> acx = coarse_.multiply(xc);
        const std::vector<double> acy = coarse_.multiply(yc);
        report_->acSymmetryError = auxiliaryRelativeBilinearMismatch(
            vectorDot(xc, acy), vectorDot(acx, yc));
        report_->acPositiveEnergy = vectorDot(xc, acx);

        const int callsBefore = report_->preconditionerCalls;
        const double secondsBefore = report_->preconditionerSeconds;
        std::vector<double> mx;
        std::vector<double> my;
        apply(x, mx);
        apply(y, my);
        report_->preconditionerSymmetryError = auxiliaryRelativeBilinearMismatch(
            vectorDot(x, my), vectorDot(mx, y));
        report_->preconditionerPositiveEnergy = vectorDot(x, mx);
        report_->preconditionerCalls = callsBefore;
        report_->preconditionerSeconds = secondsBefore;

        if (!std::isfinite(report_->a2SymmetryError)
            || report_->a2SymmetryError >= 1.0e-12) {
            fail("a2_symmetry_gate_failed");
        }
        if (!std::isfinite(report_->acSymmetryError)
            || report_->acSymmetryError >= 1.0e-12) {
            fail("ac_symmetry_gate_failed");
        }
        if (!(report_->acPositiveEnergy > 0.0)
            || !std::isfinite(report_->acPositiveEnergy)) {
            fail("ac_positive_energy_gate_failed");
        }
        if (!std::isfinite(report_->preconditionerSymmetryError)
            || report_->preconditionerSymmetryError >= 1.0e-12) {
            fail("preconditioner_symmetry_gate_failed");
        }
        if (!(report_->preconditionerPositiveEnergy > 0.0)
            || !std::isfinite(report_->preconditionerPositiveEnergy)) {
            fail("preconditioner_positive_energy_gate_failed");
        }
    }
};

static void writeP2P1AuxiliaryReport(const P2P1AuxiliaryReport& report,
                                     const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << std::setprecision(16)
        << "status,stop_reason,p2_dofs,p1_dofs,p_nnz,p_rows_one,p_rows_two,"
        << "p_build_s,constant_error,linear_x_error,linear_y_error,linear_z_error,"
        << "ac_nnz,ac_build_s,ac_estimated_factor_bytes,ac_phase11_s,ac_phase22_s,"
        << "ac_factor_bytes,smoother_setup_s,smoother_bytes,a2_symmetry_error,"
        << "ac_symmetry_error,ac_positive_energy,preconditioner_symmetry_error,"
        << "preconditioner_positive_energy,setup_s,pcg_iterations,preconditioner_calls,"
        << "preconditioner_total_s,preconditioner_average_s,solve_s,total_s,"
        << "true_residual,peak_memory_bytes\n"
        << report.status << ',' << csvEscape(report.stopReason) << ','
        << report.p2Dofs << ',' << report.p1Dofs << ',' << report.prolongationNnz << ','
        << report.prolongationRowsWithOneEntry << ','
        << report.prolongationRowsWithTwoEntries << ','
        << report.prolongationBuildSeconds << ',' << report.constantError << ','
        << report.linearXError << ',' << report.linearYError << ',' << report.linearZError << ','
        << report.coarseNnz << ',' << report.coarseAssemblySeconds << ','
        << report.coarseEstimatedFactorBytes << ',' << report.coarseSymbolicSeconds << ','
        << report.coarseNumericalSeconds << ',' << report.coarseFactorBytes << ','
        << report.smootherSetupSeconds << ',' << report.smootherBytes << ','
        << report.a2SymmetryError << ',' << report.acSymmetryError << ','
        << report.acPositiveEnergy << ',' << report.preconditionerSymmetryError << ','
        << report.preconditionerPositiveEnergy << ',' << report.setupSeconds << ','
        << report.pcgIterations << ',' << report.preconditionerCalls << ','
        << report.preconditionerSeconds << ','
        << (report.preconditionerCalls > 0
            ? report.preconditionerSeconds / report.preconditionerCalls : 0.0) << ','
        << report.solveSeconds << ',' << report.totalSeconds << ',' << report.trueResidual << ','
        << report.peakMemoryBytes << '\n';
}
