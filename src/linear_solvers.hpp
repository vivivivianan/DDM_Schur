#pragma once

// System construction, direct solvers, preconditioners, Krylov solvers, and transient/steady stepping.
// This file is intentionally included from main.cpp after the preceding SIPG modules.

static void buildSystemMatrix(const SparseMatrix& mass, const SparseMatrix& stiffness, double dt, double stiffnessWeight, SparseMatrix& system)
{
    system.appendScaledEntries(mass, 1.0 / dt);
    system.appendScaledEntries(stiffness, stiffnessWeight);
}

static std::string coordinateMatchKey(const Vec3& p, double tolerance)
{
    const long long ix = static_cast<long long>(std::llround(p.x / tolerance));
    const long long iy = static_cast<long long>(std::llround(p.y / tolerance));
    const long long iz = static_cast<long long>(std::llround(p.z / tolerance));
    return std::to_string(ix) + "," + std::to_string(iy) + "," + std::to_string(iz);
}

static int assembleMatchedNodeTieInterface(const Mesh& mesh,
                                           SparseMatrix& stiffness,
                                           double penalty,
                                           const std::filesystem::path& outputDir)
{
    std::set<int> leftDofs;
    std::set<int> rightDofs;
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        for (int i = 0; i < 10; ++i) {
            leftDofs.insert(left.dof[static_cast<size_t>(i)]);
            rightDofs.insert(right.dof[static_cast<size_t>(i)]);
        }
    }

    constexpr double tolerance = 1.0e-14;
    std::map<std::string, std::vector<int>> rightByCoordinate;
    for (int dof : rightDofs) {
        rightByCoordinate[coordinateMatchKey(mesh.nodes[static_cast<size_t>(dof)].p, tolerance)].push_back(dof);
    }

    std::ofstream out(outputDir / "rram_node_tie_interface_diagnostics.csv");
    out << "left_dof,right_dof,distance,penalty,left_x,left_y,left_z,right_x,right_y,right_z\n";
    out << std::setprecision(16);
    int pairCount = 0;
    for (int leftDof : leftDofs) {
        const Node& leftNode = mesh.nodes[static_cast<size_t>(leftDof)];
        const auto found = rightByCoordinate.find(coordinateMatchKey(leftNode.p, tolerance));
        if (found == rightByCoordinate.end()) {
            continue;
        }
        int bestRight = -1;
        double bestDistance = std::numeric_limits<double>::infinity();
        for (int rightDof : found->second) {
            const double distance = norm(leftNode.p - mesh.nodes[static_cast<size_t>(rightDof)].p);
            if (distance < bestDistance) {
                bestDistance = distance;
                bestRight = rightDof;
            }
        }
        if (bestRight < 0 || bestDistance > tolerance * 2.0) {
            continue;
        }
        stiffness.add(leftDof, leftDof, penalty);
        stiffness.add(leftDof, bestRight, -penalty);
        stiffness.add(bestRight, leftDof, -penalty);
        stiffness.add(bestRight, bestRight, penalty);
        const Node& rightNode = mesh.nodes[static_cast<size_t>(bestRight)];
        out << leftDof << ','
            << bestRight << ','
            << bestDistance << ','
            << penalty << ','
            << leftNode.p.x << ',' << leftNode.p.y << ',' << leftNode.p.z << ','
            << rightNode.p.x << ',' << rightNode.p.y << ',' << rightNode.p.z << '\n';
        ++pairCount;
    }
    return pairCount;
}

static std::vector<double> makeDirichletAdjustedSystem(const Mesh& mesh, SparseMatrix& system)
{
    std::vector<double> fixedAdjust(static_cast<size_t>(system.size()), 0.0);
    std::vector<char> isDirichlet(mesh.nodes.size(), 0);
    for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
        isDirichlet[static_cast<size_t>(i)] = mesh.nodes[static_cast<size_t>(i)].dirichlet ? 1 : 0;
    }

    std::vector<MatrixEntry> filtered;
    filtered.reserve(system.triplets.size());
    system.forEachEntry([&](int row, int col, double value) {
        if (row < 0 || row >= static_cast<int>(mesh.nodes.size())
            || col < 0 || col >= static_cast<int>(mesh.nodes.size())) {
            return;
        }
        if (isDirichlet[static_cast<size_t>(row)]) {
            return;
        }
        if (isDirichlet[static_cast<size_t>(col)]) {
            fixedAdjust[static_cast<size_t>(row)] += value * mesh.nodes[static_cast<size_t>(col)].dirichletValue;
            return;
        }
        filtered.push_back({row, col, value});
    });

    for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
        if (isDirichlet[static_cast<size_t>(i)]) {
            filtered.push_back({i, i, 1.0});
        }
    }
    system.triplets = std::move(filtered);
    system.rowPtr.clear();
    system.colInd.clear();
    system.values.clear();
    system.csrReady = false;
    return fixedAdjust;
}

static void applyDirichletRhs(const Mesh& mesh, const std::vector<double>& fixedAdjust, std::vector<double>& rhs)
{
    if (fixedAdjust.empty()) {
        return;
    }
    parallelFor(mesh.nodes.size(), [&](size_t i) {
        const Node& node = mesh.nodes[static_cast<size_t>(i)];
        if (node.dirichlet) {
            rhs[static_cast<size_t>(i)] = node.dirichletValue;
        } else {
            rhs[static_cast<size_t>(i)] -= fixedAdjust[static_cast<size_t>(i)];
        }
    });
}

static std::vector<MatrixEntry> sparseMatrixEntries(const SparseMatrix& a)
{
    if (a.csrReady) {
        std::vector<MatrixEntry> entries;
        entries.reserve(a.values.size());
        for (int row = 0; row < a.size(); ++row) {
            for (int k = a.rowPtr[static_cast<size_t>(row)]; k < a.rowPtr[static_cast<size_t>(row + 1)]; ++k) {
                entries.push_back({row, a.colInd[static_cast<size_t>(k)], a.values[static_cast<size_t>(k)]});
            }
        }
        return entries;
    }

    std::vector<MatrixEntry> entries;
    entries.reserve(a.triplets.size());
    a.forEachEntry([&](int row, int col, double value) {
        entries.push_back({row, col, value});
    });
    return entries;
}

class SubdomainDirectSolver {
public:
    SubdomainDirectSolver() = default;
    SubdomainDirectSolver(int n, const std::vector<MatrixEntry>& entries)
    {
        factor(n, entries);
    }

    SubdomainDirectSolver(const SubdomainDirectSolver&) = delete;
    SubdomainDirectSolver& operator=(const SubdomainDirectSolver&) = delete;

    SubdomainDirectSolver(SubdomainDirectSolver&& other) noexcept
    {
        moveFrom(std::move(other));
    }

    SubdomainDirectSolver& operator=(SubdomainDirectSolver&& other) noexcept
    {
        if (this != &other) {
            release();
            moveFrom(std::move(other));
        }
        return *this;
    }

    ~SubdomainDirectSolver()
    {
        release();
    }

    void solve(const std::vector<double>& rhs, std::vector<double>& x)
    {
#ifdef USE_MKL_PARDISO
        if (!numericallyFactorized_) {
            throw std::runtime_error("PARDISO solve requested before numerical factorization.");
        }
        x.assign(static_cast<size_t>(n_), 0.0);
        MKL_INT phase = 33;
        MKL_INT nrhs = 1;
        MKL_INT error = 0;
        pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                &nrhs, iparm_, &msglvl_, const_cast<double*>(rhs.data()), x.data(), &error);
        if (error != 0) {
            throw std::runtime_error("PARDISO solve failed with error " + std::to_string(error));
        }
#else
        solveDenseCholesky(rhs, x);
#endif
    }

    double symbolicAnalysisSeconds() const
    {
        return symbolicAnalysisSeconds_;
    }

    double numericalFactorizationSeconds() const
    {
        return numericalFactorizationSeconds_;
    }

    int symbolicAnalysisCalls() const
    {
        return symbolicAnalysisCalls_;
    }

    int numericalFactorizationCalls() const
    {
        return numericalFactorizationCalls_;
    }

    size_t memoryBytes() const
    {
        size_t bytes = values_.size() * sizeof(double)
            + rowPtr_.size() * sizeof(Index)
            + colInd_.size() * sizeof(Index)
            + perm_.size() * sizeof(Index);
#ifndef USE_MKL_PARDISO
        bytes += denseFactor_.size() * sizeof(double);
#else
        bytes += static_cast<size_t>(std::max<Index>(0, pardisoMemoryKb_)) * 1024;
#endif
        return bytes;
    }

private:
#ifdef USE_MKL_PARDISO
    using Index = MKL_INT;
#else
    using Index = int;
#endif

    Index n_ = 0;
    std::vector<Index> rowPtr_;
    std::vector<Index> colInd_;
    std::vector<double> values_;
    std::vector<Index> perm_;
    std::vector<double> denseFactor_;
    double symbolicAnalysisSeconds_ = 0.0;
    double numericalFactorizationSeconds_ = 0.0;
    int symbolicAnalysisCalls_ = 0;
    int numericalFactorizationCalls_ = 0;
    bool numericallyFactorized_ = false;

#ifdef USE_MKL_PARDISO
    void* pt_[64] = {};
    MKL_INT iparm_[64] = {};
    MKL_INT maxfct_ = 1;
    MKL_INT mnum_ = 1;
    MKL_INT mtype_ = 2;
    MKL_INT msglvl_ = 0;
    MKL_INT pardisoMemoryKb_ = 0;
#endif

    void factor(int n, const std::vector<MatrixEntry>& entries)
    {
        n_ = static_cast<Index>(n);
        perm_.assign(static_cast<size_t>(n_), 0);
#ifdef USE_MKL_PARDISO
        buildSymmetricUpperCsr(n, entries);
        for (int i = 0; i < 64; ++i) {
            pt_[i] = nullptr;
            iparm_[i] = 0;
        }
        iparm_[0] = 1;     // No solver default.
        iparm_[1] = 2;     // Fill-in reducing ordering: METIS.
        iparm_[7] = 2;     // Iterative refinement steps.
        iparm_[9] = 13;    // Pivot perturbation.
        iparm_[10] = 1;    // Scaling.
        iparm_[12] = 1;    // Matching.
        iparm_[17] = -1;   // Report nonzeros in factors.
        iparm_[18] = -1;   // Report MFLOPS.

        MKL_INT nrhs = 1;
        MKL_INT error = 0;
        std::vector<double> dummyRhs(static_cast<size_t>(n_), 0.0);
        std::vector<double> dummyX(static_cast<size_t>(n_), 0.0);

        // Keep the PARDISO lifecycle explicit.  Phase 11 builds and stores the
        // symbolic analysis exactly once for this immutable local matrix.
        MKL_INT phase = 11;
        const auto symbolicStart = std::chrono::steady_clock::now();
        pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                &nrhs, iparm_, &msglvl_, dummyRhs.data(), dummyX.data(), &error);
        symbolicAnalysisSeconds_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - symbolicStart).count();
        ++symbolicAnalysisCalls_;
        if (error != 0) {
            throw std::runtime_error("PARDISO symbolic analysis failed with error " + std::to_string(error));
        }

        // Phase 22 performs the numerical factorization once.  Subsequent
        // right-hand sides go directly to phase 33 in solve(); there is no
        // analysis or refactorization path while this matrix is unchanged.
        phase = 22;
        error = 0;
        const auto numericalStart = std::chrono::steady_clock::now();
        pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                &nrhs, iparm_, &msglvl_, dummyRhs.data(), dummyX.data(), &error);
        numericalFactorizationSeconds_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - numericalStart).count();
        ++numericalFactorizationCalls_;
        if (error != 0) {
            throw std::runtime_error("PARDISO numerical factorization failed with error " + std::to_string(error));
        }
        numericallyFactorized_ = true;
        pardisoMemoryKb_ = std::max<MKL_INT>(0, iparm_[14] + iparm_[15] + iparm_[16]);
#else
        const auto symbolicStart = std::chrono::steady_clock::now();
        ++symbolicAnalysisCalls_;
        symbolicAnalysisSeconds_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - symbolicStart).count();
        const auto numericalStart = std::chrono::steady_clock::now();
        buildDenseCholesky(n, entries);
        numericalFactorizationSeconds_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - numericalStart).count();
        ++numericalFactorizationCalls_;
        numericallyFactorized_ = true;
#endif
    }

    void buildSymmetricUpperCsr(int n, const std::vector<MatrixEntry>& entries)
    {
        std::vector<std::map<int, double>> rows(static_cast<size_t>(n));
        for (const MatrixEntry& e : entries) {
            if (e.col < e.row) {
                continue;
            }
            rows[static_cast<size_t>(e.row)][e.col] += e.value;
        }
        rowPtr_.assign(static_cast<size_t>(n + 1), 1);
        for (int i = 0; i < n; ++i) {
            rowPtr_[static_cast<size_t>(i + 1)] = rowPtr_[static_cast<size_t>(i)] + static_cast<Index>(rows[static_cast<size_t>(i)].size());
            for (const auto& entry : rows[static_cast<size_t>(i)]) {
                colInd_.push_back(static_cast<Index>(entry.first + 1));
                values_.push_back(entry.second);
            }
        }
    }

#ifndef USE_MKL_PARDISO
    void buildDenseCholesky(int n, const std::vector<MatrixEntry>& entries)
    {
        denseFactor_.assign(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
        for (const MatrixEntry& e : entries) {
            denseFactor_[static_cast<size_t>(e.row) * static_cast<size_t>(n) + static_cast<size_t>(e.col)] += e.value;
        }
        for (int k = 0; k < n; ++k) {
            double diag = denseFactor_[static_cast<size_t>(k) * static_cast<size_t>(n) + static_cast<size_t>(k)];
            for (int j = 0; j < k; ++j) {
                const double lkj = denseFactor_[static_cast<size_t>(k) * static_cast<size_t>(n) + static_cast<size_t>(j)];
                diag -= lkj * lkj;
            }
            if (diag <= 0.0) {
                throw std::runtime_error("Dense fallback Cholesky detected a non-SPD block.");
            }
            denseFactor_[static_cast<size_t>(k) * static_cast<size_t>(n) + static_cast<size_t>(k)] = std::sqrt(diag);
            for (int i = k + 1; i < n; ++i) {
                double value = denseFactor_[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(k)];
                for (int j = 0; j < k; ++j) {
                    value -= denseFactor_[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(j)]
                           * denseFactor_[static_cast<size_t>(k) * static_cast<size_t>(n) + static_cast<size_t>(j)];
                }
                denseFactor_[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(k)] =
                    value / denseFactor_[static_cast<size_t>(k) * static_cast<size_t>(n) + static_cast<size_t>(k)];
            }
        }
    }

    void solveDenseCholesky(const std::vector<double>& rhs, std::vector<double>& x) const
    {
        const int n = static_cast<int>(n_);
        std::vector<double> y(static_cast<size_t>(n), 0.0);
        for (int i = 0; i < n; ++i) {
            double sum = rhs[static_cast<size_t>(i)];
            for (int j = 0; j < i; ++j) {
                sum -= denseFactor_[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(j)] * y[static_cast<size_t>(j)];
            }
            y[static_cast<size_t>(i)] = sum / denseFactor_[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(i)];
        }
        x.assign(static_cast<size_t>(n), 0.0);
        for (int i = n - 1; i >= 0; --i) {
            double sum = y[static_cast<size_t>(i)];
            for (int j = i + 1; j < n; ++j) {
                sum -= denseFactor_[static_cast<size_t>(j) * static_cast<size_t>(n) + static_cast<size_t>(i)] * x[static_cast<size_t>(j)];
            }
            x[static_cast<size_t>(i)] = sum / denseFactor_[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(i)];
        }
    }
#endif

    void release()
    {
#ifdef USE_MKL_PARDISO
        if (n_ > 0) {
            MKL_INT phase = -1;
            MKL_INT nrhs = 1;
            MKL_INT error = 0;
            std::vector<double> dummyRhs(static_cast<size_t>(n_), 0.0);
            std::vector<double> dummyX(static_cast<size_t>(n_), 0.0);
            pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                    values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                    &nrhs, iparm_, &msglvl_, dummyRhs.data(), dummyX.data(), &error);
        }
#endif
        n_ = 0;
        numericallyFactorized_ = false;
    }

    void moveFrom(SubdomainDirectSolver&& other) noexcept
    {
        n_ = other.n_;
        rowPtr_ = std::move(other.rowPtr_);
        colInd_ = std::move(other.colInd_);
        values_ = std::move(other.values_);
        perm_ = std::move(other.perm_);
        denseFactor_ = std::move(other.denseFactor_);
        symbolicAnalysisSeconds_ = other.symbolicAnalysisSeconds_;
        numericalFactorizationSeconds_ = other.numericalFactorizationSeconds_;
        symbolicAnalysisCalls_ = other.symbolicAnalysisCalls_;
        numericalFactorizationCalls_ = other.numericalFactorizationCalls_;
        numericallyFactorized_ = other.numericallyFactorized_;
#ifdef USE_MKL_PARDISO
        for (int i = 0; i < 64; ++i) {
            pt_[i] = other.pt_[i];
            other.pt_[i] = nullptr;
            iparm_[i] = other.iparm_[i];
        }
        maxfct_ = other.maxfct_;
        mnum_ = other.mnum_;
        mtype_ = other.mtype_;
        msglvl_ = other.msglvl_;
        pardisoMemoryKb_ = other.pardisoMemoryKb_;
#endif
        other.n_ = 0;
        other.numericallyFactorized_ = false;
    }
};

class GeneralSparseDirectSolver {
public:
    GeneralSparseDirectSolver() = default;
    GeneralSparseDirectSolver(int n, const std::vector<MatrixEntry>& entries)
    {
        factor(n, entries);
    }

    GeneralSparseDirectSolver(const GeneralSparseDirectSolver&) = delete;
    GeneralSparseDirectSolver& operator=(const GeneralSparseDirectSolver&) = delete;

    GeneralSparseDirectSolver(GeneralSparseDirectSolver&& other) noexcept
    {
        moveFrom(std::move(other));
    }

    GeneralSparseDirectSolver& operator=(GeneralSparseDirectSolver&& other) noexcept
    {
        if (this != &other) {
            release();
            moveFrom(std::move(other));
        }
        return *this;
    }

    ~GeneralSparseDirectSolver()
    {
        release();
    }

    void solve(const std::vector<double>& rhs, std::vector<double>& x)
    {
#ifdef USE_MKL_PARDISO
        x.assign(static_cast<size_t>(n_), 0.0);
        MKL_INT phase = 33;
        MKL_INT nrhs = 1;
        MKL_INT error = 0;
        pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                &nrhs, iparm_, &msglvl_, const_cast<double*>(rhs.data()), x.data(), &error);
        if (error != 0) {
            throw std::runtime_error("General PARDISO solve failed with error " + std::to_string(error));
        }
#else
        (void)rhs;
        (void)x;
        throw std::runtime_error("General sparse direct solve requires MKL PARDISO.");
#endif
    }

    size_t memoryBytes() const
    {
        size_t bytes = values_.size() * sizeof(double)
            + rowPtr_.size() * sizeof(Index)
            + colInd_.size() * sizeof(Index)
            + perm_.size() * sizeof(Index);
#ifdef USE_MKL_PARDISO
        bytes += static_cast<size_t>(std::max<Index>(0, pardisoMemoryKb_)) * 1024;
#endif
        return bytes;
    }

private:
#ifdef USE_MKL_PARDISO
    using Index = MKL_INT;
#else
    using Index = int;
#endif

    Index n_ = 0;
    std::vector<Index> rowPtr_;
    std::vector<Index> colInd_;
    std::vector<double> values_;
    std::vector<Index> perm_;

#ifdef USE_MKL_PARDISO
    void* pt_[64] = {};
    MKL_INT iparm_[64] = {};
    MKL_INT maxfct_ = 1;
    MKL_INT mnum_ = 1;
    MKL_INT mtype_ = 11;
    MKL_INT msglvl_ = 0;
    MKL_INT pardisoMemoryKb_ = 0;
#endif

    void factor(int n, const std::vector<MatrixEntry>& entries)
    {
        n_ = static_cast<Index>(n);
        perm_.assign(static_cast<size_t>(n_), 0);
#ifdef USE_MKL_PARDISO
        buildGeneralCsr(n, entries);
        for (int i = 0; i < 64; ++i) {
            pt_[i] = nullptr;
            iparm_[i] = 0;
        }
        iparm_[0] = 1;
        iparm_[1] = 2;
        iparm_[7] = 2;
        iparm_[9] = 13;
        iparm_[10] = 1;
        iparm_[12] = 1;
        iparm_[17] = -1;
        iparm_[18] = -1;

        MKL_INT phase = 12;
        MKL_INT nrhs = 1;
        MKL_INT error = 0;
        std::vector<double> dummyRhs(static_cast<size_t>(n_), 0.0);
        std::vector<double> dummyX(static_cast<size_t>(n_), 0.0);
        pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                &nrhs, iparm_, &msglvl_, dummyRhs.data(), dummyX.data(), &error);
        if (error != 0) {
            throw std::runtime_error("General PARDISO factorization failed with error " + std::to_string(error));
        }
        pardisoMemoryKb_ = std::max<MKL_INT>(0, iparm_[14] + iparm_[15] + iparm_[16]);
#else
        (void)entries;
        throw std::runtime_error("General sparse direct factorization requires MKL PARDISO.");
#endif
    }

    void buildGeneralCsr(int n, const std::vector<MatrixEntry>& entries)
    {
        std::vector<std::map<int, double>> rows(static_cast<size_t>(n));
        for (const MatrixEntry& e : entries) {
            if (e.row < 0 || e.row >= n || e.col < 0 || e.col >= n) {
                continue;
            }
            rows[static_cast<size_t>(e.row)][e.col] += e.value;
        }
        rowPtr_.assign(static_cast<size_t>(n + 1), 1);
        for (int i = 0; i < n; ++i) {
            for (const auto& entry : rows[static_cast<size_t>(i)]) {
                if (std::abs(entry.second) == 0.0) {
                    continue;
                }
                colInd_.push_back(static_cast<Index>(entry.first + 1));
                values_.push_back(entry.second);
            }
            rowPtr_[static_cast<size_t>(i + 1)] = static_cast<Index>(colInd_.size() + 1);
        }
    }

    void release()
    {
#ifdef USE_MKL_PARDISO
        if (n_ > 0) {
            MKL_INT phase = -1;
            MKL_INT nrhs = 1;
            MKL_INT error = 0;
            std::vector<double> dummyRhs(static_cast<size_t>(n_), 0.0);
            std::vector<double> dummyX(static_cast<size_t>(n_), 0.0);
            pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                    values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                    &nrhs, iparm_, &msglvl_, dummyRhs.data(), dummyX.data(), &error);
        }
#endif
        n_ = 0;
    }

    void moveFrom(GeneralSparseDirectSolver&& other) noexcept
    {
        n_ = other.n_;
        rowPtr_ = std::move(other.rowPtr_);
        colInd_ = std::move(other.colInd_);
        values_ = std::move(other.values_);
        perm_ = std::move(other.perm_);
#ifdef USE_MKL_PARDISO
        for (int i = 0; i < 64; ++i) {
            pt_[i] = other.pt_[i];
            other.pt_[i] = nullptr;
            iparm_[i] = other.iparm_[i];
        }
        maxfct_ = other.maxfct_;
        mnum_ = other.mnum_;
        mtype_ = other.mtype_;
        msglvl_ = other.msglvl_;
        pardisoMemoryKb_ = other.pardisoMemoryKb_;
#endif
        other.n_ = 0;
    }
};

static bool confirmSpdWithPardiso(const SparseMatrix& system, std::string& errorMessage)
{
    try {
        const std::vector<MatrixEntry> entries = sparseMatrixEntries(system);
        SubdomainDirectSolver solver(system.size(), entries);
        errorMessage.clear();
        return true;
    } catch (const std::exception& ex) {
        errorMessage = ex.what();
        return false;
    }
}

static bool confirmGeneralWithPardiso(const SparseMatrix& system, std::string& errorMessage)
{
    try {
        const std::vector<MatrixEntry> entries = sparseMatrixEntries(system);
        GeneralSparseDirectSolver solver(system.size(), entries);
        errorMessage.clear();
        return true;
    } catch (const std::exception& ex) {
        errorMessage = ex.what();
        return false;
    }
}

struct PardisoInertiaDiagnostic {
    std::string status = "not_available";
    int positivePivots = -1;
    int negativePivots = -1;
    int zeroTinyPivots = -1;
    std::string message;
};

static PardisoInertiaDiagnostic computePardisoInertia(const SparseMatrix& system)
{
    PardisoInertiaDiagnostic result;
#ifdef USE_MKL_PARDISO
    using Index = MKL_INT;
    try {
        const int n = system.size();
        std::vector<std::map<int, double>> rows(static_cast<size_t>(n));
        system.forEachEntry([&](int row, int col, double value) {
            if (row < 0 || row >= n || col < 0 || col >= n || col < row) {
                return;
            }
            rows[static_cast<size_t>(row)][col] += value;
        });

        std::vector<Index> rowPtr(static_cast<size_t>(n + 1), 1);
        std::vector<Index> colInd;
        std::vector<double> values;
        for (int i = 0; i < n; ++i) {
            for (const auto& entry : rows[static_cast<size_t>(i)]) {
                if (std::abs(entry.second) == 0.0) {
                    continue;
                }
                colInd.push_back(static_cast<Index>(entry.first + 1));
                values.push_back(entry.second);
            }
            rowPtr[static_cast<size_t>(i + 1)] = static_cast<Index>(colInd.size() + 1);
        }
        std::vector<std::map<int, double>>().swap(rows);

        void* pt[64] = {};
        MKL_INT iparm[64] = {};
        MKL_INT maxfct = 1;
        MKL_INT mnum = 1;
        MKL_INT mtype = -2;
        MKL_INT phase = 12;
        MKL_INT nrhs = 1;
        MKL_INT msglvl = 0;
        MKL_INT error = 0;
        MKL_INT nMkl = static_cast<MKL_INT>(n);
        std::vector<Index> perm(static_cast<size_t>(n), 0);
        std::vector<double> dummyRhs(static_cast<size_t>(n), 0.0);
        std::vector<double> dummyX(static_cast<size_t>(n), 0.0);
        iparm[0] = 1;
        iparm[1] = 2;
        iparm[7] = 2;
        iparm[9] = 13;
        iparm[10] = 1;
        iparm[12] = 1;
        iparm[17] = -1;
        iparm[18] = -1;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nMkl,
                values.data(), rowPtr.data(), colInd.data(), perm.data(),
                &nrhs, iparm, &msglvl, dummyRhs.data(), dummyX.data(), &error);
        if (error == 0) {
            result.status = "success";
            result.positivePivots = static_cast<int>(iparm[21]);
            result.negativePivots = static_cast<int>(iparm[22]);
            result.zeroTinyPivots = n - result.positivePivots - result.negativePivots;
            result.message.clear();
        } else {
            result.status = "failure";
            result.message = "PARDISO LDLT inertia failed with error " + std::to_string(error);
        }
        phase = -1;
        MKL_INT releaseError = 0;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nMkl,
                values.data(), rowPtr.data(), colInd.data(), perm.data(),
                &nrhs, iparm, &msglvl, dummyRhs.data(), dummyX.data(), &releaseError);
    } catch (const std::exception& ex) {
        result.status = "failure";
        result.message = ex.what();
    }
#else
    result.message = "MKL PARDISO is not enabled";
#endif
    return result;
}

class SubdomainIcSolver {
public:
    struct Workspace {
        std::vector<double> scaledRhs;
        std::vector<double> y;
        std::vector<double> xPerm;
    };

    SubdomainIcSolver() = default;
    SubdomainIcSolver(int n,
                      const std::vector<MatrixEntry>& entries,
                      double requestedShift,
                      int subdomain,
                      bool diagonalScaling,
                      double diagScalingEps,
                      const std::string& shiftMode,
                      const std::string& ordering)
    {
        factor(n, entries, requestedShift, subdomain, diagonalScaling, diagScalingEps, shiftMode, ordering);
    }

    void solve(const std::vector<double>& rhs, std::vector<double>& x) const
    {
        Workspace workspace;
        solve(rhs, x, workspace);
    }

    void solve(const std::vector<double>& rhs, std::vector<double>& x, Workspace& workspace) const
    {
        workspace.scaledRhs.resize(static_cast<size_t>(n_));
        for (int i = 0; i < n_; ++i) {
            const int oldIndex = permNewToOld_[static_cast<size_t>(i)];
            workspace.scaledRhs[static_cast<size_t>(i)] =
                scaling_[static_cast<size_t>(i)] * rhs[static_cast<size_t>(oldIndex)];
        }

        workspace.y.resize(static_cast<size_t>(n_));
        for (int i = 0; i < n_; ++i) {
            double sum = workspace.scaledRhs[static_cast<size_t>(i)];
            for (const auto& entry : lowerRows_[static_cast<size_t>(i)]) {
                sum -= entry.second * workspace.y[static_cast<size_t>(entry.first)];
            }
            workspace.y[static_cast<size_t>(i)] = sum / diagonal_[static_cast<size_t>(i)];
        }

        workspace.xPerm.resize(static_cast<size_t>(n_));
        for (int i = n_ - 1; i >= 0; --i) {
            double sum = workspace.y[static_cast<size_t>(i)];
            for (const auto& entry : lowerColumns_[static_cast<size_t>(i)]) {
                sum -= entry.second * workspace.xPerm[static_cast<size_t>(entry.first)];
            }
            workspace.xPerm[static_cast<size_t>(i)] = sum / diagonal_[static_cast<size_t>(i)];
        }

        x.resize(static_cast<size_t>(n_));
        for (int i = 0; i < n_; ++i) {
            const int oldIndex = permNewToOld_[static_cast<size_t>(i)];
            x[static_cast<size_t>(oldIndex)] =
                scaling_[static_cast<size_t>(i)] * workspace.xPerm[static_cast<size_t>(i)];
        }
    }

    size_t memoryBytes() const
    {
        size_t bytes = diagonal_.size() * sizeof(double);
        bytes += scaling_.size() * sizeof(double);
        for (const auto& row : lowerRows_) {
            bytes += row.size() * (sizeof(int) + sizeof(double));
        }
        for (const auto& col : lowerColumns_) {
            bytes += col.size() * (sizeof(int) + sizeof(double));
        }
        return bytes;
    }

    const IcFactorDiagnostics& diagnostics() const { return diagnostics_; }
    const std::vector<IcFactorDiagnostics>& trialDiagnostics() const { return trialDiagnostics_; }

private:
    int n_ = 0;
    std::vector<double> diagonal_;
    std::vector<double> scaling_;
    std::vector<int> permNewToOld_;
    std::vector<int> permOldToNew_;
    std::vector<std::vector<std::pair<int, double>>> lowerRows_;
    std::vector<std::vector<std::pair<int, double>>> lowerColumns_;
    IcFactorDiagnostics diagnostics_;
    std::vector<IcFactorDiagnostics> trialDiagnostics_;

    static std::vector<double> shiftCandidates(double requestedShift, const std::string& shiftMode)
    {
        const double boundedShift = std::max(0.0, requestedShift);
        if (shiftMode == "none") {
            return {0.0};
        }
        if (shiftMode == "fixed") {
            return {boundedShift};
        }

        std::vector<double> values{
            boundedShift,
            0.0,
            1.0e-12,
            1.0e-10,
            1.0e-8,
            1.0e-6,
            1.0e-4,
            1.0e-2,
            1.0e-1,
            1.0
        };
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end(), [](double a, double b) {
            return std::abs(a - b) <= 1.0e-30 * std::max(1.0, std::max(std::abs(a), std::abs(b)));
        }), values.end());
        return values;
    }

    static std::vector<int> naturalOrdering(int n)
    {
        std::vector<int> order(static_cast<size_t>(n));
        std::iota(order.begin(), order.end(), 0);
        return order;
    }

    static std::vector<std::vector<int>> buildAdjacency(int n, const std::vector<MatrixEntry>& entries)
    {
        std::vector<std::vector<int>> adjacency(static_cast<size_t>(n));
        for (const MatrixEntry& e : entries) {
            if (e.row < 0 || e.row >= n || e.col < 0 || e.col >= n || e.row == e.col) {
                continue;
            }
            adjacency[static_cast<size_t>(e.row)].push_back(e.col);
            adjacency[static_cast<size_t>(e.col)].push_back(e.row);
        }
        for (auto& row : adjacency) {
            std::sort(row.begin(), row.end());
            row.erase(std::unique(row.begin(), row.end()), row.end());
        }
        return adjacency;
    }

    static std::vector<int> rcmOrdering(int n, const std::vector<MatrixEntry>& entries)
    {
        const auto adjacency = buildAdjacency(n, entries);
        std::vector<char> visited(static_cast<size_t>(n), 0);
        std::vector<int> order;
        order.reserve(static_cast<size_t>(n));
        while (static_cast<int>(order.size()) < n) {
            int start = -1;
            size_t bestDegree = std::numeric_limits<size_t>::max();
            for (int i = 0; i < n; ++i) {
                if (!visited[static_cast<size_t>(i)]
                    && adjacency[static_cast<size_t>(i)].size() < bestDegree) {
                    start = i;
                    bestDegree = adjacency[static_cast<size_t>(i)].size();
                }
            }
            if (start < 0) {
                break;
            }
            std::vector<int> queue{start};
            visited[static_cast<size_t>(start)] = 1;
            for (size_t head = 0; head < queue.size(); ++head) {
                const int node = queue[head];
                order.push_back(node);
                std::vector<int> neighbors = adjacency[static_cast<size_t>(node)];
                std::sort(neighbors.begin(), neighbors.end(), [&](int a, int b) {
                    return adjacency[static_cast<size_t>(a)].size() < adjacency[static_cast<size_t>(b)].size();
                });
                for (int nb : neighbors) {
                    if (!visited[static_cast<size_t>(nb)]) {
                        visited[static_cast<size_t>(nb)] = 1;
                        queue.push_back(nb);
                    }
                }
            }
        }
        std::reverse(order.begin(), order.end());
        return order;
    }

    static std::vector<int> degreeOrdering(int n, const std::vector<MatrixEntry>& entries)
    {
        const auto adjacency = buildAdjacency(n, entries);
        std::vector<int> order = naturalOrdering(n);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            return adjacency[static_cast<size_t>(a)].size() < adjacency[static_cast<size_t>(b)].size();
        });
        return order;
    }

    static std::vector<int> makeOrdering(int n, const std::vector<MatrixEntry>& entries, const std::string& ordering)
    {
        if (ordering == "rcm") {
            return rcmOrdering(n, entries);
        }
        if (ordering == "amd") {
            return degreeOrdering(n, entries);
        }
        return naturalOrdering(n);
    }

    void setPermutation(const std::vector<int>& order)
    {
        permNewToOld_ = order;
        permOldToNew_.assign(static_cast<size_t>(n_), 0);
        for (int i = 0; i < n_; ++i) {
            permOldToNew_[static_cast<size_t>(permNewToOld_[static_cast<size_t>(i)])] = i;
        }
    }

    bool factorWithShift(const std::vector<MatrixEntry>& entries,
                         const std::vector<double>& rawDiagonal,
                         double meanDiagonal,
                         double shift,
                         bool useDiagonalScaling,
                         double diagScalingEps,
                         IcFactorDiagnostics& diagnostics)
    {
        std::vector<std::map<int, double>> pattern(static_cast<size_t>(n_));
        scaling_.assign(static_cast<size_t>(n_), 1.0);
        const double absoluteFloor = std::max(1.0e-300, diagScalingEps);
        constexpr double relativeFloor = 1.0e-12;
        const double pivotFloor = std::max(absoluteFloor, relativeFloor * meanDiagonal);

        std::vector<double> shiftedDiagonal = rawDiagonal;
        for (int i = 0; i < n_; ++i) {
            const int oldIndex = permNewToOld_[static_cast<size_t>(i)];
            shiftedDiagonal[static_cast<size_t>(i)] =
                rawDiagonal[static_cast<size_t>(oldIndex)] + shift * meanDiagonal;
            if (useDiagonalScaling) {
                const double scaleDiag = std::max(pivotFloor, std::abs(shiftedDiagonal[static_cast<size_t>(i)]));
                scaling_[static_cast<size_t>(i)] = 1.0 / std::sqrt(scaleDiag);
            }
        }

        for (const MatrixEntry& entry : entries) {
            if (entry.row < 0 || entry.row >= n_ || entry.col < 0 || entry.col >= n_) {
                continue;
            }
            const int row = permOldToNew_[static_cast<size_t>(entry.row)];
            const int col = permOldToNew_[static_cast<size_t>(entry.col)];
            const double scaledValue = entry.value
                * scaling_[static_cast<size_t>(row)]
                * scaling_[static_cast<size_t>(col)];
            if (row == col) {
                pattern[static_cast<size_t>(row)][col] += shift * meanDiagonal
                    * scaling_[static_cast<size_t>(row)]
                    * scaling_[static_cast<size_t>(col)];
            }
            if (col <= row) {
                pattern[static_cast<size_t>(row)][col] += scaledValue;
            }
        }

        diagonal_.assign(static_cast<size_t>(n_), 1.0);
        lowerRows_.assign(static_cast<size_t>(n_), {});
        std::vector<std::unordered_map<int, size_t>> position(static_cast<size_t>(n_));

        for (int i = 0; i < n_; ++i) {
            auto diagIt = pattern[static_cast<size_t>(i)].find(i);
            const double diagValue = diagIt == pattern[static_cast<size_t>(i)].end() ? 1.0 : diagIt->second;
            diagonal_[static_cast<size_t>(i)] = diagValue;
            for (const auto& entry : pattern[static_cast<size_t>(i)]) {
                if (entry.first >= i) {
                    continue;
                }
                position[static_cast<size_t>(i)][entry.first] = lowerRows_[static_cast<size_t>(i)].size();
                lowerRows_[static_cast<size_t>(i)].push_back(entry);
            }
        }

        for (int i = 0; i < n_; ++i) {
            std::vector<std::pair<int, double>>& row = lowerRows_[static_cast<size_t>(i)];
            for (size_t idx = 0; idx < row.size(); ++idx) {
                const int j = row[idx].first;
                double sum = row[idx].second;
                for (size_t kIdx = 0; kIdx < idx; ++kIdx) {
                    const int k = row[kIdx].first;
                    if (k >= j) {
                        break;
                    }
                    const auto found = position[static_cast<size_t>(j)].find(k);
                    if (found != position[static_cast<size_t>(j)].end()) {
                        sum -= row[kIdx].second
                             * lowerRows_[static_cast<size_t>(j)][found->second].second;
                    }
                }
                const double pivot = std::max(pivotFloor, diagonal_[static_cast<size_t>(j)]);
                row[idx].second = sum / pivot;
                if (!std::isfinite(row[idx].second)) {
                    ++diagnostics.nonFiniteLCount;
                }
            }

            double diag = diagonal_[static_cast<size_t>(i)];
            for (const auto& entry : row) {
                diag -= entry.second * entry.second;
            }
            const double rawPivot = diag;
            if (!std::isfinite(rawPivot)) {
                ++diagnostics.nonFiniteLCount;
            }
            if (rawPivot <= 0.0 || !std::isfinite(rawPivot)) {
                ++diagnostics.pivotNonpositiveCount;
                if (diagnostics.firstBadPivotRow < 0) {
                    diagnostics.firstBadPivotRow = i;
                    diagnostics.firstBadPivotValue = rawPivot;
                }
            }
            if (!std::isfinite(diag) || diag <= pivotFloor) {
                ++diagnostics.pivotTinyCount;
                diag = pivotFloor;
            }
            diagonal_[static_cast<size_t>(i)] = std::sqrt(diag);
            diagnostics.pivotMin = std::min(diagnostics.pivotMin, diag);
            diagnostics.pivotMax = std::max(diagnostics.pivotMax, diag);
            if (!std::isfinite(diagonal_[static_cast<size_t>(i)])) {
                ++diagnostics.nonFiniteLCount;
            }
        }

        lowerColumns_.assign(static_cast<size_t>(n_), {});
        for (int row = 0; row < n_; ++row) {
            for (const auto& entry : lowerRows_[static_cast<size_t>(row)]) {
                lowerColumns_[static_cast<size_t>(entry.first)].push_back({row, entry.second});
            }
        }
        diagnostics.breakdown = diagnostics.nonFiniteLCount > 0 || diagnostics.pivotNonpositiveCount > 0;
        diagnostics.accepted = !diagnostics.breakdown;
        return diagnostics.accepted;
    }

    void factor(int n,
                const std::vector<MatrixEntry>& entries,
                double requestedShift,
                int subdomain,
                bool diagonalScaling,
                double diagScalingEps,
                const std::string& shiftMode,
                const std::string& ordering)
    {
        n_ = n;
        setPermutation(makeOrdering(n_, entries, ordering));
        std::vector<double> rawDiagonal(static_cast<size_t>(n_), 0.0);
        for (const MatrixEntry& entry : entries) {
            if (entry.row == entry.col && entry.row >= 0 && entry.row < n_) {
                rawDiagonal[static_cast<size_t>(entry.row)] += entry.value;
            }
        }
        double diagSum = 0.0;
        int diagCount = 0;
        for (double value : rawDiagonal) {
            if (std::isfinite(value) && value > 0.0) {
                diagSum += value;
                ++diagCount;
            }
        }
        const double meanDiagonal = diagCount > 0 ? diagSum / static_cast<double>(diagCount) : 1.0;

        IcFactorDiagnostics bestDiagnostics;
        bestDiagnostics.subdomain = subdomain;
        bestDiagnostics.dofs = n_;
        bestDiagnostics.ordering = ordering;
        bestDiagnostics.diagonalScaling = diagonalScaling;
        try {
            SubdomainDirectSolver exactSolver(n_, entries);
            bestDiagnostics.exactSpd = true;
        } catch (...) {
            bestDiagnostics.exactSpd = false;
        }

        bool accepted = false;
        const std::vector<double> candidates = shiftCandidates(requestedShift, shiftMode);
        for (double shift : candidates) {
            IcFactorDiagnostics trial;
            trial.subdomain = subdomain;
            trial.dofs = n_;
            trial.ordering = ordering;
            trial.appliedShift = shift;
            trial.diagonalScaling = diagonalScaling;
            trial.exactSpd = bestDiagnostics.exactSpd;
            if (factorWithShift(entries, rawDiagonal, meanDiagonal, shift, diagonalScaling, diagScalingEps, trial)) {
                diagnostics_ = trial;
                trialDiagnostics_.push_back(trial);
                accepted = true;
                break;
            }
            trialDiagnostics_.push_back(trial);
            bestDiagnostics = trial;
        }
        if (!accepted) {
            factorWithShift(entries,
                            rawDiagonal,
                            meanDiagonal,
                            candidates.empty() ? 0.0 : candidates.back(),
                            diagonalScaling,
                            diagScalingEps,
                            bestDiagnostics);
            diagnostics_ = bestDiagnostics;
        }
    }
};

class BlockJacobiPreconditioner {
public:
    BlockJacobiPreconditioner(const Mesh& mesh, const SparseMatrix& a)
    {
        build(mesh, a);
    }

    void apply(const std::vector<double>& r, std::vector<double>& z)
    {
        z.assign(r.size(), 0.0);
        parallelForCoarse(blockDofs_.size(), [&](size_t block) {
            ScopedMklSingleThread mklThreads;
            const std::vector<int>& dofs = blockDofs_[block];
            std::vector<double> rhs(dofs.size(), 0.0);
            for (size_t i = 0; i < dofs.size(); ++i) {
                rhs[i] = r[static_cast<size_t>(dofs[i])];
            }
            std::vector<double> localSolution;
            solvers_[block].solve(rhs, localSolution);
            for (size_t i = 0; i < dofs.size(); ++i) {
                z[static_cast<size_t>(dofs[i])] = localSolution[i];
            }
        });
    }

    size_t memoryBytes() const
    {
        size_t bytes = globalToBlock_.size() * sizeof(int)
            + globalToLocal_.size() * sizeof(int);
        for (const auto& block : blockDofs_) {
            bytes += block.size() * sizeof(int);
        }
        for (const auto& solver : solvers_) {
            bytes += solver.memoryBytes();
        }
        return bytes;
    }

private:
    std::vector<std::vector<int>> blockDofs_;
    std::vector<int> globalToBlock_;
    std::vector<int> globalToLocal_;
    std::vector<SubdomainDirectSolver> solvers_;

    void build(const Mesh& mesh, const SparseMatrix& a)
    {
        int maxSubdomain = 0;
        for (const Node& node : mesh.nodes) {
            maxSubdomain = std::max(maxSubdomain, node.subdomain);
        }
        blockDofs_.assign(static_cast<size_t>(maxSubdomain + 1), {});
        globalToBlock_.assign(mesh.nodes.size(), -1);
        globalToLocal_.assign(mesh.nodes.size(), -1);

        for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
            const int block = mesh.nodes[static_cast<size_t>(i)].subdomain;
            globalToBlock_[static_cast<size_t>(i)] = block;
            globalToLocal_[static_cast<size_t>(i)] = static_cast<int>(blockDofs_[static_cast<size_t>(block)].size());
            blockDofs_[static_cast<size_t>(block)].push_back(i);
        }

        solvers_.clear();
        solvers_.resize(blockDofs_.size());
        std::vector<std::string> factorErrors(blockDofs_.size());
        parallelForCoarse(blockDofs_.size(), [&](size_t blockIndex) {
            try {
                ScopedMklSingleThread mklThreads;
                const int block = static_cast<int>(blockIndex);
                std::vector<MatrixEntry> entries;
                for (int globalRow : blockDofs_[blockIndex]) {
                    const int localRow = globalToLocal_[static_cast<size_t>(globalRow)];
                    for (int k = a.rowPtr[static_cast<size_t>(globalRow)]; k < a.rowPtr[static_cast<size_t>(globalRow + 1)]; ++k) {
                        const int globalCol = a.colInd[static_cast<size_t>(k)];
                        if (globalCol < 0 || globalCol >= static_cast<int>(globalToBlock_.size())) {
                            continue;
                        }
                        if (globalToBlock_[static_cast<size_t>(globalCol)] != block) {
                            continue;
                        }
                        entries.push_back({localRow, globalToLocal_[static_cast<size_t>(globalCol)], a.values[static_cast<size_t>(k)]});
                    }
                }
                solvers_[blockIndex] = SubdomainDirectSolver(static_cast<int>(blockDofs_[blockIndex].size()), entries);
            } catch (const std::exception& err) {
                factorErrors[blockIndex] = err.what();
            }
        });
        for (size_t blockIndex = 0; blockIndex < factorErrors.size(); ++blockIndex) {
            if (!factorErrors[blockIndex].empty()) {
                throw std::runtime_error("Block " + std::to_string(blockIndex)
                    + " PARDISO factorization failed: " + factorErrors[blockIndex]);
            }
        }
    }
};

class BlockJacobiIcPreconditioner {
public:
    BlockJacobiIcPreconditioner(const Mesh& mesh,
                                const SparseMatrix& a,
                                double icShift,
                                bool diagonalScaling,
                                double diagScalingEps,
                                const std::string& shiftMode,
                                const std::string& ordering)
    {
        build(mesh, a, icShift, diagonalScaling, diagScalingEps, shiftMode, ordering);
    }

    void apply(const std::vector<double>& r, std::vector<double>& z)
    {
        z.resize(r.size());
        std::fill(z.begin(), z.end(), 0.0);
        parallelForCoarse(blockDofs_.size(), [&](size_t block) {
            const std::vector<int>& dofs = blockDofs_[block];
            std::vector<double>& rhs = rhsBuffers_[block];
            rhs.resize(dofs.size());
            for (size_t i = 0; i < dofs.size(); ++i) {
                rhs[i] = r[static_cast<size_t>(dofs[i])];
            }
            std::vector<double>& localSolution = solutionBuffers_[block];
            solvers_[block].solve(rhs, localSolution, workspaces_[block]);
            for (size_t i = 0; i < dofs.size(); ++i) {
                z[static_cast<size_t>(dofs[i])] = localSolution[i];
            }
        });
    }

    size_t memoryBytes() const
    {
        size_t bytes = globalToBlock_.size() * sizeof(int)
            + globalToLocal_.size() * sizeof(int);
        for (const auto& block : blockDofs_) {
            bytes += block.size() * sizeof(int);
        }
        for (const auto& buffer : rhsBuffers_) {
            bytes += buffer.size() * sizeof(double);
        }
        for (const auto& buffer : solutionBuffers_) {
            bytes += buffer.size() * sizeof(double);
        }
        for (const auto& workspace : workspaces_) {
            bytes += workspace.scaledRhs.size() * sizeof(double)
                + workspace.y.size() * sizeof(double)
                + workspace.xPerm.size() * sizeof(double);
        }
        for (const auto& solver : solvers_) {
            bytes += solver.memoryBytes();
        }
        return bytes;
    }

    const std::vector<IcFactorDiagnostics>& diagnostics() const { return diagnostics_; }
    const std::vector<IcFactorDiagnostics>& trialDiagnostics() const { return trialDiagnostics_; }
    bool hasBreakdown() const
    {
        for (const IcFactorDiagnostics& item : diagnostics_) {
            if (item.breakdown || item.nonFiniteLCount > 0) {
                return true;
            }
        }
        return false;
    }

private:
    std::vector<std::vector<int>> blockDofs_;
    std::vector<int> globalToBlock_;
    std::vector<int> globalToLocal_;
    std::vector<SubdomainIcSolver> solvers_;
    std::vector<std::vector<double>> rhsBuffers_;
    std::vector<std::vector<double>> solutionBuffers_;
    std::vector<SubdomainIcSolver::Workspace> workspaces_;
    std::vector<IcFactorDiagnostics> diagnostics_;
    std::vector<IcFactorDiagnostics> trialDiagnostics_;

    void build(const Mesh& mesh,
               const SparseMatrix& a,
               double icShift,
               bool diagonalScaling,
               double diagScalingEps,
               const std::string& shiftMode,
               const std::string& ordering)
    {
        int maxSubdomain = 0;
        for (const Node& node : mesh.nodes) {
            maxSubdomain = std::max(maxSubdomain, node.subdomain);
        }
        blockDofs_.assign(static_cast<size_t>(maxSubdomain + 1), {});
        globalToBlock_.assign(mesh.nodes.size(), -1);
        globalToLocal_.assign(mesh.nodes.size(), -1);

        for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
            const int block = mesh.nodes[static_cast<size_t>(i)].subdomain;
            globalToBlock_[static_cast<size_t>(i)] = block;
            globalToLocal_[static_cast<size_t>(i)] = static_cast<int>(blockDofs_[static_cast<size_t>(block)].size());
            blockDofs_[static_cast<size_t>(block)].push_back(i);
        }

        solvers_.clear();
        solvers_.resize(blockDofs_.size());
        rhsBuffers_.resize(blockDofs_.size());
        solutionBuffers_.resize(blockDofs_.size());
        workspaces_.resize(blockDofs_.size());
        for (size_t blockIndex = 0; blockIndex < blockDofs_.size(); ++blockIndex) {
            const size_t nLocal = blockDofs_[blockIndex].size();
            rhsBuffers_[blockIndex].resize(nLocal);
            solutionBuffers_[blockIndex].resize(nLocal);
            workspaces_[blockIndex].scaledRhs.resize(nLocal);
            workspaces_[blockIndex].y.resize(nLocal);
            workspaces_[blockIndex].xPerm.resize(nLocal);
        }
        diagnostics_.assign(blockDofs_.size(), {});
        trialDiagnostics_.clear();
        parallelForCoarse(blockDofs_.size(), [&](size_t blockIndex) {
            const int block = static_cast<int>(blockIndex);
            std::vector<MatrixEntry> entries;
            for (int globalRow : blockDofs_[blockIndex]) {
                const int localRow = globalToLocal_[static_cast<size_t>(globalRow)];
                for (int k = a.rowPtr[static_cast<size_t>(globalRow)]; k < a.rowPtr[static_cast<size_t>(globalRow + 1)]; ++k) {
                    const int globalCol = a.colInd[static_cast<size_t>(k)];
                    if (globalCol < 0 || globalCol >= static_cast<int>(globalToBlock_.size())) {
                        continue;
                    }
                    if (globalToBlock_[static_cast<size_t>(globalCol)] != block) {
                        continue;
                    }
                    entries.push_back({localRow, globalToLocal_[static_cast<size_t>(globalCol)], a.values[static_cast<size_t>(k)]});
                }
            }
            solvers_[blockIndex] = SubdomainIcSolver(static_cast<int>(blockDofs_[blockIndex].size()),
                                                     entries,
                                                     icShift,
                                                     block,
                                                     diagonalScaling,
                                                     diagScalingEps,
                                                     shiftMode,
                                                     ordering);
            diagnostics_[blockIndex] = solvers_[blockIndex].diagnostics();
        });
        for (const auto& solver : solvers_) {
            const auto& trials = solver.trialDiagnostics();
            trialDiagnostics_.insert(trialDiagnostics_.end(), trials.begin(), trials.end());
        }
    }
};

class SubdomainIlutSolver {
public:
    SubdomainIlutSolver() = default;
    SubdomainIlutSolver(int n,
                        const std::vector<MatrixEntry>& entries,
                        double dropTolerance,
                        int fillFactor,
                        bool diagonalScaling = false,
                        double diagScalingEps = 1.0e-30)
    {
        factor(n, entries, dropTolerance, fillFactor, diagonalScaling, diagScalingEps);
    }

    void solve(const std::vector<double>& rhs, std::vector<double>& x) const
    {
        std::vector<double> scaledRhs(static_cast<size_t>(n_), 0.0);
        for (int i = 0; i < n_; ++i) {
            scaledRhs[static_cast<size_t>(i)] = scaling_[static_cast<size_t>(i)] * rhs[static_cast<size_t>(i)];
        }

        std::vector<double> y(static_cast<size_t>(n_), 0.0);
        for (int i = 0; i < n_; ++i) {
            double sum = scaledRhs[static_cast<size_t>(i)];
            for (const auto& entry : lowerRows_[static_cast<size_t>(i)]) {
                sum -= entry.second * y[static_cast<size_t>(entry.first)];
            }
            y[static_cast<size_t>(i)] = sum;
        }

        x.assign(static_cast<size_t>(n_), 0.0);
        for (int i = n_ - 1; i >= 0; --i) {
            double sum = y[static_cast<size_t>(i)];
            for (const auto& entry : upperRows_[static_cast<size_t>(i)]) {
                sum -= entry.second * x[static_cast<size_t>(entry.first)];
            }
            x[static_cast<size_t>(i)] = sum / upperDiagonal_[static_cast<size_t>(i)];
        }
        for (int i = 0; i < n_; ++i) {
            x[static_cast<size_t>(i)] *= scaling_[static_cast<size_t>(i)];
        }
    }

    size_t memoryBytes() const
    {
        size_t bytes = upperDiagonal_.size() * sizeof(double);
        bytes += scaling_.size() * sizeof(double);
        for (const auto& row : lowerRows_) {
            bytes += row.size() * (sizeof(int) + sizeof(double));
        }
        for (const auto& row : upperRows_) {
            bytes += row.size() * (sizeof(int) + sizeof(double));
        }
        return bytes;
    }

private:
    int n_ = 0;
    std::vector<std::vector<std::pair<int, double>>> lowerRows_;
    std::vector<std::vector<std::pair<int, double>>> upperRows_;
    std::vector<double> upperDiagonal_;
    std::vector<double> scaling_;

    static void compressRows(std::vector<std::vector<std::pair<int, double>>>& rows)
    {
        for (auto& row : rows) {
            std::sort(row.begin(), row.end(), [](const auto& a, const auto& b) {
                return a.first < b.first;
            });
            size_t out = 0;
            for (size_t i = 0; i < row.size();) {
                const int col = row[i].first;
                double sum = 0.0;
                while (i < row.size() && row[i].first == col) {
                    sum += row[i].second;
                    ++i;
                }
                if (std::abs(sum) > 0.0) {
                    row[out++] = {col, sum};
                }
            }
            row.resize(out);
        }
    }

    static void dropAndLimit(std::vector<std::pair<int, double>>& row,
                             double threshold,
                             int limit)
    {
        row.erase(std::remove_if(row.begin(), row.end(), [&](const auto& entry) {
            return std::abs(entry.second) < threshold || !std::isfinite(entry.second);
        }), row.end());
        if (limit > 0 && static_cast<int>(row.size()) > limit) {
            std::nth_element(row.begin(),
                             row.begin() + limit,
                             row.end(),
                             [](const auto& a, const auto& b) {
                                 return std::abs(a.second) > std::abs(b.second);
                             });
            row.resize(static_cast<size_t>(limit));
        }
        std::sort(row.begin(), row.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
    }

    void factor(int n,
                const std::vector<MatrixEntry>& entries,
                double dropTolerance,
                int fillFactor,
                bool diagonalScaling,
                double diagScalingEps)
    {
        n_ = n;
        lowerRows_.assign(static_cast<size_t>(n_), {});
        upperRows_.assign(static_cast<size_t>(n_), {});
        upperDiagonal_.assign(static_cast<size_t>(n_), 1.0);
        scaling_.assign(static_cast<size_t>(n_), 1.0);

        std::vector<std::vector<std::pair<int, double>>> rows(static_cast<size_t>(n_));
        std::vector<int> originalNnz(static_cast<size_t>(n_), 0);
        std::vector<double> rowScale(static_cast<size_t>(n_), 0.0);
        std::vector<double> rawDiagonal(static_cast<size_t>(n_), 0.0);
        for (const MatrixEntry& entry : entries) {
            if (entry.row == entry.col && entry.row >= 0 && entry.row < n_ && std::isfinite(entry.value)) {
                rawDiagonal[static_cast<size_t>(entry.row)] += entry.value;
            }
        }
        if (diagonalScaling) {
            const double floor = std::max(1.0e-300, diagScalingEps);
            for (int i = 0; i < n_; ++i) {
                scaling_[static_cast<size_t>(i)] =
                    1.0 / std::sqrt(std::max(floor, std::abs(rawDiagonal[static_cast<size_t>(i)])));
            }
        }
        double diagonalScale = 0.0;
        int diagonalCount = 0;
        for (const MatrixEntry& entry : entries) {
            if (entry.row < 0 || entry.row >= n_ || entry.col < 0 || entry.col >= n_) {
                continue;
            }
            const double scaledValue = entry.value
                * scaling_[static_cast<size_t>(entry.row)]
                * scaling_[static_cast<size_t>(entry.col)];
            rows[static_cast<size_t>(entry.row)].push_back({entry.col, scaledValue});
            rowScale[static_cast<size_t>(entry.row)] =
                std::max(rowScale[static_cast<size_t>(entry.row)], std::abs(scaledValue));
            if (entry.row == entry.col && std::isfinite(scaledValue)) {
                diagonalScale += std::abs(scaledValue);
                ++diagonalCount;
            }
        }
        compressRows(rows);
        for (int i = 0; i < n_; ++i) {
            originalNnz[static_cast<size_t>(i)] =
                std::max<int>(1, static_cast<int>(rows[static_cast<size_t>(i)].size()));
            rowScale[static_cast<size_t>(i)] = std::max(rowScale[static_cast<size_t>(i)], 1.0);
        }
        const double pivotFloor = std::max(std::max(1.0e-300, diagScalingEps),
            1.0e-12 * (diagonalCount > 0 ? diagonalScale / static_cast<double>(diagonalCount) : 1.0));
        const int boundedFillFactor = std::max(1, fillFactor);

        for (int i = 0; i < n_; ++i) {
            std::unordered_map<int, double> work;
            work.reserve(rows[static_cast<size_t>(i)].size() * 2 + 8);
            for (const auto& entry : rows[static_cast<size_t>(i)]) {
                work[entry.first] += entry.second;
            }

            std::vector<int> lowerColumns;
            lowerColumns.reserve(work.size());
            for (const auto& entry : work) {
                if (entry.first < i) {
                    lowerColumns.push_back(entry.first);
                }
            }
            std::sort(lowerColumns.begin(), lowerColumns.end());

            const double threshold = std::max(0.0, dropTolerance) * rowScale[static_cast<size_t>(i)];
            for (int col : lowerColumns) {
                auto found = work.find(col);
                if (found == work.end()) {
                    continue;
                }
                double lij = found->second / upperDiagonal_[static_cast<size_t>(col)];
                if (!std::isfinite(lij) || std::abs(lij) < threshold) {
                    work.erase(found);
                    continue;
                }
                found->second = lij;
                for (const auto& upper : upperRows_[static_cast<size_t>(col)]) {
                    work[upper.first] -= lij * upper.second;
                }
            }

            std::vector<std::pair<int, double>> lower;
            std::vector<std::pair<int, double>> upper;
            double diagonal = 0.0;
            for (const auto& entry : work) {
                if (entry.first < i) {
                    lower.push_back(entry);
                } else if (entry.first == i) {
                    diagonal += entry.second;
                } else {
                    upper.push_back(entry);
                }
            }

            if (!std::isfinite(diagonal) || std::abs(diagonal) < pivotFloor) {
                diagonal = diagonal < 0.0 ? -pivotFloor : pivotFloor;
            }
            upperDiagonal_[static_cast<size_t>(i)] = diagonal;
            const int keep = std::max(1, boundedFillFactor * originalNnz[static_cast<size_t>(i)]);
            dropAndLimit(lower, threshold, keep);
            dropAndLimit(upper, threshold, keep);
            lowerRows_[static_cast<size_t>(i)] = std::move(lower);
            upperRows_[static_cast<size_t>(i)] = std::move(upper);
        }
    }
};

class BlockJacobiGeneralPardisoPreconditioner {
public:
    BlockJacobiGeneralPardisoPreconditioner(const Mesh& mesh, const SparseMatrix& a)
    {
        build(mesh, a);
    }

    void apply(const std::vector<double>& r, std::vector<double>& z)
    {
        z.assign(r.size(), 0.0);
        parallelForCoarse(blockDofs_.size(), [&](size_t block) {
            ScopedMklSingleThread mklThreads;
            const std::vector<int>& dofs = blockDofs_[block];
            std::vector<double> rhs(dofs.size(), 0.0);
            for (size_t i = 0; i < dofs.size(); ++i) {
                rhs[i] = r[static_cast<size_t>(dofs[i])];
            }
            std::vector<double> localSolution;
            solvers_[block].solve(rhs, localSolution);
            for (size_t i = 0; i < dofs.size(); ++i) {
                z[static_cast<size_t>(dofs[i])] = localSolution[i];
            }
        });
    }

    size_t memoryBytes() const
    {
        size_t bytes = globalToBlock_.size() * sizeof(int)
            + globalToLocal_.size() * sizeof(int);
        for (const auto& block : blockDofs_) {
            bytes += block.size() * sizeof(int);
        }
        for (const auto& solver : solvers_) {
            bytes += solver.memoryBytes();
        }
        return bytes;
    }

private:
    std::vector<std::vector<int>> blockDofs_;
    std::vector<int> globalToBlock_;
    std::vector<int> globalToLocal_;
    std::vector<GeneralSparseDirectSolver> solvers_;

    void build(const Mesh& mesh, const SparseMatrix& a)
    {
        int maxSubdomain = 0;
        for (const Node& node : mesh.nodes) {
            maxSubdomain = std::max(maxSubdomain, node.subdomain);
        }
        blockDofs_.assign(static_cast<size_t>(maxSubdomain + 1), {});
        globalToBlock_.assign(mesh.nodes.size(), -1);
        globalToLocal_.assign(mesh.nodes.size(), -1);
        for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
            const int block = mesh.nodes[static_cast<size_t>(i)].subdomain;
            globalToBlock_[static_cast<size_t>(i)] = block;
            globalToLocal_[static_cast<size_t>(i)] = static_cast<int>(blockDofs_[static_cast<size_t>(block)].size());
            blockDofs_[static_cast<size_t>(block)].push_back(i);
        }

        solvers_.resize(blockDofs_.size());
        std::vector<std::string> factorErrors(blockDofs_.size());
        parallelForCoarse(blockDofs_.size(), [&](size_t blockIndex) {
            try {
                ScopedMklSingleThread mklThreads;
                const int block = static_cast<int>(blockIndex);
                std::vector<MatrixEntry> entries;
                for (int globalRow : blockDofs_[blockIndex]) {
                    const int localRow = globalToLocal_[static_cast<size_t>(globalRow)];
                    for (int k = a.rowPtr[static_cast<size_t>(globalRow)]; k < a.rowPtr[static_cast<size_t>(globalRow + 1)]; ++k) {
                        const int globalCol = a.colInd[static_cast<size_t>(k)];
                        if (globalCol < 0 || globalCol >= static_cast<int>(globalToBlock_.size())) {
                            continue;
                        }
                        if (globalToBlock_[static_cast<size_t>(globalCol)] != block) {
                            continue;
                        }
                        entries.push_back({localRow, globalToLocal_[static_cast<size_t>(globalCol)], a.values[static_cast<size_t>(k)]});
                    }
                }
                solvers_[blockIndex] = GeneralSparseDirectSolver(static_cast<int>(blockDofs_[blockIndex].size()), entries);
            } catch (const std::exception& err) {
                factorErrors[blockIndex] = err.what();
            }
        });
        for (size_t blockIndex = 0; blockIndex < factorErrors.size(); ++blockIndex) {
            if (!factorErrors[blockIndex].empty()) {
                throw std::runtime_error("Block " + std::to_string(blockIndex)
                    + " general PARDISO factorization failed: " + factorErrors[blockIndex]);
            }
        }
    }
};

class BlockJacobiIlutPreconditioner {
public:
    BlockJacobiIlutPreconditioner(const Mesh& mesh,
                                  const SparseMatrix& a,
                                  double dropTolerance,
                                  int fillFactor,
                                  bool diagonalScaling = false,
                                  double diagScalingEps = 1.0e-30)
    {
        build(mesh, a, dropTolerance, fillFactor, diagonalScaling, diagScalingEps);
    }

    void apply(const std::vector<double>& r, std::vector<double>& z) const
    {
        z.assign(r.size(), 0.0);
        parallelForCoarse(blockDofs_.size(), [&](size_t block) {
            const std::vector<int>& dofs = blockDofs_[block];
            std::vector<double> rhs(dofs.size(), 0.0);
            for (size_t i = 0; i < dofs.size(); ++i) {
                rhs[i] = r[static_cast<size_t>(dofs[i])];
            }
            std::vector<double> localSolution;
            solvers_[block].solve(rhs, localSolution);
            for (size_t i = 0; i < dofs.size(); ++i) {
                z[static_cast<size_t>(dofs[i])] = localSolution[i];
            }
        });
    }

    size_t memoryBytes() const
    {
        size_t bytes = globalToBlock_.size() * sizeof(int)
            + globalToLocal_.size() * sizeof(int);
        for (const auto& block : blockDofs_) {
            bytes += block.size() * sizeof(int);
        }
        for (const auto& solver : solvers_) {
            bytes += solver.memoryBytes();
        }
        return bytes;
    }

private:
    std::vector<std::vector<int>> blockDofs_;
    std::vector<int> globalToBlock_;
    std::vector<int> globalToLocal_;
    std::vector<SubdomainIlutSolver> solvers_;

    void build(const Mesh& mesh,
               const SparseMatrix& a,
               double dropTolerance,
               int fillFactor,
               bool diagonalScaling,
               double diagScalingEps)
    {
        int maxSubdomain = 0;
        for (const Node& node : mesh.nodes) {
            maxSubdomain = std::max(maxSubdomain, node.subdomain);
        }
        blockDofs_.assign(static_cast<size_t>(maxSubdomain + 1), {});
        globalToBlock_.assign(mesh.nodes.size(), -1);
        globalToLocal_.assign(mesh.nodes.size(), -1);
        for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
            const int block = mesh.nodes[static_cast<size_t>(i)].subdomain;
            globalToBlock_[static_cast<size_t>(i)] = block;
            globalToLocal_[static_cast<size_t>(i)] = static_cast<int>(blockDofs_[static_cast<size_t>(block)].size());
            blockDofs_[static_cast<size_t>(block)].push_back(i);
        }

        solvers_.resize(blockDofs_.size());
        parallelForCoarse(blockDofs_.size(), [&](size_t blockIndex) {
            const int block = static_cast<int>(blockIndex);
            std::vector<MatrixEntry> entries;
            for (int globalRow : blockDofs_[blockIndex]) {
                const int localRow = globalToLocal_[static_cast<size_t>(globalRow)];
                for (int k = a.rowPtr[static_cast<size_t>(globalRow)]; k < a.rowPtr[static_cast<size_t>(globalRow + 1)]; ++k) {
                    const int globalCol = a.colInd[static_cast<size_t>(k)];
                    if (globalCol < 0 || globalCol >= static_cast<int>(globalToBlock_.size())) {
                        continue;
                    }
                    if (globalToBlock_[static_cast<size_t>(globalCol)] != block) {
                        continue;
                    }
                    entries.push_back({localRow, globalToLocal_[static_cast<size_t>(globalCol)], a.values[static_cast<size_t>(k)]});
                }
            }
            solvers_[blockIndex] = SubdomainIlutSolver(static_cast<int>(blockDofs_[blockIndex].size()),
                                                       entries,
                                                       dropTolerance,
                                                       fillFactor,
                                                       diagonalScaling,
                                                       diagScalingEps);
        });
    }
};

class RasIlutPreconditioner {
public:
    RasIlutPreconditioner(const Mesh& mesh,
                          const SparseMatrix& a,
                          int overlap,
                          double dropTolerance,
                          int fillFactor,
                          bool diagonalScaling = false,
                          double diagScalingEps = 1.0e-30)
    {
        build(mesh, a, overlap, dropTolerance, fillFactor, diagonalScaling, diagScalingEps);
    }

    void apply(const std::vector<double>& r, std::vector<double>& z) const
    {
        z.assign(r.size(), 0.0);
        parallelForCoarse(localDofs_.size(), [&](size_t block) {
            const std::vector<int>& localDofs = localDofs_[block];
            std::vector<double> rhs(localDofs.size(), 0.0);
            for (size_t i = 0; i < localDofs.size(); ++i) {
                rhs[i] = r[static_cast<size_t>(localDofs[i])];
            }
            std::vector<double> localSolution;
            solvers_[block].solve(rhs, localSolution);
            for (int globalDof : coreDofs_[block]) {
                const int local = localIndex_[block].at(globalDof);
                z[static_cast<size_t>(globalDof)] = localSolution[static_cast<size_t>(local)];
            }
        });
    }

    size_t memoryBytes() const
    {
        size_t bytes = 0;
        for (const auto& block : coreDofs_) {
            bytes += block.size() * sizeof(int);
        }
        for (const auto& block : localDofs_) {
            bytes += block.size() * sizeof(int);
        }
        for (const auto& map : localIndex_) {
            bytes += map.size() * (sizeof(int) * 2 + sizeof(size_t));
        }
        for (const auto& solver : solvers_) {
            bytes += solver.memoryBytes();
        }
        return bytes;
    }

private:
    std::vector<std::vector<int>> coreDofs_;
    std::vector<std::vector<int>> localDofs_;
    std::vector<std::unordered_map<int, int>> localIndex_;
    std::vector<SubdomainIlutSolver> solvers_;

    void build(const Mesh& mesh,
               const SparseMatrix& a,
               int overlap,
               double dropTolerance,
               int fillFactor,
               bool diagonalScaling,
               double diagScalingEps)
    {
        int maxSubdomain = 0;
        for (const Node& node : mesh.nodes) {
            maxSubdomain = std::max(maxSubdomain, node.subdomain);
        }
        coreDofs_.assign(static_cast<size_t>(maxSubdomain + 1), {});
        for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
            coreDofs_[static_cast<size_t>(mesh.nodes[static_cast<size_t>(i)].subdomain)].push_back(i);
        }

        localDofs_.assign(coreDofs_.size(), {});
        localIndex_.assign(coreDofs_.size(), {});
        solvers_.resize(coreDofs_.size());
        const int overlapLayers = std::max(0, overlap);
        parallelForCoarse(coreDofs_.size(), [&](size_t blockIndex) {
            std::vector<int> local = coreDofs_[blockIndex];
            std::vector<int> frontier = coreDofs_[blockIndex];
            std::vector<char> inLocal(mesh.nodes.size(), 0);
            for (int dof : local) {
                inLocal[static_cast<size_t>(dof)] = 1;
            }
            for (int layer = 0; layer < overlapLayers; ++layer) {
                std::vector<int> nextFrontier;
                for (int globalRow : frontier) {
                    for (int k = a.rowPtr[static_cast<size_t>(globalRow)]; k < a.rowPtr[static_cast<size_t>(globalRow + 1)]; ++k) {
                        const int globalCol = a.colInd[static_cast<size_t>(k)];
                        if (globalCol < 0 || globalCol >= static_cast<int>(mesh.nodes.size())) {
                            continue;
                        }
                        if (!inLocal[static_cast<size_t>(globalCol)]) {
                            inLocal[static_cast<size_t>(globalCol)] = 1;
                            local.push_back(globalCol);
                            nextFrontier.push_back(globalCol);
                        }
                    }
                }
                frontier = std::move(nextFrontier);
                if (frontier.empty()) {
                    break;
                }
            }
            std::sort(local.begin(), local.end());
            local.erase(std::unique(local.begin(), local.end()), local.end());
            localDofs_[blockIndex] = std::move(local);

            auto& map = localIndex_[blockIndex];
            map.reserve(localDofs_[blockIndex].size());
            for (int i = 0; i < static_cast<int>(localDofs_[blockIndex].size()); ++i) {
                map[localDofs_[blockIndex][static_cast<size_t>(i)]] = i;
            }

            std::vector<MatrixEntry> entries;
            for (int globalRow : localDofs_[blockIndex]) {
                const int localRow = map[globalRow];
                for (int k = a.rowPtr[static_cast<size_t>(globalRow)]; k < a.rowPtr[static_cast<size_t>(globalRow + 1)]; ++k) {
                    const int globalCol = a.colInd[static_cast<size_t>(k)];
                    const auto found = map.find(globalCol);
                    if (found != map.end()) {
                        entries.push_back({localRow, found->second, a.values[static_cast<size_t>(k)]});
                    }
                }
            }
            solvers_[blockIndex] = SubdomainIlutSolver(static_cast<int>(localDofs_[blockIndex].size()),
                                                       entries,
                                                       dropTolerance,
                                                       fillFactor,
                                                       diagonalScaling,
                                                       diagScalingEps);
        });
    }
};

class RasIcPreconditioner {
public:
    RasIcPreconditioner(const Mesh& mesh,
                        const SparseMatrix& a,
                        int overlap,
                        double icShift,
                        bool diagonalScaling,
                        double diagScalingEps,
                        const std::string& shiftMode,
                        const std::string& ordering)
    {
        build(mesh, a, overlap, icShift, diagonalScaling, diagScalingEps, shiftMode, ordering);
    }

    void apply(const std::vector<double>& r, std::vector<double>& z) const
    {
        z.assign(r.size(), 0.0);
        parallelForCoarse(localDofs_.size(), [&](size_t block) {
            const std::vector<int>& localDofs = localDofs_[block];
            std::vector<double> rhs(localDofs.size(), 0.0);
            for (size_t i = 0; i < localDofs.size(); ++i) {
                rhs[i] = r[static_cast<size_t>(localDofs[i])];
            }
            std::vector<double> localSolution;
            solvers_[block].solve(rhs, localSolution);
            for (int globalDof : coreDofs_[block]) {
                const int local = localIndex_[block].at(globalDof);
                z[static_cast<size_t>(globalDof)] = localSolution[static_cast<size_t>(local)];
            }
        });
    }

    size_t memoryBytes() const
    {
        size_t bytes = 0;
        for (const auto& block : coreDofs_) {
            bytes += block.size() * sizeof(int);
        }
        for (const auto& block : localDofs_) {
            bytes += block.size() * sizeof(int);
        }
        for (const auto& map : localIndex_) {
            bytes += map.size() * (sizeof(int) * 2 + sizeof(size_t));
        }
        for (const auto& solver : solvers_) {
            bytes += solver.memoryBytes();
        }
        return bytes;
    }

    bool hasBreakdown() const
    {
        for (const auto& solver : solvers_) {
            const IcFactorDiagnostics& item = solver.diagnostics();
            if (item.breakdown || item.nonFiniteLCount > 0) {
                return true;
            }
        }
        return false;
    }
    const std::vector<IcFactorDiagnostics>& diagnostics() const { return diagnostics_; }
    const std::vector<IcFactorDiagnostics>& trialDiagnostics() const { return trialDiagnostics_; }

private:
    std::vector<std::vector<int>> coreDofs_;
    std::vector<std::vector<int>> localDofs_;
    std::vector<std::unordered_map<int, int>> localIndex_;
    std::vector<SubdomainIcSolver> solvers_;
    std::vector<IcFactorDiagnostics> diagnostics_;
    std::vector<IcFactorDiagnostics> trialDiagnostics_;

    void build(const Mesh& mesh,
               const SparseMatrix& a,
               int overlap,
               double icShift,
               bool diagonalScaling,
               double diagScalingEps,
               const std::string& shiftMode,
               const std::string& ordering)
    {
        int maxSubdomain = 0;
        for (const Node& node : mesh.nodes) {
            maxSubdomain = std::max(maxSubdomain, node.subdomain);
        }
        coreDofs_.assign(static_cast<size_t>(maxSubdomain + 1), {});
        for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
            coreDofs_[static_cast<size_t>(mesh.nodes[static_cast<size_t>(i)].subdomain)].push_back(i);
        }

        localDofs_.assign(coreDofs_.size(), {});
        localIndex_.assign(coreDofs_.size(), {});
        solvers_.resize(coreDofs_.size());
        diagnostics_.assign(coreDofs_.size(), {});
        trialDiagnostics_.clear();
        const int overlapLayers = std::max(0, overlap);
        parallelForCoarse(coreDofs_.size(), [&](size_t blockIndex) {
            std::vector<int> local = coreDofs_[blockIndex];
            std::vector<int> frontier = coreDofs_[blockIndex];
            std::vector<char> inLocal(mesh.nodes.size(), 0);
            for (int dof : local) {
                inLocal[static_cast<size_t>(dof)] = 1;
            }
            for (int layer = 0; layer < overlapLayers; ++layer) {
                std::vector<int> nextFrontier;
                for (int globalRow : frontier) {
                    for (int k = a.rowPtr[static_cast<size_t>(globalRow)]; k < a.rowPtr[static_cast<size_t>(globalRow + 1)]; ++k) {
                        const int globalCol = a.colInd[static_cast<size_t>(k)];
                        if (globalCol < 0 || globalCol >= static_cast<int>(mesh.nodes.size())) {
                            continue;
                        }
                        if (!inLocal[static_cast<size_t>(globalCol)]) {
                            inLocal[static_cast<size_t>(globalCol)] = 1;
                            local.push_back(globalCol);
                            nextFrontier.push_back(globalCol);
                        }
                    }
                }
                frontier = std::move(nextFrontier);
                if (frontier.empty()) {
                    break;
                }
            }
            std::sort(local.begin(), local.end());
            local.erase(std::unique(local.begin(), local.end()), local.end());
            localDofs_[blockIndex] = std::move(local);

            auto& map = localIndex_[blockIndex];
            map.reserve(localDofs_[blockIndex].size());
            for (int i = 0; i < static_cast<int>(localDofs_[blockIndex].size()); ++i) {
                map[localDofs_[blockIndex][static_cast<size_t>(i)]] = i;
            }

            std::vector<MatrixEntry> entries;
            for (int globalRow : localDofs_[blockIndex]) {
                const int localRow = map[globalRow];
                for (int k = a.rowPtr[static_cast<size_t>(globalRow)]; k < a.rowPtr[static_cast<size_t>(globalRow + 1)]; ++k) {
                    const int globalCol = a.colInd[static_cast<size_t>(k)];
                    const auto found = map.find(globalCol);
                    if (found != map.end()) {
                        entries.push_back({localRow, found->second, a.values[static_cast<size_t>(k)]});
                    }
                }
            }
            solvers_[blockIndex] = SubdomainIcSolver(static_cast<int>(localDofs_[blockIndex].size()),
                                                     entries,
                                                     icShift,
                                                     static_cast<int>(blockIndex),
                                                     diagonalScaling,
                                                     diagScalingEps,
                                                     shiftMode,
                                                     ordering);
        });
        for (size_t blockIndex = 0; blockIndex < solvers_.size(); ++blockIndex) {
            diagnostics_[blockIndex] = solvers_[blockIndex].diagnostics();
            const auto& trials = solvers_[blockIndex].trialDiagnostics();
            trialDiagnostics_.insert(trialDiagnostics_.end(), trials.begin(), trials.end());
        }
    }
};

static bool solveSmallDenseSystem(std::vector<std::vector<double>> a,
                                  std::vector<double> b,
                                  std::vector<double>& x)
{
    const int n = static_cast<int>(b.size());
    x.assign(static_cast<size_t>(n), 0.0);
    if (n == 0) {
        return true;
    }
    double scale = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            scale = std::max(scale, std::abs(a[static_cast<size_t>(i)][static_cast<size_t>(j)]));
        }
    }
    const double pivotFloor = std::max(1.0e-30, 1.0e-14 * std::max(1.0, scale));
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double pivotAbs = std::abs(a[static_cast<size_t>(col)][static_cast<size_t>(col)]);
        for (int row = col + 1; row < n; ++row) {
            const double candidate = std::abs(a[static_cast<size_t>(row)][static_cast<size_t>(col)]);
            if (candidate > pivotAbs) {
                pivotAbs = candidate;
                pivot = row;
            }
        }
        if (pivotAbs < pivotFloor) {
            a[static_cast<size_t>(col)][static_cast<size_t>(col)] += pivotFloor;
            pivot = col;
            pivotAbs = std::abs(a[static_cast<size_t>(col)][static_cast<size_t>(col)]);
        }
        if (pivotAbs < pivotFloor) {
            return false;
        }
        if (pivot != col) {
            std::swap(a[static_cast<size_t>(pivot)], a[static_cast<size_t>(col)]);
            std::swap(b[static_cast<size_t>(pivot)], b[static_cast<size_t>(col)]);
        }
        const double diag = a[static_cast<size_t>(col)][static_cast<size_t>(col)];
        for (int row = col + 1; row < n; ++row) {
            const double factor = a[static_cast<size_t>(row)][static_cast<size_t>(col)] / diag;
            if (factor == 0.0) {
                continue;
            }
            a[static_cast<size_t>(row)][static_cast<size_t>(col)] = 0.0;
            for (int j = col + 1; j < n; ++j) {
                a[static_cast<size_t>(row)][static_cast<size_t>(j)] -=
                    factor * a[static_cast<size_t>(col)][static_cast<size_t>(j)];
            }
            b[static_cast<size_t>(row)] -= factor * b[static_cast<size_t>(col)];
        }
    }
    for (int row = n - 1; row >= 0; --row) {
        double sum = b[static_cast<size_t>(row)];
        for (int col = row + 1; col < n; ++col) {
            sum -= a[static_cast<size_t>(row)][static_cast<size_t>(col)] * x[static_cast<size_t>(col)];
        }
        const double diag = a[static_cast<size_t>(row)][static_cast<size_t>(row)];
        if (std::abs(diag) < pivotFloor) {
            return false;
        }
        x[static_cast<size_t>(row)] = sum / diag;
    }
    return true;
}

struct SubdomainConstantCoarseSpace {
    size_t globalSize = 0;
    std::vector<std::vector<int>> dofsByMode;
    std::vector<int> globalToMode;
    std::vector<double> invNormByMode;

    int dim() const { return static_cast<int>(dofsByMode.size()); }

    size_t memoryBytes() const
    {
        size_t bytes = globalToMode.size() * sizeof(int)
            + invNormByMode.size() * sizeof(double);
        for (const auto& dofs : dofsByMode) {
            bytes += dofs.size() * sizeof(int);
        }
        return bytes;
    }
};

static SubdomainConstantCoarseSpace buildSubdomainConstantCoarseSpace(const Mesh& mesh)
{
    SubdomainConstantCoarseSpace space;
    space.globalSize = mesh.nodes.size();
    space.globalToMode.assign(mesh.nodes.size(), -1);

    int maxSubdomain = 0;
    for (const Node& node : mesh.nodes) {
        maxSubdomain = std::max(maxSubdomain, node.subdomain);
    }
    std::vector<std::vector<int>> dofsBySubdomain(static_cast<size_t>(maxSubdomain + 1));
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        if (!mesh.nodes[i].dirichlet) {
            dofsBySubdomain[static_cast<size_t>(mesh.nodes[i].subdomain)].push_back(static_cast<int>(i));
        }
    }

    for (std::vector<int>& dofs : dofsBySubdomain) {
        if (dofs.empty()) {
            continue;
        }
        const int mode = static_cast<int>(space.dofsByMode.size());
        const double invNorm = 1.0 / std::sqrt(static_cast<double>(dofs.size()));
        for (int dof : dofs) {
            space.globalToMode[static_cast<size_t>(dof)] = mode;
        }
        space.invNormByMode.push_back(invNorm);
        space.dofsByMode.push_back(std::move(dofs));
    }
    return space;
}

static SubdomainConstantCoarseSpace buildSubdomainMaterialConstantCoarseSpace(const Mesh& mesh,
                                                                              const CaseConfig& physics)
{
    SubdomainConstantCoarseSpace space;
    space.globalSize = mesh.nodes.size();
    space.globalToMode.assign(mesh.nodes.size(), -1);

    std::vector<std::vector<int>> incidentTets(mesh.nodes.size());
    for (int tetId = 0; tetId < static_cast<int>(mesh.tets.size()); ++tetId) {
        const Tet& tet = mesh.tets[static_cast<size_t>(tetId)];
        for (int dof : tet.dof) {
            if (dof >= 0 && dof < static_cast<int>(incidentTets.size())) {
                incidentTets[static_cast<size_t>(dof)].push_back(tetId);
            }
        }
    }

    std::map<std::string, int> modeByKey;
    for (size_t dof = 0; dof < mesh.nodes.size(); ++dof) {
        if (mesh.nodes[dof].dirichlet) {
            continue;
        }
        const int subdomain = mesh.nodes[dof].subdomain;
        std::string materialKey = "default";
        for (int tetId : incidentTets[dof]) {
            const Tet& tet = mesh.tets[static_cast<size_t>(tetId)];
            if (tet.subdomain != subdomain) {
                continue;
            }
            const Material& material = materialForTet(physics, tet);
            materialKey = material.name + ":" + std::to_string(tet.domainEntity);
            break;
        }
        const std::string key = std::to_string(subdomain) + "|" + materialKey;
        auto found = modeByKey.find(key);
        if (found == modeByKey.end()) {
            const int mode = static_cast<int>(space.dofsByMode.size());
            found = modeByKey.emplace(key, mode).first;
            space.dofsByMode.push_back({});
            space.invNormByMode.push_back(1.0);
        }
        const int mode = found->second;
        space.globalToMode[dof] = mode;
        space.dofsByMode[static_cast<size_t>(mode)].push_back(static_cast<int>(dof));
    }

    for (size_t mode = 0; mode < space.dofsByMode.size(); ++mode) {
        const size_t count = space.dofsByMode[mode].size();
        space.invNormByMode[mode] = count > 0 ? 1.0 / std::sqrt(static_cast<double>(count)) : 1.0;
    }
    return space;
}

static SubdomainConstantCoarseSpace buildInterfaceConstantCoarseSpace(const Mesh& mesh)
{
    std::vector<std::string> keyByDof(mesh.nodes.size());
    for (size_t dof = 0; dof < mesh.nodes.size(); ++dof) {
        if (!mesh.nodes[dof].dirichlet) {
            keyByDof[dof] = "subdomain:" + std::to_string(mesh.nodes[dof].subdomain);
        }
    }

    for (const InterfaceFace& face : mesh.interfaceFaces) {
        if (face.leftTet < 0 || face.leftTet >= static_cast<int>(mesh.tets.size())
            || face.rightTet < 0 || face.rightTet >= static_cast<int>(mesh.tets.size())) {
            continue;
        }
        const Tet& leftTet = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& rightTet = mesh.tets[static_cast<size_t>(face.rightTet)];
        const int a = std::min(leftTet.subdomain, rightTet.subdomain);
        const int b = std::max(leftTet.subdomain, rightTet.subdomain);
        const std::string key = "interface:" + std::to_string(a) + "-" + std::to_string(b);
        for (int local : face.leftLocal) {
            if (local >= 0 && local < static_cast<int>(leftTet.dof.size())) {
                const int dof = leftTet.dof[static_cast<size_t>(local)];
                if (dof >= 0 && dof < static_cast<int>(mesh.nodes.size()) && !mesh.nodes[static_cast<size_t>(dof)].dirichlet) {
                    keyByDof[static_cast<size_t>(dof)] = key;
                }
            }
        }
        for (int local : face.rightLocal) {
            if (local >= 0 && local < static_cast<int>(rightTet.dof.size())) {
                const int dof = rightTet.dof[static_cast<size_t>(local)];
                if (dof >= 0 && dof < static_cast<int>(mesh.nodes.size()) && !mesh.nodes[static_cast<size_t>(dof)].dirichlet) {
                    keyByDof[static_cast<size_t>(dof)] = key;
                }
            }
        }
    }

    SubdomainConstantCoarseSpace space;
    space.globalSize = mesh.nodes.size();
    space.globalToMode.assign(mesh.nodes.size(), -1);
    std::map<std::string, int> modeByKey;
    for (size_t dof = 0; dof < keyByDof.size(); ++dof) {
        if (keyByDof[dof].empty()) {
            continue;
        }
        auto found = modeByKey.find(keyByDof[dof]);
        if (found == modeByKey.end()) {
            const int mode = static_cast<int>(space.dofsByMode.size());
            found = modeByKey.emplace(keyByDof[dof], mode).first;
            space.dofsByMode.push_back({});
            space.invNormByMode.push_back(1.0);
        }
        const int mode = found->second;
        space.globalToMode[dof] = mode;
        space.dofsByMode[static_cast<size_t>(mode)].push_back(static_cast<int>(dof));
    }
    for (size_t mode = 0; mode < space.dofsByMode.size(); ++mode) {
        const size_t count = space.dofsByMode[mode].size();
        space.invNormByMode[mode] = count > 0 ? 1.0 / std::sqrt(static_cast<double>(count)) : 1.0;
    }
    return space;
}

static SubdomainConstantCoarseSpace buildRequestedCoarseSpace(const Mesh& mesh,
                                                              const CaseConfig& physics,
                                                              const std::string& coarseSpaceName)
{
    if (coarseSpaceName == "subdomain_material_constant") {
        return buildSubdomainMaterialConstantCoarseSpace(mesh, physics);
    }
    if (coarseSpaceName == "interface_constant") {
        return buildInterfaceConstantCoarseSpace(mesh);
    }
    return buildSubdomainConstantCoarseSpace(mesh);
}

static std::vector<std::vector<double>> buildSubdomainConstantCoarseVectors(const Mesh& mesh)
{
    const SubdomainConstantCoarseSpace space = buildSubdomainConstantCoarseSpace(mesh);
    std::vector<std::vector<double>> vectors(static_cast<size_t>(space.dim()),
                                             std::vector<double>(space.globalSize, 0.0));
    for (size_t mode = 0; mode < space.dofsByMode.size(); ++mode) {
        const double invNorm = space.invNormByMode[mode];
        for (int dof : space.dofsByMode[mode]) {
            vectors[mode][static_cast<size_t>(dof)] = invNorm;
        }
    }
    return vectors;
}

class CoarseCorrectedRasIlutPreconditioner {
public:
    CoarseCorrectedRasIlutPreconditioner(const Mesh& mesh,
                                         const SparseMatrix& a,
                                         int overlap,
                                         double dropTolerance,
                                         int fillFactor,
                                         std::vector<std::vector<double>> coarseVectors,
                                         bool diagonalScaling = false,
                                         double diagScalingEps = 1.0e-30)
        : local_(mesh, a, overlap, dropTolerance, fillFactor, diagonalScaling, diagScalingEps),
          coarseVectors_(std::move(coarseVectors))
    {
        const auto setupStart = std::chrono::steady_clock::now();
        buildCoarseMatrix(a);
        coarseSetupSeconds_ =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - setupStart).count();
    }

    void apply(const std::vector<double>& r, std::vector<double>& z) const
    {
        local_.apply(r, z);
        if (coarseVectors_.empty()) {
            return;
        }
        const auto coarseSolveStart = std::chrono::steady_clock::now();
        std::vector<double> rhs(coarseVectors_.size(), 0.0);
        for (size_t mode = 0; mode < coarseVectors_.size(); ++mode) {
            const std::vector<double>& vector = coarseVectors_[mode];
            double dot = 0.0;
            for (size_t i = 0; i < r.size(); ++i) {
                dot += vector[i] * r[i];
            }
            rhs[mode] = dot;
        }
        std::vector<double> coeffs;
        if (!solveSmallDenseSystem(coarseMatrix_, rhs, coeffs)) {
            coarseResidualNorm_ = std::numeric_limits<double>::infinity();
            coarseSolveSeconds_ +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - coarseSolveStart).count();
            return;
        }
        for (size_t mode = 0; mode < coarseVectors_.size(); ++mode) {
            const double alpha = coeffs[mode];
            const std::vector<double>& vector = coarseVectors_[mode];
            for (size_t i = 0; i < z.size(); ++i) {
                z[i] += alpha * vector[i];
            }
        }
        coarseResidualNorm_ = denseResidualNorm(coarseMatrix_, coeffs, rhs);
        coarseSolveSeconds_ +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - coarseSolveStart).count();
    }

    size_t memoryBytes() const
    {
        size_t bytes = local_.memoryBytes();
        for (const auto& vector : coarseVectors_) {
            bytes += vector.size() * sizeof(double);
        }
        for (const auto& row : coarseMatrix_) {
            bytes += row.size() * sizeof(double);
        }
        return bytes;
    }

    int coarseDim() const { return static_cast<int>(coarseVectors_.size()); }
    double coarseSetupSeconds() const { return coarseSetupSeconds_; }
    double coarseSolveSeconds() const { return coarseSolveSeconds_; }
    size_t coarseMatrixNnz() const { return coarseMatrixNnz_; }
    double coarseResidualNorm() const { return coarseResidualNorm_; }

private:
    RasIlutPreconditioner local_;
    std::vector<std::vector<double>> coarseVectors_;
    std::vector<std::vector<double>> coarseMatrix_;
    size_t coarseMatrixNnz_ = 0;
    double coarseSetupSeconds_ = 0.0;
    mutable double coarseSolveSeconds_ = 0.0;
    mutable double coarseResidualNorm_ = std::numeric_limits<double>::quiet_NaN();

    void buildCoarseMatrix(const SparseMatrix& a)
    {
        coarseMatrix_.assign(coarseVectors_.size(), std::vector<double>(coarseVectors_.size(), 0.0));
        coarseMatrixNnz_ = 0;
        for (size_t col = 0; col < coarseVectors_.size(); ++col) {
            const std::vector<double> az = a.multiply(coarseVectors_[col]);
            for (size_t row = 0; row < coarseVectors_.size(); ++row) {
                double dot = 0.0;
                for (size_t i = 0; i < az.size(); ++i) {
                    dot += coarseVectors_[row][i] * az[i];
                }
                coarseMatrix_[row][col] = dot;
                if (std::abs(dot) > 0.0) {
                    ++coarseMatrixNnz_;
                }
            }
        }
    }

    static double denseResidualNorm(const std::vector<std::vector<double>>& matrix,
                                    const std::vector<double>& x,
                                    const std::vector<double>& rhs)
    {
        double sum = 0.0;
        for (size_t row = 0; row < matrix.size(); ++row) {
            double ax = 0.0;
            for (size_t col = 0; col < x.size(); ++col) {
                ax += matrix[row][col] * x[col];
            }
            const double residual = ax - rhs[row];
            sum += residual * residual;
        }
        return std::sqrt(sum);
    }
};

template <typename LocalPreconditioner>
class CoarseCorrectedPreconditioner {
public:
    CoarseCorrectedPreconditioner(LocalPreconditioner local,
                                  const SparseMatrix& a,
                                  std::vector<std::vector<double>> coarseVectors)
        : local_(std::move(local)),
          coarseVectors_(std::move(coarseVectors))
    {
        const auto setupStart = std::chrono::steady_clock::now();
        buildCoarseMatrix(a);
        coarseSetupSeconds_ =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - setupStart).count();
    }

    void apply(const std::vector<double>& r, std::vector<double>& z)
    {
        local_.apply(r, z);
        if (coarseVectors_.empty()) {
            return;
        }
        const auto coarseSolveStart = std::chrono::steady_clock::now();
        std::vector<double> rhs(coarseVectors_.size(), 0.0);
        for (size_t mode = 0; mode < coarseVectors_.size(); ++mode) {
            const std::vector<double>& vector = coarseVectors_[mode];
            double dot = 0.0;
            for (size_t i = 0; i < r.size(); ++i) {
                dot += vector[i] * r[i];
            }
            rhs[mode] = dot;
        }
        std::vector<double> coeffs;
        if (!solveSmallDenseSystem(coarseMatrix_, rhs, coeffs)) {
            coarseResidualNorm_ = std::numeric_limits<double>::infinity();
            coarseSolveSeconds_ +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - coarseSolveStart).count();
            return;
        }
        for (size_t mode = 0; mode < coarseVectors_.size(); ++mode) {
            const double alpha = coeffs[mode];
            const std::vector<double>& vector = coarseVectors_[mode];
            for (size_t i = 0; i < z.size(); ++i) {
                z[i] += alpha * vector[i];
            }
        }
        coarseResidualNorm_ = denseResidualNorm(coarseMatrix_, coeffs, rhs);
        coarseSolveSeconds_ +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - coarseSolveStart).count();
    }

    size_t memoryBytes() const
    {
        size_t bytes = local_.memoryBytes();
        for (const auto& vector : coarseVectors_) {
            bytes += vector.size() * sizeof(double);
        }
        for (const auto& row : coarseMatrix_) {
            bytes += row.size() * sizeof(double);
        }
        return bytes;
    }

    int coarseDim() const { return static_cast<int>(coarseVectors_.size()); }
    double coarseSetupSeconds() const { return coarseSetupSeconds_; }
    double coarseSolveSeconds() const { return coarseSolveSeconds_; }
    size_t coarseMatrixNnz() const { return coarseMatrixNnz_; }
    double coarseResidualNorm() const { return coarseResidualNorm_; }

private:
    LocalPreconditioner local_;
    std::vector<std::vector<double>> coarseVectors_;
    std::vector<std::vector<double>> coarseMatrix_;
    size_t coarseMatrixNnz_ = 0;
    double coarseSetupSeconds_ = 0.0;
    double coarseSolveSeconds_ = 0.0;
    double coarseResidualNorm_ = std::numeric_limits<double>::quiet_NaN();

    void buildCoarseMatrix(const SparseMatrix& a)
    {
        coarseMatrix_.assign(coarseVectors_.size(), std::vector<double>(coarseVectors_.size(), 0.0));
        coarseMatrixNnz_ = 0;
        for (size_t col = 0; col < coarseVectors_.size(); ++col) {
            const std::vector<double> az = a.multiply(coarseVectors_[col]);
            for (size_t row = 0; row < coarseVectors_.size(); ++row) {
                double dot = 0.0;
                for (size_t i = 0; i < az.size(); ++i) {
                    dot += coarseVectors_[row][i] * az[i];
                }
                coarseMatrix_[row][col] = dot;
                if (std::abs(dot) > 0.0) {
                    ++coarseMatrixNnz_;
                }
            }
        }
    }

    static double denseResidualNorm(const std::vector<std::vector<double>>& matrix,
                                    const std::vector<double>& x,
                                    const std::vector<double>& rhs)
    {
        double sum = 0.0;
        for (size_t row = 0; row < matrix.size(); ++row) {
            double ax = 0.0;
            for (size_t col = 0; col < x.size(); ++col) {
                ax += matrix[row][col] * x[col];
            }
            const double residual = ax - rhs[row];
            sum += residual * residual;
        }
        return std::sqrt(sum);
    }
};

template <typename LocalPreconditioner>
class SubdomainConstantCoarseCorrectedPreconditioner {
public:
    SubdomainConstantCoarseCorrectedPreconditioner(LocalPreconditioner local,
                                                   const SparseMatrix& a,
                                                   SubdomainConstantCoarseSpace coarseSpace)
        : local_(std::move(local)),
          coarseSpace_(std::move(coarseSpace))
    {
        const auto setupStart = std::chrono::steady_clock::now();
        buildCoarseMatrix(a);
        rhs_.resize(coarseSpace_.dofsByMode.size());
        coeffs_.resize(coarseSpace_.dofsByMode.size());
        coarseSetupSeconds_ =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - setupStart).count();
    }

    void apply(const std::vector<double>& r, std::vector<double>& z)
    {
        local_.apply(r, z);
        if (coarseSpace_.dofsByMode.empty()) {
            return;
        }

        const auto coarseSolveStart = std::chrono::steady_clock::now();
        std::fill(rhs_.begin(), rhs_.end(), 0.0);
        for (size_t mode = 0; mode < coarseSpace_.dofsByMode.size(); ++mode) {
            double sum = 0.0;
            for (int dof : coarseSpace_.dofsByMode[mode]) {
                sum += r[static_cast<size_t>(dof)];
            }
            rhs_[mode] = coarseSpace_.invNormByMode[mode] * sum;
        }

        if (!solveSmallDenseSystem(coarseMatrix_, rhs_, coeffs_)) {
            coarseResidualNorm_ = std::numeric_limits<double>::infinity();
            coarseSolveSeconds_ +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - coarseSolveStart).count();
            return;
        }

        for (size_t mode = 0; mode < coarseSpace_.dofsByMode.size(); ++mode) {
            const double alpha = coeffs_[mode] * coarseSpace_.invNormByMode[mode];
            for (int dof : coarseSpace_.dofsByMode[mode]) {
                z[static_cast<size_t>(dof)] += alpha;
            }
        }

        coarseResidualNorm_ = denseResidualNorm(coarseMatrix_, coeffs_, rhs_);
        coarseSolveSeconds_ +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - coarseSolveStart).count();
    }

    size_t memoryBytes() const
    {
        size_t bytes = local_.memoryBytes() + coarseSpace_.memoryBytes()
            + rhs_.size() * sizeof(double)
            + coeffs_.size() * sizeof(double);
        for (const auto& row : coarseMatrix_) {
            bytes += row.size() * sizeof(double);
        }
        return bytes;
    }

    int coarseDim() const { return coarseSpace_.dim(); }
    double coarseSetupSeconds() const { return coarseSetupSeconds_; }
    double coarseSolveSeconds() const { return coarseSolveSeconds_; }
    size_t coarseMatrixNnz() const { return coarseMatrixNnz_; }
    double coarseResidualNorm() const { return coarseResidualNorm_; }

private:
    LocalPreconditioner local_;
    SubdomainConstantCoarseSpace coarseSpace_;
    std::vector<std::vector<double>> coarseMatrix_;
    std::vector<double> rhs_;
    std::vector<double> coeffs_;
    size_t coarseMatrixNnz_ = 0;
    double coarseSetupSeconds_ = 0.0;
    double coarseSolveSeconds_ = 0.0;
    double coarseResidualNorm_ = std::numeric_limits<double>::quiet_NaN();

    void buildCoarseMatrix(const SparseMatrix& a)
    {
        const size_t dim = coarseSpace_.dofsByMode.size();
        coarseMatrix_.assign(dim, std::vector<double>(dim, 0.0));
        if (dim == 0) {
            return;
        }
        for (int row = 0; row < a.n; ++row) {
            const int rowMode = coarseSpace_.globalToMode[static_cast<size_t>(row)];
            if (rowMode < 0) {
                continue;
            }
            const double rowScale = coarseSpace_.invNormByMode[static_cast<size_t>(rowMode)];
            for (int k = a.rowPtr[static_cast<size_t>(row)];
                 k < a.rowPtr[static_cast<size_t>(row + 1)];
                 ++k) {
                const int col = a.colInd[static_cast<size_t>(k)];
                if (col < 0 || col >= static_cast<int>(coarseSpace_.globalToMode.size())) {
                    continue;
                }
                const int colMode = coarseSpace_.globalToMode[static_cast<size_t>(col)];
                if (colMode < 0) {
                    continue;
                }
                const double colScale = coarseSpace_.invNormByMode[static_cast<size_t>(colMode)];
                coarseMatrix_[static_cast<size_t>(rowMode)][static_cast<size_t>(colMode)] +=
                    rowScale * a.values[static_cast<size_t>(k)] * colScale;
            }
        }
        coarseMatrixNnz_ = 0;
        for (const auto& row : coarseMatrix_) {
            for (double value : row) {
                if (std::abs(value) > 0.0) {
                    ++coarseMatrixNnz_;
                }
            }
        }
    }

    static double denseResidualNorm(const std::vector<std::vector<double>>& matrix,
                                    const std::vector<double>& x,
                                    const std::vector<double>& rhs)
    {
        double sum = 0.0;
        for (size_t row = 0; row < matrix.size(); ++row) {
            double ax = 0.0;
            for (size_t col = 0; col < x.size(); ++col) {
                ax += matrix[row][col] * x[col];
            }
            const double residual = ax - rhs[row];
            sum += residual * residual;
        }
        return std::sqrt(sum);
    }
};

static double vectorDot(const std::vector<double>& a, const std::vector<double>& b)
{
    double result = 0.0;
#ifdef _OPENMP
    const int threadCount = static_cast<int>(solverParallelWorkers());
#pragma omp parallel for schedule(static) reduction(+:result) num_threads(threadCount)
    for (long long i = 0; i < static_cast<long long>(a.size()); ++i) {
        result += a[static_cast<size_t>(i)] * b[static_cast<size_t>(i)];
    }
#else
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
#endif
    return result;
}

static bool vectorHasNonFinite(const std::vector<double>& x)
{
    for (double value : x) {
        if (!std::isfinite(value)) {
            return true;
        }
    }
    return false;
}

static double relativeResidualNorm(const SparseMatrix& a,
                                   const std::vector<double>& x,
                                   const std::vector<double>& b)
{
    const std::vector<double> ax = a.multiply(x);
    double residualSquared = 0.0;
    double rhsSquared = 0.0;
    for (size_t i = 0; i < b.size(); ++i) {
        const double r = b[i] - ax[i];
        residualSquared += r * r;
        rhsSquared += b[i] * b[i];
    }
    return std::sqrt(residualSquared) / std::sqrt(std::max(1.0e-300, rhsSquared));
}

static void recordRhsDiagnostics(SolverStatistics* stats,
                                 const std::vector<double>& rhs,
                                 const CaseConfig& physics,
                                 const std::vector<double>& fixedAdjust)
{
    if (stats == nullptr) {
        return;
    }
    stats->rhsNorm = std::sqrt(std::max(0.0, vectorDot(rhs, rhs)));
    stats->rhsMin = rhs.empty() ? std::numeric_limits<double>::quiet_NaN()
                                : *std::min_element(rhs.begin(), rhs.end());
    stats->rhsMax = rhs.empty() ? std::numeric_limits<double>::quiet_NaN()
                                : *std::max_element(rhs.begin(), rhs.end());
    stats->rhsNonzeroCount = 0;
    stats->rhsNearZeroCount = 0;
    const double nearZero = 1.0e-14 * std::max(1.0, stats->rhsNorm);
    for (double value : rhs) {
        if (std::abs(value) > 0.0) {
            ++stats->rhsNonzeroCount;
        }
        if (std::abs(value) <= nearZero) {
            ++stats->rhsNearZeroCount;
        }
    }
    stats->heatSourceTotal = 0.0;
    for (const HeatSource& source : physics.heatSources) {
        stats->heatSourceTotal += source.heatRateW;
    }
    stats->dirichletRhsContributionNorm = fixedAdjust.empty()
        ? 0.0
        : std::sqrt(std::max(0.0, vectorDot(fixedAdjust, fixedAdjust)));
}

static void recordInitialResidualDiagnostics(SolverStatistics* stats,
                                             const SparseMatrix& a,
                                             const std::vector<double>& b,
                                             const std::vector<double>& x,
                                             double relativeTolerance)
{
    if (stats == nullptr) {
        return;
    }
    const std::vector<double> ax = a.multiply(x);
    std::vector<double> r(b.size(), 0.0);
    parallelFor(b.size(), [&](size_t i) {
        r[i] = b[i] - ax[i];
    });
    stats->rhsNorm = std::sqrt(std::max(0.0, vectorDot(b, b)));
    stats->initialXNorm = std::sqrt(std::max(0.0, vectorDot(x, x)));
    stats->initialAxNorm = std::sqrt(std::max(0.0, vectorDot(ax, ax)));
    stats->initialResidualNorm = std::sqrt(std::max(0.0, vectorDot(r, r)));
    stats->initialRelativeResidual =
        stats->initialResidualNorm / std::sqrt(std::max(1.0e-300, vectorDot(b, b)));
    stats->solverTolerance = relativeTolerance;
    stats->zeroIterationDueToInitialResidual = stats->initialRelativeResidual < relativeTolerance;
    std::cout << "Initial residual diagnostics [" << stats->name << "]:"
              << " initial_guess_type=" << stats->initialGuessType
              << " ||b||2=" << stats->rhsNorm
              << " ||x0||2=" << stats->initialXNorm
              << " ||A*x0||2=" << stats->initialAxNorm
              << " ||r0||2=" << stats->initialResidualNorm
              << " rel=" << stats->initialRelativeResidual
              << " tol=" << relativeTolerance << "\n";
    if (stats->zeroIterationDueToInitialResidual) {
        std::cout << "WARNING: initial relative residual is already below solver tolerance; zero iterations may be expected.\n";
    }
}

static std::vector<double> diagonalPreconditioner(const SparseMatrix& a)
{
    std::vector<double> invDiag(static_cast<size_t>(a.size()), 1.0);
    parallelFor(invDiag.size(), [&](size_t i) {
        const double diag = a.diagonal(static_cast<int>(i));
        if (std::abs(diag) > 1.0e-30) {
            invDiag[static_cast<size_t>(i)] = 1.0 / diag;
        }
    });
    return invDiag;
}

static double updateSolutionAndResidual(std::vector<double>& x,
                                        std::vector<double>& r,
                                        const std::vector<double>& p,
                                        const std::vector<double>& ap,
                                        double alpha)
{
    double rr = 0.0;
#ifdef _OPENMP
    const int threadCount = static_cast<int>(solverParallelWorkers());
#pragma omp parallel for schedule(static) reduction(+:rr) num_threads(threadCount)
    for (long long i = 0; i < static_cast<long long>(x.size()); ++i) {
        const size_t idx = static_cast<size_t>(i);
        x[idx] += alpha * p[idx];
        r[idx] -= alpha * ap[idx];
        rr += r[idx] * r[idx];
    }
#else
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += alpha * p[i];
        r[i] -= alpha * ap[i];
        rr += r[i] * r[i];
    }
#endif
    return rr;
}

static double applyDiagonalPreconditionerAndDot(const std::vector<double>& invDiag,
                                                const std::vector<double>& r,
                                                std::vector<double>& z)
{
    double rz = 0.0;
    z.resize(r.size());
#ifdef _OPENMP
    const int threadCount = static_cast<int>(solverParallelWorkers());
#pragma omp parallel for schedule(static) reduction(+:rz) num_threads(threadCount)
    for (long long i = 0; i < static_cast<long long>(r.size()); ++i) {
        const size_t idx = static_cast<size_t>(i);
        z[idx] = invDiag[idx] * r[idx];
        rz += r[idx] * z[idx];
    }
#else
    for (size_t i = 0; i < r.size(); ++i) {
        z[i] = invDiag[i] * r[i];
        rz += r[i] * z[i];
    }
#endif
    return rz;
}

static std::vector<double> preconditionedConjugateGradient(const SparseMatrix& a,
                                                           const std::vector<double>& b,
                                                           std::vector<double> x,
                                                           int& iterations,
                                                           SolverStatistics* stats = nullptr,
                                                           int maxIterations = 5000,
                                                           double relativeTolerance = 1.0e-10)
{
    const std::vector<double> invDiag = diagonalPreconditioner(a);
    std::vector<double> ax = a.multiply(x);
    std::vector<double> r(b.size(), 0.0);
    std::vector<double> z(b.size(), 0.0);
    parallelFor(b.size(), [&](size_t i) {
        r[i] = b[i] - ax[i];
    });
    if (stats != nullptr) {
        recordInitialResidualDiagnostics(stats, a, b, x, relativeTolerance);
    }
    double rzOld = applyDiagonalPreconditionerAndDot(invDiag, r, z);
    if (stats != nullptr) {
        stats->pcgRTr = vectorDot(r, r);
        stats->pcgRMz = rzOld;
        stats->pcgResidualNorm = std::sqrt(std::max(0.0, stats->pcgRTr));
    }
    if (!std::isfinite(rzOld) || rzOld <= 0.0) {
        if (stats != nullptr) {
            stats->status = "failed";
            stats->failureReason = "Jacobi-PCG breakdown: initial rMz is nonpositive or nonfinite";
            stats->breakdownIteration = 0;
            stats->finalRelativeResidual = relativeResidualNorm(a, x, b);
        }
        iterations = 0;
        return x;
    }
    std::vector<double> p = z;
    const double rhsNorm = std::sqrt(std::max(1.0e-60, vectorDot(b, b)));
    const double tolerance = relativeTolerance * std::max(1.0, rhsNorm);
    const double toleranceSquared = tolerance * tolerance;

    iterations = 0;
    for (; iterations < maxIterations; ++iterations) {
        const std::vector<double> ap = a.multiply(p);
        const double denom = vectorDot(p, ap);
        if (stats != nullptr) {
            stats->breakdownIteration = iterations;
            stats->pcgPAp = denom;
            stats->pcgRTr = vectorDot(r, r);
            stats->pcgRMz = rzOld;
            stats->pcgResidualNorm = std::sqrt(std::max(0.0, stats->pcgRTr));
            stats->pcgPHasNonFinite = vectorHasNonFinite(p);
            stats->pcgApHasNonFinite = vectorHasNonFinite(ap);
        }
        if (!std::isfinite(denom) || denom <= 0.0) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = "Jacobi-PCG breakdown: pAp is nonpositive or nonfinite";
            }
            break;
        }
        const double alpha = rzOld / denom;
        if (stats != nullptr) {
            stats->pcgAlpha = alpha;
        }
        if (!std::isfinite(alpha)) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = "Jacobi-PCG breakdown: alpha is nonfinite";
            }
            break;
        }
        const double residualSquared = updateSolutionAndResidual(x, r, p, ap, alpha);
        if (residualSquared < toleranceSquared) {
            ++iterations;
            if (stats != nullptr) {
                stats->status = "success";
            }
            break;
        }
        const double rzNew = applyDiagonalPreconditionerAndDot(invDiag, r, z);
        if (!std::isfinite(rzNew) || rzNew <= 0.0) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = "Jacobi-PCG breakdown: rMz is nonpositive or nonfinite";
                stats->pcgRTr = residualSquared;
                stats->pcgRMz = rzNew;
                stats->pcgResidualNorm = std::sqrt(std::max(0.0, residualSquared));
            }
            break;
        }
        const double beta = rzNew / rzOld;
        if (stats != nullptr) {
            stats->pcgBeta = beta;
        }
        if (!std::isfinite(beta)) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = "Jacobi-PCG breakdown: beta is nonfinite";
            }
            break;
        }
        parallelFor(p.size(), [&](size_t i) {
            p[i] = z[i] + beta * p[i];
        });
        rzOld = rzNew;
    }
    if (stats != nullptr) {
        stats->finalRelativeResidual = relativeResidualNorm(a, x, b);
        if (stats->status == "not_run") {
            stats->status = stats->finalRelativeResidual <= relativeTolerance ? "success" : "failed";
            if (stats->status == "failed" && stats->failureReason.empty()) {
                stats->failureReason = "Jacobi-PCG did not converge";
            }
        }
    }
    return x;
}

template <typename Preconditioner>
static std::vector<double> blockJacobiPcg(const SparseMatrix& a,
                                          const std::vector<double>& b,
                                          std::vector<double> x,
                                          Preconditioner& preconditioner,
                                          int& iterations,
                                          SolverStatistics* stats = nullptr,
                                          const std::string& solverLabel = "BJ-PCG",
                                          int maxIterations = 5000,
                                          double relativeTolerance = 1.0e-10,
                                          const std::filesystem::path& residualHistoryPath = {},
                                          bool detailedIterationDiagnostics = true)
{
    std::vector<double> ax = a.multiply(x);
    std::vector<double> r(b.size(), 0.0);
    parallelFor(b.size(), [&](size_t i) {
        r[i] = b[i] - ax[i];
    });
    if (stats != nullptr) {
        recordInitialResidualDiagnostics(stats, a, b, x, relativeTolerance);
    }

    std::vector<double> z;
    auto preconditionerApplyStart = std::chrono::steady_clock::now();
    preconditioner.apply(r, z);
    auto preconditionerApplyEnd = std::chrono::steady_clock::now();
    if (stats != nullptr) {
        stats->preconditionerApplySeconds +=
            std::chrono::duration<double>(preconditionerApplyEnd - preconditionerApplyStart).count();
        ++stats->preconditionerApplyCalls;
    }
    std::vector<double> p = z;
    double rzOld = vectorDot(r, z);
    if (stats != nullptr) {
        stats->pcgRTr = vectorDot(r, r);
        stats->pcgRMz = rzOld;
        stats->pcgResidualNorm = std::sqrt(std::max(0.0, stats->pcgRTr));
    }
    if (!std::isfinite(rzOld) || rzOld <= 0.0
        || (detailedIterationDiagnostics && vectorHasNonFinite(z))) {
        if (stats != nullptr) {
            stats->status = "failed";
            stats->failureReason = solverLabel + " breakdown: initial rMz is nonpositive or nonfinite";
            stats->breakdownIteration = 0;
            stats->finalRelativeResidual = relativeResidualNorm(a, x, b);
        }
        iterations = 0;
        return x;
    }
    const double rhsNorm = std::sqrt(std::max(1.0e-60, vectorDot(b, b)));
    const double tolerance = relativeTolerance * std::max(1.0, rhsNorm);
    const double toleranceSquared = tolerance * tolerance;
    // ========== ���� 1���������� CSV �ļ� ==========
    std::ofstream residualHistory;
    if (!residualHistoryPath.empty()) {
        residualHistory.open(residualHistoryPath);
        if (residualHistory.is_open()) {
            residualHistory << "iteration,relative_residual,residual_norm,rhs_norm\n";
            const double initialResidual = std::sqrt(std::max(0.0, vectorDot(r, r)));
            residualHistory << 0 << ','
                            << std::scientific << std::setprecision(16)
                            << initialResidual / rhsNorm << ','
                            << initialResidual << ','
                            << rhsNorm << '\n';
        }
    }
    // ==
    iterations = 0;
    for (; iterations < maxIterations; ++iterations) {
        const std::vector<double> ap = a.multiply(p);
        const double denom = vectorDot(p, ap);
        if (stats != nullptr) {
            stats->breakdownIteration = iterations;
            stats->pcgPAp = denom;
            if (detailedIterationDiagnostics) {
                stats->pcgRTr = vectorDot(r, r);
                stats->pcgResidualNorm = std::sqrt(std::max(0.0, stats->pcgRTr));
                stats->pcgPHasNonFinite = vectorHasNonFinite(p);
                stats->pcgApHasNonFinite = vectorHasNonFinite(ap);
            }
            stats->pcgRMz = rzOld;
        }
        if (!std::isfinite(denom) || denom <= 0.0) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = solverLabel + " breakdown: pAp is nonpositive or nonfinite";
            }
            break;
        }
        const double alpha = rzOld / denom;
        if (stats != nullptr) {
            stats->pcgAlpha = alpha;
        }
        if (!std::isfinite(alpha)) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = solverLabel + " breakdown: alpha is nonfinite";
            }
            break;
        }
        const double residualSquared = updateSolutionAndResidual(x, r, p, ap, alpha);

        // ========== ���� 2������ǰ�в�д�� CSV �ļ� ==========
        if (residualHistory.is_open()) {
            const double residualNorm = std::sqrt(std::max(0.0, residualSquared));
            residualHistory << (iterations + 1) << ','
                            << std::scientific << std::setprecision(16)
                            << residualNorm / rhsNorm << ','
                            << residualNorm << ','
                            << rhsNorm << '\n';
        }
        // ===================================================

        if (!std::isfinite(residualSquared)
            || (detailedIterationDiagnostics && (vectorHasNonFinite(x) || vectorHasNonFinite(r)))) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = solverLabel + " breakdown: residual or solution became nonfinite";
            }
            break;
        }
        if (residualSquared < toleranceSquared) {
            ++iterations;
            if (stats != nullptr) {
                stats->status = "success";
            }
            break;
        }
        preconditionerApplyStart = std::chrono::steady_clock::now();
        preconditioner.apply(r, z);
        preconditionerApplyEnd = std::chrono::steady_clock::now();
        if (stats != nullptr) {
            stats->preconditionerApplySeconds +=
                std::chrono::duration<double>(preconditionerApplyEnd - preconditionerApplyStart).count();
            ++stats->preconditionerApplyCalls;
        }
        const double rzNew = vectorDot(r, z);
        if (!std::isfinite(rzNew) || rzNew <= 0.0
            || (detailedIterationDiagnostics && vectorHasNonFinite(z))) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = solverLabel + " breakdown: rMz is nonpositive or nonfinite";
                stats->pcgRTr = residualSquared;
                stats->pcgRMz = rzNew;
                stats->pcgResidualNorm = std::sqrt(std::max(0.0, residualSquared));
            }
            break;
        }
        const double beta = rzNew / rzOld;
        if (stats != nullptr) {
            stats->pcgBeta = beta;
        }
        if (!std::isfinite(beta)) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = solverLabel + " breakdown: beta is nonfinite";
            }
            break;
        }
        parallelFor(p.size(), [&](size_t i) {
            p[i] = z[i] + beta * p[i];
        });
        rzOld = rzNew;
    }
    if (stats != nullptr) {
        stats->finalRelativeResidual = relativeResidualNorm(a, x, b);
        if (stats->status == "not_run") {
            stats->status = stats->finalRelativeResidual <= 1.0e-10 ? "success" : "failed";
            if (stats->status == "failed" && stats->failureReason.empty()) {
                stats->failureReason = solverLabel + " did not converge";
            }
        }
    }
    return x;
}

static double vectorNorm(const std::vector<double>& x)
{
    return std::sqrt(std::max(0.0, vectorDot(x, x)));
}

static void axpy(double alpha, const std::vector<double>& x, std::vector<double>& y)
{
    parallelFor(y.size(), [&](size_t i) {
        y[i] += alpha * x[i];
    });
}

template <typename Preconditioner>
static std::vector<double> flexibleGmres(const SparseMatrix& a,
                                         const std::vector<double>& b,
                                         std::vector<double> x,
                                         Preconditioner& preconditioner,
                                         int& iterations,
                                         SolverStatistics* stats = nullptr,
                                         const std::string& solverLabel = "FGMRES",
                                         int maxIterations = 50,
                                         double relativeTolerance = 1.0e-10,
                                         int restart = 30)
{
    const int n = a.size();
    const int restartLength = std::max(1, restart);
    const double rhsNorm = std::max(1.0e-300, vectorNorm(b));
    const double tolerance = relativeTolerance * std::max(1.0, rhsNorm);
    iterations = 0;

    auto computeResidual = [&]() {
        const std::vector<double> ax = a.multiply(x);
        std::vector<double> r(b.size(), 0.0);
        parallelFor(b.size(), [&](size_t i) {
            r[i] = b[i] - ax[i];
        });
        return r;
    };

    std::vector<double> r = computeResidual();
    double beta = vectorNorm(r);
    if (stats != nullptr) {
        recordInitialResidualDiagnostics(stats, a, b, x, relativeTolerance);
        stats->pcgResidualNorm = beta;
        stats->finalRelativeResidual = beta / rhsNorm;
    }
    if (beta <= tolerance) {
        if (stats != nullptr) {
            stats->status = "success";
        }
        return x;
    }

    while (iterations < maxIterations) {
        std::vector<std::vector<double>> v(static_cast<size_t>(restartLength + 1),
                                           std::vector<double>(static_cast<size_t>(n), 0.0));
        std::vector<std::vector<double>> z(static_cast<size_t>(restartLength),
                                           std::vector<double>(static_cast<size_t>(n), 0.0));
        std::vector<std::vector<double>> h(static_cast<size_t>(restartLength + 1),
                                           std::vector<double>(static_cast<size_t>(restartLength), 0.0));
        std::vector<double> cs(static_cast<size_t>(restartLength), 0.0);
        std::vector<double> sn(static_cast<size_t>(restartLength), 0.0);
        std::vector<double> g(static_cast<size_t>(restartLength + 1), 0.0);

        parallelFor(v[0].size(), [&](size_t i) {
            v[0][i] = r[i] / beta;
        });
        g[0] = beta;

        int innerCount = 0;
        bool converged = false;
        for (; innerCount < restartLength && iterations < maxIterations; ++innerCount) {
            auto preconditionerApplyStart = std::chrono::steady_clock::now();
            preconditioner.apply(v[static_cast<size_t>(innerCount)], z[static_cast<size_t>(innerCount)]);
            auto preconditionerApplyEnd = std::chrono::steady_clock::now();
            if (stats != nullptr) {
                stats->preconditionerApplySeconds +=
                    std::chrono::duration<double>(preconditionerApplyEnd - preconditionerApplyStart).count();
                ++stats->preconditionerApplyCalls;
            }
            if (vectorHasNonFinite(z[static_cast<size_t>(innerCount)])) {
                if (stats != nullptr) {
                    stats->status = "failed";
                    stats->failureReason = solverLabel + " breakdown: preconditioner produced NaN/Inf";
                    stats->breakdownIteration = iterations;
                }
                return x;
            }

            std::vector<double> w = a.multiply(z[static_cast<size_t>(innerCount)]);
            if (vectorHasNonFinite(w)) {
                if (stats != nullptr) {
                    stats->status = "failed";
                    stats->failureReason = solverLabel + " breakdown: A*z produced NaN/Inf";
                    stats->breakdownIteration = iterations;
                }
                return x;
            }

            for (int i = 0; i <= innerCount; ++i) {
                h[static_cast<size_t>(i)][static_cast<size_t>(innerCount)] =
                    vectorDot(w, v[static_cast<size_t>(i)]);
                axpy(-h[static_cast<size_t>(i)][static_cast<size_t>(innerCount)],
                     v[static_cast<size_t>(i)],
                     w);
            }
            h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)] = vectorNorm(w);
            if (h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)] > 0.0) {
                parallelFor(w.size(), [&](size_t i) {
                    v[static_cast<size_t>(innerCount + 1)][i] =
                        w[i] / h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)];
                });
            }

            for (int i = 0; i < innerCount; ++i) {
                const double hi = h[static_cast<size_t>(i)][static_cast<size_t>(innerCount)];
                const double hip1 = h[static_cast<size_t>(i + 1)][static_cast<size_t>(innerCount)];
                h[static_cast<size_t>(i)][static_cast<size_t>(innerCount)] =
                    cs[static_cast<size_t>(i)] * hi + sn[static_cast<size_t>(i)] * hip1;
                h[static_cast<size_t>(i + 1)][static_cast<size_t>(innerCount)] =
                    -sn[static_cast<size_t>(i)] * hi + cs[static_cast<size_t>(i)] * hip1;
            }

            const double hii = h[static_cast<size_t>(innerCount)][static_cast<size_t>(innerCount)];
            const double hnext = h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)];
            const double denom = std::hypot(hii, hnext);
            if (!(denom > 0.0) || !std::isfinite(denom)) {
                if (stats != nullptr) {
                    stats->status = "failed";
                    stats->failureReason = solverLabel + " breakdown: Arnoldi Hessenberg pivot is zero/nonfinite";
                    stats->breakdownIteration = iterations;
                }
                return x;
            }
            cs[static_cast<size_t>(innerCount)] = hii / denom;
            sn[static_cast<size_t>(innerCount)] = hnext / denom;
            h[static_cast<size_t>(innerCount)][static_cast<size_t>(innerCount)] = denom;
            h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)] = 0.0;

            const double gi = g[static_cast<size_t>(innerCount)];
            g[static_cast<size_t>(innerCount)] = cs[static_cast<size_t>(innerCount)] * gi;
            g[static_cast<size_t>(innerCount + 1)] = -sn[static_cast<size_t>(innerCount)] * gi;

            ++iterations;
            const double relativeResidual = std::abs(g[static_cast<size_t>(innerCount + 1)]) / rhsNorm;
            if (stats != nullptr) {
                stats->breakdownIteration = iterations;
                stats->pcgResidualNorm = std::abs(g[static_cast<size_t>(innerCount + 1)]);
                stats->finalRelativeResidual = relativeResidual;
            }
            if (relativeResidual <= relativeTolerance) {
                ++innerCount;
                converged = true;
                break;
            }
        }

        std::vector<double> y(static_cast<size_t>(innerCount), 0.0);
        for (int i = innerCount - 1; i >= 0; --i) {
            double sum = g[static_cast<size_t>(i)];
            for (int j = i + 1; j < innerCount; ++j) {
                sum -= h[static_cast<size_t>(i)][static_cast<size_t>(j)] * y[static_cast<size_t>(j)];
            }
            const double diag = h[static_cast<size_t>(i)][static_cast<size_t>(i)];
            if (!(std::abs(diag) > 0.0) || !std::isfinite(diag)) {
                if (stats != nullptr) {
                    stats->status = "failed";
                    stats->failureReason = solverLabel + " breakdown: singular least-squares system";
                    stats->breakdownIteration = iterations;
                }
                return x;
            }
            y[static_cast<size_t>(i)] = sum / diag;
        }
        for (int i = 0; i < innerCount; ++i) {
            axpy(y[static_cast<size_t>(i)], z[static_cast<size_t>(i)], x);
        }
        if (vectorHasNonFinite(x)) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = solverLabel + " breakdown: solution became NaN/Inf";
                stats->breakdownIteration = iterations;
                stats->finalRelativeResidual = relativeResidualNorm(a, x, b);
            }
            return x;
        }

        r = computeResidual();
        beta = vectorNorm(r);
        if (stats != nullptr) {
            stats->pcgResidualNorm = beta;
            stats->finalRelativeResidual = beta / rhsNorm;
        }
        if (beta <= tolerance) {
            if (stats != nullptr) {
                stats->status = "success";
                stats->finalRelativeResidual = beta / rhsNorm;
            }
            return x;
        }
        (void)converged;
        if (!(beta > 0.0) || !std::isfinite(beta)) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = solverLabel + " breakdown: restarted residual is zero/nonfinite";
                stats->breakdownIteration = iterations;
            }
            return x;
        }
    }

    if (stats != nullptr) {
        stats->finalRelativeResidual = relativeResidualNorm(a, x, b);
        if (stats->status == "not_run") {
            stats->status = stats->finalRelativeResidual <= relativeTolerance ? "success" : "failed";
            if (stats->status == "failed") {
                stats->failureReason = solverLabel + " did not converge";
            }
        }
    }
    return x;
}

static double interfaceAverageJump(const Mesh& mesh, const std::vector<double>& temperature);

static std::vector<double> initialTemperatureVector(const Mesh& mesh, const CaseConfig& physics)
{
    std::string type = physics.initialGuessType;
    if (type == "random_small") {
        type = "random-small";
    }
    std::vector<double> temperature(mesh.nodes.size(), physics.initialTemperature);
    if (type == "zero") {
        std::fill(temperature.begin(), temperature.end(), 0.0);
    } else if (type == "constant") {
        std::fill(temperature.begin(), temperature.end(), physics.initialTemperature);
    } else if (type == "random-small") {
        std::mt19937 rng(1234567u);
        std::uniform_real_distribution<double> dist(-1.0e-3, 1.0e-3);
        for (double& value : temperature) {
            value = dist(rng);
        }
    } else {
        parallelFor(mesh.nodes.size(), [&](size_t i) {
            if (mesh.nodes[static_cast<size_t>(i)].dirichlet) {
                temperature[static_cast<size_t>(i)] = mesh.nodes[static_cast<size_t>(i)].dirichletValue;
            }
        });
    }
    return temperature;
}

static std::vector<double> makeTransientRhs(const SparseMatrix& mass,
                                            const SparseMatrix& system,
                                            const std::vector<double>& temperature,
                                            const std::vector<double>& source,
                                            const std::vector<double>& fixedAdjust,
                                            const CaseConfig& physics)
{
    std::vector<double> rhs = mass.multiply(temperature);
    if (physics.timeIntegrator == "crank_nicolson") {
        const std::vector<double> implicitPart = system.multiply(temperature);
        parallelFor(rhs.size(), [&](size_t i) {
            rhs[i] = 2.0 * rhs[i] / physics.timeStep - implicitPart[i] + source[i];
            // `system` has already had strong Dirichlet columns eliminated.  For
            // Crank--Nicolson the old-time operator contributes a second A_fd*g
            // term.  applyDirichletRhs() below removes the new-time one; remove
            // the old-time one here.  Backward Euler has no analogous explicit
            // stiffness contribution.
            if (!fixedAdjust.empty()) {
                rhs[i] -= fixedAdjust[i];
            }
        });
    } else {
        parallelFor(rhs.size(), [&](size_t i) {
            rhs[i] = rhs[i] / physics.timeStep + source[i];
        });
    }
    return rhs;
}

template <typename StepSolve>
static std::vector<double> runTransientSolver(const std::string& solverName,
                                              const Mesh& mesh,
                                              const CaseConfig& physics,
                                              const SparseMatrix& mass,
                                              const SparseMatrix& system,
                                              const std::vector<double>& source,
                                              const std::vector<double>& fixedAdjust,
                                              SolverStatistics& stats,
                                              StepSolve&& stepSolve,
                                              std::vector<double>* finalRhs = nullptr)
{
    stats.name = solverName;
    stats.totalIterations = 0;
    stats.maxIterations = 0;
    stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
    stats.initialGuessType = physics.initialGuessType;
    if (stats.workingSetBeforeBytes == 0) {
        stats.workingSetBeforeBytes = currentWorkingSetBytes();
    }
    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
    std::vector<double> temperature = initialTemperatureVector(mesh, physics);
    std::vector<double> lastRhs;
    const auto solveStart = std::chrono::steady_clock::now();

    for (int step = 1; step <= physics.timeSteps; ++step) {
        std::vector<double> rhs = makeTransientRhs(mass, system, temperature, source, fixedAdjust, physics);
        applyDirichletRhs(mesh, fixedAdjust, rhs);
        lastRhs = rhs;
        if (finalRhs != nullptr) {
            *finalRhs = rhs;
        }

        recordRhsDiagnostics(&stats, rhs, physics, fixedAdjust);
        int iterations = 0;
        std::vector<double> x0 = physics.disableWarmStart
            ? initialTemperatureVector(mesh, physics)
            : temperature;
        temperature = stepSolve(system, rhs, std::move(x0), iterations);
        stats.totalIterations += iterations;
        stats.maxIterations = std::max(stats.maxIterations, iterations);

        if (step == 1 || step % 10 == 0 || step == physics.timeSteps) {
            const auto minmax = std::minmax_element(temperature.begin(), temperature.end());
            const double avg = std::accumulate(temperature.begin(), temperature.end(), 0.0) / static_cast<double>(temperature.size());
            std::cout << "[" << solverName << "] step " << std::setw(3) << step
                      << "  time=" << std::setw(10) << step * physics.timeStep
                      << " s  Tmin=" << std::setw(12) << *minmax.first
                      << "  Tmax=" << std::setw(12) << *minmax.second
                      << "  Tavg=" << std::setw(12) << avg
                      << "  interface_avg_jump=" << interfaceAverageJump(mesh, temperature)
                      << "  iterations=" << iterations << "\n";
        }
    }

    const auto solveEnd = std::chrono::steady_clock::now();
    stats.solveSeconds = std::chrono::duration<double>(solveEnd - solveStart).count();
    stats.workingSetAfterBytes = currentWorkingSetBytes();
    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
    if (!lastRhs.empty()) {
        stats.trueRelativeResidual = relativeResidualNorm(system, temperature, lastRhs);
    }
    if (stats.status == "not_run") {
        stats.status = vectorHasNonFinite(temperature) ? "failed" : "success";
        if (stats.status == "failed") {
            stats.failureReason = "solution contains NaN/Inf";
        }
    }
    return temperature;
}

template <typename StepSolve>
static std::vector<double> runSteadySolver(const std::string& solverName,
                                           const Mesh& mesh,
                                           const CaseConfig& physics,
                                           const SparseMatrix& system,
                                           const std::vector<double>& source,
                                           const std::vector<double>& fixedAdjust,
                                           SolverStatistics& stats,
                                           StepSolve&& stepSolve,
                                           std::vector<double>* finalRhs = nullptr)
{
    stats.name = solverName;
    stats.totalIterations = 0;
    stats.maxIterations = 0;
    stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
    stats.initialGuessType = physics.initialGuessType;
    if (stats.workingSetBeforeBytes == 0) {
        stats.workingSetBeforeBytes = currentWorkingSetBytes();
    }
    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());

    std::vector<double> rhs = source;
    applyDirichletRhs(mesh, fixedAdjust, rhs);
    if (finalRhs != nullptr) {
        *finalRhs = rhs;
    }
    recordRhsDiagnostics(&stats, rhs, physics, fixedAdjust);
    std::vector<double> temperature = initialTemperatureVector(mesh, physics);

    int iterations = 0;
    const auto solveStart = std::chrono::steady_clock::now();
    temperature = stepSolve(system, rhs, std::move(temperature), iterations);
    const auto solveEnd = std::chrono::steady_clock::now();

    stats.totalIterations = iterations;
    stats.maxIterations = iterations;
    stats.solveSeconds = std::chrono::duration<double>(solveEnd - solveStart).count();
    stats.workingSetAfterBytes = currentWorkingSetBytes();
    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
    if (std::isnan(stats.finalRelativeResidual)) {
        stats.finalRelativeResidual = relativeResidualNorm(system, temperature, rhs);
    }
    stats.trueRelativeResidual = relativeResidualNorm(system, temperature, rhs);
    if (stats.status == "not_run") {
        stats.status = vectorHasNonFinite(temperature) ? "failed" : "success";
        if (stats.status == "failed") {
            stats.failureReason = "solution contains NaN/Inf";
        }
    }

    const auto minmax = std::minmax_element(temperature.begin(), temperature.end());
    const double avg = std::accumulate(temperature.begin(), temperature.end(), 0.0) / static_cast<double>(temperature.size());
    std::cout << "[" << solverName << "] steady"
              << "  Tmin=" << std::setw(12) << *minmax.first
              << "  Tmax=" << std::setw(12) << *minmax.second
              << "  Tavg=" << std::setw(12) << avg
              << "  interface_avg_jump=" << interfaceAverageJump(mesh, temperature)
              << "  iterations=" << iterations << "\n";
    return temperature;
}

template <typename StepSolve>
static std::vector<double> runAnalysisSolver(bool steady,
                                             const std::string& solverName,
                                             const Mesh& mesh,
                                             const CaseConfig& physics,
                                             const SparseMatrix& mass,
                                             const SparseMatrix& system,
                                             const std::vector<double>& source,
                                             const std::vector<double>& fixedAdjust,
                                             SolverStatistics& stats,
                                             StepSolve&& stepSolve,
                                             std::vector<double>* finalRhs = nullptr)
{
    if (steady) {
        return runSteadySolver(solverName, mesh, physics, system, source, fixedAdjust, stats,
                               std::forward<StepSolve>(stepSolve), finalRhs);
    }
    return runTransientSolver(solverName, mesh, physics, mass, system, source, fixedAdjust, stats,
                              std::forward<StepSolve>(stepSolve), finalRhs);
}
