#pragma once
#include <testutil/data_util.hpp>
#include <testutil/thread.hpp>
#include <testutil/base_fixture.hpp>
#include <Eigen/Core>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>

namespace glmnetpp {

struct GaussianPack
{
    const double thr = 1e-14;
    const int maxit, nx, ne, nlam;
    const double alpha, flmin;

    // will be modified
    Eigen::MatrixXd ao;
    Eigen::VectorXi ia;
    Eigen::VectorXi kin;
    Eigen::VectorXd rsqo; 
    Eigen::VectorXd almo;
    int nlp = 0, jerr = 0, lmu = 0;
    
    const Eigen::VectorXd& xv;
    const Eigen::VectorXd& ulam;
    const Eigen::VectorXd& vp;
    const Eigen::MatrixXd& cl;
    const Eigen::VectorXi& ju;

    GaussianPack(int _maxit,
                 int _nx,
                 int _ne,
                 int _nlam,
                 double _alpha,
                 double _flmin,
                 const Eigen::VectorXd& _xv,
                 const Eigen::VectorXd& _ulam,
                 const Eigen::VectorXd& _vp,
                 const Eigen::MatrixXd& _cl,
                 const Eigen::VectorXi& _ju)
        : maxit(_maxit), nx(_nx), ne(_ne), nlam(_nlam)
        , alpha(_alpha), flmin(_flmin)
        , ao(_nx, _nlam)
        , ia(_ju.size())
        , kin(_nlam)
        , rsqo(_nlam)
        , almo(_nlam)
        , xv(_xv), ulam(_ulam), vp(_vp), cl(_cl), ju(_ju)
    {
        ao.setZero();
        ia.setZero();
        kin.setZero();
        rsqo.setZero();
        almo.setZero();
    }

    virtual void fit() =0;
    virtual void fit_old() =0;
};

struct gaussian_fixture
    : base_fixture
    , testing::WithParamInterface<
        std::tuple<int, int, int, int, int, double, double> 
      >
{
    void SetUp() override
    {
        size_t seed, n, p;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X(n, p);
        auto beta = dgen.make_beta(p);
        y = dgen.make_y(X, beta);
        ju = dgen.make_ju(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
        xv.setOnes(p);
        ia.setZero(p);
    }

protected:
    Eigen::MatrixXd X, cl;
    Eigen::VectorXd y, xv, ulam, vp;
    Eigen::VectorXi ju, ia;
    int nx, ne, maxit, nlam;
    double alpha, flmin;

    void check_pack(const GaussianPack& actual,
                    const GaussianPack& expected)
    {
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);

        expect_eq_vec(actual.kin, expected.kin);

        // Legacy is 1-indexed, so ia should be shifted by 1.
        // Only applies up to the indicies corresponding to active variables.
        // The best I can think of testing this is that the absolute distance is off by at most 1.
        expect_near_vec(actual.ia, expected.ia, 1);

        expect_near_mat(actual.ao, expected.ao, 1e-15);
        expect_near_vec(actual.rsqo, expected.rsqo, 1e-15);

        // This check loosens expect_near_vec.
        // When almo is large (>= 1), check for exact equality.
        // Otherwise, if too small, put a tolerance of 1e-15.
        EXPECT_EQ(actual.almo.size(), expected.almo.size());
        for (int i = 0; i < actual.almo.size(); ++i) {
            if (actual.almo[i] < 1) {
                EXPECT_NEAR(actual.almo[i], expected.almo[i], 1e-15);
            } else {
                EXPECT_DOUBLE_EQ(actual.almo[i], expected.almo[i]);
            }
        }
    }
};

struct sp_gaussian_fixture
    : gaussian_fixture
{
    void SetUp() override
    {
        size_t seed, n, p;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X_sparse(n, p);
        auto beta = dgen.make_beta(p);
        y = dgen.make_y(X, beta);
        w = dgen.make_w(n);
        ju = dgen.make_ju(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
        xv.setOnes(p);
        ia.setZero(p);

        // compute xm and xs
        Eigen::MatrixXd X_dense = X;
        xm = X_dense.transpose() * w;
        xs = ((X_dense.array().square().matrix().transpose() * w).array()
            - xm.array().square()).sqrt().matrix();
    }

protected:
    Eigen::SparseMatrix<double> X;
    Eigen::VectorXd xm, xs, w;
};

} // namespace glmnetpp
