#pragma once
#include <testutil/base_fixture.hpp>
#include <testutil/data_util.hpp>
#include <glmnetpp_bits/elnet_driver/standardize.hpp>    // separately unit-tested

namespace glmnetpp {

struct ElnetBasePack
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

    ElnetBasePack(
            int _maxit,
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
        , alpha(_alpha), flmin(_flmin), xv(_xv), ulam(_ulam)
        , vp(_vp), cl(_cl), ju(_ju)
    {
        int p = xv.size();
        ao.setZero(nx, nlam);
        ia.setZero(p);
        kin.setZero(nlam);
        rsqo.setZero(nlam); 
        almo.setZero(nlam);
    }

    virtual void fit() =0;
    virtual void fit_old() =0;
};

struct elnet_base_fixture: 
    base_fixture,
    testing::WithParamInterface<
        std::tuple<size_t, size_t, size_t, size_t, size_t, double, double> >
{
    void SetUp() override
    {
        size_t seed, n, p;
        std::tie(seed, n, p, maxit, nlam, alpha, flmin) = GetParam();
        DataGen dgen(seed); // hehe
        X = dgen.make_X(n, p);
        auto beta = dgen.make_beta(p);
        y = dgen.make_y(X, beta);
        auto w = dgen.make_w(n);
        ju = dgen.make_ju(p);
        vp = dgen.make_vp(p);
        cl = dgen.make_cl(p);
        nx = dgen.make_nx(p);
        ne = dgen.make_ne(p);
        ulam = dgen.make_ulam(nlam);
        xv.setOnes(p);
        ia.setZero(p);

        Eigen::VectorXd xm(p), xs(p);
        double ym = 0, ys = 0;
        Standardize1::eval(X, y, w, 1, 1, ju, xm, xs, ym, ys, xv);
    }

protected:
    Eigen::MatrixXd X, cl;
    Eigen::VectorXd y, xv, ulam, vp;
    Eigen::VectorXi ju, ia;
    int nx, ne, maxit, nlam;
    double alpha, flmin;

    void check_pack(const ElnetBasePack& actual,
                    const ElnetBasePack& expected) const
    {
        EXPECT_EQ(actual.lmu, expected.lmu);
        EXPECT_EQ(actual.nlp, expected.nlp);
        EXPECT_EQ(actual.jerr, expected.jerr);

        expect_eq_vec(actual.kin, expected.kin);

        // Legacy is 1-indexed, so ia should be shifted by 1.
        // Only applies up to the indicies corresponding to active variables.
        // The best I can think of testing this is that the absolute distance is off by at most 1.
        expect_near_vec(actual.ia, expected.ia, 1);

        expect_float_eq_mat(actual.ao, expected.ao);
        expect_float_eq_vec(actual.rsqo, expected.rsqo);
        expect_float_eq_vec(actual.almo, expected.almo);
    }
};

struct spelnet_base_fixture
    : elnet_base_fixture
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

        Eigen::VectorXd xv(p);
        xm.setZero(p); xs.setZero(p);
        double ym, ys;
        SpStandardize1::eval(X, y, w, true, true, ju, xm, xs, ym, ys, xv);
    }

protected: 
    Eigen::SparseMatrix<double> X;
    Eigen::VectorXd xm, xs, w;
};

} // namespace glmnetpp
