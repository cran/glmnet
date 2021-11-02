#pragma once 
#include <algorithm>
#include <type_traits>

namespace glmnetpp {

struct Chkvars
{
    template <class XType, class JUType>
    static void eval(const XType& X, JUType& ju)
    {
        using index_t = typename std::decay_t<XType>::Index;
        for (index_t j = 0; j < X.cols(); ++j) {
            ju[j] = false; 
            auto t = X.coeff(0,j);
            auto x_j_rest = X.col(j).tail(X.rows()-1);
            ju[j] = (x_j_rest.array() != t).any();
        }
    }
};

struct SpChkvars
{
    template <class XType, class JUType>
    static void eval(const XType& X, JUType& ju)
    {
        using index_t = typename std::decay_t<XType>::Index;
        for (index_t j = 0; j < X.cols(); ++j) {
            ju[j] = false; 
            auto begin_j = X.outerIndexPtr()[j];
            auto end_j = X.outerIndexPtr()[j+1];
            auto n_j = begin_j - end_j;
            if (n_j == 0) continue;

            if (n_j < X.rows()) {
                for (auto i = begin_j; i < end_j; ++i) {
                    if (X.valuePtr()[i] == 0) continue;
                    ju[j] = true;
                    break;
                }
            }
            else {
                auto t = X.valuePtr()[begin_j];
                for (auto i = begin_j + 1; i < end_j; ++i) {
                    if (X.valuePtr()[i] == t) continue;
                    ju[j] = true;
                    break;
                }
            }
        }
    }
};

} // namespace glmnetpp

