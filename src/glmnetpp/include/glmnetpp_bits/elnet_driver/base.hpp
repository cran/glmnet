#pragma once
#include <algorithm>
#include <glmnetpp_bits/util/exceptions.hpp>

namespace glmnetpp {

struct ElnetDriverBase
{
    template <class VType>
    void normalize_penalty(VType&& vq) const {
        if (vq.maxCoeff() <= 0) throw util::non_positive_penalty_error();
        vq.array() = vq.array().max(0.0);
        vq *= vq.size() / vq.sum();
    }

    template <class JDType, class JUType>
    void init_inclusion(const JDType& jd, JUType&& ju) const {
        if (jd(0) > 0) {
            for (int i = 1; i < jd(0) + 1; ++i) {
                ju[jd(i)-1] = false;
            }
        }
        // can't find true value in ju
        if (std::find_if(ju.begin(), ju.end(), [](auto x) { return x;}) == ju.end()) {
            throw util::all_excluded_error();
        } 
    }
};

} // namespace glmnetpp
