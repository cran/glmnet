#pragma once
#include <glmnetpp_bits/util/exceptions.hpp>
#include <glmnetpp_bits/elnet_path/decl.hpp>

namespace glmnetpp {

template <class ElnetPathDerived>
struct ElnetPathBase
{
private:
    using derived_t = ElnetPathDerived;

    derived_t& self() { return static_cast<derived_t&>(*this); }
    const derived_t& self() const { return static_cast<const derived_t&>(*this); }

protected:

    template <class PackType>
    void fit(PackType&& pack) const
    {
        using pack_t = std::decay_t<PackType>;
        using int_t = typename pack_t::int_t;

        auto path_config_pack = self().initialize_path(pack);
        const auto& alm = path_config_pack.alm;
        const auto& mnl = path_config_pack.mnl;
        const auto& sml0 = path_config_pack.sml0;
        const auto& rsqmax0 = path_config_pack.rsqmax0;

        const auto& flmin = pack.flmin;
        const auto& ne = pack.ne;
        auto& lmu = pack.lmu;
        auto& ao = pack.ao;
        const auto& ia = pack.ia;
        auto& kin = pack.kin;
        auto& rsqo = pack.rsqo;
        auto& almo = pack.almo;
        auto& jerr = pack.jerr;

        auto&& elnet_point = self().get_elnet_point(pack);

        for (int_t m = 0; m < pack.nlam; ++m) {

            auto point_config_pack = 
                self().initialize_point(m, pack, path_config_pack, elnet_point);
            
            auto rsq0 = elnet_point.rsq(); 

            try {
                elnet_point.fit(point_config_pack);
            } 
            catch (const util::elnet_error& e) {
                jerr = e.err_code(m);
                return;
            } 

            auto n_active = elnet_point.n_active();
            if (n_active > 0) { 
                for (int_t j = 0; j < n_active; ++j) {
                    ao(j, m) = elnet_point.beta(ia(j)-1);
                }
            } 
            auto rsq = elnet_point.rsq();
            kin(m) = n_active;
            rsqo(m) = rsq; 
            almo(m) = alm; 
            lmu = m + 1;
            if (lmu < mnl || flmin >= 1.0) continue; 
            int_t me = 0; 
            for (int_t j = 0; j < n_active; ++j) {
                 if (ao(j, m)) ++me;
            }
            if ((me > ne) ||
                (rsq - rsq0 < sml0 * rsq) ||
                (rsq > rsqmax0)) break; 
        }
    }
};

} // namespace glmnetpp
