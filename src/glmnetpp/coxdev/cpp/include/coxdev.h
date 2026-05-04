#include <cstddef>
#include <limits>
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>

#ifdef PY_INTERFACE
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#endif

#ifdef R_INTERFACE
#include <RcppEigen.h>
#endif

#define EIGEN_REF Eigen::Ref
#define MAKE_MAP_Xd(y) Eigen::Map<Eigen::VectorXd>((y).data(), (y).size())
#define MAKE_MAP_Xi(y) Eigen::Map<Eigen::VectorXi>((y).data(), (y).size())

#if defined(PY_INTERFACE)
#define ERROR_MSG(x) throw std::runtime_error(x)
#elif defined(R_INTERFACE)
#define ERROR_MSG(x) Rcpp::stop(x)
#else
#define ERROR_MSG(x) throw std::runtime_error(x)
#endif

struct CoxContext {
    Eigen::Ref<const Eigen::VectorXi> event_order;
    Eigen::Ref<const Eigen::VectorXi> start_order;
    Eigen::Ref<const Eigen::VectorXi> status;
    Eigen::Ref<const Eigen::VectorXi> first;
    Eigen::Ref<const Eigen::VectorXi> last;
    Eigen::Ref<const Eigen::VectorXd> scaling;
    Eigen::Ref<const Eigen::VectorXi> event_map;
    Eigen::Ref<const Eigen::VectorXi> start_map;
    bool have_start_times;
    bool efron;
};

class CoxDeviance {
public:
    CoxDeviance(const Eigen::VectorXd& start,
                const Eigen::VectorXd& event,
                const Eigen::VectorXi& status,
                const Eigen::VectorXd& weight,
                bool efron);

    double compute_deviance(const Eigen::VectorXd& eta,
                            const Eigen::VectorXd& sample_weight);

    void compute_hessian_matvec(const Eigen::VectorXd& arg,
                                Eigen::VectorXd& out);

    // Getters for buffers (for result construction or debugging)
    const Eigen::VectorXd& get_gradient() const { return grad_buffer; }
    const Eigen::VectorXd& get_diag_hessian() const { return diag_hessian_buffer; }
    const Eigen::VectorXd& get_linear_predictor() const { return linear_predictor; }
    const Eigen::VectorXd& get_sample_weight() const { return sample_weight; }
    double get_loglik_sat() const { return loglik_sat; }

    // Getters for preprocessed data
    const Eigen::VectorXi& get_event_order() const { return event_order; }
    const Eigen::VectorXi& get_start_order() const { return start_order; }
    const Eigen::VectorXi& get_first() const { return _first; }
    const Eigen::VectorXi& get_last() const { return _last; }
    const Eigen::VectorXi& get_start_map() const { return _start_map; }
    const Eigen::VectorXi& get_event_map() const { return _event_map; }
    const Eigen::VectorXd& get_scaling() const { return _scaling; }
    const Eigen::VectorXi& get_status() const { return _status; }
    const Eigen::VectorXd& get_event() const { return _event; }
    const Eigen::VectorXd& get_start() const { return _start; }

private:
    // Preprocessing results
    Eigen::VectorXd _start, _event, _scaling;
    Eigen::VectorXi _status, _first, _last, _start_map, _event_map;
    Eigen::VectorXi event_order, start_order;
    
    // Buffers
    Eigen::VectorXd T_1_term, T_2_term, forward_scratch_buffer, w_avg_buffer, exp_w_buffer;
    Eigen::VectorXd grad_buffer, diag_hessian_buffer, diag_part_buffer, hess_matvec_buffer;
    Eigen::VectorXd linear_predictor, sample_weight;
    
    std::vector<Eigen::VectorXd> event_reorder_buffers;
    std::vector<Eigen::VectorXd> forward_cumsum_buffers;
    std::vector<Eigen::VectorXd> reverse_cumsum_buffers;
    std::vector<Eigen::VectorXd> risk_sum_buffers;

    double loglik_sat;
    bool have_start_times;
    bool _efron;

    void setup_buffers(int n);

    CoxContext get_context() const {
        return {event_order, start_order, _status, _first, _last, _scaling, _event_map, _start_map, have_start_times, _efron};
    }
};

class StratifiedCoxDeviance {
public:
    StratifiedCoxDeviance(const Eigen::VectorXd& start,
                          const Eigen::VectorXd& event,
                          const Eigen::VectorXi& status,
                          const Eigen::VectorXi& strata,
                          const Eigen::VectorXd& weight,
                          bool efron);

    double compute_deviance(const Eigen::VectorXd& eta,
                            const Eigen::VectorXd& sample_weight);

    void compute_hessian_matvec(const Eigen::VectorXd& arg,
                                Eigen::VectorXd& out);

    const Eigen::VectorXd& get_gradient() const { return grad_buffer; }
    const Eigen::VectorXd& get_diag_hessian() const { return diag_hessian_buffer; }
    const Eigen::VectorXd& get_linear_predictor() const { return linear_predictor; }
    const Eigen::VectorXd& get_sample_weight() const { return sample_weight; }
    double get_loglik_sat() const { return loglik_sat; }

    const std::vector<int>& get_unique_strata() const { return unique_strata; }
    const std::vector<std::vector<int>>& get_stratum_indices() const { return stratum_indices; }
    const std::vector<std::shared_ptr<CoxDeviance>>& get_cox_devs() const { return cox_devs; }

private:
    std::vector<int> unique_strata;
    std::vector<std::vector<int>> stratum_indices;
    std::vector<std::shared_ptr<CoxDeviance>> cox_devs;

    Eigen::VectorXd grad_buffer;
    Eigen::VectorXd diag_hessian_buffer;
    Eigen::VectorXd linear_predictor;
    Eigen::VectorXd sample_weight;
    double loglik_sat;
};

// Standalone function declarations
void forward_cumsum(const Eigen::Ref<const Eigen::VectorXd> sequence,
		    Eigen::Ref<Eigen::VectorXd> output);

void reverse_cumsums(const Eigen::Ref<const Eigen::VectorXd> sequence,
                     Eigen::Ref<Eigen::VectorXd> event_buffer,
                     Eigen::Ref<Eigen::VectorXd> start_buffer,
                     const Eigen::Ref<const Eigen::VectorXi> event_order,
                     const Eigen::Ref<const Eigen::VectorXi> start_order,
		     bool do_event = false,
		     bool do_start = false);

void to_native_from_event(Eigen::Ref<Eigen::VectorXd> arg,
			  const Eigen::Ref<const Eigen::VectorXi> event_order,
			  Eigen::Ref<Eigen::VectorXd> reorder_buffer);

void to_event_from_native(const Eigen::Ref<const Eigen::VectorXd> arg,
                          const Eigen::Ref<const Eigen::VectorXi> event_order,
                          Eigen::Ref<Eigen::VectorXd> reorder_buffer);

void forward_prework(const Eigen::Ref<const Eigen::VectorXi> status,
                     const Eigen::Ref<const Eigen::VectorXd> w_avg,
                     const Eigen::Ref<const Eigen::VectorXd> scaling,
                     const Eigen::Ref<const Eigen::VectorXd> risk_sums,
                     int i,
                     int j,
                     Eigen::Ref<Eigen::VectorXd> moment_buffer,
		     const Eigen::Ref<const Eigen::VectorXd> arg,		     
                     bool use_w_avg = true);

double compute_sat_loglik(const Eigen::Ref<const Eigen::VectorXi> first,
			  const Eigen::Ref<const Eigen::VectorXi> last,
			  const Eigen::Ref<const Eigen::VectorXd> weight,
			  const Eigen::Ref<const Eigen::VectorXi> event_order,
			  const Eigen::Ref<const Eigen::VectorXi> status,
			  const Eigen::Ref<const Eigen::VectorXd> scaling,
			  Eigen::Ref<Eigen::VectorXd> W_status,
			  bool efron);

void sum_over_events(const CoxContext& ctx,
                     Eigen::Ref<Eigen::VectorXd> C_arg,
                     Eigen::Ref<Eigen::VectorXd> C_arg_scale,
                     Eigen::Ref<Eigen::VectorXd> forward_scratch_buffer,
                     Eigen::Ref<Eigen::VectorXd> value_buffer);

void sum_over_risk_set(const Eigen::Ref<const Eigen::VectorXd> arg,
                       const CoxContext& ctx,
                       Eigen::Ref<Eigen::VectorXd> risk_sum_buffer,
                       Eigen::Ref<Eigen::VectorXd> event_cumsum,
                       Eigen::Ref<Eigen::VectorXd> start_cumsum);

double cox_dev(const Eigen::Ref<const Eigen::VectorXd> eta,
	       const Eigen::Ref<const Eigen::VectorXd> sample_weight,
	       const Eigen::Ref<const Eigen::VectorXd> exp_w,
               const CoxContext& ctx,
	       double loglik_sat,
	       Eigen::Ref<Eigen::VectorXd> T_1_term,
	       Eigen::Ref<Eigen::VectorXd> T_2_term,
	       Eigen::Ref<Eigen::VectorXd> grad_buffer,
	       Eigen::Ref<Eigen::VectorXd> diag_hessian_buffer,
	       Eigen::Ref<Eigen::VectorXd> diag_part_buffer,
	       Eigen::Ref<Eigen::VectorXd> w_avg_buffer,
               Eigen::Ref<Eigen::VectorXd> eta_event,
               Eigen::Ref<Eigen::VectorXd> w_event,
               Eigen::Ref<Eigen::VectorXd> exp_eta_w_event,
               Eigen::Ref<Eigen::VectorXd> risk_sums,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum0,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum1,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum2,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum3,
               Eigen::Ref<Eigen::VectorXd> forward_cumsum4,
	       Eigen::Ref<Eigen::VectorXd> forward_scratch_buffer,
               Eigen::Ref<Eigen::VectorXd> event_cumsum,
               Eigen::Ref<Eigen::VectorXd> start_cumsum);

void hessian_matvec(const Eigen::Ref<const Eigen::VectorXd> arg,
                    const CoxContext& ctx,
                    Eigen::Ref<Eigen::VectorXd> risk_sums,
                    Eigen::Ref<const Eigen::VectorXd> diag_part,
                    Eigen::Ref<const Eigen::VectorXd> w_avg,
                    Eigen::Ref<const Eigen::VectorXd> exp_w,
                    Eigen::Ref<Eigen::VectorXd> risk_sum_arg,
                    Eigen::Ref<Eigen::VectorXd> forward_cumsum0,
                    Eigen::Ref<Eigen::VectorXd> forward_cumsum1,
                    Eigen::Ref<Eigen::VectorXd> forward_scratch_buffer,
                    Eigen::Ref<Eigen::VectorXd> event_cumsum,
                    Eigen::Ref<Eigen::VectorXd> start_cumsum,
                    Eigen::Ref<Eigen::VectorXd> hess_matvec_buffer);

void setup_preprocess(const Eigen::Ref<const Eigen::VectorXd> start,
                      const Eigen::Ref<const Eigen::VectorXd> event,
                      const Eigen::Ref<const Eigen::VectorXi> status,
                      const Eigen::Ref<const Eigen::VectorXd> weight,
                      Eigen::VectorXd& _start,
                      Eigen::VectorXd& _event,
                      Eigen::VectorXi& _status,
                      Eigen::VectorXi& _first,
                      Eigen::VectorXi& _last,
                      Eigen::VectorXd& _scaling,
                      Eigen::VectorXi& _start_map,
                      Eigen::VectorXi& _event_map,
                      Eigen::VectorXi& event_order,
                      Eigen::VectorXi& start_order);
