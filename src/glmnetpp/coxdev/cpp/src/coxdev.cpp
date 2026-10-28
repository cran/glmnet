#include "coxdev.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>

#ifndef DEBUG_PRINT
#define DEBUG_PRINT(x) std::cout << x << std::endl
#endif

void forward_cumsum(const Eigen::Ref<const Eigen::VectorXd> sequence,
		    Eigen::Ref<Eigen::VectorXd> output)
{
  if (sequence.size() + 1 != output.size()) {
    ERROR_MSG("forward_cumsum: output size must be one longer than input's.");
  }
      
  double sum = 0.0;
  output(0) = sum;
  for (int i = 1; i < output.size(); ++i) {
    sum = sum + sequence(i - 1);
    output(i) = sum;
  }
}

void reverse_cumsums(const Eigen::Ref<const Eigen::VectorXd> sequence,
                     Eigen::Ref<Eigen::VectorXd> event_buffer,
                     Eigen::Ref<Eigen::VectorXd> start_buffer,
                     const Eigen::Ref<const Eigen::VectorXi> event_order,
                     const Eigen::Ref<const Eigen::VectorXi> start_order,
		     bool do_event,
		     bool do_start)
{
  double sum = 0.0;
  
  int n = sequence.size();
  if (do_event) {
    if (sequence.size() + 1 != event_buffer.size()) {
      ERROR_MSG("reverse_cumsums: event_buffer size must be one more than input's.");
    }
    event_buffer(n) = sum;
    for (int i = n - 1; i >= 0;  --i) {
      sum = sum + sequence(event_order(i));
      event_buffer(i) = sum;
    }
  }

  if (do_start) {
    if (sequence.size() + 1 != start_buffer.size()) {
      ERROR_MSG("reverse_cumsums: event_buffer size must be one more than input's.");
    }
    sum = 0.0;
    start_buffer(n) = sum;
    for (int i = n - 1; i >= 0;  --i) {
      sum = sum + sequence(start_order(i));
      start_buffer(i) = sum;
    }
  }
}

void to_native_from_event(Eigen::Ref<Eigen::VectorXd> arg,
			  const Eigen::Ref<const Eigen::VectorXi> event_order,
			  Eigen::Ref<Eigen::VectorXd> reorder_buffer)
{
  reorder_buffer = arg;
  for (int i = 0; i < event_order.size(); ++i) {
    arg(event_order(i)) = reorder_buffer(i);
  }
}

void to_event_from_native(const Eigen::Ref<const Eigen::VectorXd> arg,
                          const Eigen::Ref<const Eigen::VectorXi> event_order,
                          Eigen::Ref<Eigen::VectorXd> reorder_buffer)
{
  for (int i = 0; i < event_order.size(); ++i) {
    reorder_buffer(i) = arg(event_order(i));
  }
}

void forward_prework(const Eigen::Ref<const Eigen::VectorXi> status,
                     const Eigen::Ref<const Eigen::VectorXd> w_avg,
                     const Eigen::Ref<const Eigen::VectorXd> scaling,
                     const Eigen::Ref<const Eigen::VectorXd> risk_sums,
                     int i,
                     int j,
                     Eigen::Ref<Eigen::VectorXd> moment_buffer,
		     const Eigen::Ref<const Eigen::VectorXd> arg,		     
                     bool use_w_avg)
{
  bool has_arg = arg.size() > 0;
  for (int k = 0; k < status.size(); ++k) {
    if (use_w_avg) {
      if (w_avg(k) > 0) {
        moment_buffer(k) = status(k) * w_avg(k) * std::pow(scaling(k), i) / std::pow(risk_sums(k), j);
      } else {
        moment_buffer(k) = 0.0;
      }
    } else {
      if (status(k) > 0) {
        moment_buffer(k) = status(k) * std::pow(scaling(k), i) / std::pow(risk_sums(k), j);
      } else {
        moment_buffer(k) = 0.0;
      }
    }
    if (has_arg) {
      moment_buffer(k) *= arg(k);
    }
  }
}

double compute_sat_loglik(const Eigen::Ref<const Eigen::VectorXi> first,
			  const Eigen::Ref<const Eigen::VectorXi> last,
			  const Eigen::Ref<const Eigen::VectorXd> weight,
			  const Eigen::Ref<const Eigen::VectorXi> event_order,
			  const Eigen::Ref<const Eigen::VectorXi> status,
			  const Eigen::Ref<const Eigen::VectorXd> scaling,
			  Eigen::Ref<Eigen::VectorXd> W_status,
			  bool efron)
{
  
  Eigen::VectorXd weight_event_order_times_status(event_order.size());
  for (int i = 0; i < event_order.size(); ++i) {
    weight_event_order_times_status(i) = weight(event_order(i)) * status(i);
  }
  forward_cumsum(MAKE_MAP_Xd(weight_event_order_times_status), W_status);

  Eigen::VectorXd sums(last.size());
  for (int i = 0; i < last.size(); ++i) {
    sums(i) = W_status(last(i) + 1) - W_status(first(i));
  }
  double loglik_sat = 0.0;
  int prev_first = -1;

 for (int i = 0; i < first.size(); ++i) {
    int f = first(i);
    int l = last(i);
    
    // Process each unique cluster of tied failure times only once
    if (f != prev_first) {
      double W_C = sums(i); // W_C: Total valid weight in this cluster

      if ((W_C > 0) && (weight(event_order(i)) > 0)) {
        // 1. Structural Term (Identical to Breslow)
        loglik_sat -= W_C * std::log(W_C);
	
        // 2. Efron Penalty Term
        int K_C_plus = last(i) + 1 - first(i);
        double sum_log_penalty = 0.0;

        // Iterate through the cluster to count K_C^+ and sum the scaling penalty
        for (int j = f; j <= l; ++j) {
          // Only include individuals who actually failed and have positive weight
          // These should all be failing
          if (weight(event_order(j)) > 0 && status(j) == 1) {
            // Using the scaling vector (\sigma) passed as an argument
            // Note: scaling(j) should be m / K_C^+ where m is 0, 1, ..., K_C^+ - 1
            sum_log_penalty += std::log(1.0 - efron * scaling(j));
          }
        }

        // Apply the penalty, weighted by the average weight per valid observation
        if (K_C_plus > 0) {
          loglik_sat -= (W_C / K_C_plus) * sum_log_penalty;
        }
      }
      prev_first = f;
    }
  }
  
  return loglik_sat;
}

void sum_over_events(const CoxContext& ctx,
                     Eigen::Ref<Eigen::VectorXd> C_arg,
                     Eigen::Ref<Eigen::VectorXd> C_arg_scale,
                     Eigen::Ref<Eigen::VectorXd> forward_scratch_buffer,
                     Eigen::Ref<Eigen::VectorXd> value_buffer)
{
  forward_cumsum(forward_scratch_buffer, C_arg); //length=n+1

  if (ctx.have_start_times) {
    for (int i = 0; i < ctx.last.size(); ++i) {
      value_buffer(i) = C_arg(ctx.last(i) + 1) - C_arg(ctx.start_map(i));
    }
  } else {
    for (int i = 0; i < ctx.last.size(); ++i) {
      value_buffer(i) = C_arg(ctx.last(i) + 1);
    }
  }
  if (ctx.efron) {
    forward_scratch_buffer = forward_scratch_buffer.array() * ctx.scaling.array();
    forward_cumsum(forward_scratch_buffer, C_arg_scale); // length=n+1
    for (int i = 0; i < ctx.last.size(); ++i) {
      value_buffer(i) -= (C_arg_scale(ctx.last(i) + 1) - C_arg_scale(ctx.first(i)));
    }
  }
}

void sum_over_risk_set(const Eigen::Ref<const Eigen::VectorXd> arg,
                       const CoxContext& ctx,
                       Eigen::Ref<Eigen::VectorXd> risk_sum_buffer,
                       Eigen::Ref<Eigen::VectorXd> event_cumsum,
                       Eigen::Ref<Eigen::VectorXd> start_cumsum)
{
  reverse_cumsums(arg,
		  event_cumsum,
		  start_cumsum,
		  ctx.event_order,
		  ctx.start_order,
		  true, // do_event 
		  ctx.have_start_times); // do_start
    
  if (ctx.have_start_times) {
    for (int i = 0; i < ctx.first.size(); ++i) {
      risk_sum_buffer(i) = event_cumsum(ctx.first(i)) - start_cumsum(ctx.event_map(i));
    }
  } else {
    for (int i = 0; i < ctx.first.size(); ++i) {
      risk_sum_buffer(i) = event_cumsum(ctx.first(i));
    }
  }
        
  if (ctx.efron) {
    for (int i = 0; i < ctx.first.size(); ++i) {
      risk_sum_buffer(i) = risk_sum_buffer(i) - ( event_cumsum(ctx.first(i)) - event_cumsum(ctx.last(i) + 1) ) * ctx.scaling(i);
    }
  }
}

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
               Eigen::Ref<Eigen::VectorXd> start_cumsum)
{
  to_event_from_native(eta, ctx.event_order, eta_event);
  to_event_from_native(sample_weight, ctx.event_order, w_event);
  to_event_from_native(exp_w, ctx.event_order, exp_eta_w_event);

  sum_over_risk_set(exp_w, ctx, risk_sums, event_cumsum, start_cumsum);

  for (int i = 0; i < w_avg_buffer.size(); ++i) {
    if (ctx.status(i) == 1) {
      w_avg_buffer(i) = (forward_cumsum0(ctx.last(i) + 1) - forward_cumsum0(ctx.first(i))) / ((double) (ctx.last(i) + 1 - ctx.first(i))); 
      } else {
      w_avg_buffer(i) = 0;
    }  
  }

  double loglik_penalty = 0.0;
  for (int i = 0; i < risk_sums.size(); ++i) {
    if (w_avg_buffer(i) > 0) {
      loglik_penalty += std::log(risk_sums(i)) * w_avg_buffer(i) * ctx.status(i);
    }
  }

  double loglik = ( w_event.array() * eta_event.array() * ctx.status.cast<double>().array() ).sum() - loglik_penalty;
    
  Eigen::VectorXd dummy;
  Eigen::Map<Eigen::VectorXd> dummy_map(dummy.data(), dummy.size());  

  forward_prework(ctx.status, w_avg_buffer, ctx.scaling, risk_sums, 0, 1, forward_scratch_buffer, dummy_map, true);
  forward_cumsum(forward_scratch_buffer, forward_cumsum0);
  
  forward_prework(ctx.status, w_avg_buffer, ctx.scaling, risk_sums, 0, 2, forward_scratch_buffer, dummy_map, true);
  forward_cumsum(forward_scratch_buffer, forward_cumsum1);
  
  if (!ctx.efron) {
    if (ctx.have_start_times) {
      for (int i = 0; i < ctx.last.size(); ++i) {
	T_1_term(i) = forward_cumsum0(ctx.last(i) + 1) - forward_cumsum0(ctx.start_map(i));
	T_2_term(i) = forward_cumsum1(ctx.last(i) + 1) - forward_cumsum1(ctx.start_map(i));
      }
    } else {
      for (int i = 0; i < ctx.last.size(); ++i) {
	T_1_term(i) = forward_cumsum0(ctx.last(i) + 1);
	T_2_term(i) = forward_cumsum1(ctx.last(i) + 1);
      }
    }
  } else {
    forward_prework(ctx.status, w_avg_buffer, ctx.scaling, risk_sums, 1, 1, forward_scratch_buffer, dummy_map, true);
    forward_cumsum(forward_scratch_buffer, forward_cumsum2);

    forward_prework(ctx.status, w_avg_buffer, ctx.scaling, risk_sums, 1, 2, forward_scratch_buffer, dummy_map, true);
    forward_cumsum(forward_scratch_buffer, forward_cumsum3);

    forward_prework(ctx.status, w_avg_buffer, ctx.scaling, risk_sums, 2, 2, forward_scratch_buffer, dummy_map, true);
    forward_cumsum(forward_scratch_buffer, forward_cumsum4);

    for (int i = 0; i < ctx.last.size(); ++i) {
      T_1_term(i) = (forward_cumsum0(ctx.last(i) + 1) - 
		     (forward_cumsum2(ctx.last(i) + 1) - forward_cumsum2(ctx.first(i))));
      T_2_term(i) = ((forward_cumsum4(ctx.last(i) + 1) - forward_cumsum4(ctx.first(i))) 
		      - 2 * (forward_cumsum3(ctx.last(i) + 1) - forward_cumsum3(ctx.first(i))) + 
		      forward_cumsum1(ctx.last(i) + 1));
    }
    if (ctx.have_start_times) {
      for (int i = 0; i < ctx.start_map.size(); ++i) {
	T_1_term(i) -= forward_cumsum0(ctx.start_map(i));
      }
      for (int i = 0; i < ctx.first.size(); ++i) {      
	T_2_term(i) -= forward_cumsum1(ctx.start_map(i));
      }
    }
  }
  
  diag_part_buffer = exp_eta_w_event.array() * T_1_term.array();
  grad_buffer = w_event.array() * ctx.status.cast<double>().array() - diag_part_buffer.array();
  grad_buffer.array() *= -2.0;
  
  diag_hessian_buffer = exp_eta_w_event.array().pow(2) * T_2_term.array() - diag_part_buffer.array();
  diag_hessian_buffer.array() *= -2.0;
  
  to_native_from_event(grad_buffer, ctx.event_order, forward_scratch_buffer);
  to_native_from_event(diag_hessian_buffer, ctx.event_order, forward_scratch_buffer);
  to_native_from_event(diag_part_buffer, ctx.event_order, forward_scratch_buffer);
  
  double deviance = 2.0 * (loglik_sat - loglik);
  return(deviance);
}

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
                    Eigen::Ref<Eigen::VectorXd> hess_matvec_buffer)
{
  Eigen::VectorXd exp_w_times_arg = exp_w.array() * arg.array();

  sum_over_risk_set(MAKE_MAP_Xd(exp_w_times_arg),
                    ctx,
                    risk_sum_arg,
                    event_cumsum,
                    start_cumsum);

  for (int i = 0; i < ctx.status.size(); ++i) {
    if (w_avg(i) > 0) {
      forward_scratch_buffer(i) = (ctx.status(i) * w_avg(i) * risk_sum_arg(i)) / (risk_sums(i) * risk_sums(i));
    } else {
      forward_scratch_buffer(i) = 0.0;
    }
  }

  sum_over_events(ctx,
                  forward_cumsum0,
                  forward_cumsum1,
                  forward_scratch_buffer,
                  hess_matvec_buffer);
  
  to_native_from_event(hess_matvec_buffer, ctx.event_order, forward_scratch_buffer);
  hess_matvec_buffer = hess_matvec_buffer.array() * exp_w.array() - (diag_part.array() * arg.array());
}

/**
 * Among observations with tied times, the sort order is:
 * 1. Events (status == 1) occur first.
 * 2. Within events, those with positive weights occur first.
 * 
 * This ensures every other observation is a candidate for being in the risk set. 
 * Those with weight == 0 are ultimately not included in the risk set.
 */
std::vector<int> lexsort(const Eigen::VectorXi & a,  // is_start
                         const Eigen::VectorXd & b,  // weight
			 const Eigen::VectorXi & c,  // status
			 const Eigen::VectorXd & d)  // event
{
  std::vector<int> idx(a.size());
  std::iota(idx.begin(), idx.end(), 0); 
  
  auto comparator = [&](int i, int j) {
    if (d[i] != d[j]) return d[i] < d[j];
    if (c[i] != c[j]) return c[i] < c[j];
    if (b[i] != b[j]) return b[i] < b[j];
    if (a[i] != a[j]) return a[i] < a[j];
    // a total order, so that the order of ties (and with it the draws of
    // sample_score) does not depend on the standard library's std::sort
    return i < j;
  };
  
  std::sort(idx.begin(), idx.end(), comparator);
  
  return idx;
}

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
                      Eigen::VectorXi& start_order)
{
  int nevent = status.size();
  Eigen::VectorXi ones = Eigen::VectorXi::Ones(nevent);
  Eigen::VectorXi zeros = Eigen::VectorXi::Zero(nevent);

  Eigen::VectorXd stacked_time(nevent + nevent);
  stacked_time.segment(0, nevent) = start;
  stacked_time.segment(nevent, nevent) = event;

  Eigen::VectorXi stacked_status_c(nevent + nevent);
  stacked_status_c.segment(0, nevent) = ones;
  stacked_status_c.segment(nevent, nevent) = ones - status; 

  Eigen::VectorXi stacked_is_start(nevent + nevent);
  stacked_is_start.segment(0, nevent) = ones;
  stacked_is_start.segment(nevent, nevent) = zeros;

  Eigen::VectorXd stacked_weight(nevent + nevent);
  stacked_weight.segment(0, nevent) = -weight;      
  stacked_weight.segment(nevent, nevent) = -weight;

  Eigen::VectorXi stacked_index(nevent + nevent);
  stacked_index.segment(0, nevent) = Eigen::VectorXi::LinSpaced(nevent, 0, nevent - 1);
  stacked_index.segment(nevent, nevent) =  Eigen::VectorXi::LinSpaced(nevent, 0, nevent - 1);

  std::vector<int> sort_order = lexsort(stacked_is_start, stacked_weight, stacked_status_c, stacked_time);
  Eigen::VectorXi argsort = Eigen::Map<const Eigen::VectorXi>(sort_order.data(), sort_order.size());

  Eigen::VectorXd sorted_time(stacked_time.size()), sorted_status(stacked_status_c.size()),
    sorted_is_start(stacked_is_start.size()), sorted_index(stacked_index.size()), sorted_weight(stacked_weight.size());
  for (int i = 0; i < sorted_time.size(); ++i) {
    int j = argsort(i);
    sorted_time(i) = stacked_time(j);
    sorted_status(i) = 1 - stacked_status_c(j);
    sorted_is_start(i) = stacked_is_start(j);
    sorted_index(i) = stacked_index(j);    
    sorted_weight(i) = -stacked_weight(j);     
  }

  int event_count = 0, start_count = 0;
  std::vector<int> event_order_vec, start_order_vec, start_map_vec, event_map_vec, first_vec;
  int first_event = -1, num_successive_event = 1;
  double last_row_time;
  bool last_row_time_set = false;

  for (int i = 0; i < sorted_time.size(); ++i) {
    double _time = sorted_time(i); 
    int _status_val = (int)sorted_status(i);
    int _is_start = (int)sorted_is_start(i);
    int _index = (int)sorted_index(i);
    double _weight = sorted_weight(i);
    if (_is_start == 1) { 
      start_order_vec.push_back(_index);
      start_map_vec.push_back(event_count);
      start_count++;
    } else { 
      if (_status_val == 1)
	{
	  if ((last_row_time_set  && _time > last_row_time) || (_weight == 0)) {
	  first_event += num_successive_event;
	  num_successive_event = 1;
	} else {
	  num_successive_event++;
	}
	first_vec.push_back(first_event); 
      } 
      else {
	first_event += num_successive_event;
	num_successive_event = 1;
	first_vec.push_back(first_event); 
      }
      event_map_vec.push_back(start_count);
      event_order_vec.push_back(_index);
      event_count++;
    }
    last_row_time = _time;
    last_row_time_set = true;
  }

  _first = Eigen::Map<Eigen::VectorXi>(first_vec.data(), first_vec.size());
  start_order = Eigen::Map<Eigen::VectorXi>(start_order_vec.data(), start_order_vec.size());
  event_order = Eigen::Map<Eigen::VectorXi>(event_order_vec.data(), event_order_vec.size());
  Eigen::VectorXi start_map_tmp = Eigen::Map<Eigen::VectorXi>(start_map_vec.data(), start_map_vec.size());
  _event_map = Eigen::Map<Eigen::VectorXi>(event_map_vec.data(), event_map_vec.size());

  _start_map.resize(nevent);
  for (int i = 0; i < nevent; ++i) {
    _start_map(start_order(i)) = start_map_tmp(i);
  }

  _status.resize(nevent);
  for (int i = 0; i < nevent; ++i) {
    _status(i) = status(event_order(i));
  }
  
  Eigen::VectorXi _start_map_cp = _start_map;
  _start_map.resize(nevent);
  for (int i = 0; i < nevent; ++i) {
    _start_map(i) = _start_map_cp(event_order(i));
  }

  _event.resize(nevent);
  for (int i = 0; i < nevent; ++i) {
    _event(i) = event(event_order(i));
  }

  _start.resize(nevent);
  for (int i = 0; i < nevent; ++i) {
    _start(i) = event(start_order(i));
  }

  std::vector<int> last_vec;
  int last_event = nevent - 1, first_size = _first.size();
  for (int i = 0; i < first_size; ++i) {
    int f = _first(first_size - i - 1);
    last_vec.push_back(last_event);
    if (f - (nevent - 1 - i) == 0) {
      last_event = f - 1;
    }
  }

  _last.resize(nevent);
  for (int i = 0; i < nevent; ++i) {
    _last(i) = last_vec[nevent - i - 1];
  }

  _scaling.resize(nevent);
  for (int i = 0; i < nevent; ) {
    int f = _first(i);
    int l = _last(i);
    for (int j = f; j <= l; ++j) {
      _scaling(j) = (j - f) / (l + 1. - f);
    }
    i = l + 1;
  }

  bool check_ok = true;
  for (int i = 0; (i < _first.size()) && (check_ok); ++i) {
    check_ok = (_first[_start_map[i]] == _start_map[i]);
  }
  if (!check_ok) {
    ERROR_MSG("first_start disagrees with start_map");
  }
}

void CoxDeviance::setup_buffers(int n) {
    T_1_term.setZero(n);
    T_2_term.setZero(n);
    forward_scratch_buffer.setZero(n);
    w_avg_buffer.setZero(n);
    exp_w_buffer.setZero(n);
    grad_buffer.setZero(n);
    diag_hessian_buffer.setZero(n);
    diag_part_buffer.setZero(n);
    hess_matvec_buffer.setZero(n);
    
    event_reorder_buffers.clear();
    for (int i = 0; i < 3; ++i) event_reorder_buffers.push_back(Eigen::VectorXd::Zero(n));
    
    forward_cumsum_buffers.clear();
    for (int i = 0; i < 5; ++i) forward_cumsum_buffers.push_back(Eigen::VectorXd::Zero(n + 1));
    
    reverse_cumsum_buffers.clear();
    for (int i = 0; i < 4; ++i) reverse_cumsum_buffers.push_back(Eigen::VectorXd::Zero(n + 1));
    
    risk_sum_buffers.clear();
    for (int i = 0; i < 2; ++i) risk_sum_buffers.push_back(Eigen::VectorXd::Zero(n));
}

CoxDeviance::CoxDeviance(const Eigen::VectorXd& start,
                         const Eigen::VectorXd& event,
                         const Eigen::VectorXi& status,
                         const Eigen::VectorXd& weight,
                         bool efron) : _efron(efron) {
    setup_preprocess(start, event, status, weight, _start, _event, _status, _first, _last, _scaling, _start_map, _event_map, event_order, start_order);

    int n = _status.size();
    setup_buffers(n);
    
    have_start_times = (start.array() > -std::numeric_limits<double>::infinity()).any();
    sample_weight = weight;
    raw_start = start;
    raw_event = event;
    raw_status = status;
}

double CoxDeviance::compute_deviance(const Eigen::VectorXd& eta,
                                     const Eigen::VectorXd& sw) {
    linear_predictor = eta;
    sample_weight = sw;
    
    Eigen::VectorXd eta_centered = eta.array() - eta.mean();
    exp_w_buffer = sw.array() * (eta_centered.array().min(30.0)).exp();

    loglik_sat = compute_sat_loglik(_first, _last, sample_weight, event_order, _status, _scaling, forward_cumsum_buffers[0], _efron);

    return cox_dev(eta_centered,
                   sample_weight,
                   exp_w_buffer,
                   get_context(),
                   loglik_sat,
                   T_1_term,
                   T_2_term,
                   grad_buffer,
                   diag_hessian_buffer,
                   diag_part_buffer,
                   w_avg_buffer,
                   event_reorder_buffers[0],
                   event_reorder_buffers[1],
                   event_reorder_buffers[2],
                   risk_sum_buffers[0],
                   forward_cumsum_buffers[0],
                   forward_cumsum_buffers[1],
                   forward_cumsum_buffers[2],
                   forward_cumsum_buffers[3],
                   forward_cumsum_buffers[4],
                   forward_scratch_buffer,
                   reverse_cumsum_buffers[0],
                   reverse_cumsum_buffers[1]);
}

void CoxDeviance::compute_hessian_matvec(const Eigen::VectorXd& arg,
                                         Eigen::VectorXd& out) {
    Eigen::VectorXd neg_arg = -arg;

    hessian_matvec(neg_arg,
                   get_context(),
                   risk_sum_buffers[0],
                   diag_part_buffer,
                   w_avg_buffer,
                   exp_w_buffer,
                   risk_sum_buffers[1],
                   forward_cumsum_buffers[0],
                   forward_cumsum_buffers[1],
                   forward_scratch_buffer,
                   reverse_cumsum_buffers[2],
                   reverse_cumsum_buffers[3],
                   hess_matvec_buffer);
    out = hess_matvec_buffer;
}

StratifiedCoxDeviance::StratifiedCoxDeviance(const Eigen::VectorXd& start,
                                           const Eigen::VectorXd& event,
                                           const Eigen::VectorXi& status,
                                           const Eigen::VectorXi& strata,
                                           const Eigen::VectorXd& weight,
                                           bool efron) {
    int n = status.size();
    grad_buffer.setZero(n);
    diag_hessian_buffer.setZero(n);
    
    std::vector<int> strata_std(strata.data(), strata.data() + strata.size());
    unique_strata = strata_std;
    std::sort(unique_strata.begin(), unique_strata.end());
    auto last = std::unique(unique_strata.begin(), unique_strata.end());
    unique_strata.erase(last, unique_strata.end());

    stratum_indices.resize(unique_strata.size());
    for (int i = 0; i < n; ++i) {
        int s = strata(i);
        auto it = std::lower_bound(unique_strata.begin(), unique_strata.end(), s);
        int idx = std::distance(unique_strata.begin(), it);
        stratum_indices[idx].push_back(i);
    }

    for (size_t i = 0; i < unique_strata.size(); ++i) {
        const auto& indices = stratum_indices[i];
        int n_stratum = indices.size();
        Eigen::VectorXd s_start(n_stratum);
        Eigen::VectorXd s_event(n_stratum);
        Eigen::VectorXi s_status(n_stratum);
        Eigen::VectorXd s_weight(n_stratum);

        for (int j = 0; j < n_stratum; ++j) {
            int orig_idx = indices[j];
            s_start(j) = start(orig_idx);
            s_event(j) = event(orig_idx);
            s_status(j) = status(orig_idx);
            s_weight(j) = weight(orig_idx);
        }

        cox_devs.push_back(std::make_shared<CoxDeviance>(s_start, s_event, s_status, s_weight, efron));
    }
}

double StratifiedCoxDeviance::compute_deviance(const Eigen::VectorXd& eta,
                                               const Eigen::VectorXd& sw) {
    linear_predictor = eta;
    sample_weight = sw;
    double total_deviance = 0.0;
    double total_loglik_sat = 0.0;
    
    grad_buffer.setZero(eta.size());
    diag_hessian_buffer.setZero(eta.size());

    for (size_t i = 0; i < unique_strata.size(); ++i) {
        const auto& indices = stratum_indices[i];
        int n_stratum = indices.size();
        Eigen::VectorXd s_eta(n_stratum);
        Eigen::VectorXd s_sw(n_stratum);

        for (int j = 0; j < n_stratum; ++j) {
            int orig_idx = indices[j];
            s_eta(j) = eta(orig_idx);
            s_sw(j) = sw(orig_idx);
        }

        double dev = cox_devs[i]->compute_deviance(s_eta, s_sw);
        total_deviance += dev;
        total_loglik_sat += cox_devs[i]->get_loglik_sat();

        const auto& s_grad = cox_devs[i]->get_gradient();
        const auto& s_diag_hessian = cox_devs[i]->get_diag_hessian();

        for (int j = 0; j < n_stratum; ++j) {
            int orig_idx = indices[j];
            grad_buffer(orig_idx) = s_grad(j);
            diag_hessian_buffer(orig_idx) = s_diag_hessian(j);
        }
    }

    loglik_sat = total_loglik_sat;
    return total_deviance;
}

void StratifiedCoxDeviance::compute_hessian_matvec(const Eigen::VectorXd& arg,
                                                   Eigen::VectorXd& out) {
    out.setZero(arg.size());

    for (size_t i = 0; i < unique_strata.size(); ++i) {
        const auto& indices = stratum_indices[i];
        int n_stratum = indices.size();
        Eigen::VectorXd s_arg(n_stratum);
        Eigen::VectorXd s_out(n_stratum);

        for (int j = 0; j < n_stratum; ++j) {
            int orig_idx = indices[j];
            s_arg(j) = arg(orig_idx);
        }

        cox_devs[i]->compute_hessian_matvec(s_arg, s_out);

        for (int j = 0; j < n_stratum; ++j) {
            int orig_idx = indices[j];
            out(orig_idx) = s_out(j);
        }
    }
}


// ---------------------------------------------------------------------------
// Sampling the score

namespace {

// Random numbers from the bits of std::mt19937_64, whose output the C++
// standard specifies exactly. The distributions of <random> differ between
// standard libraries, so they are not used: the same seed gives the same
// draws on every platform, from Python or R.
class PortableRNG {
public:
  explicit PortableRNG(std::uint64_t seed) : engine(seed) {}

  double uniform() {            // in (0, 1)
    return ((engine() >> 11) + 0.5) * (1.0 / 9007199254740992.0);
  }

  double normal() {             // Box-Muller
    if (have_spare) {
      have_spare = false;
      return spare;
    }
    double r = std::sqrt(-2.0 * std::log(uniform()));
    double theta = 2.0 * 3.14159265358979323846 * uniform();
    spare = r * std::sin(theta);
    have_spare = true;
    return r * std::cos(theta);
  }

private:
  std::mt19937_64 engine;
  bool have_spare = false;
  double spare = 0.0;
};

// splitmix64: the seeds of strata, derived from one seed
std::uint64_t mix_seed(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// A Fenwick tree of nonnegative weights: draws an index >= first with
// probability proportional to its weight, and changes a weight, in O(log n).
class WeightTree {
public:
  explicit WeightTree(const Eigen::VectorXd& w)
    : n(w.size()), weight(w), tree(w.size() + 1, 0.0), top(1) {
    for (int i = 1; i <= n; ++i) {
      tree[i] += w(i - 1);
      int parent = i + (i & -i);
      if (parent <= n) tree[parent] += tree[i];
    }
    while (2 * top <= n) top *= 2;
  }

  double get(int i) const { return weight(i); }

  void set(int i, double value) {
    double delta = value - weight(i);
    weight(i) = value;
    for (int j = i + 1; j <= n; j += j & -j) tree[j] += delta;
  }

  // sum of the weights of the indices < i
  double prefix(int i) const {
    double sum = 0.0;
    for (int j = i; j > 0; j -= j & -j) sum += tree[j];
    return sum;
  }

  // an index >= first, with probability proportional to its weight; -1 if
  // their total is not positive. Rounding can return an index of weight 0,
  // which the caller rejects.
  int draw(int first, PortableRNG& rng) const {
    double base = prefix(first);
    double total = prefix(n) - base;
    if (!(total > 0)) return -1;
    double target = base + rng.uniform() * total;
    int pos = 0;
    for (int step = top; step > 0 && n > 0; step >>= 1) {
      if (pos + step <= n && tree[pos + step] <= target) {
        pos += step;
        target -= tree[pos];
      }
    }
    return std::min(std::max(pos, first), n - 1);
  }

private:
  int n;
  Eigen::VectorXd weight;
  std::vector<double> tree;
  int top;
};

// A distinct failure time: the failures with positive weight in event-order
// positions first, ..., last, d of them, at time t.
struct FailureTime {
  int first;
  int last;
  int d;
  double time;
};

std::vector<FailureTime> failure_times(const CoxContext& ctx,
                                       const Eigen::VectorXd& w_event,
                                       const Eigen::VectorXd& event_time) {
  std::vector<FailureTime> times;
  int n = ctx.status.size();
  for (int i = 0; i < n; ) {
    int l = std::max(ctx.last(i), i);
    if (ctx.status(i) == 1 && ctx.first(i) == i && w_event(i) > 0) {
      int d = 0;
      for (int j = i; j <= l; ++j) {
        if (ctx.status(j) == 1 && w_event(j) > 0) ++d;
      }
      if (d > 0) times.push_back({i, l, d, event_time(i)});
    }
    i = l + 1;
  }
  return times;
}

// For values c at the failure positions t (event order), the sum over t whose
// risk set contains j of c_t, for each position j; with efron, the tied
// failures at t count with weight 1 - scaling(t), as in the Efron risk sums.
// This is the computation of T_1 in cox_dev, with c in place of w_avg / risk_sums.
Eigen::VectorXd risk_set_accumulate(const CoxContext& ctx,
                                    const Eigen::VectorXd& c,
                                    bool efron) {
  int n = c.size();
  Eigen::VectorXd C(n + 1), out(n);
  forward_cumsum(c, C);
  for (int j = 0; j < n; ++j) {
    out(j) = C(ctx.last(j) + 1) - (ctx.have_start_times ? C(ctx.start_map(j)) : 0.0);
  }
  if (efron) {
    Eigen::VectorXd c_scaled = c.array() * ctx.scaling.array();
    forward_cumsum(c_scaled, C);
    for (int j = 0; j < n; ++j) {
      out(j) -= C(ctx.last(j) + 1) - C(ctx.first(j));
    }
  }
  return out;
}

// the common positive weight (the model-based draws need equal weights)
double common_weight(const Eigen::VectorXd& w) {
  double w0 = 0.0;
  for (int i = 0; i < w.size(); ++i) {
    if (w(i) > 0) {
      if (w0 == 0.0) {
        w0 = w(i);
      } else if (std::abs(w(i) - w0) > 1e-12 * w0) {
        ERROR_MSG("sample_score: the model-based methods need equal positive weights; use method 'wild'");
      }
    }
  }
  return w0;
}

}  // namespace

// wild_score and score_residuals reuse what compute_deviance leaves in its
// buffers, which depend only on the linear predictor and weights:
// risk_sum_buffers[0] (S, the risk sums of w exp(eta), with Efron's weights
// if efron), w_avg_buffer (wbar, the mean weight of the failures tied at each
// time), event_reorder_buffers[1] and [2] (w and w exp(eta) in event order)
// and grad_buffer. A call then costs one risk-set sum of its argument.

// the failure terms t: a failure with positive weight and risk sum
inline bool failure_term(const CoxContext& ctx, const Eigen::VectorXd& w_ev,
                         const Eigen::VectorXd& S, int t) {
  return ctx.status(t) == 1 && w_ev(t) > 0 && S(t) > 0;
}

Eigen::VectorXd CoxDeviance::wild_score(const Eigen::VectorXd& multipliers) {
  int n = _status.size();
  if (linear_predictor.size() != n) {
    ERROR_MSG("wild_score: call compute_deviance first");
  }
  if (multipliers.size() != n) {
    ERROR_MSG("wild_score: multipliers must have one entry per observation");
  }
  CoxContext ctx = get_context();
  const Eigen::VectorXd& S = risk_sum_buffers[0];
  const Eigen::VectorXd& w_ev = event_reorder_buffers[1];
  const Eigen::VectorXd& e_ev = event_reorder_buffers[2];

  // the risk sums of G w exp(eta), with Efron's weights if efron
  Eigen::VectorXd SG(n), event_cumsum(n + 1), start_cumsum(n + 1);
  Eigen::VectorXd Ge = multipliers.array() * exp_w_buffer.array();
  sum_over_risk_set(Ge, ctx, SG, event_cumsum, start_cumsum);

  // in event order: the multipliers, and the mean of G w over the failures
  Eigen::VectorXd G_ev(n), status_Gw(n), CGw(n + 1);
  for (int i = 0; i < n; ++i) {
    G_ev(i) = multipliers(event_order(i));
    status_Gw(i) = ctx.status(i) * w_ev(i) * G_ev(i);
  }
  forward_cumsum(status_Gw, CGw);

  // For the failure term t at a time with failures D and average weight wbar:
  //   c_t = sum_{i in D} G_i w_i / d - wbar * SG_t / S_t,
  // and the draw is r_j = G_j M_j - e_j sum_{t containing j} c_t / S_t.
  Eigen::VectorXd c = Eigen::VectorXd::Zero(n);
  for (int t = 0; t < n; ++t) {
    if (failure_term(ctx, w_ev, S, t)) {
      int f = ctx.first(t), l = ctx.last(t);
      double Gwbar = (CGw(l + 1) - CGw(f)) / (l + 1 - f);
      c(t) = (Gwbar - w_avg_buffer(t) * SG(t) / S(t)) / S(t);
    }
  }
  Eigen::VectorXd T = risk_set_accumulate(ctx, c, ctx.efron);

  Eigen::VectorXd out(n);
  for (int i = 0; i < n; ++i) {
    int k = event_order(i);
    out(k) = G_ev(i) * (-0.5 * grad_buffer(k)) - e_ev(i) * T(i);
  }
  return out;
}

Eigen::VectorXd CoxDeviance::score_residuals(const Eigen::VectorXd& x) const {
  int n = _status.size();
  if (linear_predictor.size() != n) {
    ERROR_MSG("score_residuals: call compute_deviance first");
  }
  if (x.size() != n) {
    ERROR_MSG("score_residuals: x must have one entry per observation");
  }
  CoxContext ctx = get_context();
  const Eigen::VectorXd& S = risk_sum_buffers[0];
  const Eigen::VectorXd& w_ev = event_reorder_buffers[1];
  const Eigen::VectorXd& e_ev = event_reorder_buffers[2];

  // the risk sums of w exp(eta) x, with Efron's weights if efron
  Eigen::VectorXd SX(n), event_cumsum(n + 1), start_cumsum(n + 1);
  Eigen::VectorXd ex = exp_w_buffer.array() * x.array();
  sum_over_risk_set(ex, ctx, SX, event_cumsum, start_cumsum);

  // For the failure terms t, with their average weight wbar_t and xbar_t the
  // mean of x over the risk set of t (with weights w exp(eta), and Efron's
  // for ties):
  //   u_i = x_i M_i - (w_i delta_i / d) sum_{t at i's time} xbar_t
  //         + e_i sum_{t containing i} wbar_t xbar_t / S_t,
  // with M = -gradient / 2.
  Eigen::VectorXd xbar = Eigen::VectorXd::Zero(n), c = Eigen::VectorXd::Zero(n);
  for (int t = 0; t < n; ++t) {
    if (failure_term(ctx, w_ev, S, t)) {
      xbar(t) = SX(t) / S(t);
      c(t) = w_avg_buffer(t) * xbar(t) / S(t);
    }
  }
  Eigen::VectorXd T = risk_set_accumulate(ctx, c, ctx.efron);
  Eigen::VectorXd Cx(n + 1);
  forward_cumsum(xbar, Cx);

  Eigen::VectorXd out(n);
  for (int i = 0; i < n; ++i) {
    int k = event_order(i);
    double own = 0.0;
    if (ctx.status(i) == 1 && w_ev(i) > 0) {
      int f = ctx.first(i), l = ctx.last(i);
      own = w_ev(i) * (Cx(l + 1) - Cx(f)) / (l + 1 - f);
    }
    out(k) = x(k) * (-0.5 * grad_buffer(k)) - own + e_ev(i) * T(i);
  }
  return out;
}

BaselineHazard CoxDeviance::baseline_hazard() const {
  int n = _status.size();
  if (linear_predictor.size() != n) {
    ERROR_MSG("baseline_hazard: call compute_deviance first");
  }
  CoxContext ctx = get_context();

  Eigen::VectorXd S(n), event_cumsum(n + 1), start_cumsum(n + 1);
  sum_over_risk_set(exp_w_buffer, ctx, S, event_cumsum, start_cumsum);
  Eigen::VectorXd w_ev(n);
  for (int i = 0; i < n; ++i) w_ev(i) = sample_weight(event_order(i));
  Eigen::VectorXd status_w = ctx.status.cast<double>().array() * w_ev.array();
  Eigen::VectorXd Cw(n + 1);
  forward_cumsum(status_w, Cw);

  // the increment at a failure time is the sum over its failure terms t of
  // wbar_t / S_t: d wbar / S for Breslow's ties, Efron's estimate for his.
  // exp_w is w exp(eta - mean(eta)), so the increments for eta as given are
  // these times exp(-mean(eta)).
  double scale = std::exp(-linear_predictor.mean());
  std::vector<double> time, hazard;
  int current = -1;
  for (int t = 0; t < n; ++t) {
    if (ctx.status(t) == 1 && w_ev(t) > 0 && S(t) > 0) {
      int f = ctx.first(t), l = ctx.last(t);
      double wbar = (Cw(l + 1) - Cw(f)) / (l + 1 - f);
      if (f != current) {
        current = f;
        time.push_back(_event(t));
        hazard.push_back(0.0);
      }
      hazard.back() += scale * wbar / S(t);
    }
  }

  BaselineHazard out;
  int m = time.size();
  out.time = Eigen::Map<Eigen::VectorXd>(time.data(), m);
  out.hazard = Eigen::Map<Eigen::VectorXd>(hazard.data(), m);
  out.cumulative_hazard.resize(m);
  double total = 0.0;
  for (int k = 0; k < m; ++k) {
    total += hazard[k];
    out.cumulative_hazard(k) = total;
  }
  out.stratum = Eigen::VectorXi::Zero(m);
  return out;
}

Eigen::VectorXd CoxDeviance::schoenfeld_residuals(const Eigen::VectorXd& x) const {
  int n = _status.size();
  if (linear_predictor.size() != n) {
    ERROR_MSG("schoenfeld_residuals: call compute_deviance first");
  }
  if (x.size() != n) {
    ERROR_MSG("schoenfeld_residuals: x must have one entry per observation");
  }
  CoxContext ctx = get_context();
  const Eigen::VectorXd& S = risk_sum_buffers[0];
  const Eigen::VectorXd& w_ev = event_reorder_buffers[1];

  // xbar_t, the mean of x over the risk set of each failure term t (with
  // Efron's weights if efron), as in score_residuals
  Eigen::VectorXd SX(n), event_cumsum(n + 1), start_cumsum(n + 1);
  Eigen::VectorXd ex = exp_w_buffer.array() * x.array();
  sum_over_risk_set(ex, ctx, SX, event_cumsum, start_cumsum);
  Eigen::VectorXd xbar = Eigen::VectorXd::Zero(n);
  for (int t = 0; t < n; ++t) {
    if (failure_term(ctx, w_ev, S, t)) xbar(t) = SX(t) / S(t);
  }
  Eigen::VectorXd Cx(n + 1);
  forward_cumsum(xbar, Cx);

  // each failure less the average of xbar over the failures tied with it
  Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; ++i) {
    if (failure_term(ctx, w_ev, S, i)) {
      int f = ctx.first(i), l = ctx.last(i);
      int k = event_order(i);
      out(k) = x(k) - (Cx(l + 1) - Cx(f)) / (l + 1 - f);
    }
  }
  return out;
}

Eigen::VectorXd CoxDeviance::risk_set_variance(const Eigen::VectorXd& z) const {
  int n = _status.size();
  if (linear_predictor.size() != n) {
    ERROR_MSG("risk_set_variance: call compute_deviance first");
  }
  if (z.size() != n) {
    ERROR_MSG("risk_set_variance: z must have one entry per observation");
  }
  CoxContext ctx = get_context();
  const Eigen::VectorXd& S = risk_sum_buffers[0];
  const Eigen::VectorXd& w_ev = event_reorder_buffers[1];

  // the variance is unchanged by a shift of z: centering it first keeps the
  // mean of the squares less the square of the mean accurate
  Eigen::VectorXd zc = z.array() - z.mean();
  Eigen::VectorXd SZ(n), SZZ(n), event_cumsum(n + 1), start_cumsum(n + 1);
  Eigen::VectorXd ez = exp_w_buffer.array() * zc.array();
  sum_over_risk_set(ez, ctx, SZ, event_cumsum, start_cumsum);
  Eigen::VectorXd ezz = ez.array() * zc.array();
  sum_over_risk_set(ezz, ctx, SZZ, event_cumsum, start_cumsum);

  // the variance of each failure term t (with Efron's weights if efron),
  // then its average over the failures tied with i
  Eigen::VectorXd var = Eigen::VectorXd::Zero(n);
  for (int t = 0; t < n; ++t) {
    if (failure_term(ctx, w_ev, S, t)) {
      double mean = SZ(t) / S(t);
      var(t) = std::max(SZZ(t) / S(t) - mean * mean, 0.0);
    }
  }
  Eigen::VectorXd Cv(n + 1);
  forward_cumsum(var, Cv);

  Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; ++i) {
    if (failure_term(ctx, w_ev, S, i)) {
      int f = ctx.first(i), l = ctx.last(i);
      out(event_order(i)) = (Cv(l + 1) - Cv(f)) / (l + 1 - f);
    }
  }
  return out;
}

Eigen::VectorXd CoxDeviance::expected_events() const {
  int n = _status.size();
  if (linear_predictor.size() != n) {
    ERROR_MSG("expected_events: call compute_deviance first");
  }
  CoxContext ctx = get_context();
  const Eigen::VectorXd& S = risk_sum_buffers[0];
  const Eigen::VectorXd& w_ev = event_reorder_buffers[1];

  // the hazard increment of each failure term t is wbar_t / S_t (see
  // baseline_hazard); summed over the terms whose risk set contains i, with
  // Efron's weights for i's own failure time
  Eigen::VectorXd c = Eigen::VectorXd::Zero(n);
  for (int t = 0; t < n; ++t) {
    if (failure_term(ctx, w_ev, S, t)) c(t) = w_avg_buffer(t) / S(t);
  }
  Eigen::VectorXd T = risk_set_accumulate(ctx, c, ctx.efron);

  // S is a risk sum of w exp(eta - mean(eta)), so exp(eta - mean(eta)) T is
  // the expected number of events for eta as given
  double center = linear_predictor.mean();
  Eigen::VectorXd out(n);
  for (int i = 0; i < n; ++i) {
    int k = event_order(i);
    out(k) = std::exp(linear_predictor(k) - center) * T(i);
  }
  return out;
}

Eigen::MatrixXd CoxDeviance::sample_score(const std::string& method,
                                            int n_draws,
                                            std::uint64_t seed) {
  int n = _status.size();
  if (linear_predictor.size() != n) {
    ERROR_MSG("sample_score: call compute_deviance first");
  }
  if (n_draws < 0) {
    ERROR_MSG("sample_score: n_draws must be nonnegative");
  }
  PortableRNG rng(seed);
  Eigen::MatrixXd out(n, n_draws);

  if (method == "wild") {
    Eigen::VectorXd G(n);
    for (int b = 0; b < n_draws; ++b) {
      for (int i = 0; i < n; ++i) G(i) = rng.normal();
      out.col(b) = wild_score(G);
    }
    return out;
  }

  bool multinomial = (method == "multinomial");
  bool bernoulli = (method == "bernoulli");
  bool fisher = (method == "fisher");
  if (!(multinomial || bernoulli || fisher)) {
    ERROR_MSG("sample_score: method must be one of 'multinomial', 'bernoulli', "
              "'fisher' or 'wild'");
  }

  CoxContext ctx = get_context();
  double w0 = common_weight(sample_weight);

  // in event order
  Eigen::VectorXd e_ev(n), w_ev(n), start_ev(n);
  for (int i = 0; i < n; ++i) {
    int k = event_order(i);
    e_ev(i) = exp_w_buffer(k);
    w_ev(i) = sample_weight(k);
    start_ev(i) = raw_start(k);
  }
  std::vector<FailureTime> times = failure_times(ctx, w_ev, _event);

  // the risk sums without Efron's weights: the draws are from the whole risk set
  CoxContext breslow = ctx;
  breslow.efron = false;
  Eigen::VectorXd S(n), event_cumsum(n + 1), start_cumsum(n + 1);
  sum_over_risk_set(exp_w_buffer, breslow, S, event_cumsum, start_cumsum);

  auto at_risk = [&](int j, double t) {
    return e_ev(j) > 0 && (!have_start_times || start_ev(j) < t);
  };

  WeightTree tree(e_ev);

  // a subject at risk at ft, with probability proportional to its weight in tree
  auto draw_at_risk = [&](const FailureTime& ft) -> int {
    // the tree holds every subject with event time >= t: reject those who
    // entered at t or later
    for (int attempt = 0; attempt < 64; ++attempt) {
      int j = tree.draw(ft.first, rng);
      if (j < 0) break;
      if (tree.get(j) > 0 && at_risk(j, ft.time)) return j;
    }
    // many late entries: draw directly among those at risk
    double total = 0.0;
    for (int j = ft.first; j < n; ++j) {
      if (tree.get(j) > 0 && at_risk(j, ft.time)) total += tree.get(j);
    }
    if (!(total > 0)) {
      ERROR_MSG("sample_score: a failure time with an empty risk set");
    }
    double target = rng.uniform() * total;
    int last_eligible = -1;
    for (int j = ft.first; j < n; ++j) {
      if (tree.get(j) > 0 && at_risk(j, ft.time)) {
        last_eligible = j;
        target -= tree.get(j);
        if (target < 0) return j;
      }
    }
    return last_eligible;
  };

  // the failure positions of ft (they carry the per-term values for risk_set_accumulate)
  auto spread = [&](const FailureTime& ft, double value, Eigen::VectorXd& c) {
    for (int t = ft.first; t <= ft.last; ++t) {
      if (ctx.status(t) == 1 && w_ev(t) > 0) c(t) = value;
    }
  };

  if (fisher) {
    // At each failure time the set of d failures is drawn from Fisher's
    // noncentral hypergeometric distribution, P(A) proportional to the product
    // of e_j over A (|A| = d), the conditional law of independent Bernoullis
    // given d failures; r is its score under the exact (discrete) partial
    // likelihood, r_j = w0 sum_k (delta_kj - E_k delta_j). Both use the
    // elementary symmetric polynomials E_k of the risk weights, computed on the
    // log scale by the recursion of Gail, Lubin and Rubinstein (1981).
    const double neg_inf = -std::numeric_limits<double>::infinity();
    auto log_add = [neg_inf](double a, double b) {
      if (a == neg_inf) return b;
      if (b == neg_inf) return a;
      return std::max(a, b) + std::log1p(std::exp(-std::abs(a - b)));
    };
    Eigen::MatrixXd r_ev = Eigen::MatrixXd::Zero(n, n_draws);
    std::vector<int> risk;
    std::vector<double> loge, suffix, prefix;
    for (const FailureTime& ft : times) {
      risk.clear();
      for (int j = ft.first; j < n; ++j) {
        if (at_risk(j, ft.time)) risk.push_back(j);
      }
      int m = risk.size();
      int d = ft.d;
      if (m < d) {
        ERROR_MSG("sample_score: more failures than subjects at risk");
      }
      loge.resize(m);
      for (int i = 0; i < m; ++i) loge[i] = std::log(e_ev(risk[i]));

      // suffix(i, k) = log E_k(risk[i..m-1]), prefix(i, k) = log E_k(risk[0..i-1])
      int D = d + 1;
      suffix.assign((m + 1) * D, neg_inf);
      prefix.assign((m + 1) * D, neg_inf);
      suffix[m * D] = 0.0;
      for (int i = m - 1; i >= 0; --i) {
        suffix[i * D] = 0.0;
        for (int k = 1; k <= d; ++k) {
          suffix[i * D + k] = log_add(suffix[(i + 1) * D + k],
                                      loge[i] + suffix[(i + 1) * D + k - 1]);
        }
      }
      prefix[0] = 0.0;
      for (int i = 1; i <= m; ++i) {
        prefix[i * D] = 0.0;
        for (int k = 1; k <= d; ++k) {
          prefix[i * D + k] = log_add(prefix[(i - 1) * D + k],
                                      loge[i - 1] + prefix[(i - 1) * D + k - 1]);
        }
      }
      double log_total = suffix[d];

      // the exact expectation: E delta_i = e_i E_{d-1}(risk without i) / E_d(risk)
      for (int i = 0; i < m; ++i) {
        double log_without = neg_inf;
        for (int a = 0; a <= d - 1; ++a) {
          log_without = log_add(log_without,
                                prefix[i * D + a] + suffix[(i + 1) * D + d - 1 - a]);
        }
        double p = std::exp(loge[i] + log_without - log_total);
        r_ev.row(risk[i]).array() -= w0 * p;
      }

      // sequential draws: include risk[i] with probability e_i E_{k-1}(after i) / E_k(from i)
      for (int b = 0; b < n_draws; ++b) {
        int k = d;
        for (int i = 0; i < m && k > 0; ++i) {
          double p = std::exp(loge[i] + suffix[(i + 1) * D + k - 1] - suffix[i * D + k]);
          if (rng.uniform() < p) {
            r_ev(risk[i], b) += w0;
            --k;
          }
        }
      }
    }
    for (int i = 0; i < n; ++i) out.row(event_order(i)) = r_ev.row(i);
    return out;
  }

  // each draw is the (Breslow) score of failures drawn from the model: at
  // each event time, sum over the failures of x_j minus the number of
  // failures times the risk-set mean, i.e. r_j = w0 (count_j - e_j sum_k D_k / S_k)
  // with D_k the number of failures drawn at time k
  Eigen::VectorXd r(n), c(n);
  for (int b = 0; b < n_draws; ++b) {
    r.setZero();
    c.setZero();
    if (bernoulli) {
      // each subject at risk fails independently with probability min(1, d e_j / S)
      for (const FailureTime& ft : times) {
        double S_k = S(ft.first);
        double D = 0.0;
        for (int j = ft.first; j < n; ++j) {
          if (!at_risk(j, ft.time)) continue;
          double pi = std::min(1.0, ft.d * e_ev(j) / S_k);
          if (rng.uniform() < pi) {
            r(j) += 1.0;
            D += 1.0;
          }
        }
        // spread puts the value on each of the d failure positions
        spread(ft, D / (ft.d * S_k), c);
      }
    } else {
      // d failures at each time, drawn with replacement
      for (const FailureTime& ft : times) {
        for (int m = 0; m < ft.d; ++m) r(draw_at_risk(ft)) += 1.0;
        spread(ft, 1.0 / S(ft.first), c);
      }
    }
    Eigen::VectorXd T = risk_set_accumulate(ctx, c, false);
    r.array() -= e_ev.array() * T.array();
    r *= w0;
    for (int i = 0; i < n; ++i) out(event_order(i), b) = r(i);
  }
  return out;
}

Eigen::MatrixXd StratifiedCoxDeviance::sample_score(const std::string& method,
                                                      int n_draws,
                                                      std::uint64_t seed) {
  int n = linear_predictor.size();
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, std::max(n_draws, 0));
  for (size_t s = 0; s < unique_strata.size(); ++s) {
    const auto& indices = stratum_indices[s];
    Eigen::MatrixXd draws = cox_devs[s]->sample_score(method, n_draws,
                                                        mix_seed(seed + mix_seed(s + 1)));
    for (size_t j = 0; j < indices.size(); ++j) {
      out.row(indices[j]) = draws.row(j);
    }
  }
  return out;
}

Eigen::VectorXd StratifiedCoxDeviance::wild_score(const Eigen::VectorXd& multipliers) {
  int n = linear_predictor.size();
  if (multipliers.size() != n) {
    ERROR_MSG("wild_score: multipliers must have one entry per observation");
  }
  Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
  for (size_t s = 0; s < unique_strata.size(); ++s) {
    const auto& indices = stratum_indices[s];
    Eigen::VectorXd G(indices.size());
    for (size_t j = 0; j < indices.size(); ++j) G(j) = multipliers(indices[j]);
    Eigen::VectorXd draw = cox_devs[s]->wild_score(G);
    for (size_t j = 0; j < indices.size(); ++j) out(indices[j]) = draw(j);
  }
  return out;
}

Eigen::VectorXd StratifiedCoxDeviance::score_residuals(const Eigen::VectorXd& x) const {
  int n = linear_predictor.size();
  if (x.size() != n) {
    ERROR_MSG("score_residuals: x must have one entry per observation");
  }
  Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
  for (size_t s = 0; s < unique_strata.size(); ++s) {
    const auto& indices = stratum_indices[s];
    Eigen::VectorXd xs(indices.size());
    for (size_t j = 0; j < indices.size(); ++j) xs(j) = x(indices[j]);
    Eigen::VectorXd u = cox_devs[s]->score_residuals(xs);
    for (size_t j = 0; j < indices.size(); ++j) out(indices[j]) = u(j);
  }
  return out;
}

Eigen::VectorXd StratifiedCoxDeviance::schoenfeld_residuals(const Eigen::VectorXd& x) const {
  int n = linear_predictor.size();
  if (x.size() != n) {
    ERROR_MSG("schoenfeld_residuals: x must have one entry per observation");
  }
  Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
  for (size_t s = 0; s < unique_strata.size(); ++s) {
    const auto& indices = stratum_indices[s];
    Eigen::VectorXd xs(indices.size());
    for (size_t j = 0; j < indices.size(); ++j) xs(j) = x(indices[j]);
    Eigen::VectorXd r = cox_devs[s]->schoenfeld_residuals(xs);
    for (size_t j = 0; j < indices.size(); ++j) out(indices[j]) = r(j);
  }
  return out;
}

Eigen::VectorXd StratifiedCoxDeviance::risk_set_variance(const Eigen::VectorXd& z) const {
  int n = linear_predictor.size();
  if (z.size() != n) {
    ERROR_MSG("risk_set_variance: z must have one entry per observation");
  }
  Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
  for (size_t s = 0; s < unique_strata.size(); ++s) {
    const auto& indices = stratum_indices[s];
    Eigen::VectorXd zs(indices.size());
    for (size_t j = 0; j < indices.size(); ++j) zs(j) = z(indices[j]);
    Eigen::VectorXd v = cox_devs[s]->risk_set_variance(zs);
    for (size_t j = 0; j < indices.size(); ++j) out(indices[j]) = v(j);
  }
  return out;
}

Eigen::VectorXd StratifiedCoxDeviance::expected_events() const {
  int n = linear_predictor.size();
  Eigen::VectorXd out = Eigen::VectorXd::Zero(n);
  for (size_t s = 0; s < unique_strata.size(); ++s) {
    const auto& indices = stratum_indices[s];
    Eigen::VectorXd e = cox_devs[s]->expected_events();
    for (size_t j = 0; j < indices.size(); ++j) out(indices[j]) = e(j);
  }
  return out;
}

BaselineHazard StratifiedCoxDeviance::baseline_hazard() const {
  // the strata in turn: each stratum's hazard is already for its linear
  // predictor as given
  std::vector<BaselineHazard> parts;
  int m = 0;
  for (size_t s = 0; s < unique_strata.size(); ++s) {
    parts.push_back(cox_devs[s]->baseline_hazard());
    m += parts.back().time.size();
  }
  BaselineHazard out;
  out.time.resize(m);
  out.hazard.resize(m);
  out.cumulative_hazard.resize(m);
  out.stratum.resize(m);
  int k = 0;
  for (size_t s = 0; s < parts.size(); ++s) {
    int ms = parts[s].time.size();
    out.time.segment(k, ms) = parts[s].time;
    out.hazard.segment(k, ms) = parts[s].hazard;
    out.cumulative_hazard.segment(k, ms) = parts[s].cumulative_hazard;
    out.stratum.segment(k, ms).setConstant(unique_strata[s]);
    k += ms;
  }
  return out;
}

// Harrell's concordance index; see coxdev.h for the definition and the
// references. Ported from glmstar (src/concordance.cpp), for one predictor.

namespace {

// prefix sums of weights indexed 0, ..., n - 1, under additions, in O(log n)
// (Fenwick, 1994)
class FenwickSum {
public:
  explicit FenwickSum(std::size_t n) : tree(n + 1, 0.0) {}
  void add(std::size_t i, double w) {
    for (++i; i < tree.size(); i += i & (~i + 1)) tree[i] += w;
  }
  // the sum of the weights of the indices < i
  double prefix(std::size_t i) const {
    double s = 0.0;
    for (; i > 0; i -= i & (~i + 1)) s += tree[i];
    return s;
  }
private:
  std::vector<double> tree;
};

}  // namespace

Concordance concordance(const Eigen::VectorXd& pred,
                        const Eigen::VectorXd& start,
                        const Eigen::VectorXd& event,
                        const Eigen::VectorXi& status,
                        const Eigen::VectorXd& weight) {
  const int m = pred.size();
  if (start.size() != m || event.size() != m || status.size() != m || weight.size() != m) {
    ERROR_MSG("concordance: pred, start, event, status and weight must have the same length");
  }
  Concordance out{0.0, 0.0};
  if (m == 0) return out;

  std::vector<int> by_stop(m), by_start(m), order(m);
  std::iota(by_stop.begin(), by_stop.end(), 0);
  std::iota(by_start.begin(), by_start.end(), 0);
  std::sort(by_stop.begin(), by_stop.end(),
            [&](int a, int b) { return event(a) > event(b); });
  std::sort(by_start.begin(), by_start.end(),
            [&](int a, int b) { return start(a) > start(b); });

  // dense ranks of pred, ties sharing a rank
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int a, int b) { return pred(a) < pred(b); });
  std::vector<std::size_t> rank(m);
  std::size_t r = 0;
  for (int k = 0; k < m; ++k) {
    if (k > 0 && pred(order[k]) > pred(order[k - 1])) ++r;
    rank[order[k]] = r;
  }

  // the subjects in decreasing time; the tree holds the weights of those
  // at risk at the current time
  FenwickSum tree(r + 1);
  double total = 0.0;
  int p_stop = 0, p_start = 0;
  while (p_stop < m) {
    // up to the next failure time, everyone with a later stop time
    while (p_stop < m && !status(by_stop[p_stop])) {
      int j = by_stop[p_stop++];
      tree.add(rank[j], weight(j));
      total += weight(j);
    }
    if (p_stop == m) break;
    const double t = event(by_stop[p_stop]);

    // those stopping at t: censored ones are at risk, failures are compared
    int group_end = p_stop;
    while (group_end < m && event(by_stop[group_end]) == t) {
      int j = by_stop[group_end++];
      if (!status(j)) {
        tree.add(rank[j], weight(j));
        total += weight(j);
      }
    }
    // those whose interval starts at or after t are not at risk at t
    while (p_start < m && start(by_start[p_start]) >= t) {
      int j = by_start[p_start++];
      tree.add(rank[j], -weight(j));
      total -= weight(j);
    }
    for (int k = p_stop; k < group_end; ++k) {
      int i = by_stop[k];
      if (!status(i)) continue;
      double less = tree.prefix(rank[i]);
      double equal = tree.prefix(rank[i] + 1) - less;
      out.numerator += weight(i) * (less + 0.5 * equal);
      out.denominator += weight(i) * total;
    }
    // the failures at t are at risk at earlier times
    for (int k = p_stop; k < group_end; ++k) {
      int i = by_stop[k];
      if (status(i)) {
        tree.add(rank[i], weight(i));
        total += weight(i);
      }
    }
    p_stop = group_end;
  }
  return out;
}

Concordance CoxDeviance::concordance(const Eigen::VectorXd& pred,
                                     const Eigen::VectorXd& weight) const {
  if (pred.size() != raw_event.size() || weight.size() != raw_event.size()) {
    ERROR_MSG("concordance: pred and weight must have one entry per observation");
  }
  return ::concordance(pred, raw_start, raw_event, raw_status, weight);
}

Concordance StratifiedCoxDeviance::concordance(const Eigen::VectorXd& pred,
                                               const Eigen::VectorXd& weight) const {
  int n = pred.size();
  int total = 0;
  for (const auto& indices : stratum_indices) total += indices.size();
  if (n != total || weight.size() != n) {
    ERROR_MSG("concordance: pred and weight must have one entry per observation");
  }
  Concordance out{0.0, 0.0};
  for (size_t s = 0; s < unique_strata.size(); ++s) {
    const auto& indices = stratum_indices[s];
    Eigen::VectorXd ps(indices.size()), ws(indices.size());
    for (size_t j = 0; j < indices.size(); ++j) {
      ps(j) = pred(indices[j]);
      ws(j) = weight(indices[j]);
    }
    Concordance c = cox_devs[s]->concordance(ps, ws);
    out.numerator += c.numerator;
    out.denominator += c.denominator;
  }
  return out;
}
