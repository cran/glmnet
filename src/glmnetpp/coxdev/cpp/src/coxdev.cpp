#include "coxdev.h"
#include <cmath>
#include <iostream>
#include <numeric>

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

    forward_prework(ctx.status, w_avg_buffer, ctx.scaling, risk_sums, 2, 1, forward_scratch_buffer, dummy_map, true);
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
	T_2_term(i) -= forward_cumsum1(ctx.first(i));
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
    return a[i] < a[j];
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

