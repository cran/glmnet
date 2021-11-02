#include <Rcpp.h>

using namespace Rcpp;

// This is the new C++ way of simplifying the progress bar.
// We directly call the setTxtProgressBar function with the given value.
void setpb_cpp(SEXP pb, int val) {
    Environment utils = Environment::namespace_env("utils");
    Function f = utils["setTxtProgressBar"];
    f(Named("pb", pb), Named("value", val));
}
