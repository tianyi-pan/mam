// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_mam_TMBExports
#include <TMB.hpp>
#include "conditional_model.hpp"
#include "marginal_mean.hpp"
#include "marginal_mean2.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "conditional_model") {
    return conditional_model(this);
  } else if(model == "marginal_mean") {
    return marginal_mean(this);
  } else if(model == "marginal_mean2") {
    return marginal_mean2(this);
  } else {
    Rf_error("Unknown model.");
  }
  return 0;
}
