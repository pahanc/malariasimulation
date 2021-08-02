/*
 * mosquito_emergence.cpp
 *
 *  Created on: 4 Aug 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "mosquito_emergence.h"
#include "mosquito_ode.h"
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;



//' @title Mosquito emergence process
//' @description Move mosquitos from NonExistent to Sm in line with the number of
//' pupals in the ODE models
//'
//' @param odes a list of odes for each species of mosquito
//' @param state the variable for the mosquito state
//' @param species the variable for the mosquito species
//' @param species_names a vector of category names for the species variable
////' @param mosq_suppression a vector of proportional suppresssion of mosquito emergence (due to GM mossies)
//' @param dpl the delay for pupal growth (in timesteps)
//[[Rcpp::export]]
Rcpp::XPtr<process_t> create_mosquito_emergence_process_cpp(
    Rcpp::List odes,
    Rcpp::XPtr<CategoricalVariable> state,
    Rcpp::XPtr<CategoricalVariable> species,
    std::vector<std::string> species_names,
    std::vector<double> mosq_suppression,
    double dpl) {
    //double mosq_suppression=0.5;
    auto rate = .5 * 1./dpl;
    return Rcpp::XPtr<process_t>(
        new process_t([=] (size_t t) {
	    auto n = 0u;
	    //Rcpp::Rcout << "mosq_suppression[t] " << mosq_suppression[t] << endl;
	    for (Rcpp::XPtr<MosquitoModel> ode : odes) {
                n += ode->get_state()[get_idx(ODEState::P)] * rate* mosq_suppression[t];
            }
            auto source = state->get_index_of(std::vector<std::string>{"NonExistent"});
            if (source.size() < n) {
                std::stringstream m;
                m << "Not enough mosquitos (short by ";
                m << n - source.size();
                m << "). Please raise the value of parameters$mosquito_limit. ";
                m << "If you have used `parameterise_mosquito_equilibrium`, ";
                m << "your seasonality parameters may have lead to more mosquitoes than expected.";
                Rcpp::stop(m.str());
            }
            if (n > 0) {
                auto species_i = 0u;
                auto sourceit = source.begin();
                for (Rcpp::XPtr<MosquitoModel> ode : odes) {
                    auto to_set = static_cast<size_t>(ode->get_state()[get_idx(ODEState::P)] * rate  * mosq_suppression[t]);
                    auto target = individual_index_t(state->size);
                    for (auto i = 0u; i < to_set; ++i) {
                        target.insert(*sourceit);
                        ++sourceit;
                    }
                    state->queue_update("Sm", target);
                    species->queue_update(species_names[species_i], target);
                    ++species_i;
                }
            }
        }),
        true
    );
}
