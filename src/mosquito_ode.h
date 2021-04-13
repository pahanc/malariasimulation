/*
 * mosquito_ode.h
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#ifndef SRC_MOSQUITO_ODE_H_
#define SRC_MOSQUITO_ODE_H_

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
// [[Rcpp::depends(BH)]]
#include <boost/numeric/odeint.hpp>
#pragma GCC diagnostic pop

#include <array>
#include <queue>
#include <cmath>
#include <functional>
#include <type_traits>
#include "mosquito_biology.h"

/*
 * The states are:
 * 0 - F  - Larval Food Supply
 * 1 - L  - Larval Stage
 * 2 - P  - Pupal Stage
 * 3 - T -  Larval Development Time
 */
enum class ODEState : size_t {Fl,L,P,Tl};

// Provide a convenience function for getting at the index of an enumeration
template<typename T>
constexpr auto get_idx(T value)
{
    return static_cast<std::underlying_type_t<T>>(value);
}

using state_t = std::array<double, 4>;
using integration_function_t = std::function<void (const state_t&, state_t&, double)>;

struct MosquitoModel {
    /* Parameters */
    const double beta; //egg laying rate
    const double de; //delay for early larval growth
    const double mue; //death rate for early larvae
    const double K0; //baseline carrying capacity
    const double gamma; //carrying capacity parameter for late larvae
    const double dl; //delay for late larval growth
    const double mul; //death rate for late larvae
    const double dp; //delay for for pupal growth
    const double mup; //death rate for pupae
    size_t total_M; //the number of adult female mosquitos in the model
    const bool model_seasonality; //whether to model seasonality
    const double days_per_timestep; //scale of the fourier model for seasonality
    const double g0; //fourier shape parameter
    const std::vector<double> g; //fourier shape parameters
    const std::vector<double> h; //fourier shape parameters
    std::vector<double> history_f; 
    std::vector<double> history_m; 
    const double R_bar; //average rainfall
    const double G0;
    const double KF;
    const double Amax;
    std::queue<double> lagged_incubating; //last tau values for incubating mosquitos

    MosquitoModel(
        std::vector<double> init,
        double beta,
        double de,
        double mue,
        double K0,
        double gamma,
        double dl,
        double mul,
        double dp,
        double mup,
        size_t total_M,
        bool model_seasonality,
        double days_per_timestep,
        double g0,
        std::vector<double> g,
        std::vector<double> h,
	std::vector<double> history_f, 
        std::vector<double> history_m, 
        double R_bar,
	double G0,
	double KF,
	double Amax 
    );
    virtual void step(size_t);
    virtual state_t get_state();
    virtual ~MosquitoModel() {};
private:
    //solver fields
    boost::numeric::odeint::dense_output_runge_kutta<
        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_dopri5<state_t>
        >
    >rk;
    boost::numeric::odeint::runge_kutta4< state_t > stepper;
    const double r_tolerance = 1.0e-6;
    const double a_tolerance = 1.0e-6;
    integration_function_t ode;

    double t = 0.;
    const double dt = 1.;
    state_t state;
};

integration_function_t create_ode(MosquitoModel& model);

#endif /* SRC_MOSQUITO_ODE_H_ */
