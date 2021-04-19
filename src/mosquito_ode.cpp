/*
 * mosquito_ode.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include <individual.h>
#include "mosquito_ode.h"
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

integration_function_t create_ode(MosquitoModel& model) {
    return [&model](const state_t& x , state_t& dxdt , double t) {
        auto K = carrying_capacity(
            t,
            model.model_seasonality,
            model.days_per_timestep,
            model.g0,
            model.g,
            model.h,
            model.K0,
            model.R_bar
        );
	
	double g_t,g_mTL;
	int t_dlay,T_eq;
	std::ofstream histories_m,histories_f;
	T_eq=(1/model.mue)*log(0.5*model.beta/model.mum);


	t_dlay=t-x[get_idx(ODEState::Tl)];

	g_t=x[get_idx(ODEState::Fl)]/(model.KF+x[get_idx(ODEState::Fl)]);//larval development rate function
	g_mTL=model.history_f[t_dlay+1.5*T_eq]/(model.KF + model.history_f[t_dlay+1.5*T_eq]);//larval development rate at 
	//time delay TL

	dxdt[get_idx(ODEState::Fl)]=model.G0-model.Amax*x[get_idx(ODEState::L)]*x[get_idx(ODEState::Fl)]/(model.KF+x[get_idx(ODEState::Fl)]);//larval food supply

        dxdt[get_idx(ODEState::L)] = model.beta * (model.total_M) //new eggs
            - model.beta*model.history_m[t_dlay+1.5*T_eq]*exp(-model.mue*x[get_idx(ODEState::Tl)])*g_t/g_mTL //growth to pupal stage
            - x[get_idx(ODEState::L)] * model.mue; //larval deaths * (1 + (x[get_idx(ODEState::L)] + x[get_idx(ODEState::L)]) / K); 
        
	dxdt[get_idx(ODEState::P)] = model.beta*model.history_m[t_dlay+1.5*T_eq]*exp(-model.mue*x[get_idx(ODEState::Tl)])*g_t/g_mTL
	       	//growth to pupae
            - x[get_idx(ODEState::P)] / model.dp //growth to adult
            - x[get_idx(ODEState::P)] * model.mup; // death of pupae
	
	dxdt[get_idx(ODEState::Tl)]=1-g_t/g_mTL;//time varying larval development rate
	
	if (/*t>60 && */ t<2 /*&& int(100*t) % 5==0*/ ) Rcpp::Rcout  << "dxdt[0] " << dxdt[0]  << " t " << t << " t_dlay " << t-x[get_idx(ODEState::Tl)] <<  endl;
	if (/*t>60 && */  t<2) Rcpp::Rcout << "Fl " << x[get_idx(ODEState::Fl)] << " L " << x[get_idx(ODEState::L)] << " P " << x[get_idx(ODEState::P)]<< " Tl " << x[get_idx(ODEState::Tl)] << " t  " << t << endl;
	if (/*t>60 && */  t<2) Rcpp::Rcout << " beta " <<  model.beta  << " total_M " << model.total_M << endl;
	if (/*t>60 && */  t<2) Rcpp::Rcout << "g_t " << g_t << " g_mTL " << g_mTL << endl;


	if (t==1000) {
		histories_m.open("histories_m.txt");
		histories_f.open("histories_f.txt");
		Rcpp::Rcout << " history_f " << model.history_f[1] << endl;
		for (std::vector<double>::iterator i=model.history_f.begin(); i<model.history_f.end(); i++){
			histories_f << *i << endl;
			//histories_m << model.history_m[i] << endl;
		}
	}
	//Rcpp:Rcout << " KF " << model.KF << " mue " << model.mue << endl;
	//Rcpp::Rcout << " G0 " << model.G0 << endl;
    };
}

MosquitoModel::MosquitoModel(
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
    double G0,
    double KF,
    double Amax, 
    double mum,
    double R_bar
    ):
    beta(beta),
    de(de),
    mue(mue),
    K0(K0),
    gamma(gamma),
    dl(dl),
    mul(mul),
    dp(dp),
    mup(mup),
    total_M(total_M),
    model_seasonality(model_seasonality),
    days_per_timestep(days_per_timestep),
    g0(g0),
    g(g),
    h(h),
    history_f(history_f), 
    history_m(history_m),
    G0(G0),
    KF(KF),
    Amax(Amax),
    mum(mum),
    R_bar(R_bar)
    {
    for (auto i = 0u; i < state.size(); ++i) {
        state[i] = init[i];
    }
    ode = create_ode(*this);
    rk = boost::numeric::odeint::make_dense_output(
        a_tolerance,
        r_tolerance,
        boost::numeric::odeint::runge_kutta_dopri5<state_t>()
    );
}

void MosquitoModel::step(size_t new_total_M) {
    total_M = new_total_M;
    boost::numeric::odeint::integrate_const(stepper, ode, state, t, t + dt, 0.05);
    history_f.push_back(state[get_idx(ODEState::Fl)]);
    history_m.push_back(total_M);
    //Rcpp::Rcout << " t " << t << " dt " << dt << endl;
    //Rcpp::Rcout << "history_f " << history_f[t] << " t " << t <<endl;
    ++t;
}

state_t MosquitoModel::get_state() {
    return state;
}

//[[Rcpp::export]]
Rcpp::XPtr<MosquitoModel> create_mosquito_model(
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
    double G0,
    double KF,
    double Amax, 
    double mum,
    double R_bar
    ) {
    auto model = new MosquitoModel(
        init,
        beta,
        de,
        mue,
        K0,
        gamma,
        dl,
        mul,
        dp,
        mup,
        total_M,
        model_seasonality,
        days_per_timestep,
        g0,
        g,
        h,
	history_f,
	history_m,
	G0,
	KF,
	Amax,
	mum,
        R_bar	
    );
    return Rcpp::XPtr<MosquitoModel>(model, true);
}

//[[Rcpp::export]]
std::vector<double> mosquito_model_get_states(Rcpp::XPtr<MosquitoModel> model) {
    auto state = model->get_state();
    return std::vector<double>(state.cbegin(), state.cend());
}

//[[Rcpp::export]]
void mosquito_model_step(Rcpp::XPtr<MosquitoModel> model, size_t total_M) {
    model->step(total_M);
}
