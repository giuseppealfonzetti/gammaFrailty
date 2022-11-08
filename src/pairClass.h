#ifndef pairClass_H
#define pairClass_H

#include <boost/math/special_functions.hpp>
#include "utils.h"

class pair_class{

public:

    // Inputs
    unsigned int j, jp;        // periods indices
    unsigned int n_j, n_jp;    // count values
    double alpha_j, alpha_jp;  // intercepts
    Eigen::VectorXd x;         // external regressors
    Eigen::VectorXd beta;      // regression coefficients
    double xi;                 // inverse of gamma parameetrs
    double rho;                // correlation adjacent periods
    double lxi;
    double artanhrho;

    // Derived helpers
    unsigned int p;            // number of periods
    unsigned int r;            // number of regressors
    unsigned int d;            // number of parameters
    unsigned int adiff;        // absolute difference indices
    double drho_jjp;           // derivative of rho^{|j-j'|}

    // Intermediate quantities
    double m_1, m_2;
    double u_j, u_jp;
    double log_u_j, log_u_jp;
    double rho_jjp;
    double D_j, D_jp;
    double Delta;
    double f;
    double S;

    //////////////////////////
    //////////////////////////
    /// INTERNALFUNCTIONS  ///
    //////////////////////////
    //////////////////////////

    // SETUP input quantities
    void setup(
            unsigned int J,
            unsigned int JP,
            unsigned int N_J,
            unsigned int N_JP,
            Eigen::VectorXd X,
            Eigen::VectorXd BETA,
            double ALPHA_J,
            double ALPHA_JP,
            double LXI,
            double artanhRHO,
            unsigned int P
    ){
        j = J;
        jp = JP;
        n_j = N_J;
        n_jp = N_JP;
        x = X;
        beta = BETA;
        alpha_j = ALPHA_J;
        alpha_jp = ALPHA_JP;
        // xi = exp(LXI);
        xi = exp(-LXI);
        lxi = LXI;
        //xi = rofz_cpp(LXI) + 1.1;
        //xi = pow(LXI, 2);
        artanhrho = artanhRHO;
        rho = rofz_cpp(artanhRHO);
        p = P;
        r = BETA.size();
        d = p + r + 2;
    }

    // public intermediate computing (only for debugging)
    double compute_Q_(unsigned int ind, bool verboseFLAG = false){
        double Q = 1;
        for(unsigned int ind2 = m_2; ind2 <= m_1 + m_2 - ind - 1; ind2++){
            Q *= 1 + ind2*xi;
        }

        return Q;
    }
    double compute_Sind_(unsigned int ind, bool verboseFLAG = false){
        double Q = compute_Q_(ind);
        double S = pow(-xi*f, ind) * Rf_choose(m_1, ind) * Rf_choose(m_2, ind) * boost::math::factorial<double>(ind) * Q;

        if(verboseFLAG)Rcpp::Rcout << "pow(-xi*f, ind):"<<pow(-xi*f, ind)<<"\nRf_choose(m_1, ind):"<<Rf_choose(m_1, ind)
                                   <<"\nRf_choose(m_2, ind):"<<Rf_choose(m_2, ind)<<"\nboost::math::factorial<double>(ind):"
                                   << boost::math::factorial<double>(ind) << "\nQ:"<<Q<<"\n";
        return S;
    }

    // compute intermediate quantities and derivatives
    void compute_intermediate_(){
        compute_m_();
        compute_rho_();
        compute_u_();
        compute_Delta_();
        compute_D_();
        compute_f_();
        //compute_S_();
    }
    void compute_dintermediate_(bool verboseS = false){

        compute_du_();
        compute_dDelta_();
        compute_dD_();
        compute_df_();
        compute_dlogS_(verboseS);
    }

    // return loglikelihood and gradient contribution
    double compute_ll_(){

        double out = 0;
        double tmp0 = 0;
        for(unsigned int ind = 0; ind < m_2; ind++){
            tmp0 += log(1 + ind*xi);
        }
        double tmp1 = n_j*log_u_j + n_jp*log_u_jp +
            n_j*log(D_j) + n_jp*log(D_jp) - (n_j + n_jp + 1/xi)*log(Delta);

        out = tmp0 + tmp1 + log(S);

        return out;
    }
    Eigen::VectorXd compute_gradient_(){
        Eigen::VectorXd out = Eigen::VectorXd::Zero(d);

        double tmp_xi= 0;
        for(unsigned int ind = 0; ind < m_2; ind++){
            tmp_xi += ind/(1 + ind*xi);
        }

        Eigen::VectorXd tmp0 = -pow(xi*Delta,-1)*dDelta; tmp0(0) += pow(xi, -2)*log(Delta);
        Eigen::VectorXd tmp1 = (n_j/u_j)*du_j + (n_jp/u_jp)*du_jp + (n_j/D_j)*dD_j + (n_jp/D_jp)*dD_jp - ((n_j+n_jp)/Delta)*dDelta;
        out = tmp0 + tmp1 + dlogS;
        out[0] += tmp_xi;

        // out[0] *= xi;
        out[0] *= -xi;
        //out[0] *= 2*lxi;
        //out[0] *= drofz_cpp(lxi);
        out[1] *= drofz_cpp(artanhrho);

        return out;
    }

    // return lists with intermediate quantities (only for debugging)
    Rcpp::List get_intermediate_(){

        Rcpp::List output =
            Rcpp::List::create(
                Rcpp::Named("m_1"    ) = m_1,
                Rcpp::Named("m_2"    ) = m_2,
                Rcpp::Named("u_j"    ) = u_j,
                Rcpp::Named("u_jp"   ) = u_jp,
                Rcpp::Named("rho_jjp") = rho_jjp,
                Rcpp::Named("D_j"    ) = D_j,
                Rcpp::Named("D_jp"   ) = D_jp,
                Rcpp::Named("Delta"  ) = Delta,
                Rcpp::Named("f"      ) = f,
                Rcpp::Named("S"      ) = S
            );

        return output;
    }
    Rcpp::List get_dintermediate_(){

        Rcpp::List output =
            Rcpp::List::create(
                Rcpp::Named("du_j"    ) = du_j,
                Rcpp::Named("du_jp"   ) = du_jp,
                Rcpp::Named("drho_jjp") = drho_jjp,
                Rcpp::Named("dD_j"    ) = dD_j,
                Rcpp::Named("dD_jp"   ) = dD_jp,
                Rcpp::Named("dDelta"  ) = dDelta,
                Rcpp::Named("df"      ) = df,
                Rcpp::Named("dlogS"   ) = dlogS
            );

        return output;
    }

private:

    // Compute intermediate quantities
    void   compute_m_(){
        m_1 = std::min(n_j, n_jp);
        m_2 = std::max(n_j, n_jp);

    }
    void   compute_rho_(){
        int diff = j-jp;
        adiff = abs(diff);
        rho_jjp = pow(rho, adiff);
        drho_jjp = pow(rho, adiff-1)*adiff;
    }
    void   compute_u_(){
        log_u_j = alpha_j + beta.dot(x);
        u_j = exp(log_u_j);
        log_u_jp = alpha_jp + beta.dot(x);
        u_jp = exp(log_u_jp);
    }
    void   compute_Delta_(){
        Delta = 1+ xi*(u_j+u_jp) + pow(xi,2)*u_j*u_jp*(1-rho_jjp);
    }
    void   compute_D_(){
        D_j = 1 + xi*u_jp*(1-rho_jjp);
        D_jp = 1 + xi*u_j*(1-rho_jjp);
    }
    void   compute_f_(){
        f = (Delta*(1-rho_jjp))/(D_j*D_jp);
    }
    void compute_S_(bool verboseFLAG = false){
        double tmp = 0;
        for(unsigned int ind = 0; ind <= m_1; ind++){
            tmp += compute_Sind_(ind);
        }

        S = tmp;
    }

    /////////////////
    // DERIVATIVES //
    /////////////////
    // (xi, rho,   beta   , alpha_1,..., alpha_p)
    // (0 , 1  , 2,...,r+1, r+2    ,..., r+p+1  )

    Eigen::VectorXd du_j; Eigen::VectorXd du_jp;
    void compute_du_(){
        du_j.resize(d); du_jp.resize(d);

        // wrt xi and rho
        du_j.setZero(); du_jp.setZero();

        // wrt beta
        du_j.segment(2,r) = u_j*x; du_jp.segment(2,r) = u_jp*x;

        // wrt alpha
        du_j(r+2+j) = u_j; du_jp(r+2+jp) = u_jp;
    }

    double dDelta_u_j, dDelta_u_jp; Eigen::VectorXd dDelta;
    void compute_dDelta_(){
        dDelta.resize(d);
        dDelta.setZero();

        // partial: wrt u_j and u_jp
        dDelta_u_j = xi*D_j; dDelta_u_jp = xi*D_jp;

        // wrt xi
        dDelta(0) = u_j + u_jp + 2*xi*u_j*u_jp*(1-rho_jjp);

        // wrt rho
        dDelta(1) = -pow(xi,2)*u_j*u_jp*drho_jjp;

        // wrt beta
        dDelta.segment(2,r) = dDelta_u_j * du_j.segment(2,r) + dDelta_u_jp * du_jp.segment(2,r);

        // wrt alpha
        dDelta(r+2+j ) = dDelta_u_j  * du_j( r+2+j);
        dDelta(r+2+jp) = dDelta_u_jp * du_jp(r+2+jp);
    }

    double dD_u; Eigen::VectorXd dD_j, dD_jp;
    void compute_dD_(){
        dD_j.resize(d); dD_jp.resize(d);
        dD_j.setZero(); dD_jp.setZero();

        // partial: D_j wrt u_jp and D_jp wrt u_j
        dD_u = xi*(1-rho_jjp);

        // wrt xi
        dD_j(0)  = u_jp*(1-rho_jjp);
        dD_jp(0) = u_j *(1-rho_jjp);

        // wrt rho
        dD_j(1)  = -xi*u_jp*drho_jjp;
        dD_jp(1) = -xi*u_j *drho_jjp;

        // wrt beta
        dD_j.segment(2,r)  = dD_u*du_jp.segment(2,r);
        dD_jp.segment(2,r) = dD_u* du_j.segment(2,r);

        // wrt alpha
        dD_j(r+2+jp) = dD_u*du_jp(r+2+jp);
        dD_jp(r+2+j) = dD_u*du_j( r+2+j );
    }

    double df_D_j, df_D_jp; double df_Delta; Eigen::VectorXd df;
    void compute_df_(){
        df.resize(d);
        df.setZero();

        // partial: wrt D_j and D_jp
        df_D_j  = -Delta * ( 1 - rho_jjp )/( pow(D_j, 2) * D_jp);
        df_D_jp = -Delta * ( 1 - rho_jjp )/( pow(D_jp,2) * D_j );

        // partial: wrt Delta
        df_Delta = (1 - rho_jjp)/(D_j*D_jp);

        // wrt xi
        df(0) =  df_Delta*dDelta(0) + df_D_j*dD_j(0) + df_D_jp*dD_jp(0);

        // wrt rho
        df(1) = df_Delta*dDelta(1) + df_D_j*dD_j(1) + df_D_jp*dD_jp(1) - (Delta/(D_j*D_jp))*drho_jjp;

        // wrt beta
        df.segment(2,r) = df_Delta*dDelta.segment(2,r) + df_D_j*dD_j.segment(2,r) + df_D_jp*dD_jp.segment(2,r);

        // wrt alpha
        df(r+2+j ) = df_Delta*dDelta(r+2+j ) + df_D_jp*dD_jp(r+2+j);
        df(r+2+jp) = df_Delta*dDelta(r+2+jp) + df_D_j *dD_j(r+2+jp);
    }

    double compute_dQ_xi_dividedbyQ_(unsigned int ind, bool verboseFLAG = false){
        double tmp_sum = 0;
        for(unsigned int ind2 = m_2; ind2 <= m_1 + m_2 - ind - 1; ind2++){
            double tmp = 1 + ind2*xi;
            tmp_sum += ind2/tmp;
        }

        return tmp_sum;
    }

    Eigen::VectorXd compute_dSind_(unsigned int ind, bool verboseFLAG = false){
        Eigen::VectorXd out = Eigen::VectorXd::Zero(d);

        double Sind = compute_Sind_(ind);

        // wrt xi
        out(0) = Sind*(compute_dQ_xi_dividedbyQ_(ind) + ind*pow(xi,-1) + ind*pow(f,-1)*df(0));

        // wrt rho
        out(1) = Sind*(ind/f)*df(1);

        // wrt beta
        out.segment(2,r) = Sind*(ind/f)*df.segment(2,r);

        // wrt alpha
        out(r+2+jp) = Sind*(ind/f)*df(r+2+jp);
        out(r+2+j)  = Sind*(ind/f)*df(r+2+j );

        return out;
    }

    Eigen::VectorXd dlogS;
    void compute_dlogS_(bool verboseFLAG = false){
        dlogS.resize(d);
        dlogS.setZero();

        double tmp_sum = 0;
        Eigen::VectorXd tmp_dSum = Eigen::VectorXd::Zero(d);

        for(unsigned int ind = 0; ind <= m_1; ind++){
            double Sind = compute_Sind_(ind);
            if(verboseFLAG) Rcpp::Rcout << "index:" << ind << " Sind:" << Sind << "\n";
            tmp_sum  += Sind;

            Eigen::VectorXd out = Eigen::VectorXd::Zero(d);
            {
                // wrt xi
                out(0) = Sind*(compute_dQ_xi_dividedbyQ_(ind) + ind*pow(xi,-1) + ind*pow(f,-1)*df(0));

                // wrt rho
                out(1) = Sind*(ind/f)*df(1);

                // wrt beta
                out.segment(2,r) = Sind*(ind/f)*df.segment(2,r);

                // wrt alpha
                out(r+2+jp) = Sind*(ind/f)*df(r+2+jp);
                out(r+2+j)  = Sind*(ind/f)*df(r+2+j );
            }

            tmp_dSum += out;
        }


        S = tmp_sum;
        dlogS = tmp_dSum/S;
    }

};

#endif
