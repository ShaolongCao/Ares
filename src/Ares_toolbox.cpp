#include <cmath>
#include <algorithm>  
#include <vector>
#include <limits>
#include <stdlib.h>
#include <math.h>
using namespace Rcpp;
//using namespace arma;

// [[Rcpp::export]]
double Mixed_Binomial_LogLikelihood_repara(NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr, double omega, double mu){
	NumericVector prob=purity*CNv*omega/(2.0*(1.0-purity)/mu+purity*(CNv*omega+CNr));
	NumericVector theta = V/R;
	NumericVector res = Pi_nk*R*(theta-prob)*(theta-prob)/(theta*(1.0-theta));
	return(sum(res));
}

// [[Rcpp::export]]
double Min_Mu_Binomial_repara(double lower, double upper, double tol, NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr, double omega) {
   /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = std::numeric_limits<double>::epsilon();
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = lower;
    b = upper;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
	
	fx = Mixed_Binomial_LogLikelihood_repara(Pi_nk, V, R, purity, CNv, CNr, omega, x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;
        
        /* check stopping criterion */
        
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1) { /* fit parabola */
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
        }
        
        if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        }
        else { /* a parabolic-interpolation step */
            
            d = p / q;
            u = x + d;
            
            /* f must not be evaluated too close to ax or bx */
            
            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }
        
        /* f must not be evaluated too close to x */
        
        if (fabs(d) >= tol1)
            u = x + d;
        else if (d > 0.)
            u = x + tol1;
        else
            u = x - tol1;
        
		fu = Mixed_Binomial_LogLikelihood_repara(Pi_nk, V, R, purity, CNv, CNr, omega, u);
        
        /*  update  a, b, v, w, and x */
        
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    /* end of main loop */
    
    return x;
}

// [[Rcpp::export]]
double Min_Omega_Binomial_repara(double lower, double upper, double tol, NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr, double mu) {
   /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = std::numeric_limits<double>::epsilon();
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = lower;
    b = upper;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
	fx = Mixed_Binomial_LogLikelihood_repara(Pi_nk, V, R, purity, CNv, CNr, x, mu);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;
        
        /* check stopping criterion */
        
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1) { /* fit parabola */
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
        }
        
        if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        }
        else { /* a parabolic-interpolation step */
            
            d = p / q;
            u = x + d;
            
            /* f must not be evaluated too close to ax or bx */
            
            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }
        
        /* f must not be evaluated too close to x */
        
        if (fabs(d) >= tol1)
            u = x + d;
        else if (d > 0.)
            u = x + tol1;
        else
            u = x - tol1;
        fu = Mixed_Binomial_LogLikelihood_repara(Pi_nk, V, R, purity, CNv, CNr, u, mu);
        
        /*  update  a, b, v, w, and x */
        
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    /* end of main loop */
    
    return x;
}

// [[Rcpp::export]]
List Solve_mixed_Binomial_repara(NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr){
  double omega=1.0;
  double mu=1.0;
  double err=1.0;
  int iter=1;
  double obj=0.0;
  double obj_old=0.0;
  obj_old=Mixed_Binomial_LogLikelihood_repara(Pi_nk, V, R, purity, CNv, CNr, omega, mu);
  while(err>1e-5 && iter<100){
	omega=Min_Omega_Binomial_repara(0.01, 100.0, 1e-6, Pi_nk, V, R, purity, CNv, CNr, mu);
	mu=Min_Mu_Binomial_repara(0.01, 100.0, 1e-6, Pi_nk, V, R, purity, CNv, CNr, omega);
	obj=Mixed_Binomial_LogLikelihood_repara(Pi_nk, V, R, purity, CNv, CNr, omega, mu);
    iter++;
    err=fabs(obj_old-obj)/obj_old;
    obj_old=obj;
  }
  return(List::create(Named("omega")=omega,Named("mu")=mu));
}


// [[Rcpp::export]]
double Mixed_Binomial_LogLikelihood_newton(NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr, double omega, double mu){
	NumericVector theta=purity*CNv*omega*omega/(2*(1-purity)/mu/mu+purity*(CNv*omega*omega+CNr));
	NumericVector res=Pi_nk*(V*log(theta)+(R-V)*log(1-theta));
	return(-1.0*sum(res));
}

// [[Rcpp::export]]
double Mixed_Binomial_LogLikelihood(NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr, double omega, double mu){
	NumericVector theta=purity*CNv*omega/(2*(1-purity)/mu+purity*(CNv*omega+CNr));
	NumericVector res=Pi_nk*(V*log(theta)+(R-V)*log(1-theta));
	return(-1.0*sum(res));
}

// [[Rcpp::export]]
double Min_Mu_Binomial(double lower, double upper, double tol, NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr, double omega) {
   /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = std::numeric_limits<double>::epsilon();
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = lower;
    b = upper;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
	
	fx = Mixed_Binomial_LogLikelihood(Pi_nk, V, R, purity, CNv, CNr, omega, x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;
        
        /* check stopping criterion */
        
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1) { /* fit parabola */
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
        }
        
        if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        }
        else { /* a parabolic-interpolation step */
            
            d = p / q;
            u = x + d;
            
            /* f must not be evaluated too close to ax or bx */
            
            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }
        
        /* f must not be evaluated too close to x */
        
        if (fabs(d) >= tol1)
            u = x + d;
        else if (d > 0.)
            u = x + tol1;
        else
            u = x - tol1;
        
		fu = Mixed_Binomial_LogLikelihood(Pi_nk, V, R, purity, CNv, CNr, omega, u);
        
        /*  update  a, b, v, w, and x */
        
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    /* end of main loop */
    
    return x;
}

// [[Rcpp::export]]
double Min_Omega_Binomial(double lower, double upper, double tol, NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr, double mu) {
   /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = std::numeric_limits<double>::epsilon();
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = lower;
    b = upper;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
	fx = Mixed_Binomial_LogLikelihood(Pi_nk, V, R, purity, CNv, CNr, x, mu);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;
        
        /* check stopping criterion */
        
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1) { /* fit parabola */
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
        }
        
        if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        }
        else { /* a parabolic-interpolation step */
            
            d = p / q;
            u = x + d;
            
            /* f must not be evaluated too close to ax or bx */
            
            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }
        
        /* f must not be evaluated too close to x */
        
        if (fabs(d) >= tol1)
            u = x + d;
        else if (d > 0.)
            u = x + tol1;
        else
            u = x - tol1;
        fu = Mixed_Binomial_LogLikelihood(Pi_nk, V, R, purity, CNv, CNr, u, mu);
        
        /*  update  a, b, v, w, and x */
        
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    /* end of main loop */
    
    return x;
}

// [[Rcpp::export]]
List Solve_mixed_Binomial_Goldensection(NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr){
  double omega=1.0;
  double mu=1.0;
  double err=1.0;
  int iter=1;
  double obj=0.0;
  double obj_old=0.0;
  obj_old=Mixed_Binomial_LogLikelihood(Pi_nk, V, R, purity, CNv, CNr, omega, mu);
  while(err>1e-5 && iter<100){
	omega=Min_Omega_Binomial(0.01, 100.0, 1e-6, Pi_nk, V, R, purity, CNv, CNr, mu);
	mu=Min_Mu_Binomial(0.01, 100.0, 1e-6, Pi_nk, V, R, purity, CNv, CNr, omega);
	obj=Mixed_Binomial_LogLikelihood(Pi_nk, V, R, purity, CNv, CNr, omega, mu);
    iter++;
    err=fabs(obj_old-obj)/obj_old;
    obj_old=obj;
  }
  return(List::create(Named("omega")=omega,Named("mu")=mu));
}

// [[Rcpp::export]]
List Solve_mixed_Binomial_newton_Goldensection(NumericVector Pi_nk, NumericVector V, NumericVector R, NumericVector purity, NumericVector CNv, NumericVector CNr){
	int n=V.size();
  double omega=1.0;
  double mu=1.0;
  double err=1.0;
  int iter=1;
  double obj=0.0;
  double obj_old=0.0;
	obj_old=Mixed_Binomial_LogLikelihood(Pi_nk, V, R, purity, CNv, CNr, omega, mu);
  while(err>1e-3 && iter<20){
	omega=Min_Omega_Binomial(0.01, 100.0, 1e-4, Pi_nk, V, R, purity, CNv, CNr, mu);
	mu=Min_Mu_Binomial(0.01, 100.0, 1e-4, Pi_nk, V, R, purity, CNv, CNr, omega);
	obj=Mixed_Binomial_LogLikelihood(Pi_nk, V, R, purity, CNv, CNr, omega, mu);
    iter++;
    err=fabs(obj_old-obj)/obj_old;
    obj_old=obj;
  }
  
	// Start Newton method
  omega=sqrt(omega);
  mu=sqrt(mu);
  err=1.0;
  iter=1;
  obj=0.0;
  obj_old=0.0;
  double lrate=0.5;
  NumericVector a(n);
  NumericVector b(n);
  NumericVector c(n);
  NumericVector theta(n);
  NumericVector FF(n);
  NumericVector FF2(n);
  
  NumericVector dThetadOmega(n);
  NumericVector dThetadMu(n);
  NumericVector d2ThetadOmega2(n);
  NumericVector d2ThetadMu2(n);
  NumericVector d2ThetadOmegaMu(n);
  
  a=purity*CNv;
  b=2.0*(1.0-purity);
  c=purity*CNr;
  theta=a*omega*omega/(b/mu/mu+a*omega*omega+c);
  FF= -Pi_nk*(V/theta+(R-V)/(theta-1));
  FF2=Pi_nk*(V/theta/theta+(R-V)/(theta-1)/(theta-1));
  
  dThetadOmega=2.0*a*(b/mu/mu+c)*omega/pow((a*omega*omega+b/mu/mu+c),2);
  dThetadMu=2.0*a*b*omega*omega*mu/pow((a*omega*omega+c)*mu*mu+b,2);
  d2ThetadOmega2 = -2.0*a*(b/mu/mu+c)*(3*a*omega*omega-b/mu/mu-c)/pow(a*omega*omega+b/mu/mu+c,2);
  d2ThetadMu2 = -2.0*a*b*omega*omega*(3*(a*omega*omega+c)*mu*mu-b)/pow((a*omega*omega+c)*mu*mu+b,3);
  d2ThetadOmegaMu = -4.0*a*b*omega*mu*(a*omega*omega*mu*mu-c*mu*mu-b)/pow((a*omega*omega+c)*mu*mu+b,3);
  
  double dFdOmega=sum(FF*dThetadOmega);
  double dFdMu=sum(FF*dThetadMu);
  
  double H11=sum(FF2*dThetadOmega*dThetadOmega+FF*d2ThetadOmega2);
  double H12=sum(FF2*dThetadOmega*dThetadMu+FF*d2ThetadOmegaMu);
  double H22=sum(FF2*dThetadMu*dThetadMu+FF*d2ThetadMu2);
  
  obj_old=Mixed_Binomial_LogLikelihood_newton(Pi_nk, V, R, purity, CNv, CNr, omega, mu);
  while(err>1e-6 && iter<50){
	theta=a*omega*omega/(b/mu/mu+a*omega*omega+c);
	FF= -Pi_nk*(V/theta+(R-V)/(theta-1));
	FF2=Pi_nk*(V/theta/theta+(R-V)/(theta-1)/(theta-1));
	
	dThetadOmega=2.0*a*(b/mu/mu+c)*omega/pow((a*omega*omega+b/mu/mu+c),2);
	dThetadMu=2.0*a*b*omega*omega*mu/pow((a*omega*omega+c)*mu*mu+b,2);
	d2ThetadOmega2 = -2.0*a*(b/mu/mu+c)*(3*a*omega*omega-b/mu/mu-c)/pow(a*omega*omega+b/mu/mu+c,2);
	d2ThetadMu2 = -2.0*a*b*omega*omega*(3*(a*omega*omega+c)*mu*mu-b)/pow((a*omega*omega+c)*mu*mu+b,3);
	d2ThetadOmegaMu = -4.0*a*b*omega*mu*(a*omega*omega*mu*mu-c*mu*mu-b)/pow((a*omega*omega+c)*mu*mu+b,3);
	dFdOmega=sum(FF*dThetadOmega);
	dFdMu=sum(FF*dThetadMu);
	H11=sum(FF2*dThetadOmega*dThetadOmega+FF*d2ThetadOmega2);
	H12=sum(FF2*dThetadOmega*dThetadMu+FF*d2ThetadOmegaMu);
	H22=sum(FF2*dThetadMu*dThetadMu+FF*d2ThetadMu2);
	
	omega = omega - lrate*(H22*dFdOmega-H12*dFdMu)/(H11*H22-H12*H12);
	mu = mu + lrate*(H12*dFdOmega-H11*dFdMu)/(H11*H22-H12*H12);
	obj=Mixed_Binomial_LogLikelihood_newton(Pi_nk, V, R, purity, CNv, CNr, omega, mu);
    iter++;
    err=fabs(obj_old-obj)/obj_old;
    obj_old=obj;
  }
  omega=omega*omega;
  mu=mu*mu;
  return(List::create(Named("omega")=omega, Named("mu")=mu));
}



/* SNP model */

// [[Rcpp::export]]
double Ares_mix_SNP_LogLikelihood(NumericVector Pi_nk, NumericVector B, NumericVector R, NumericVector purity, NumericVector CNA, NumericVector CNB, double alpha, double omega, double mu){
	/* Objective function*/
	NumericVector theta = 1/ (1 + ( 1-purity + purity*CNA*mu)/((1-purity)*alpha + purity*CNB*omega*mu)) ;
	NumericVector prob = B/R;
	NumericVector res = Pi_nk*R*(theta-prob)*(theta-prob)/(prob*(1.0-prob));
	return(sum(res));
}

// [[Rcpp::export]]
double Min_Alpha_Ares_mix_SNP(double lower, double upper, double tol, NumericVector Pi_nk, NumericVector B, NumericVector R, NumericVector purity, NumericVector CNA, NumericVector CNB, double omega, double mu) {
   /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = std::numeric_limits<double>::epsilon();
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = lower;
    b = upper;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
	fx = Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, x, omega, mu);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;
        
        /* check stopping criterion */
        
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1) { /* fit parabola */
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
        }
        
        if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        }
        else { /* a parabolic-interpolation step */
            
            d = p / q;
            u = x + d;
            
            /* f must not be evaluated too close to ax or bx */
            
            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }
        
        /* f must not be evaluated too close to x */
        
        if (fabs(d) >= tol1)
            u = x + d;
        else if (d > 0.)
            u = x + tol1;
        else
            u = x - tol1;
        fu = Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, u, omega, mu);
        
        /*  update  a, b, v, w, and x */
        
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    /* end of main loop */
    
    return x;
}

// [[Rcpp::export]]
double Min_Omega_Ares_mix_SNP(double lower, double upper, double tol, NumericVector Pi_nk, NumericVector B, NumericVector R, NumericVector purity, NumericVector CNA, NumericVector CNB, double alpha, double mu) {
   /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = std::numeric_limits<double>::epsilon();
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = lower;
    b = upper;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
	
	fx = Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha, x, mu);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;
        
        /* check stopping criterion */
        
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1) { /* fit parabola */
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
        }
        
        if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        }
        else { /* a parabolic-interpolation step */
            
            d = p / q;
            u = x + d;
            
            /* f must not be evaluated too close to ax or bx */
            
            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }
        
        /* f must not be evaluated too close to x */
        
        if (fabs(d) >= tol1)
            u = x + d;
        else if (d > 0.)
            u = x + tol1;
        else
            u = x - tol1;
        
		fu = Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha, u, mu);
        
        /*  update  a, b, v, w, and x */
        
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    /* end of main loop */
    
    return x;
}

// [[Rcpp::export]]
double Min_Mu_Ares_mix_SNP(double lower, double upper, double tol, NumericVector Pi_nk, NumericVector B, NumericVector R, NumericVector purity, NumericVector CNA, NumericVector CNB, double alpha, double omega) {
   /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = std::numeric_limits<double>::epsilon();
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = lower;
    b = upper;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
	
	fx = Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha, omega, x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;
        
        /* check stopping criterion */
        
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1) { /* fit parabola */
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
        }
        
        if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        }
        else { /* a parabolic-interpolation step */
            
            d = p / q;
            u = x + d;
            
            /* f must not be evaluated too close to ax or bx */
            
            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }
        
        /* f must not be evaluated too close to x */
        
        if (fabs(d) >= tol1)
            u = x + d;
        else if (d > 0.)
            u = x + tol1;
        else
            u = x - tol1;
        
		fu = Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha, omega, u);
        
        /*  update  a, b, v, w, and x */
        
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    /* end of main loop */
    
    return x;
}


// [[Rcpp::export]]
List Solve_Ares_mix_SNP(NumericVector Pi_nk, NumericVector B, NumericVector R, NumericVector purity, NumericVector CNA, NumericVector CNB){
  double alpha=1.0;
  double omega=1.0;
  double mu=1.0;
  double err=1.0;
  int iter=1;
  double obj=0.0;
  double obj_old=0.0;
  obj_old=Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha, omega, mu);
  while(err>1e-5 && iter<100){
	omega=Min_Omega_Ares_mix_SNP(0.01, 100.0, 1e-5, Pi_nk, B, R, purity, CNA, CNB, alpha, mu);
	mu=Min_Mu_Ares_mix_SNP(0.01, 100.0, 1e-5, Pi_nk, B, R, purity, CNA, CNB, alpha, omega);
	alpha=Min_Alpha_Ares_mix_SNP(0.01, 100.0, 1e-5, Pi_nk, B, R, purity, CNA, CNB, omega, mu);
	obj=Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha, omega, mu);
    iter++;
    err=fabs((obj_old-obj)/obj_old);
    obj_old=obj;
	//printf("%f\n",obj);
  }
  return(List::create(Named("alpha")=alpha,Named("omega")=omega,Named("mu")=mu));
}


// [[Rcpp::export]]
List Ares_mix_SNP_Mstep(NumericMatrix Pi, NumericVector B, NumericVector R, NumericVector purity, NumericVector CNA, NumericVector CNB, int K_cluster){
	NumericVector Alpha(K_cluster);
	NumericVector Omega(K_cluster);
	NumericVector Mu(K_cluster);
	int n=R.size();
	NumericVector Pi_nk(n);
	double alpha=1.0;
	double omega=1.0;
	double mu=1.0;
	double err=1.0;
	int iter=1;
	double obj=0.0;
	double obj_old=0.0;
	for(int i =0; i<K_cluster; i++){
		alpha=1.0;
		omega=1.0;
		mu=1.0;
		err=1.0;
		iter=1;
		obj=0.0;
		obj_old=0.0;
		Pi_nk=Pi(_,i);
		obj_old=Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha, omega, mu);
		while(err>1e-4 && iter<50){
			//Golden section search for each parameter of interest
			omega=Min_Omega_Ares_mix_SNP(0.01, 100.0, 1e-4, Pi_nk, B, R, purity, CNA, CNB, alpha, mu);
			mu=Min_Mu_Ares_mix_SNP(0.01, 100.0, 1e-4, Pi_nk, B, R, purity, CNA, CNB, alpha, omega);
			alpha=Min_Alpha_Ares_mix_SNP(0.01, 100.0, 1e-4, Pi_nk, B, R, purity, CNA, CNB, omega, mu);
			obj=Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha, omega, mu);
			iter++;
			err=fabs((obj_old-obj)/obj_old);
			obj_old=obj;
		}
		Alpha[i]=alpha;
		Omega[i]=omega;
		Mu[i]=mu;
	}
	return(List::create(Named("Alpha")=Alpha,Named("Omega")=Omega,Named("Mu")=Mu));
}

// [[Rcpp::export]]
List Ares_mix_SNP_Gradient(NumericMatrix Pi, NumericVector B, NumericVector R, NumericVector purity, NumericVector CNA, NumericVector CNB, int K_cluster){
	NumericVector Alpha(K_cluster);
	NumericVector Omega(K_cluster);
	NumericVector Mu(K_cluster);
	NumericVector Derivative(3);
	NumericVector parameter(3);
	double lrate0=0.00001;
	double lrate;
	double decay_rate=0.5;
	int n=R.size();
	NumericVector Pi_nk(n);
	NumericVector prob = B/R;
	NumericVector theta(n);
	NumericVector dTheta(n);
	NumericVector D1(n);
	NumericVector dThetadalpha(n);
	NumericVector dThetadomega(n);
	NumericVector dThetadmu(n);
	double alpha=1.0;
	double omega=1.0;
	double mu=1.0;
	double alpha2=1.0;
	double omega2=1.0;
	double mu2=1.0;
	double err=1.0;
	int iter=1;
	double obj=0.0;
	double obj_old=0.0;
	
	for(int i =0; i<K_cluster; i++){
		alpha2=1.0;
		omega2=1.0;
		mu2=1.0;
		err=1.0;
		iter=1;
		obj=0.0;
		obj_old=0.0;
		Pi_nk=Pi(_,i);
		
		alpha=1.0;
		omega=1.0;
		mu=1.0;
		
		//Start Gradient descent
		err=1.0;
		iter=1;
		obj_old=Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha2, omega2, mu2);
		parameter[0]=alpha;
		parameter[1]=omega;
		parameter[2]=mu;
		while(err>1e-5 && iter<100){
			lrate=lrate0/(1+decay_rate*iter);
			theta = 1/ (1 + ( 1-purity + purity*CNA*mu2)/((1-purity)*alpha2 + purity*CNB*omega2*mu2));
			dTheta=2*Pi_nk*R/prob/(1-prob)*(theta-prob);
			
			D1 = ((1-purity)*alpha2 + 1-purity + purity*mu2*(CNB*omega2+CNA))*((1-purity)*alpha2 + 1-purity + purity*mu2*(CNB*omega2+CNA));
			dThetadalpha = 2*alpha*(1- purity)*(CNA*purity*mu2 - purity + 1)/D1;
			dThetadomega = 2*CNB*mu2*omega*purity*(CNA*purity*mu2 - purity + 1)/D1;
			dThetadmu = 2*mu*purity*(CNA*alpha2 - CNB*omega2)*(purity - 1)/D1;
			Derivative[0]=sum(dThetadalpha*dTheta);
			Derivative[1]=sum(dThetadomega*dTheta);
			Derivative[2]=sum(dThetadmu*dTheta);
			
			parameter = parameter - lrate*Derivative;
			alpha=parameter[0];
			omega=parameter[1];
			mu=parameter[2];
			
			alpha2=alpha*alpha;
			omega2=omega*omega;
			mu2=mu*mu;
			obj=Ares_mix_SNP_LogLikelihood(Pi_nk, B, R, purity, CNA, CNB, alpha2, omega2, mu2);
			//printf("obj: %f in iter: %d \n", obj, iter);
			err=fabs((obj_old-obj)/obj_old);
			obj_old=obj;
			iter++;
		}
		Alpha[i]=alpha*alpha;
		Omega[i]=omega*omega;
		Mu[i]=mu*mu;
	}
	return(List::create(Named("Alpha")=Alpha,Named("Omega")=Omega,Named("Mu")=Mu));
}


/* Combined model with shared mu */

// [[Rcpp::export]]
double Ares_mix_Combined_LogLikelihood(NumericVector Pi_nk_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericVector Pi_nk_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr, double alpha, double omega_SNP, double omega_SNV, double mu){
	/* Objective function*/
	NumericVector theta_SNP = 1/ (1 + ( 1-purity_SNP + purity_SNP*CNA*mu)/((1-purity_SNP)*alpha + purity_SNP*CNB*omega_SNP*mu)) ;
	NumericVector prob_SNP = B/R_SNP;
	NumericVector res_SNP = Pi_nk_SNP*R_SNP*(theta_SNP-prob_SNP)*(theta_SNP-prob_SNP)/(prob_SNP*(1.0-prob_SNP));
	NumericVector theta_SNV = purity_SNV*CNv*omega_SNV/(2.0*(1.0-purity_SNV)/mu+purity_SNV*(CNv*omega_SNV+CNr));
	NumericVector prob_SNV = V/R_SNV;
	NumericVector res_SNV = Pi_nk_SNV*R_SNV*(theta_SNV-prob_SNV)*(theta_SNV-prob_SNV)/(prob_SNV*(1.0-prob_SNV));
	return(sum(res_SNP)+sum(res_SNV));
}

// [[Rcpp::export]]
double Min_Mu_Ares_mix_Combined(double lower, double upper, double tol, NumericVector Pi_nk_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericVector Pi_nk_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr, double alpha, double omega_SNP, double omega_SNV) {
   /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = std::numeric_limits<double>::epsilon();
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = lower;
    b = upper;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
	
	fx = Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
        xm = (a + b) * .5;
        tol1 = eps * fabs(x) + tol3;
        t2 = tol1 * 2.;
        
        /* check stopping criterion */
        
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
        p = 0.;
        q = 0.;
        r = 0.;
        if (fabs(e) > tol1) { /* fit parabola */
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
        }
        
        if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
        }
        else { /* a parabolic-interpolation step */
            
            d = p / q;
            u = x + d;
            
            /* f must not be evaluated too close to ax or bx */
            
            if (u - a < t2 || b - u < t2) {
                d = tol1;
                if (x >= xm) d = -d;
            }
        }
        
        /* f must not be evaluated too close to x */
        
        if (fabs(d) >= tol1)
            u = x + d;
        else if (d > 0.)
            u = x + tol1;
        else
            u = x - tol1;
        
		fu = Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, u);
        
        /*  update  a, b, v, w, and x */
        
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    /* end of main loop */
    
    return x;
}

// [[Rcpp::export]]
List Solve_Ares_mix_Combined_Mu(NumericVector Pi_nk_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericVector Pi_nk_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr){
	double alpha=1.0;
	double omega_SNP=1.0;
	double omega_SNV=1.0;
	double mu=1.0;
	double err=1.0;
	int iter=1;
	double obj=0.0;
	double obj_old=0.0;
	while(err>1e-5 && iter<50){
		//Golden section search for each parameter of interest
		mu=Min_Mu_Ares_mix_Combined(0.01, 200.0, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV);
		obj=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
		iter++;
		err=fabs((obj_old-obj)/obj_old);
		obj_old=obj;
	}
	return(List::create(Named("Alpha")=alpha,Named("Omega_SNP")=omega_SNP,Named("Omega_SNV")=omega_SNV,Named("Mu")=mu));
}

// [[Rcpp::export]]
List Solve_Ares_mix_Combined(NumericVector Pi_nk_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericVector Pi_nk_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr){
	double alpha=1.0;
	double omega_SNP=1.0;
	double omega_SNV=1.0;
	double mu=1.0;
	double err=1.0;
	int iter=1;
	double obj=0.0;
	double obj_old=0.0;
	while(err>1e-5 && iter<50){
		//Golden section search for each parameter of interest
		omega_SNV = Min_Omega_Binomial_repara(0.01, 100, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, mu);
		alpha = Min_Alpha_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, omega_SNP, mu);
		omega_SNP = Min_Omega_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, mu);
		mu = Min_Mu_Ares_mix_Combined(0.01, 100.0, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV);
		obj=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
		iter++;
		err=fabs((obj_old-obj)/obj_old);
		obj_old=obj;
	}
	return(List::create(Named("Alpha")=alpha,Named("Omega_SNP")=omega_SNP,Named("Omega_SNV")=omega_SNV,Named("Mu")=mu));
}


/* Combined model */

// [[Rcpp::export]]
List Ares_mix_Combined_Mstep_Mu(NumericMatrix Pi_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericMatrix Pi_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr, int K_cluster){
	NumericVector Alpha(K_cluster);
	NumericVector Omega_SNP(K_cluster);
	NumericVector Omega_SNV(K_cluster);
	NumericVector Mu(K_cluster);
	int n_SNP=R_SNP.size();
	int n_SNV=R_SNV.size();
	NumericVector Pi_nk_SNP(n_SNP);
	NumericVector Pi_nk_SNV(n_SNV);
	double alpha=1.0;
	double omega_SNP=1.0;
	double omega_SNV=1.0;
	double mu=1.0;
	double err=1.0;
	int iter=1;
	double obj=0.0;
	double obj_old=0.0;
	for(int i =0; i<K_cluster; i++){
		alpha=1.0;
		omega_SNP=1.0;
		omega_SNV=1.0;
		mu=1.0;
		err=1.0;
		iter=1;
		obj=0.0;
		obj_old=0.0;
		Pi_nk_SNP=Pi_SNP(_,i);
		Pi_nk_SNV=Pi_SNV(_,i);
		obj_old=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
		while(err>1e-4 && iter<50){
			//Golden section search for each parameter of interest
			mu=Min_Mu_Ares_mix_Combined(0.01, 200.0, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV);
			obj=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
			iter++;
			err=fabs((obj_old-obj)/obj_old);
			obj_old=obj;
		}
		Alpha[i]=alpha;
		Omega_SNP[i]=omega_SNP;
		Omega_SNV[i]=omega_SNV;
		Mu[i]=mu;
	}
	return(List::create(Named("Alpha")=Alpha,Named("Omega_SNP")=Omega_SNP,Named("Omega_SNV")=Omega_SNV,Named("Mu")=Mu));
}


// [[Rcpp::export]]
List Ares_mix_Combined_Mstep(NumericMatrix Pi_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericMatrix Pi_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr, int K_cluster){
	NumericVector Alpha(K_cluster);
	NumericVector Omega_SNP(K_cluster);
	NumericVector Omega_SNV(K_cluster);
	NumericVector Mu(K_cluster);
	int n_SNP=R_SNP.size();
	int n_SNV=R_SNV.size();
	NumericVector Pi_nk_SNP(n_SNP);
	NumericVector Pi_nk_SNV(n_SNV);
	double alpha=1.0;
	double omega_SNP=1.0;
	double omega_SNV=1.0;
	double mu=1.0;
	double err=1.0;
	int iter=1;
	double obj=0.0;
	double obj_old=0.0;
	for(int i =0; i<K_cluster; i++){
		alpha=1.0;
		omega_SNP=1.0;
		omega_SNV=1.0;
		mu=1.0;
		err=1.0;
		iter=1;
		obj=0.0;
		obj_old=0.0;
		Pi_nk_SNP=Pi_SNP(_,i);
		Pi_nk_SNV=Pi_SNV(_,i);
		obj_old=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
		if(i==0){
			while(err>1e-4 && iter<50){
				//Golden section search for each parameter of interest
				mu=Min_Mu_Ares_mix_Combined(0.01, 100.0, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV);
				obj=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
				iter++;
				err=fabs((obj_old-obj)/obj_old);
				obj_old=obj;
			}
			Alpha[i]=alpha;
			Omega_SNP[i]=omega_SNP;
			Omega_SNV[i]=omega_SNV;
			Mu[i]=mu;
		}else{
			while(err>1e-4 && iter<50){
				//Golden section search for each parameter of interest
				omega_SNV = Min_Omega_Binomial_repara(0.01, 100, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, mu);
				alpha = Min_Alpha_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, omega_SNP, mu);
				omega_SNP = Min_Omega_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, mu);
				mu=Min_Mu_Ares_mix_Combined(0.01, 100.0, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV);
				obj=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
				iter++;
				err=fabs((obj_old-obj)/obj_old);
				obj_old=obj;
			}
			Alpha[i]=alpha;
			Omega_SNP[i]=omega_SNP;
			Omega_SNV[i]=omega_SNV;
			Mu[i]=mu;
		}
		
	}
	return(List::create(Named("Alpha")=Alpha,Named("Omega_SNP")=Omega_SNP,Named("Omega_SNV")=Omega_SNV,Named("Mu")=Mu));
}

// [[Rcpp::export]]
List Ares_mix_Combined_Mstep2(NumericMatrix Pi_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericMatrix Pi_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr, int K_cluster){
	NumericVector Alpha(K_cluster);
	NumericVector Omega_SNP(K_cluster);
	NumericVector Omega_SNV(K_cluster);
	NumericVector Mu(K_cluster);
	int n_SNP=R_SNP.size();
	int n_SNV=R_SNV.size();
	NumericVector Pi_nk_SNP(n_SNP);
	NumericVector Pi_nk_SNV(n_SNV);
	double alpha=1.0;
	double omega_SNP=1.0;
	double omega_SNV=1.0;
	double mu=1.0;
	double err=1.0;
	int iter=1;
	double obj=0.0;
	double obj_old=0.0;
	for(int i =0; i<K_cluster; i++){
		alpha=1.0;
		omega_SNP=1.0;
		omega_SNV=1.0;
		mu=1.0;
		err=1.0;
		iter=1;
		obj=0.0;
		obj_old=0.0;
		Pi_nk_SNP=Pi_SNP(_,i);
		Pi_nk_SNV=Pi_SNV(_,i);
		obj_old=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
		if(i==0){
			while(err>1e-4 && iter<50){
				//Golden section search for each parameter of interest
				omega_SNV = Min_Omega_Binomial_repara(0.01, 100, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, mu);
				mu=Min_Mu_Ares_mix_Combined(0.01, 100.0, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV);
				obj=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
				iter++;
				err=fabs((obj_old-obj)/obj_old);
				obj_old=obj;
			}
			Alpha[i]=alpha;
			Omega_SNP[i]=omega_SNP;
			Omega_SNV[i]=omega_SNV;
			Mu[i]=mu;
		}else if(i==1){
			while(err>1e-4 && iter<50){
				//Golden section search for each parameter of interest
				alpha = Min_Alpha_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, omega_SNP, mu);
				omega_SNP = Min_Omega_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, mu);
				mu=Min_Mu_Ares_mix_Combined(0.01, 100.0, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV);
				obj=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
				iter++;
				err=fabs((obj_old-obj)/obj_old);
				obj_old=obj;
			}
			Alpha[i]=alpha;
			Omega_SNP[i]=omega_SNP;
			Omega_SNV[i]=omega_SNV;
			Mu[i]=mu;
		}else{
			while(err>1e-4 && iter<50){
				//Golden section search for each parameter of interest
				omega_SNV = Min_Omega_Binomial_repara(0.01, 100, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, mu);
				alpha = Min_Alpha_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, omega_SNP, mu);
				omega_SNP = Min_Omega_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, mu);
				mu=Min_Mu_Ares_mix_Combined(0.01, 100.0, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV);
				obj=Ares_mix_Combined_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu);
				iter++;
				err=fabs((obj_old-obj)/obj_old);
				obj_old=obj;
			}
			Alpha[i]=alpha;
			Omega_SNP[i]=omega_SNP;
			Omega_SNV[i]=omega_SNV;
			Mu[i]=mu;
		}
		
	}
	return(List::create(Named("Alpha")=Alpha,Named("Omega_SNP")=Omega_SNP,Named("Omega_SNV")=Omega_SNV,Named("Mu")=Mu));
}


/* Combined model SNP and SNv independent*/

// [[Rcpp::export]]
List Solve_Ares_mix_Combined_Independent(NumericVector Pi_nk_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericVector Pi_nk_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr){
	double alpha=1.0;
	double omega_SNP=1.0;
	double omega_SNV=1.0;
	double mu_SNP=1.0;
	double mu_SNV=1.0;
	double err_SNP=1.0;
	double err_SNV=1.0;
	int iter_SNP=1;
	int iter_SNV=1;
	double obj_SNP=0.0;
	double obj_SNV=0.0;
	double obj_SNP_old=Ares_mix_SNP_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP, mu_SNP);
	double obj_SNV_old=Mixed_Binomial_LogLikelihood_repara(Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV, mu_SNV);
	while(err_SNP>1e-4 && iter_SNP<50){
		//Golden section search for each parameter of interest
		mu_SNP = Min_Mu_Ares_mix_SNP(0.01, 100.0, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP);
		alpha = Min_Alpha_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, omega_SNP, mu_SNV);
		omega_SNP = Min_Omega_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, mu_SNV);
		obj_SNP=Ares_mix_SNP_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP, mu_SNP);
		iter_SNP++;
		err_SNP=fabs((obj_SNP_old-obj_SNP)/obj_SNP_old);
		obj_SNP_old=obj_SNP;
	}
	while(err_SNV>1e-4 && iter_SNV<50){
		//Golden section search for each parameter of interest
		mu_SNV = Min_Mu_Binomial_repara(0.01, 100, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV);
		omega_SNV = Min_Omega_Binomial_repara(0.01, 100, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, mu_SNV);
		obj_SNV=Mixed_Binomial_LogLikelihood_repara(Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV, mu_SNV);
		iter_SNV++;
		err_SNV=fabs((obj_SNV_old-obj_SNV)/obj_SNV_old);
		obj_SNV_old=obj_SNV;
	}
	
	return(List::create(Named("Alpha")=alpha,Named("Omega_SNP")=omega_SNP,Named("Omega_SNV")=omega_SNV,Named("Mu_SNP")=mu_SNP,Named("Mu_SNV")=mu_SNV));
}


// [[Rcpp::export]]
List Ares_mix_Combined_Independent(NumericMatrix Pi_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericMatrix Pi_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr, int K_cluster){
	NumericVector Alpha(K_cluster);
	NumericVector Omega_SNP(K_cluster);
	NumericVector Omega_SNV(K_cluster);
	NumericVector Mu_SNP(K_cluster);
	NumericVector Mu_SNV(K_cluster);
	int n_SNP=R_SNP.size();
	int n_SNV=R_SNV.size();
	NumericVector Pi_nk_SNP(n_SNP);
	NumericVector Pi_nk_SNV(n_SNV);
	double alpha=1.0;
	double omega_SNP=1.0;
	double omega_SNV=1.0;
	double mu_SNP=1.0;
	double mu_SNV=1.0;
	double err_SNP=1.0;
	double err_SNV=1.0;
	int iter_SNP=1;
	int iter_SNV=1;
	double obj_SNP=0.0;
	double obj_SNV=0.0;
	double obj_SNP_old=0.0;
	double obj_SNV_old=0.0;
	for(int i =0; i<K_cluster; i++){
		alpha=1.0;
		omega_SNP=1.0;
		omega_SNV=1.0;
		mu_SNP=1.0;
		mu_SNV=1.0;
		err_SNP=1.0;
		err_SNV=1.0;
		iter_SNP=1;
		iter_SNV=1;
		obj_SNP=0.0;
		obj_SNV=0.0;
		Pi_nk_SNP=Pi_SNP(_,i);
		Pi_nk_SNV=Pi_SNV(_,i);
		obj_SNP_old=Ares_mix_SNP_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP, mu_SNP);
		obj_SNV_old=Mixed_Binomial_LogLikelihood_repara(Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV, mu_SNV);
		if(i==0){
			while(err_SNP>1e-4 && iter_SNP<50){
				//Golden section search for each parameter of interest
				mu_SNP=Min_Mu_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP);
				obj_SNP=Ares_mix_SNP_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP, mu_SNP);
				iter_SNP++;
				err_SNP=fabs((obj_SNP_old-obj_SNP)/obj_SNP_old);
				obj_SNP_old=obj_SNP;
			}
			while(err_SNV>1e-4 && iter_SNV<50){
				//Golden section search for each parameter of interest
				mu_SNV=Min_Mu_Binomial_repara(0.01, 100.0, 1e-6, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV);
				obj_SNV=Mixed_Binomial_LogLikelihood_repara(Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV, mu_SNV);
				iter_SNV++;
				err_SNV=fabs((obj_SNV_old-obj_SNV)/obj_SNV_old);
				obj_SNV_old=obj_SNV;
			}
			
		}else{
			while(err_SNP>1e-4 && iter_SNP<50){
				//Golden section search for each parameter of interest
				mu_SNP=Min_Mu_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP);
				omega_SNP=Min_Omega_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, mu_SNP);
				alpha=Min_Alpha_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, omega_SNP, mu_SNP);
				obj_SNP=Ares_mix_SNP_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP, mu_SNP);
				iter_SNP++;
				err_SNP=fabs((obj_SNP_old-obj_SNP)/obj_SNP_old);
				obj_SNP_old=obj_SNP;
			}
			while(err_SNV>1e-4 && iter_SNV<50){
				//Golden section search for each parameter of interest
				mu_SNV=Min_Mu_Binomial_repara(0.01, 100.0, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV);
				omega_SNV=Min_Omega_Binomial_repara(0.01, 100.0, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, mu_SNV);
				obj_SNV=Mixed_Binomial_LogLikelihood_repara(Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV, mu_SNV);
				iter_SNV++;
				err_SNV=fabs((obj_SNV_old-obj_SNV)/obj_SNV_old);
				obj_SNV_old=obj_SNV;
			}
		}
		Alpha[i]=alpha;
		Omega_SNP[i]=omega_SNP;
		Omega_SNV[i]=omega_SNV;
		Mu_SNP[i]=mu_SNP;
		Mu_SNV[i]=mu_SNV;
	}
	return(List::create(Named("Alpha")=Alpha,Named("Omega_SNP")=Omega_SNP,Named("Omega_SNV")=Omega_SNV,Named("Mu_SNP")=Mu_SNP,Named("Mu_SNV")=Mu_SNV));
}

// [[Rcpp::export]]
List Ares_mix_Joint(NumericMatrix Pi_SNP, NumericVector B, NumericVector R_SNP, NumericVector purity_SNP, NumericVector CNA, NumericVector CNB, NumericMatrix Pi_SNV, NumericVector V, NumericVector R_SNV, NumericVector purity_SNV, NumericVector CNv, NumericVector CNr, int K_cluster){
	NumericVector Alpha(K_cluster);
	NumericVector Omega_SNP(K_cluster);
	NumericVector Omega_SNV(K_cluster);
	NumericVector Mu_SNP(K_cluster);
	NumericVector Mu_SNV(K_cluster);
	int n_SNP=R_SNP.size();
	int n_SNV=R_SNV.size();
	NumericVector Pi_nk_SNP(n_SNP);
	NumericVector Pi_nk_SNV(n_SNV);
	double alpha=1.0;
	double omega_SNP=1.0;
	double omega_SNV=1.0;
	double mu_SNP=1.0;
	double mu_SNV=1.0;
	double err_SNP=1.0;
	double err_SNV=1.0;
	int iter_SNP=1;
	int iter_SNV=1;
	double obj_SNP=0.0;
	double obj_SNV=0.0;
	double obj_SNP_old=0.0;
	double obj_SNV_old=0.0;
	for(int i =0; i<K_cluster; i++){
		alpha=1.0;
		omega_SNP=1.0;
		omega_SNV=1.0;
		mu_SNP=1.0;
		mu_SNV=1.0;
		err_SNP=1.0;
		err_SNV=1.0;
		iter_SNP=1;
		iter_SNV=1;
		obj_SNP=0.0;
		obj_SNV=0.0;
		Pi_nk_SNP=Pi_SNP(_,i);
		Pi_nk_SNV=Pi_SNV(_,i);
		obj_SNP_old=Ares_mix_SNP_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP, mu_SNP);
		obj_SNV_old=Mixed_Binomial_LogLikelihood_repara(Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV, mu_SNV);
		while(err_SNP>1e-4 && iter_SNP<50){
			//Golden section search for each parameter of interest
			mu_SNP=Min_Mu_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP);
			omega_SNP=Min_Omega_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, mu_SNP);
			alpha=Min_Alpha_Ares_mix_SNP(0.01, 100, 1e-4, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, omega_SNP, mu_SNP);
			obj_SNP=Ares_mix_SNP_LogLikelihood(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, alpha, omega_SNP, mu_SNP);
			iter_SNP++;
			err_SNP=fabs((obj_SNP_old-obj_SNP)/obj_SNP_old);
			obj_SNP_old=obj_SNP;
		}
		while(err_SNV>1e-4 && iter_SNV<50){
			//Golden section search for each parameter of interest
			mu_SNV=Min_Mu_Binomial_repara(0.01, 100.0, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV);
			omega_SNV=Min_Omega_Binomial_repara(0.01, 100.0, 1e-4, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, mu_SNV);
			obj_SNV=Mixed_Binomial_LogLikelihood_repara(Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, omega_SNV, mu_SNV);
			iter_SNV++;
			err_SNV=fabs((obj_SNV_old-obj_SNV)/obj_SNV_old);
			obj_SNV_old=obj_SNV;
		}
		Alpha[i]=alpha;
		Omega_SNP[i]=omega_SNP;
		Omega_SNV[i]=omega_SNV;
		Mu_SNP[i]=mu_SNP;
		Mu_SNV[i]=mu_SNV;
	}
	return(List::create(Named("Alpha")=Alpha,Named("Omega_SNP")=Omega_SNP,Named("Omega_SNV")=Omega_SNV,Named("Mu_SNP")=Mu_SNP,Named("Mu_SNV")=Mu_SNV));
}
