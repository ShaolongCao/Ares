# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Mixed_Binomial_LogLikelihood_repara <- function(Pi_nk, V, R, purity, CNv, CNr, omega, mu) {
    .Call(`_Ares_Mixed_Binomial_LogLikelihood_repara`, Pi_nk, V, R, purity, CNv, CNr, omega, mu)
}

Min_Mu_Binomial_repara <- function(lower, upper, tol, Pi_nk, V, R, purity, CNv, CNr, omega) {
    .Call(`_Ares_Min_Mu_Binomial_repara`, lower, upper, tol, Pi_nk, V, R, purity, CNv, CNr, omega)
}

Min_Omega_Binomial_repara <- function(lower, upper, tol, Pi_nk, V, R, purity, CNv, CNr, mu) {
    .Call(`_Ares_Min_Omega_Binomial_repara`, lower, upper, tol, Pi_nk, V, R, purity, CNv, CNr, mu)
}

Solve_mixed_Binomial_repara <- function(Pi_nk, V, R, purity, CNv, CNr) {
    .Call(`_Ares_Solve_mixed_Binomial_repara`, Pi_nk, V, R, purity, CNv, CNr)
}

Mixed_Binomial_LogLikelihood_newton <- function(Pi_nk, V, R, purity, CNv, CNr, omega, mu) {
    .Call(`_Ares_Mixed_Binomial_LogLikelihood_newton`, Pi_nk, V, R, purity, CNv, CNr, omega, mu)
}

Mixed_Binomial_LogLikelihood <- function(Pi_nk, V, R, purity, CNv, CNr, omega, mu) {
    .Call(`_Ares_Mixed_Binomial_LogLikelihood`, Pi_nk, V, R, purity, CNv, CNr, omega, mu)
}

Min_Mu_Binomial <- function(lower, upper, tol, Pi_nk, V, R, purity, CNv, CNr, omega) {
    .Call(`_Ares_Min_Mu_Binomial`, lower, upper, tol, Pi_nk, V, R, purity, CNv, CNr, omega)
}

Min_Omega_Binomial <- function(lower, upper, tol, Pi_nk, V, R, purity, CNv, CNr, mu) {
    .Call(`_Ares_Min_Omega_Binomial`, lower, upper, tol, Pi_nk, V, R, purity, CNv, CNr, mu)
}

Solve_mixed_Binomial_Goldensection <- function(Pi_nk, V, R, purity, CNv, CNr) {
    .Call(`_Ares_Solve_mixed_Binomial_Goldensection`, Pi_nk, V, R, purity, CNv, CNr)
}

Solve_mixed_Binomial_newton_Goldensection <- function(Pi_nk, V, R, purity, CNv, CNr) {
    .Call(`_Ares_Solve_mixed_Binomial_newton_Goldensection`, Pi_nk, V, R, purity, CNv, CNr)
}

Ares_mix_SNP_LogLikelihood <- function(Pi_nk, B, R, purity, CNA, CNB, alpha, omega, mu) {
    .Call(`_Ares_Ares_mix_SNP_LogLikelihood`, Pi_nk, B, R, purity, CNA, CNB, alpha, omega, mu)
}

Min_Alpha_Ares_mix_SNP <- function(lower, upper, tol, Pi_nk, B, R, purity, CNA, CNB, omega, mu) {
    .Call(`_Ares_Min_Alpha_Ares_mix_SNP`, lower, upper, tol, Pi_nk, B, R, purity, CNA, CNB, omega, mu)
}

Min_Omega_Ares_mix_SNP <- function(lower, upper, tol, Pi_nk, B, R, purity, CNA, CNB, alpha, mu) {
    .Call(`_Ares_Min_Omega_Ares_mix_SNP`, lower, upper, tol, Pi_nk, B, R, purity, CNA, CNB, alpha, mu)
}

Min_Mu_Ares_mix_SNP <- function(lower, upper, tol, Pi_nk, B, R, purity, CNA, CNB, alpha, omega) {
    .Call(`_Ares_Min_Mu_Ares_mix_SNP`, lower, upper, tol, Pi_nk, B, R, purity, CNA, CNB, alpha, omega)
}

Solve_Ares_mix_SNP <- function(Pi_nk, B, R, purity, CNA, CNB) {
    .Call(`_Ares_Solve_Ares_mix_SNP`, Pi_nk, B, R, purity, CNA, CNB)
}

Ares_mix_SNP_Mstep <- function(Pi, B, R, purity, CNA, CNB, K_cluster) {
    .Call(`_Ares_Ares_mix_SNP_Mstep`, Pi, B, R, purity, CNA, CNB, K_cluster)
}

Ares_mix_SNP_Gradient <- function(Pi, B, R, purity, CNA, CNB, K_cluster) {
    .Call(`_Ares_Ares_mix_SNP_Gradient`, Pi, B, R, purity, CNA, CNB, K_cluster)
}

Ares_mix_Combined_LogLikelihood <- function(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu) {
    .Call(`_Ares_Ares_mix_Combined_LogLikelihood`, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV, mu)
}

Min_Mu_Ares_mix_Combined <- function(lower, upper, tol, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV) {
    .Call(`_Ares_Min_Mu_Ares_mix_Combined`, lower, upper, tol, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr, alpha, omega_SNP, omega_SNV)
}

Solve_Ares_mix_Combined_Mu <- function(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr) {
    .Call(`_Ares_Solve_Ares_mix_Combined_Mu`, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr)
}

Solve_Ares_mix_Combined <- function(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr) {
    .Call(`_Ares_Solve_Ares_mix_Combined`, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr)
}

Ares_mix_Combined_Mstep_Mu <- function(Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster) {
    .Call(`_Ares_Ares_mix_Combined_Mstep_Mu`, Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster)
}

Ares_mix_Combined_Mstep <- function(Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster) {
    .Call(`_Ares_Ares_mix_Combined_Mstep`, Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster)
}

Ares_mix_Combined_Mstep2 <- function(Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster) {
    .Call(`_Ares_Ares_mix_Combined_Mstep2`, Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster)
}

Solve_Ares_mix_Combined_Independent <- function(Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr) {
    .Call(`_Ares_Solve_Ares_mix_Combined_Independent`, Pi_nk_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_nk_SNV, V, R_SNV, purity_SNV, CNv, CNr)
}

Ares_mix_Combined_Independent <- function(Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster) {
    .Call(`_Ares_Ares_mix_Combined_Independent`, Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster)
}

Ares_mix_Joint <- function(Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster) {
    .Call(`_Ares_Ares_mix_Joint`, Pi_SNP, B, R_SNP, purity_SNP, CNA, CNB, Pi_SNV, V, R_SNV, purity_SNV, CNv, CNr, K_cluster)
}

