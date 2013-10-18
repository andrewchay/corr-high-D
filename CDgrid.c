#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#define EPSILON (1e-9)
static double *vector(int n)
{
    double *v;
    v = Calloc(n,double);
    return(v);
}
static int *vector_int(int n)
{
    int *v;
    v = Calloc(n,int);
    return(v);
}
static void free_vector(double *v)
{
  Free(v);
}
static void free_vector_int(int *v)
{
    Free(v);
}
static int sgn(double z)
{
    if (z > 0) return(1);
    else if (z < 0) return(-1);
    else return(0);
}
static double positive(double z)
{
    if (z > 0) return(z); else return(0);
}
static double soft(double lambda, double z)
{
    if (z > lambda) return(z - lambda);
    else
    if (z < -lambda) return(z + lambda);
    else return(0);
}
static double fmcp(double lambda, double kappa, double z)
{
    if (fabs(z) > 1 / kappa * lambda) return(z); else
    {
        if (z > lambda) return((z - lambda) / (1 - kappa));
        else if (z < -lambda) return ((z + lambda) / (1 - kappa));
        else return(0);
    }
}
static double fscad(double lambda, double kappa, double z)
{
    if (fabs(z) > 1 / kappa * lambda) return(z); else
    {
        if (fabs(z) > 2 * lambda) return(((1 / kappa - 1) * z - sgn(z) * lambda / kappa) / (1 / kappa - 2));
        else return(sgn(z) * positive(fabs(z) - lambda));
    }
}

static int equalZero(double num)
{
    if (fabs(num) < EPSILON) return(1); else return(0);
}

static int checkConvergence(double *beta, double *beta_old, double eps,
                            int len)
{
    int i, converged = 1;
    for (i = 0; i <= len; i++)
    {
        if (!equalZero(beta[i]) && !equalZero(beta_old[i]))
        {
            if (fabs((beta[i]-beta_old[i]) / beta_old[i]) > eps)
            {
                converged = 0;
                break;
            }
        }
        else if (!equalZero(beta[i]) && equalZero(beta_old[i]))
        {converged = 0;break;}
        else if (equalZero(beta[i]) && !equalZero(beta_old[i]))
        {converged = 0;break;}
    }
    return(converged);
}

static void lasso(double *x, double *y, double *init, double lambda,
                  double eps, int max, int n,
                  int p, int *iter)
{
    double *beta_old, *beta, *u;
    double  z;
    int i, j, itetime = 0, flag = 1;
    beta = vector(p);
    beta_old = vector(p);
    for (i = 0; i < p; i++)
    {
        beta_old[i] = init[i];
        beta[i] = init[i];
    }
    u = vector(n);
    for (i = 0; i < n; i++)
        {u[i] = y[i];
        for (j = 0; j < p; j++)
            u[i] = u[i] - x[i + j * n] * beta_old[j];
        }
    while (flag && itetime <= max)
    {
        if (itetime > 0) for (i = 0; i < p; i++) beta_old[i] = beta[i];
        for (i = 0; i < p; i++)
        {
            z = 0;
            for (j = 0; j < n; j++) z = x[i * n + j] * u[j] + z;
            z = z + beta_old[i];
            beta[i] = soft(lambda, z);
            if (!equalZero(beta[i] - beta_old[i]))
                for (j = 0; j < n; j++)
                    u[j] = u[j] - (beta[i] - beta_old[i]) * x[i * n + j];
        }
        itetime++;
        flag = !checkConvergence(beta, beta_old, eps, p);
    }


    for (i = 0; i < p; i++) init[i] = beta[i];
    iter[0] = itetime;
    free_vector(beta);
    free_vector(beta_old);
    free_vector(u);
}

static void scaled_lasso(double *x, double *y, double *init, double *sigma,
                         double lambda, double eps, int max, int n,
                         int p, int *iter)
{
    double *beta_old, *u;
    double z, sigma_tmp, l;
    int i, j, itetime = 0, flag = 1;
    beta_old = vector(p);
    u = vector(n);
    while (flag && itetime <= max)
    {
        for (i = 0; i < p; i++) beta_old[i] = init[i];
        if (itetime == 0) l = lambda; else l = lambda * sigma[0];
        lasso(x, y, init, l, eps, max, n, p, iter);
        z = 0;
        for (i = 0; i < n; i++)
        {
            u[i] = y[i];
            for (j = 0; j < p; j++)
                u[i] = u[i] - x[i + j * n] * init[j];
            z = z + pow(u[i], 2);
        }
        sigma_tmp = sqrt(z / n);
        itetime++;
        flag = !(checkConvergence(init, beta_old, eps, p) &&
                 fabs((sigma[0]-sigma_tmp) / sigma[0]) <= eps);
        sigma[0] = sigma_tmp;
    }
    iter[0] = itetime;
    free_vector(beta_old);
    free_vector(u);
}

SEXP scaled(SEXP Rx, SEXP Ry, SEXP Rlambda, SEXP Reps, SEXP Rmax, SEXP Rpenalty)
{
    const char *penalty;
    int  i, j, k, n, p, max, iter[1];
    double eps, sigma[1];
    double lambda, *y, *x, *init;
    void lasso(), free_vec(), scaled_lasso();
    SEXP betahat, sigmahat, Riter, return_list;
    Rx = coerceVector(Rx, REALSXP);
    Ry = coerceVector(Ry, REALSXP);
    Rlambda = coerceVector(Rlambda, REALSXP);
    Reps = coerceVector(Reps, REALSXP);
    Rmax = coerceVector(Rmax, INTSXP);
    p = Rf_ncols(Rx);
    n = Rf_nrows(Rx);
    lambda = REAL(Rlambda)[0];
    penalty = CHAR(STRING_ELT(Rpenalty, 0));
    eps = REAL(Reps)[0];
    max = INTEGER(Rmax)[0];
    x = REAL(Rx);
    y = REAL(Ry);
    init = vector(p);
    PROTECT(betahat = Rf_allocVector(REALSXP, p));
    PROTECT(sigmahat = Rf_allocVector(REALSXP, 1));
    PROTECT(Riter = Rf_allocVector(INTSXP, 1));
    PROTECT(return_list = Rf_allocVector(VECSXP, 4));
    for (i = 0;i < p; i++) init[i] = 0;
    if (strcmp(penalty, "LASSO") == 0) {
        lasso(x, y, init, lambda, eps, max, n, p, iter);
        for (j = 0; j < p; j++) REAL(betahat)[j] = init[j];
        REAL(sigmahat)[0] = 0;
        INTEGER(Riter)[0] = iter[0];}
    if (strcmp(penalty, "SCALED") == 0) {
        scaled_lasso(x, y, init, sigma, lambda, eps, 500, n, p, iter);
        for (j = 0; j < p; j++) REAL(betahat)[j] = init[j];
        REAL(sigmahat)[0] = sigma[0];
        INTEGER(Riter)[0] = iter[0];}
    SET_VECTOR_ELT(return_list, 0, betahat);
    SET_VECTOR_ELT(return_list, 1, sigmahat);
    SET_VECTOR_ELT(return_list, 2, Rlambda);
    SET_VECTOR_ELT(return_list, 3, Riter);
    UNPROTECT(4);
    free_vector(init);
    return(return_list);
}

static R_CMethodDef DotCEntries[] = {
  {"scaled", (DL_FUNC) &scaled, 9},
  {NULL}
};

void R_init_pareto(DllInfo *info)
{
  R_registerRoutines(info, DotCEntries, NULL, NULL, NULL);
}
