/* dqag.c -- modified version of QUADPACK routine DQAG.
 * (C)1999, C. Bond. All right reserved.
 *
 * There are no changes to the basic computational method. Only
 * the temporary storage strategy is changed to utilize the
 * local stack at the appropriate level. This reduces the
 * need for memory allocation of arrays at higher levels and
 * the resulting passing of memory pointers down the line.
 *
 */
#include "cquadpack.h"

/* DQAG - Approximation to definite integral. (From QUADPACK)
 *
 *  Calls DQAGE with appropriate parameters assigned.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 *
 */
double ldqag(dq_function_type f,double a,double b,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
    double result;
    int last;

    double area,area1,area2,area12,a1,a2,b1,b2,c,defabs;
    double defab1,defab2,errbnd,errmax,error1,error2;
    double erro12,errsum,resabs;
    double alist[LIMIT],blist[LIMIT],rlist[LIMIT],elist[LIMIT];
    int iroff1,iroff2,k,keyf,maxerr,nrmax,iord[LIMIT],limit;

    limit = LIMIT - 1;
    *ier = 0;
    *neval = 0;
    last = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    defabs = 0.0;
    resabs = 0.0;
    if ((epsabs < 0.0) && (epsrel < 0.0))
        *ier = 6;
    if (*ier == 6) return result;

    const double lepsabs = log(epsabs);
    const double lepsrel = log(epsrel);

/* First approximation to the integral. */
    *neval = 0;
    result = LG_K21(f,a,b,abserr,&defabs,&resabs, user_data);
    last = 0;
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;

/* Test on accuracy. */
    errbnd = GSL_MAX_DBL(lepsabs, lepsrel + fabs(result));
    if ((*abserr <= log(50 * epmach) + defabs) && (*abserr > errbnd))
        *ier = 2;
    if (limit == 0) *ier = 1;
    if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) ||
        (gsl_isinf(*abserr) == -1)) goto _60;

/* Initialization. */
    errmax = *abserr;
    maxerr = 0;
    area = result;
    errsum = *abserr;
    nrmax = 0;
    iroff1 = 0;
    iroff2 = 0;

/* Main Loop. */
    for (last = 1; last <= limit; (last)++) {
/* Bisect the subinterval with the largest error estimate. */
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        area1 = LG_K21(f,a1,b1,&error1,&resabs,&defab1, user_data);
        area2 = LG_K21(f,a2,b2,&error2,&resabs,&defab2, user_data);

/* Improve previous approximations to integral and error,
        and test for accuracy. */
        (*neval) += 1;
        area12 = logsumexp(area1, area2);
        erro12 = logsumexp(error1, error2);
        errsum = logsubexp(logsumexp(errsum, erro12), errmax);
        area = logsubexp(logsumexp(area, area12), rlist[maxerr]);
        if ((defab1 != error1) && (defab2 != error2)) {
            double delta = LOGDIFF(rlist[maxerr], area12);
            if ((delta <= log(1.0e-5) + fabs(area12)) &&
                (erro12 >= log(0.99) + errmax))
                    iroff1++;
            if ((last > 9) && (erro12 > errmax))
                iroff2++;
        }
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = GSL_MAX_DBL(lepsabs,lepsrel + fabs(area));
        if (errsum > errbnd)  {

/* Test for roundoff error and eventually set error flag. */
            if ((iroff1 > 6) || (iroff2 > 20))
                *ier = 2;

/* Set error flag in the case that the number of subintervals
    equals the limit. */
            if (last == limit)
                *ier = 1;

/* Set error flag in the case of bad integrand behavior at a
    point of the integration range. */
            if (max(fabs(a1),fabs(b2)) <= (1.0 + c * 1000.0 * epmach) *
                (fabs(a2)+1.0e4 * uflow))
            *ier = 3;
        }
/* Append the newly-created intervals to the list. */

        if (error2 <= error1) {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = error1;
            elist[last] = error2;
        }
        else {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = area2;
            rlist[last] = area1;
            elist[maxerr] = error2;
            elist[last] = error1;
        }

/* Call DQSORT to maintain the descending ordering in the list of
    error estimates and select the subinterval with the
    largest error estimate (to be bisected next). */

        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);
        if ((*ier != 0) || (errsum <= errbnd)) break;
    }

/* Compute final result. */

    result = -INFINITY;
    for (k = 0; k <= last; k++) {
        result = logsumexp(result, rlist[k]);
    }
    *abserr = errsum;
_60:
    *neval = 21 * (2 * (*neval) + 1);
    return result;
}
