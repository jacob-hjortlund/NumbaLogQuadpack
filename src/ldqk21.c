#include "cquadpack.h"

double LG_K21(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
    static double XGK21[11] = {
        0.99565716302580808074,
        0.97390652851717172008,
        0.93015749135570822600,
        0.86506336668898451073,
        0.78081772658641689706,
        0.67940956829902440623,
        0.56275713466860468334,
        0.43339539412924719080,
        0.29439286270146019813,
        0.14887433898163121088,
        0.00000000000000000000};
    static double WGK21[11] = {
        0.01169463886737187428,
        0.03255816230796472748,
        0.05475589657435199603,
        0.07503967481091995277,
        0.09312545458369760554,
        0.10938715880229764190,
        0.12349197626206585108,
        0.13470921731147332593,
        0.14277593857706008080,
        0.14773910490133849137,
        0.14944555400291690566};
    static double WG10[5] = {
        0.06667134430868813759,
        0.14945134915058059315,
        0.21908636251598204400,
        0.26926671930999635509,
        0.29552422471475287017};
    double fv1[10],fv2[10];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = -INFINITY;
    fc=(*f)(centr, user_data);
    resk = fc + log(WGK21[10]);
    *resabs = resk;
    for (j = 0; j < 5; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK21[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = logsumexp(fval1, fval2);
        resg = logsumexp(resg, log(WG10[j]) + fsum);
        resk = logsumexp(resk, log(WGK21[jtw]) + fsum);
        *resabs = logsumexp(*resabs, log(WGK21[jtw]) + fsum);
    }
    for (j = 0; j < 5; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK21[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = logsumexp(fval1, fval2);
        resk = logsumexp(resk, log(WGK21[jtwm1]) + fsum);
        *resabs = logsumexp(*resabs, log(WGK21[jtwm1]) + fsum);
    }
    reskh = resk + log(0.5);
    *resasc = log(WGK21[10]) + LOGDIFF(fc, reskh);
    for (j = 0; j < 10; j++ )
        *resasc = logsumexp(*resasc, logsumexp(log(WGK21[j]) + LOGDIFF(fv1[j], reskh),
            LOGDIFF(fv2[j], reskh)));
    result = resk + log(hlgth);
    *resabs = (*resabs) + log(dhlgth);
    *resasc = (*resasc) + log(dhlgth);
    *abserr = LOGDIFF(resk, resg) + log(dhlgth);
    *abserr = rescale_error(*abserr, *resabs , *resasc);
    return result;
}
