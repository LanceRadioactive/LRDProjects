#include "proto.h"
int main()
{
    data d;
    printf("Please provide all the starting data - ellipsoid coefficients (A, B, C), line vector (M, N, L) and starting point (Xo, Yo, Zo).\n");

    do{
        d.ae = proof('A');
        if(d.ae==0) printf("Invalid value!\n");
    }while(d.ae==0);
    do{
        d.be = proof('B');
        if(d.be==0) printf("Invalid value!\n");
    }while(d.be==0);
    do{
        d.ce = proof('C');
        if(d.ce==0) printf("Invalid value!\n");
    }while(d.ce==0);
    do{
        d.ml = proof('M');
        d.nl = proof('N');
        d.ll = proof('L');
    if((d.ml==0) && (d.nl==0) && (d.ll==0)) printf("Please insert a non-zero vector!\n");
    }while((d.ml==0) && (d.nl==0) && (d.ll==0));
    d.xl = proof('X');
    d.yl = proof('Y');
    d.zl = proof('Z');

    printf("Performing intersection check.\n");
    float t = (d.xl*d.ml+d.yl*d.nl+d.zl*d.ll)/(d.ae*d.ae+d.be*d.be+d.ce*d.ce);
    float x = d.xl+d.ml*t;
    float y = d.yl+d.nl*t;
    float z = d.zl+d.ll*t;
    if(x*x/(d.ae*d.ae)+y*y/(d.be*d.be)+z*z/(d.ce*d.ce)<=1) {printf("Intersection detected!"); exit(0);};

    vect result = {0, 0, 0}; vect prevresult;
    printf("Data received. Please provide a margin of gradient error. Lower value means more precision but longer calculation time.\n");
    float err;
    do {
        err = proof('E');
        if (err<=0) printf("Margin must be above 0!\n");
    }while(err<=0);
    printf("Please input the step divisor. Lower value means more precision but longer calculation time. (common choice: 2):\n");
    float k;
    do {
        k = proof('K');
        if (k<=1) printf("Divisor must be above 1!\n");
    }while (k<=1);
    int i = 0; int f; float dist, prevdist, n;
    vect grad, prevgrad; pts points, prevpoints;
    printf("Initializing. Initial assumption: t=0, Theta=0, Alpha=0\n");
    calccoord(result, d, &points);
    calcgrad(result, d, points, &grad);
    dist = calcdist(points); t = 0;

    while(grad.a*grad.a+grad.f*grad.f+grad.t*grad.t>err*err) {
        i++; n=1; f=0;
        printf("Step %i: dT = %f, dTheta = %f, dAlpha = %f\n", i, grad.t, grad.f, grad.a);

        prevdist = dist; prevgrad = grad; prevpoints = points; prevresult = result;

        do {
            result.a = result.a - grad.a/n;
            result.f = result.f - grad.f/n;
            result.t = result.t - grad.t/n;

            calccoord(result, d, &points);
            calcgrad(result, d, points, &grad);
            dist = calcdist(points);
            if (dist>prevdist) {
                    n = n*k;
                    //n++;
                    grad = prevgrad;
                    points = prevpoints;
                    result = prevresult;
                    f++;
                    }
        } while (dist>prevdist);

        printf("Step divisor power used: %i\n", f);
        printf("Current assumption: t = %f, Theta = %f, Alpha = %f\n", result.t, result.f, result.a);
        printf("Current distance: %f\n", dist);
        if (dist==prevdist) t++; else t=0;
        if(t==5){
            printf("Recurrence detected, quitting. Warning: gradient is above error margin!\n");
            break;
        }
        //system("PAUSE");
    }

    printf("Calculations finished after %i steps. Result: %f; dT = %f, dTheta = %f, dAlpha = %f\nGradient length = %f\n", i, dist, grad.t, grad.f, grad.a, sqrt(grad.t*grad.t+grad.f*grad.f+grad.a*grad.a));
    system("PAUSE");

}
