/**************************************************************************
 * PCA3D.c
 * 
 * Eigen decomposition of 3D symmetric matrices
 *
 * Author : Alexis Mignon (c)
 * e-mail : alexis.mignon@gmail.com
 *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
/**
 * Solution of the cardan form of 3rd degree equation
 *
 */
int solve_cardan(double p, double q, double* z0, double* z1, double* z2){
    if (p == 0.0) {
        return -1;
    }
    else {
        double q2 = q*q;
        double p3 = p*p*p;
        double D = q2 + 4.0*p3/27.0;
        
        if (fabs(D) < 1e-10) D = 0.0;
        if (D>0.0) {
            return -1;
        }
		else if  (D < 0.0) {
            double acosq = acos(-q/2.0 * sqrt(27/-p3));
            double _2sqrt = 2 * sqrt(-p/3.0) ;
			*z0 =  _2sqrt * cos(1.0/3.0 * acosq );
			*z1 =  _2sqrt * cos(1.0/3.0 * acosq + 2*M_PI/3);
			*z2 =  _2sqrt * cos(1.0/3.0 * acosq + 4*M_PI/3);
        }
		else {
			*z0 = 3*q/p;
			*z1 = *z2 = -3.0*q/(2.0*p);
        }
    }
	return 0;
}

/**
 * Solution of 3rd degree equation
 *
 */
int solve_3rd(double a, double b, double c, double d, double* x0, double* x1, double* x2){
    int res;
    if (a==0.0 && b == 0.0 && c == 0.0 && d == 0.0) return -1;
	if (b == 0.0) {
        res = solve_cardan(c/a,d/a,x0,x1,x2);
        return res;
    }
    double p, q, s;
    double a2 = a*a;
    double b2 = b*b;

    
	p = -b2/(3*a2) + c/a;
	q = b/(27*a)*(2*b2/a2 - 9.0*c/a) + d/a;
	res = solve_cardan(p,q,x0,x1,x2);
	s = -b/(3*a);
    *x0 += s;
    *x1 += s;
    *x2 += s;
    return res;
}

/**
 * Find an eigen vector of matrix A knowing an eigenvalue
 *
 */
void find_eigen_vector(double* A, double eigen_value, double*u){
	double a11 = A[0];
    double a12 = A[1];
    double a13 = A[2];
	double a22 = A[4];
    double a23 = A[5];
    
	u[0] = a12 * a23 - a13 * (a22 - eigen_value);
    u[1] = a12 * a13 - a23 * (a11 - eigen_value);
    u[2] = (a11 - eigen_value) * (a22 - eigen_value) - a12*a12;
    double norm = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
	u[0] /= norm;
	u[1] /= norm;
	u[2] /= norm;
}


void cross_product(double* u, double* v, double* w){
	w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
}

/**
 * Performs eigen decomposition of a 3x3 symmetric matrix
 */
 
int eigen_sym_3D(double* A, double* w, double* U){
	double a11 = A[0];
    double a12 = A[1];
    double a13 = A[2];
	double a22 = A[4];
    double a23 = A[5];
    double a33 = A[8];
    
	double a = -1.0;
	double b = a11+a22+a33;
	double c = (-a22-a11)*a33+a23*a23-a11*a22+a13*a13+a12*a12;
	double d = a11*a22*a33-a12*a12*a33-a11*a23*a23+2*a12*a13*a23-a13*a13*a22;
	
    double x0,x1,x2;
    int res;
    
    res = solve_3rd(a,b,c,d,&x0,&x1,&x2);
    if (res != 0) return -1;

    double tmp;
    if (x0<x1) {
        tmp = x0;
        x0 = x1;
        x1 = tmp;
    }

    if (x0<x2) {
        tmp = x0;
        x0 = x2;
        x2 = tmp;
    }
    
    if (x1<x2) {
        tmp = x1;
        x1 = x2;
        x2 = tmp;
    }
    w[0] = x0;
    w[1] = x1;
    w[2] = x2;
    
    double u0[3], u1[3], u2[3];
    
	if (x1 == x2) {
        if (x0 == x1) {
            U[0] = 1.0; U[1] = 0.0; U[2] = 0.0;
            U[3] = 0.0; U[4] = 1.0; U[5] = 0.0;
            U[6] = 0.0; U[7] = 0.0; U[7] = 1.0;
            return 0;
        }
	 	else {
			find_eigen_vector(A,x0,u0);
			find_eigen_vector(A,x1,u1);
			cross_product(u0,u1,u2);
        }
    }
	else {
        find_eigen_vector(A,x0,u0);
		find_eigen_vector(A,x1,u1);
		find_eigen_vector(A,x2,u2);
    }

    U[0] = u0[0];
    U[3] = u0[1];
    U[6] = u0[2];
    
    U[1] = u1[0];
    U[4] = u1[1];
    U[7] = u1[2];
    
    U[2] = u2[0];
    U[5] = u2[1];
    U[8] = u2[2];
    return 0;
}

//int main() {
    //double A[9] = {
        //0.09660759,  0.08685899,  0.08720136 ,
        //0.08685899,  0.08838769,  0.0867068 ,
        //0.08720136,  0.0867068 ,  0.08728963
    //};
    
    //double w[3];
    //double U[9];
    //eigen_sym_3D(A,w,U);
    
    //int i,j;
    //for (i=0;i<3;i++) printf("%f ",w[i]); printf("\n");
    //printf("\n");
    //for (i=0;i<3;i++){
        //for (j=0;j<3;j++) printf("%f ",U[3*i+j]); printf("\n");
    //}
    
    //return 0;
//}
