#include <iostream>
#include "final - 120315.h"

using namespace std;

int main() {

//--------------------------------------------------test of matrix-------------------------------------------------------
	
	Matrix A(3,3);
	A(0,0)=-11; A(0,1)=2; A(0,2)=3;
	A(1,0)=1; A(1,1)=0; A(1,2)=3;
	A(2,0)=2; A(2,1)=2; A(2,2)=4;
	cout <<"Matrix A: " << A << endl;
	cout << "\n"<< endl;
	
	//test the getitem & setitem
	double d = A.getitem(2,1);  
	cout << "get A[2,1]: " << d << endl;
	A.setitem(0,0,-7);  
	cout << "set A[0,0]= -7, the new matrix: "<< A << endl;
	cout <<"New Matrix A: " << A << endl;
	cout << "\n"<< endl;

	//test the Transposed, Inverse, Absolute & negative value of Matrix; swap_rows & swap_cols
	
	Matrix B(3,3);
	Matrix C(3,3);
	B = t(A);   cout << "Transposed: " << B << endl;
	C = inv(A); cout << "Inverse: " << C << endl;
	C = abs(A); cout << "Absolute value: " << C << endl;
	C = neg(A); cout << "Negative value: " << C << endl;
	C = swap_rows(A,0,1); cout << "swap row 1,3: " << C << endl;
	C = swap_cols(A,0,2); cout << "swap column 1,3: " << C << endl;
	cout << "\n"<< endl;
	
	//test the operations
	cout << "Matrix A: " << A << endl;
	cout << "Matrix B: " << B << endl;
	C = A+B; cout << "A+B: " << C << endl;
	C = A-B; cout << "A-B: " << C << endl;
	C = A*B; cout << "A*B: " << C << endl;
	C = A/B; cout << "A/B: " << C << endl;
	double a = 5;
	C = 5/A; cout << "5/A: " << C << endl;
	C = A/5; cout << "A/5: " << C << endl;
	cout << "\n"<< endl;

	//test the norm
	double b = 0, c = 0; 
	b = norm_1(A);
	c = norm_inf(A);
	cout << "1-norm of A: " << b << endl
		 << "infinity-norm of A: " << c << endl;

	Matrix V(3,1);
	V(0,0)=-11; V(1,0)=2; V(2,0)=3;
	double v = 0;
	v = norm_2(V);
	cout << "Matrix V: " << V << endl
		 << "2-norm of V: " << v << endl;
	cout << "\n"<< endl;
	
	//test the "!=" & "=="
	if (A==A)
		cout << "A,A EQUAL" << endl;
	if (A!=B)
		cout << "A,B NotEQUAL " << endl;
	
	//test the diagonal
	B = diagonal(4,1); //construct an I matix
	C = diagonal_e(3);
	cout << "An I matrix(): " << B << endl
		 << "A diagonal matrix(e): " << C << endl;
	cout << "\n"<< endl;
	
	//test the sym, markoviz
	Matrix D(3,3);
	D(0,0)=0.04; D(1,1)=0.09; D(2,2)=0.16;
	D(0,1)= D(1,0)=0.006; D(0,2)=D(2,0)=0.02; D(1,2)=D(2,1)=0.06;
    cout << "Matrix D: " << D << endl;
	if(sym(D)) 
		cout << "D is symetric" << endl;
	else cout <<"D is not symetric" << endl;
	cout << "\n"<< endl;

	Matrix mu(3,1);
	mu(0,0)=0.10; mu(1,0)=0.12; mu(2,0)=0.15;
	double r_free=0.05;
	cout << "mu: " << mu << endl;
	cout << "r_free: " << r_free << endl;
	Matrix x=markoviz(mu,D,r_free);
	cout << "markovize(mu,D,r_free): "<< x << endl;
	cout << "\n"<< endl;
	
//--------------------------------------------------test of integrate-------------------------------------------------------
	
	IntegrateTestFunction m;
	cout << "integrateTestFunction: sin(x)+x*x+1 in range(0,2)" << endl;
	double integrate1,integrate2,integrate3;
	integrate1 = m.integrate(0,2);
	integrate2 = m.integrate_quadrature(0,2);
	integrate3 = m.integrate_adaptative(0,2);
	cout << "integrate: " << integrate1 << endl
		 << "integrate_quadrature: " << integrate2 << endl
		 << "integrate_adaptative: " << integrate3 << endl;
	cout << "\n"<< endl;
		 
//--------------------------------------------------test of solvers-------------------------------------------------------
	
	Portfolio2 p;
    p.r1 = 0.05;
    p.r2 = 0.12;
    p.sigma1 = 0.01;
    p.sigma2 = 0.03;
    p.rho = 0.5;
    p.r_free = 0.02;
    cout << p.optimize_newton(0.5) << endl;
    MyBondFunction bond;
    bond.A = 5000;
    bond.n = 10;
    bond.t = 365;
    bond.p = 40000;
    cout << bond.solve_newton(0)*360 << endl;

	system("pause");

	return 0;
}
