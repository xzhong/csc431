#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

const double pi=3.1415;
const double e=2.718282;
const double ap=1e-6;
const double rp=1e-4;
const double h=1e-5;

typedef string Exception;

class Matrix {
public:
	int rows;
	int cols;
	vector<double> data;
	Matrix(int rows, int cols) {
		if(rows<=0 || cols<=0) 
			throw Exception("InvalidMatrix");
		this->rows=rows;
		this->cols=cols;
		this->data.resize(rows*cols);
		for(int i=0; i<rows; i++)
			for(int j=0; j<cols; j++)
				this->data[i*cols+j]= 0;
	}
	double operator()(int r, int c) const {
		if(r<0 || r>=rows || c<0 || c>=cols) 
			throw Exception("OutOfBounds");
		return data[r*cols+c];
	}
	double &operator()(int r, int c) {
		if(r<0 || r>=rows || c<0 || c>=cols) 
			throw Exception("OutOfBounds");
		return data[r*cols+c];
	}
	double getitem(int r, int c) {
		if(r<0 || r>=rows || c<0 || c>=cols) 
			throw Exception("OutOfBounds");
		return data[r*cols+c];
	}
	double setitem(int r, int c, double v) {
		data[r*cols+c] = v;
		return 0;
	}
};

ostream &operator<<(ostream &out, const Matrix& A) {
	out << "[";
	for(int r=0; r<A.rows; r++) {
		if(r>0) out << ",";
		out << "[";
		for(int c=0; c<A.cols; c++) {
			if(c>0) out << ",";
			out << A(r,c);    
		}
		out << "]";
	}
	out << "]";    
	return out;
}

Matrix operator+(const Matrix &A, const Matrix &B) {
	if(A.rows!=B.rows || A.cols!=B.cols)
		throw Exception("WrongMatrixSize");
	Matrix C(A.rows,A.cols);
	for(int r=0; r<A.rows; r++)
		for(int c=0; c<A.cols; c++)
			C(r,c)=A(r,c)+B(r,c);
	return C;
}

Matrix operator-(const Matrix &A, const Matrix &B) {
	if(A.rows!=B.rows || A.cols!=B.cols)
    	throw Exception("WrongMatrixSize");
	Matrix C(A.rows,A.cols);
	for(int r=0; r<A.rows; r++)
		for(int c=0; c<A.cols; c++)
			C(r,c)=A(r,c)-B(r,c);
	return C;
}

Matrix operator*(double a, const Matrix &B) {// scalar*matrix
	Matrix C(B.rows,B.cols);
	for(int r=0; r<B.rows; r++)
		for(int c=0; c<B.cols; c++)
			C(r,c)=a*B(r,c);
	return C;
}

Matrix operator*(const Matrix &A, const Matrix &B) {// matrix*matrix
	if(A.cols != B.rows)
		throw Exception("WrongMatrixSize");
	Matrix C(A.rows,B.cols);
    if(A.cols == B.rows){
		for(int r=0; r<A.rows; r++){
			for(int c=0; c<B.cols; c++){
				for(int k=0; k<A.cols; k++)
					C(r,c)+=A(r,k)*B(k,c);
			}
		}
	}
	return C;
}

void swap(double &a, double &b) {//assistant function
	double c=a; a=b; b=c;
}

Matrix inv(Matrix A) {// Inverse of A
	if(A.cols!=A.rows)
		throw Exception("NotSquared");
	Matrix C(A.rows,A.cols);
	double p;
	double q;
	int m;
	for(int r=0; r<C.cols;r++) 
		C(r,r)=1;
	for(int c=0; c<A.cols;c++) {    
		m=c; p=A(c,c);
		for(int i=c+1; i<A.rows; i++)
			if(abs(A(i,c)) > abs(p)){
				m=i; 
				p=A(i,c);
			}
			for(int i=0; i<A.cols; i++) {
				swap(A(m,i),A(c,i));
				swap(C(m,i),C(c,i));
			}
			for(int i=0; i<A.cols; i++) {
				A(c,i) /= p; 
				C(c,i) /= p;
			}
			for(int r=0; r<A.rows; r++) 
				if(r!=c) {
					q = A(r,c);
					for(int i=0; i<A.cols; i++) {
						A(r,i)-=q*A(c,i);
						C(r,i)-=q*C(c,i);
					}
				}
	}
	return C;
}

Matrix t(const Matrix A) {// Transposed of A
	Matrix C(A.rows, A.cols);
	for(int r=0; r<A.rows; r++)
		for(int c=0; c<A.cols; c++)
			C(r, c) = A(c, r);
	return C;
}

Matrix operator/(double a, const Matrix &B) {// scalar/matrix
	Matrix C(B.rows, B.cols);
	C = inv(B);
	return a*C;
}

Matrix operator/(const Matrix &A, double a) {// matrix/scalar
	Matrix C(A.rows, A.cols);
	return C = (1/a)*A;
}

Matrix operator/(const Matrix &A, const Matrix &B) {// matrix/matrix
	if(A.cols != B.rows)
		throw Exception("WrongMatrixSize");
	Matrix C(B.rows, B.cols);
	C = inv(B);
	return A*C;
}

Matrix neg(const Matrix A){ //negative value of A
	Matrix C(A.rows, A.cols);
	for(int i=0; i<A.rows; i++)
		for(int j=0; j<A.cols; j++)
			C(i, j) = (-1)*A(i, j);
	return C;
}

Matrix abs(const Matrix A){ //absolute value of A
	Matrix C(A.rows, A.cols);
	for(int i=0; i<A.rows; i++)
		for(int j=0; j<A.cols; j++)
			C(i, j) = abs(A(i, j));
	return C;
}

Matrix swap_rows(Matrix A, int i, int j){// swap rows i,j
	if(i<0 || i>=A.rows || j<0 || j>=A.rows) 
		throw Exception("OutOfBounds");
	Matrix C(A.rows, A.cols);
	for(int c=0; c<A.cols; c++){
		C(i,c) = A(i,c);
		A(i,c) = A(j,c);
		A(j,c) = C(i,c);
	}
	return A;
}

Matrix swap_cols(Matrix A, int i, int j){// swap cols i,j
	if(i<0 || i>=A.cols || j<0 || j>=A.cols) 
		throw Exception("OutOfBounds");
	Matrix C(A.rows, A.cols);
	for(int r=0; r<A.rows; r++){
		C(r,i) = A(r,i);
		A(r,i) = A(r,j);
		A(r,j) = C(r,i);
	}
	return A;
}

Matrix diagonal(int r,int v){// Constuct a diagonal matrix
	// v:the value in diagonal is v
	// r:the integer number of rows (also number of columns)
	Matrix C(r, r);
	for(int i=0; i<r; i++)
		C(i,i)= v ;
	return C;
}

Matrix diagonal_e(int r){// Constuct a diagonal matrix which the value in diagonal is e - six decimal
	Matrix C(r, r);
	for(int i=0; i<r; i++)
		C(i,i)= e;
	return C;
}

double norm_1(Matrix A) {// The 1-norm of A
	Matrix C(A.rows, A.cols);
	C = abs(A);
	double max=0;
	double sum=0;
	for(int c=0; c<C.cols; c++){
		sum=0;
		for(int r=0; r<C.rows;r++)
			sum += C(r,c);
		if(sum>max) 
			max = sum;
		}
	return max;
}

double norm_inf(Matrix A) {// The infinity-norm of A
	Matrix C(A.rows, A.cols);
	C = abs(A);
	double max=0;
	double sum=0;
	for(int r=0; r<C.rows; r++){
		sum=0;
		for(int c=0; c<C.cols; c++)
			sum += C(r,c);
		if(sum>max) 
			max = sum;
		}
	return max;
}

bool operator!=(const Matrix &A, const Matrix &B) {
  return (norm_1(A-B)>ap);
}

bool operator==(const Matrix &A, const Matrix &B) {
  return (norm_1(A-B)<ap);
}

double cond_1(const Matrix A) {
  return norm_1(A)*norm_1(inv(A));
}

double cond_inf(const Matrix A) {
  return norm_inf(A)*norm_inf(inv(A));
}

double max(float a, float b) {//assistant function
	if(a-b>0) return a;
	else return b;
}

bool almost_sym(Matrix A){//almost symmetric detection
	if (A.rows !=A.cols) return false;
	else {
		for(int r=0; r<A.rows; r++){
		    for(int c=0; c<=r; c++) {
			    double delta = abs(A(r,c)-A(c,r));
			    if( delta>ap && (delta/max(abs(A(r,c)),abs(A(c,r))))>rp)
				    return false;
			}
		}
		return true;
	}
}

bool sym(Matrix A){//symmetric detection
	if(A == t(A)) return true;
	else return false;
}

Matrix cholesky(const Matrix A){
	if(sym(A))
		throw Exception("NotSymmetric");
	Matrix B=A;
	for(int c=0; c<B.cols; c++){
		if(B(c,c) <=0)
			throw Exception("NotPositiveDefinitive");
		B(c,c) = sqrt(B(c,c));
		for(int r=c+1; r<B.rows; r++){
			B(r,c) /= B(c,c);
		}
		for(int i=c+1; i<B.rows; i++){
			for(int j=c+1; j<B.rows; j++)
				B(j,i) -= B(j,c)*B(i,c);
		}
	}
	for(int r=0; r<B.rows; r++){
		for(int c=r+1; c<B.cols; c++)
			B(r,c)=0;
	}
	return B;
}

Matrix markoviz(Matrix mu, Matrix A, double r_free) {
	Matrix x(A.rows, 1);
	for(int r=0; r<mu.rows; r++) 
		mu(r,0) -= r_free;
	x=inv(A)*mu;
	double y=0;
	for(int r=0; r<mu.rows; r++) 
		y += x(r,0);
	for(int r=0; r<mu.rows; r++) 
		x(r,0) /=y;
	return x;
}

Matrix fit(Matrix K, int n) {
	Matrix b(K.rows,1);
	Matrix A(K.rows,n);
	for(int i=0; i<K.rows; i++) {
		for(int j=0; j<n; j++)
			A(i,j)=pow(K(i,0),j);
		b(i,1)=K(i,1);
	}
	Matrix L=cholesky(t(A)*A);
	Matrix x=inv(t(L))*inv(L)*t(A)*b;
  return x;
}

//---------------------------------------------------------solver-------------------------------------------------------

class Function{
public:
  virtual double f(double x)=0;
  
  double Df(double x) {
	  return (f(x+h)-f(x-h))/(2.0*h);
  }

  double DDf(double x) {
	  return (f(x+h)-2.0*f(x)+f(x-h))/(h*h);
  }

  virtual double g(double x){
	  return f(x)+x;
  }

   double Dg(double x) {
	  return (g(x+h)-g(x-h))/(2.0*h);
  }

  double solve_fixed_point(double x_guess, int ns=100){
	  double x_old;
	  double x = x_guess;
	  double Dgx;
	  for(int k=0; k<ns; k++) {
		  Dgx = Dg(x);
		  if( abs(Dgx)>=1) 
			  throw Exception("Error:D(g)(x)>=1");
		  x_old = x;
		  x = g(x);
		  if( k>2 && abs(x-x_old)<max(ap,rp*abs(x))) 
			  return x;
	  }
	  throw Exception("NoConvergence");
  }  
  
  double solve_bisection(double a, double b, int ns=100){
	  double fa = f(a), fb = f(b);
	  if( fa==0) return a;
	  if( fb==0) return b;
	  if( fa*fb>0)
		  throw Exception("f(a) and f(b) must have opposite sign");
	  double x;
	  double fx;
	  for(int k=0; k<ns; k++){
		  x = (a+b)/2;
		  fx = f(x);
		  if( fx==0 || abs(b-a)<max(ap,rp*abs(x)))
			  return x;
		  else{
			  if( fx*fa<0){
				  b = x;
				  fb = fx;}
			  else{
				  a = x;
				  fa = fx;}
		  }
	  }
	  throw Exception("NoConvergence");
  }  
				 
  double solve_newton(double x_guess, int ns=100) {
	  double x_old;
	  double x = x_guess;
	  double fx, Dfx;
	  for(int k=0; k<ns; k++) {
		  fx = f(x);
		  Dfx = Df(x);
		  if( abs(Dfx)<ap) 
			  throw Exception("UnstableSolution");
		  x_old = x;
		  x = x - fx/Dfx;
		  if( k>2 && abs(x-x_old)<max(ap,rp*abs(x))) 
			  return x;
	  }
	  throw Exception("NoConvergence");
  }  

  double solve_secant(double x_guess, int ns=20){
	  double x_old, fx_old;
	  double x = x_guess, fx = f(x), Dfx = Df(x);
	  for(int k=0; k<ns; k++) {
		  if( abs(Dfx)<ap) 
			  throw Exception("UnstableSolution");
		  x_old = x;
		  fx_old = fx;
		  x = x-fx/Dfx;
		  if( k>2 && abs(x-x_old)<max(ap,rp*abs(x))) 
			  return x;
		  fx = f(x);
		  Dfx = (fx-fx_old)/(x-x_old);
	  }
	  throw Exception("NoConvergence");
  }  

  double solve_newton_stabilized(double a, double b, int ns=20){
	  double fa = f(a), fb = f(b);
	  if( fa==0) return a;
	  if( fb==0) return b;
	  if( fa*fb>0)
		  throw Exception("f(a) and f(b) must have opposite sign");
	  double x = (a+b)/2;
	  double x_old, fx_old;
	  double fx = f(x), Dfx = Df(x);
	  for(int k=0; k<ns; k++){
		  x_old = x;
		  fx_old = fx;
		  if( abs(Dfx)>ap)
			  x = x - fx/Dfx;
		  if( x==x_old || x<a || x>b)
			  x = (a+b)/2;
		  fx = f(x);
		  if( fx==0 || abs(x-x_old)<max(ap,rp*abs(x)))
			  return x;
		  Dfx = (fx-fx_old)/(x-x_old);
		  if( fx*fa<0){
			  b = x;
			  fb = fx;
		  }else{
			  a = x;
			  fa = fx;
		  }
	  }
	  throw Exception("NoConvergence");
  }
		
  double optimize_bisection(double a, double b, int ns=100){
	  double Dfa = Df(a), Dfb = Df(b);
	  if( Dfa==0) return a;
	  if( Dfb==0) return b;
	  if( Dfa*Dfb>0)
		  throw Exception("D(f)(a) and D(f)(b) must have opposite sign");
	  double x, Dfx;
	  for(int k=0; k<ns; k++){
		x = (a+b)/2;
		Dfx = Df(x);
		if (Dfx==0 || abs(b-a)<max(ap,rp*abs(x)))
			return x;
		else{
			if(Dfx*Dfa<0){
				b = x;
				Dfb = Dfx;}
			else{
				a = x;
				Dfa = Dfx;}
		}
	  }
	  throw Exception("NoConvergence");
  }

  double optimize_newton(double x_guess, int ns=20) {
	  double x_old;
	  double x = x_guess;
	  double Dfx, DDfx;
	  for( int k=0; k<ns; k++) {
		  Dfx = Df(x);
		  if( Dfx==0) 
			  return x;
		  DDfx = DDf(x);
		  if( abs(DDfx)<ap)
			  throw Exception("UnstableSolution");
		  x_old = x;
		  x = x - Dfx/DDfx;
		  if( abs(x-x_old)<max(ap,rp*abs(x))) 
			  return x;
	  }
	  throw Exception("NoConvergence");
  }

  double optimize_secant(double x_guess, int ns=100){
	  double x = x_guess; 
	  double fx=f(x), Dfx=Df(x), DDfx=DDf(x);
	  double x_old, Dfx_old;
	  for( int k=0; k<ns; k++){
		  if( Dfx==0)
			  return x;
		  if( abs(DDfx)<ap)
			  throw Exception("UnstableSolution");
		  x_old = x;
		  Dfx_old = Dfx;
		  x = x - Dfx/DDfx;
		  if( abs(x-x_old)<max(ap,rp*abs(x))) 
			  return x;
		  fx = f(x);
		  Dfx = Df(x);
		  DDfx = (Dfx-Dfx_old)/(x-x_old);
	  }
	  throw Exception("NoConvergence");
  }

  double optimize_newton_stabilized(double a, double b, int ns=20){
	  double Dfa=Df(a), Dfb=Df(b);
	  if( Dfa==0) return a; 
	  if( Dfb==0) return b;
	  if( Dfa*Dfb>0)
		  throw Exception("D(f)(a) and D(f)(b) must have opposite sign");
	  double x, fx, Dfx, DDfx;
	  x = (a+b)/2;
	  fx = f(x);
	  Dfx = Df(x);
	  DDfx = DDf(x);
	  double x_old, fx_old, Dfx_old;
	  for( int k=0; k<ns; k++){
		  if( Dfx==0)
			  return x;
		  x_old = x;
		  fx_old = fx;
		  Dfx_old = Dfx;
		  if( abs(DDfx)>ap)
			  x = x - Dfx/DDfx;
		  if( x=x_old || x<a || x>b)
			  x = (a+b)/2;
		  if( abs(x-x_old)<max(ap,rp*abs(x))) 
			  return x;
		  fx = f(x);
		  Dfx = (fx-fx_old)/(x-x_old);
		  DDfx = (Dfx-Dfx_old)/(x-x_old);
		  if( Dfx*Dfa <0){
			  b = x;
			  Dfb = Dfx;}
		  else{
			  a = x;
			  Dfa = Dfx;}
	  }
	  throw Exception("NoConvergence");
  }

  double optimize_golden_search(double a, double b, int ns=100){
	  double tau = (sqrt(5.0)-1.0)/2.0;
	  double x1 = a+(1.0-tau)*(b-a);
	  double x2 = a+tau*(b-a);
	  double fa=f(a), fb=f(b), f1=f(x1), f2=f(x2);
	  for( int k=0; k<ns; k++){
		  if( f1>f2){
			  a = x1;
			  fa = f1;
			  x1 = x2;
			  f1 = f2;
			  x2 = a+tau*(b-a);
			  f2 = f(x2);}
		  else{
			  b = x2;
			  fb = f2;
			  x2 = x1;
			  f2 = f1;
			  x1 = a+(1.0-tau)*(b-a);
			  f2 = f(x2);
		  }
		  if( k>2 && abs(b-a)<max(ap,rp*abs(b)))
			  return b;
	  }
	  throw Exception("NoConvergence");
  }

  //-------------------------------------------------------Integral---------------------------------------------------------------

  double integrate(double a, double b){
	  double I=0;
	  double h=0;
	  double I_old=0;
	  for(int n=2;; n*=2) {
		  h=(b-a)/n;
		  I_old = I;
		  I = 0;
		  for(int i=0; i<n; i++) 
			  I += h*f(a+i*h);
		  if((n>2) && abs(I-I_old)<max(ap,rp*abs(I))) 
			  return I;
	  }
	  return I;
  }

  double  integrate_quadrature(double a, double b, int n=15) {
	double h=(b-a)/n;
	Matrix A(n,n), c(n,1), w(n,1);
	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++)
			A(i,j)= pow(a+(j+1)*h,i);
		c(i,0)=(pow(b,i+1)-pow(a,i+1))/(i+1);
	}
	w=inv(A)*c;
	double I = 0;
	for(int i=0; i<n; i++) {
		I+=w(i,0)*f(a+(i+1)*h);
	}
	return I;	
  }

  double integral_adaptative(double a, double b, int n1=3, int n2=5){
	  double I1=integrate_quadrature(a,b,n1);
	  double I2=integrate_quadrature(a,b,n2);
	  if(abs(I1-I2)<ap) 
		  return I2;
	  else{
		  double c = (a+b)/2;
		  double I_left = integral_adaptative(a,c,n1,n2);
		  double I_right = integral_adaptative(c,b,n1,n2);
		  return I_left+I_right;
	  }
  }

};

//-----------------------------------------------assistant functions for test------------------------------------------------------

 class IntegrateTestFunction : public Function {
	 double f(double x) {
		 return sin(x)+x*x+1;
	 }
 };

class SimpleTransaction {
	//continuous compounding
public:
	double A, t;
	//A: final value
	//t: total time in years
	SimpleTransaction(double A, double t) {
		this->A = A;
		this->t = t;
	}
	double present_value(float r) {
		//r: annual nominal interest rate (as a decimal, not in percentage)
		return A*exp(-r*t);
	}
};

class ComplexTransaction {
public:
	vector<SimpleTransaction> transactions;
	double present_value(float r) {
		double total = 0;
		for(int k=0; k<transactions.size(); k++) {
			total += transactions[k].present_value(r);
		}
	return total;
	}
};

class MyBondFunction : public Function {
public:
  double A, t, p;
  int n;
  double f(double r) {
    ComplexTransaction bond;
    for(int i=1; i<=n; i++)
      bond.transactions.push_back(SimpleTransaction(A,t*i));
    return bond.present_value(r)-p;
  }
};

class Portfolio2 : public Function {
	//A portfolio contained two risky assets
public:
	double r1, r2, sigma1, sigma2, rho, r_free;
	//r1:expeted return of asset1
	//r2:expeted return of asset2
	//sigma1:SD of asset1
	//sigma2:SD of asset2
	//rho:correlation coefficient of the tow assets
	//r_free:risk-free interest
	double R(double x) {
		//portfolio expected return
		//x: the proportion of asset1
		return x*r1+(1.0-x)*r2;
	}
	double sigma(double x) {
		//portfolio SD
		return sqrt(x*x*sigma1*sigma1+2.0*x*(1.0-x)*sigma1*sigma2*rho+(1.0-x)*(1.0-x)*sigma2*sigma2);
	}
	double sharpe(double x) {
		//sharpe ratio
		return (R(x)-r_free)/sigma(x);
	}
	double f(double x) {
		return sharpe(x);
	}
};

