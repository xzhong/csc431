#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

const double pi=3.1415;
const double e=2.718282;
const double ap=1e-6;
const double h=1e-5;
const double rp=1e-4;

typedef string Exception;

class Function {
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
			  throw Exception("Error_D(g)(x)>=1");
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
  double optimaize_golden_search(double a, double b, int ns=100){
	  double tau = (sqrt(5.0)-1.0)/2.0;
	  double x1 = a+(1.0-tau)*(b-a);
	  double x2 = a+tau)*(b-a);
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
};