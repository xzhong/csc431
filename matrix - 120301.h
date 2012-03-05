#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

const double pi=3.1415;
const double e=2.718282;
const double ap=1e-6;
const double rp=1e-4;

typedef string Exception;

class Matrix {
public:
	int rows;
	int cols;
	vector<double> data;
	Matrix(int rows, int cols) {
		if(rows<=0 || cols<=0) throw Exception("InvalidMatrix");
		this->rows=rows;
		this->cols=cols;
		this->data.resize(rows*cols);
		for(int i=0; i<rows; i++)
			for(int j=0; j<cols; j++)
				this->data[i*cols+j]= 0;
	}
	double operator()(int r, int c) const {
		if(r<0 || r>=rows || c<0 || c>=cols) throw Exception("OutOfBounds");
		return data[r*cols+c];
	}
	double &operator()(int r, int c) {
		if(r<0 || r>=rows || c<0 || c>=cols) throw Exception("OutOfBounds");
		return data[r*cols+c];
	}
	double getitem(int r, int c) {
		if(r<0 || r>=rows || c<0 || c>=cols) throw Exception("OutOfBounds");
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
		throw Exception("MatrixNotSquared");
	Matrix C(A.rows,A.cols);
	float p;
	float q;
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
	if( A.cols != B.cols || A.rows != B.rows )
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
	if(i<0 || i>=A.rows || j<0 || j>=A.rows) throw Exception("OutOfBounds");
	Matrix C(A.rows, A.cols);
	for(int c=0; c<A.cols; c++){
		C(i,c) = A(i,c);
		A(i,c) = A(j,c);
		A(j,c) = C(i,c);
	}
	return A;
}

Matrix swap_cols(Matrix A, int i, int j){// swap cols i,j
	if(i<0 || i>=A.cols || j<0 || j>=A.cols) throw Exception("OutOfBounds");
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

bool sym(Matrix A){//symmetry detection
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

