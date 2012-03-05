#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;

class Matrix {
public:
	int rows;
	int cols;
	vector<float> data;
	Matrix(int rows, int cols) {
		this->rows=rows;
		this->cols=cols;
		this->data.resize(rows*cols);
		for(int i=0; i<rows; i++)
			for(int j=0; j<cols; j++)
				this->data[i*cols+j]= 0;
	}
	float operator()(int i, int j) const {
		return data[i*cols+j];
	}
	float &operator()(int i, int j) {
		return data[i*cols+j];
	}
	float getitem(int i, int j) {
		return data[i*cols+j];
	}
	float setitem(int i, int j, float v) {
		data[i*cols+j] = v;
		return 0;
	}
};

ostream &operator<<(ostream &out, const Matrix& A) {
	out << "[";
	for(int i=0; i<A.rows; i++) {
		if(i>0) out << ",";
		out << "[";
		for(int j=0; j<A.cols; j++) {
			if(j>0) out << ",";
			out << A(i,j);
		}
		out << "]";
	}
	out << "]";
	return out;
}

Matrix operator+(const Matrix &A, const Matrix &B) {
	if(A.rows!=B.rows || A.cols!=B.cols)
		cout << "BAD\n";
	Matrix C(A.rows,A.cols);
	for(int i=0; i<A.rows; i++)
		for(int j=0; j<A.cols; j++)
			C(i,j)=A(i,j)+B(i,j);
	return C;
}

Matrix operator-(const Matrix &A, const Matrix &B) {
	if(A.rows!=B.rows || A.cols!=B.cols)
		cout << "BAD\n";
	Matrix C(A.rows,A.cols);
	for(int i=0; i<A.rows; i++)
		for(int j=0; j<A.cols; j++)
			C(i,j)=A(i,j)-B(i,j);
	return C;
}

Matrix operator*(float a, const Matrix &B) {
	Matrix C(B.rows,B.cols);
	for(int r=0; r<B.rows; r++)
		for(int c=0; c<B.cols; c++)
			C(r,c)=a*B(r,c);
	return C;
}

Matrix operator*(const Matrix &A, const Matrix &B) {
	if(A.cols != B.rows)
		cout << "BAD\n";
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

void swap(float&a, float &b) {
	float c=a; a=b; b=c;
}

Matrix inv(Matrix A) {// Inverse of A
	if(A.cols!=A.rows)
		cout << "BAD\n";
	Matrix C(A.rows,A.cols);
	float p;
	float q;
	int m;
	for(int r=0; r<C.cols;r++)
		C(r,r)=1;
	for(int c=0; c<A.cols;c++) {
		m=c; p=A(c,c);
		for(int i=c+1; i<A.rows; i++)
			if(abs(A(i,c)) > abs(p)) {m=i; p=A(i,c);}
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

Matrix operator/(float a, const Matrix &B) {// scalar/matrix
	Matrix C(B.rows, B.cols);
	C = inv(B);
	return a*C;
}

Matrix operator/(const Matrix &A, float a) {// matrix/scalar
	Matrix C(A.rows, A.cols);
	return C = (1/a)*A;
}

Matrix operator/(const Matrix &A, const Matrix &B) {//matrix/matrix
	if( A.cols != B.cols || A.rows != B.rows )
		cout << "BAD\n";
	for(int r=0; r<A.rows; r++){
		for(int c=0; c<B.cols; c++)
			if(B(r, c) ==0){
				cout << "BAD\n";
				break;
			}
	}
	Matrix C(B.rows, B.cols);
	C = inv(B);
	return A*C;
}

Matrix t(const Matrix A) {// Transposed of A
	Matrix C(A.rows, A.cols);
	for(int i=0; i<A.rows; i++)
		for(int j=0; j<A.cols; j++)
			C(i, j) = A(j, i);
	return C;
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

Matrix swap_rows(Matrix A, int i, int j){
	Matrix C(A.rows, A.cols);
	for(int c=0; c<A.cols; c++){
		C(i,c) = A(i,c);
		A(i,c) = A(j,c);
		A(j,c) = C(i,c);
	}
	return A;
}

Matrix swap_cols(Matrix A, int i, int j){
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
		C(i,i)= 2.718282;
	return C;
}

float norm_1(Matrix A) {// The 1-norm of A
	Matrix C(A.rows, A.cols);
	C = abs(A);
	float max=0;
	if(C.rows==1 || C.cols==1) {
		for(int c=0; c<C.cols; c++)
			for(int r=0; r<C.cols; r++)
				max += C(r,c);
	}
	else{
		float sum=0;
		for(int r=0; r<C.rows; r++)
			max += C(r, 0);
		for(int c=1; c<C.cols; c++){
			for(int r=0; r<C.rows;r++)
				sum += C(r,c);
			if(sum>max)
				max = sum;
			sum = 0;
		}
	}
	return max;
}

float norm_inf(Matrix A) {// The infinity-norm of A
	Matrix C(A.rows, A.cols);
	C = abs(A);
	float max=0;
	if(C.rows==1 || C.cols==1) {
		for(int r=0; r<C.rows; r++)
			for(int c=0; c<C.cols; c++)
				max += C(r,c);
	}
	else{
		float sum=0;
		for(int c=0; c<C.cols; c++)
			max += C(0, c);
		for(int r=1; r<C.rows; r++){
			for(int c=0; c<C.cols;c++)
				sum += C(r,c);
			if(sum>max)
				max = sum;
			sum = 0;
		}
	}
	return max;
}