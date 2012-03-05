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
    for(int r=0; r<rows; r++)
      for(int c=0; c<cols; c++)
	this->data[r*cols+c]=0;
  }
  float operator()(int i, int j) const {
    return data[i*cols+j];
  }
  float &operator()(int i, int j) {
    return data[i*cols+j];
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
    cout << "BAD\n";
  Matrix C(A.rows,A.cols);
  for(int r=0; r<A.rows; r++)
    for(int c=0; c<A.cols; c++)
      C(r,c)=A(r,c)+B(r,c);
  return C;
}

Matrix operator-(const Matrix &A, const Matrix &B) {
  if(A.rows!=B.rows || A.cols!=B.cols)
    cout << "BAD\n";
  Matrix C(A.rows,A.cols);
  for(int r=0; r<A.rows; r++)
    for(int c=0; c<A.cols; c++)
      C(r,c)=A(r,c)-B(r,c);
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

Matrix operator/(const Matrix &A, const Matrix &B) {
	if( A.cols != B.cols || A.rows != B.rows )
		cout << "BAD\n";
	Matrix C(A.rows, B.cols);
	for(int r=0; r<A.rows; r++){
		for(int c=0; c<B.cols; c++)
			C(r, c) = A(r, c)/B(r, c);
	}
	return C;
}

void swap(float&a, float &b) {
  float c=a; a=b; b=c;
}

float max(float &a, float &b) {
	if(a-b>0) return a;
	else return b;
}

Matrix identity(int &a){
	Matrix C(a, a);
	for(int r=0; r<a; r++)
		C(r,r)=1;
	return C;
}

Matrix diagonal(int &a){
	Matrix C(a, a);
	for(int r=0; r<a; r++)
		C(r,r)=2.7;
	return C;
}

Matrix t(const Matrix &A) {
	Matrix C(A.rows, A.cols);
	for(int r=0; r<A.rows; r++)
		for(int c=0; c<A.cols; c++)
			C(r, c) = A(c, r);
	return C;
}

bool symmetric(Matrix &A){
	float ap=0.000001, rp=0.0001;
	if (A.rows !=A.cols) return false;
	else {
	    for(int r=0; r<A.rows; r++)
		    for(int c=0; c<A.cols; c++) {
			    float delta = abs(A(r, c)-A(c, r));
			    if(delta>ap && delta-max(abs(A(r, c)), abs(A(c, r)))*rp)
				    return false;
			    else return true;
		}
}
}

bool zero(Matrix &A){
	float ap=0.000001, rp=0.0001;
	for(int r=0; r<A.rows; r++){
		for(int c=0; c<A.cols; c++) {
			float delta = abs(A(r, c)-A(c, r));
			if(delta>ap && delta-max(abs(A(r, c)), abs(A(c, r)))*rp)
				return false;
			else return true;
		}
	}
}

float norm(Matrix &A) {
	float sum=0, max=0;
	if(A.rows==1 || A.cols==1) {
		for(int r=0; r<A.rows; r++)
			for(int c=0; c<A.cols; c++)
				max += A(r, c);
	}
	else {
		for(int c=0; c<A.cols; c++)
			max += A(0, c);
		for(int r=1; r<A.rows; r++){
			for(int c=0; c<A.cols; c++)
				sum += A(r, c);
			if( (sum-max)>0 )
				max=sum;
	return sum;
}
	}
}

Matrix inv(Matrix A) {
  if(A.cols!=A.rows)
    cout << "BAD\n";
  Matrix B(A.rows,A.cols);
  float p;
  float q;
  int m;
  for(int r=0; r<B.cols;r++) B(r,r)=1;
  for(int c=0; c<A.cols;c++) {    
    m=c; p=A(c,c);
    for(int i=c+1; i<A.rows; i++)
      if(abs(A(i,c)) > abs(p)) {m=i; p=A(i,c);}
    for(int i=0; i<A.cols; i++) {
      swap(A(m,i),A(c,i));
      swap(B(m,i),B(c,i));
    }
    for(int i=0; i<A.cols; i++) {
      A(c,i) /= p; 
      B(c,i) /= p;
    }
    for(int r=0; r<A.rows; r++) 
      if(r!=c) {
	q = A(r,c);
	for(int i=0; i<A.cols; i++) {
	  A(r,i)-=q*A(c,i);
	  B(r,i)-=q*B(c,i);
	}
      }
  }
  return B;
}

float condition_number(Matrix &A){
	Matrix C(A.rows, A.cols);
	float cn;
	C = inv(A);
	cn = norm(A)*norm(C);
	return cn;
}

int main() {
  Matrix A(3,3);
  Matrix B(3,3);
  A(0,0)=1; A(0,1)=2; A(0,2)=3;
  A(1,0)=1; A(1,1)=0; A(1,2)=3;
  A(2,0)=2; A(2,1)=2; A(2,2)=4;
  B = inv(A);
  cout << B*A << endl;
  return 0;
}


