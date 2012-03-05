
#include "matrix - 120301.h"
#include <iostream>
using namespace std;

int main() {
	Matrix A(3,3);
	Matrix B(3,3);
	Matrix C(3,3);
	A(0,0)=-11; A(0,1)=2; A(0,2)=3;
	A(1,0)=1; A(1,1)=0; A(1,2)=3;
	A(2,0)=2; A(2,1)=2; A(2,2)=4;
	cout <<" The Matrix: " << A << endl;
	cout << "\n"<< endl;
	
	//test the getitem & setitem
	float d = A.getitem(2,1);  
	cout << "get A[2,1]: " << d << endl;
	A.setitem(0,0,-7);  
	cout << "set A[0,0]= -7, the new matrix: "<< A << endl;
	cout << "\n"<< endl;

	//test the Transposed & Inverse of Matrix
	B = t(A);   cout << "Transposed: " << B << endl;
	C = inv(A); cout << "Inverse: " << C << endl;
	cout << A << endl;

	//The absolute & negative value of Matrix; swap_rows & swap_cols
	C = abs(A); cout << "Absolute value: " << C << endl;
	C = neg(A); cout << "Negative value: " << C << endl;
	C = swap_rows(A,0,1); cout << "swap row 1,3: " << C << endl;
	C = swap_cols(A,0,2); cout << "swap column 1,3: " << C << endl;
	cout << A << endl;
	cout << "\n"<< endl;
	
	//test the operations
	C = A+B; cout << "A+B: " << C << endl;
	C = A-B; cout << "A-B: " << C << endl;
	C = A*B; cout << "A*B: " << C << endl;
	C = A/B; cout << "A/B: " << C << endl;
	double a = 5;
	C = 5/A; cout << "5/A: " << C << endl;
	C = A/5; cout << "A/5: " << C << endl;
	cout << "\n"<< endl;

	//test the norm_1 & norm_inf
	float b = 0, c = 0;
	b = norm_1(A);
	c = norm_inf(A);
	cout << "1-norm: " << b << endl
		 << "infinity-norm: " << c << endl;
	cout << A << endl;
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

	//test the sym
	Matrix D(3,3);
	D(0,0)=1; D(0,1)=0.00009; D(0,2)=3;
	D(1,0)=0.000001; D(1,1)=0; D(1,2)=3; 
	D(2,0)=3; D(2,1)=3; D(2,2)=4;
	
	if (sym(D)) cout << "ture" << endl;
	else cout << "false" << endl;
	
	system("pause");

	return 0;
}

