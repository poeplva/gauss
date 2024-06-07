#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void usage() {
	cout << endl << "Usage:" << endl;
	cout << "This program uses Gaussian elimination with partial pivoting to solve linear equations." << endl;
	cout << "The program accepts two arguments:" << endl;
	cout << "The first argument is the location of the file that contains an nxn matrix, A." << endl;
	cout << "The second argument is the location of the file that contains an nx1 vector, b." << endl;
	cout << "The program will then output the solution of Ax=b to a file named x.txt in the current working directory." << endl;
}

double max(double x, double y) { // returns the maximum of two numbers
	if (x > y) return x;
	return y;
}

int number_of_lines(string filename) { // number of lines in the files gives the dimensions
	ifstream fh;
	fh.open(filename);
	string line;
	int c = 0;
	while (getline(fh, line)) { // loop through lines until there are non left, and count
		c++;
	}
	fh.close();
	return c;
}

void read_matrix(string A_file, int n, double A[]) { // reads the matrix given its dimensions
	ifstream fh;
	fh.open(A_file);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fh >> A[i*n+j]; // read until the space character
		}
	}
	fh.close();
}

void read_vector(string b_file, int n, double b[]) { // reads the vector given its dimensions
	ifstream fh;
	fh.open(b_file);
	for (int i = 0; i < n; i++) {
		fh >> b[i]; // read until the space character
	}
	fh.close();
}

void interchange_rows(double A[], int n, int m, int i1, int i2) { // interchanges rows i_1 and i_2 in a nxm matrix
	for (int j = 0; j < m; j++) { // loop through the columns and interchange the entries A_(i_1)j and A_(i_2)j one by one
		double first = A[i1*m+j];
		A[i1*m+j] = A[i2*m+j];
		A[i2*m+j] = first;
	}
}

void gaussian_elimination(double A[], double b[], int n, int k) { // does gaussian elimination in the kth column, no need to check if A_kk==0 since we are doing partial pivoting
	for (int i = k + 1; i < n; i++) { // loop through each row below the kth row
		double multiplier = -A[i*n+k]/A[k*n+k]; // mutliplying the kth row with -A_ik/A_kk and adding it to the ith row, we get int he ij entry
		for (int j = 0; j < n; j++) {
			A[i*n+j] = A[i*n+j] + multiplier * A[k*n+j];
		}
		b[i] = b[i] + multiplier * b[k];
	}
}

void partial_pivoting(double A[], double b[], int n, int j) { // partial pivoting in the jth column
	double max = abs(A[j*n+j]);
	int max_i = j;
	for (int i = j + 1; i < n; i++) { // determine the absolute maximum entry in the jth column below the jth entry
		if (abs(A[i*n+j]) > max) {
			max = abs(A[i*n+j]);
			max_i = i;
		}
	}
	if (max_i != j) { // interchange the jth row with the row below jth that has absolute maximum entry in the jth column
		interchange_rows(A, n, n, j, max_i);
		interchange_rows(b, n, 1, j, max_i);
	}
	//for (int i = 0; i = n; i++) {
	//	for (int j = 0; j <= n; j++) {
	//		if (j==n) {
	//			cout << b[i] << endl;
	//		} else {
	//			cout << A[i*n+j] << " ";
	//		}
	//	}
	//}
}

void gaussian_elimination_with_partial_pivoting(double A[], double b[], int n) { // do gaussian elimination with partial pivoting
	for (int j = 0; j < n; j++) { // loop through each column
		partial_pivoting(A, b, n, j);
		gaussian_elimination(A, b, n, j);
	}
}

bool is_singular(double A[], int n) { // check if A is singular, since we are going to pass upper triangular matrices, it is enough to check that the diagonal is nonzero
	for (int i = 0; i < n; i++) {
		if (A[i*n+i] == 0) { // if an entry in the diagonal is 0, then the matrix is singular
			return true;
		}
	}
	return false; // non-singular otherwise
}

// we have x_(n-i)=(b_(n-i)-A_(n-i)(n-1)*x_(n-1)-...-A_(n-i)(n-i+1)*x_(n-i+1))/A_(n-i)(n-i) in backward substitution, below code implements this

void back_substitution(double A[], double b[], double x[], int n) { // backward substitution
	for (int i = 1; i <= n; i++) {
		x[n-i] = b[n-i];
		for (int k = 1; k < i; k++) {
			x[n-i] = x[n-i] - A[(n-i)*n+n-k]*x[n-k];
		}
		x[n-i] = x[n-i] / A[(n-i)*n+n-i];
	}
}

// 1-norm is the maximum column sum
// inf-norm is the maximum row sum
// if our matrix is
// a	b
// c	d
// its inverse is
// d/e	-b/e
// -c/e	a/e
// where e = ad-bc is the determinant of the matrix
// then the condition number by 1-norm becomes norm_1(A)*norm_1(A^-1) = max(a+b,c+d)*max(d-b,a-c)/e
// and the condition number by inf-norm becomes norm_inf(A)*norm_inf(A^-1) = max(a+c,b+d)*max(d-c,a-b)/e
// below two functions implement these

double condition_number_one_2(double A[]) { // condition number by 1-norm of 2x2 matrix
	double det = A[0*2+0]*A[1*2+1]-A[0*2+1]*A[1*2+0];
	double max_column_sum = max(A[0*2+0]+A[1*2+0],A[0*2+1]+A[1*2+1]);
	double max_inv_column_sum = max(A[1*2+1]-A[1*2+0],A[0*2+0]-A[0*2+1]) / det;
	return (max_column_sum * max_inv_column_sum);
}

double condition_number_inf_2(double A[]) { // condition number by inf-norm of 2x2 matrix
	double det = A[0*2+0]*A[1*2+1]-A[0*2+1]*A[1*2+0];
	double max_row_sum = max(A[0*2+0]+A[0*2+1],A[1*2+0]+A[1*2+1]);
	double max_inv_row_sum = max(A[1*2+1]-A[0*2+1],A[0*2+0]-A[1*2+0]) / det;
	return (max_row_sum * max_inv_row_sum);
}

int main(int argc, char* argv[]) {
	if (argc != 3) { // if a wrong number of arguments are passed, exit
		cout << "You should pass exactly 2 arguments." << endl;
		usage();
		return 1;
	}
	string A_file, b_file;
	A_file = argv[1];
	b_file = argv[2];

	int n_A = number_of_lines(A_file);
	int n_b = number_of_lines(b_file);

	if (n_A != n_b) { // if the dimensions don't match, exit
		cout << "Dimensions do not match." << endl;
		usage();
		return 1;
	}

	int n = n_A;

	double* A = new double[n*n];
	read_matrix(A_file, n, A);
	double* b = new double[n];
	read_vector(b_file, n, b);

	double cond_1, cond_inf;
	if (n == 2) {
		cond_1 = condition_number_one_2(A);
		cond_inf = condition_number_inf_2(A);
	}

	gaussian_elimination_with_partial_pivoting(A, b, n);

	if (is_singular(A, n)) { // if the matrix is singular, exit with code 2
		cout << "Matrix is singular." << endl;
		return 2;
	}

	if (n == 2) {
		cout << "Condition number by 1-norm is " << cond_1 << endl;
		cout << "Condition number by inf-norm is " << cond_inf << endl;
		cout << endl;
	}

	double* x = new double[n];

	back_substitution(A, b, x, n);

	cout << "The solution is" << endl;
	ofstream fh;
	fh.open("x.txt");
	for (int i = 0; i < n; i++) { // print out the solution to stdout and to the file x.txt
		cout << x[i] << endl;
		fh << x[i] << endl;
	}
	fh.close();

	return 0;
}
