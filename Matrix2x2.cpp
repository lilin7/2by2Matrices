/*
* @ClassName: Matrix2x2.cpp
* @Function: Class using methods and overloaded operators to implement the Mat2x2 ADT described in assignment 3
* @author: Lin Li
* @version: 1.0
* @date: Jul 14, 2018
*/

#include <string>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include "Matrix2x2.h"
using namespace std;

namespace {
	/*
	* get the length of a double number, to be used to print the Matrix2x2 object
	* @param const double
	* @return int
	*/
	int getLength(const double d)
	{
		auto leng = 3;
		if (d == 0) { return 4; }
		if (d < 0) {
			auto x = static_cast<int>(d);
			if (x == 0) {
				return 5;
			}
			else {
				while (x) {
					x /= 10;
					leng++;
				}
				return leng + 1;
			}
		}
		else {
			auto x = static_cast<int>(d);
			if (x == 0) {
				return 4;
			}
			else {
				while (x) {
					x /= 10;
					leng++;
				}
				return leng;
			}
		}
	}
}

/*
* Default constructor: sets all four entries to zero.
* @param 
* @return void
*/
Matrix2x2::Matrix2x2() {
	values.insert(values.begin(), 4, 0.0);
}

/*
* Normal constructor: accepts initial values for all four entries as parameters.
* @param double a
* @param double b
* @param double c
* @param double d
* @return void
*/
Matrix2x2::Matrix2x2(double a, double b, double c, double d) {
	values.reserve(4);
	values.push_back(a);
	values.push_back(b);
	values.push_back(c);
	values.push_back(d);
}

//Default copy constructor
Matrix2x2::Matrix2x2(const Matrix2x2& m) {
	values.insert(values.begin(), 4, 0.0);
	(*this) = m;
}

/*
* default assignment operator
* @param const Matrix2x2& n
* @return Matrix2x2
*/
Matrix2x2 & Matrix2x2::operator =(const Matrix2x2& n) {
	for (int i = 0; i < 4; ++i) {
		values[i] = n.values[i];
	}
	return *this;
}

/*
* default destructor
* @param 
* @return void
*/
Matrix2x2::~Matrix2x2() {
}


/*
* over loaded operator + for +M, positive M is M
* @param
* @return Matrix2x2
*/
Matrix2x2 Matrix2x2::operator + (){
	return *this;
}

/*
* over loaded operator - for -M, negative M
* @param
* @return Matrix2x2
*/
Matrix2x2 Matrix2x2::operator - (){
	Matrix2x2 ret;
	for (int i = 0; i < 4; i++) {
		ret.values[i] = -values[i];
	}
	return ret;
}

/*
* Compound Assignment Operator Overload M+=M'
* @param const Matrix2x2& n
* @return Matrix2x2&
*/
Matrix2x2& Matrix2x2::operator += (const Matrix2x2& n){
	for (int i = 0; i<4; i++){
		values[i] = values[i] + n.values[i];
	}
	return *this;
}

/*
* Compound Assignment Operator Overload M+=d
* @param const double& d
* @return Matrix2x2&
*/
Matrix2x2& Matrix2x2::operator += (const double& d){
	for (int i = 0; i<4; i++) {
		values[i] = values[i] + d;
	}
	return *this;
}

/*
* Compound Assignment Operator Overload M-=M'
* @param const Matrix2x2& n
* @return Matrix2x2&
*/
Matrix2x2& Matrix2x2::operator-=(const Matrix2x2& n){
	for (int i = 0; i<4; i++){
		values[i] = values[i] - n.values[i];
	}
	return *this;
}

/*
* Compound Assignment Operator Overload M-=d
* @param const double& d
* @return Matrix2x2&
*/
Matrix2x2& Matrix2x2::operator -= (const double& d){
	for (auto& value : values){
		value = value - d;
	}
	return *this;
}

/*
* Compound Assignment Operator Overload M*=M'
* @param const Matrix2x2& n
* @return Matrix2x2&
*/
Matrix2x2& Matrix2x2::operator*=(const Matrix2x2& n){
	double a = values[0] * n.values[0] + values[1] * n.values[2];
	double b = values[0] * n.values[1] + values[1] * n.values[3];
	double c = values[2] * n.values[0] + values[3] * n.values[2];
	double d = values[2] * n.values[1] + values[3] * n.values[3];
	values[0] = a;
	values[1] = b;
	values[2] = c;
	values[3] = d;
	return *this;
}

/*
* Compound Assignment Operator Overload M*=d
* @param const double & d
* @return Matrix2x2&
*/
Matrix2x2& Matrix2x2::operator *= (const double& d){
	for (int i = 0; i < 4; ++i){
		values[i] = values[i] * d;
	}
	return *this;
}

/*
* Compound Assignment Operator Overload M/=M'
* @param Matrix2x2& m
* @return Matrix2x2
*/
Matrix2x2& Matrix2x2::operator/=(Matrix2x2& n){
	Matrix2x2 temp = *this * (n.inverse());
	values[0] = temp.values[0];
	values[1] = temp.values[1];
	values[2] = temp.values[2];
	values[3] = temp.values[3];
	return *this;
}

/*
* Compound Assignment Operator Overload M/=d
* @param const double& a
* @return Matrix2x2&
*/
Matrix2x2& Matrix2x2::operator /= (const double& d){
	for (int i = 0; i < 4; ++i) {
		values[i] = values[i] / d;
	}
	return *this;
}

/*
* overloaded operator ++ for M++ post-increment
* @param
* @return Matrix2x2
*/
Matrix2x2 Matrix2x2::operator ++ (int){
	Matrix2x2 ret = *this;
	for (int i = 0; i < 4; ++i) {
		values[i] +=1 ;
	}
	return ret;
}

/*
* overloaded operator ++ for ++M preincrement
* @param
* @return Matrix2x2&
*/
Matrix2x2& Matrix2x2::operator ++ (){
	for (int i = 0; i < 4; ++i) {
		values[i] += 1;
	}
	return *this;
}

/*
* overloaded operator -- for M-- post-decrement
* @param
* @return Matrix2x2
*/
Matrix2x2 Matrix2x2::operator -- (int){
	Matrix2x2 ret = *this;
	for (int i = 0; i < 4; ++i) {
		values[i] -= 1;
	}
	return ret;
}

/*
* overloaded operator -- for --M pre-decrement
* @param
* @return Matrix2x2&
*/
Matrix2x2& Matrix2x2::operator -- (){
	for (int i = 0; i < 4; ++i) {
		values[i] -= 1;
	}
	return *this;
}


/*
* function Inversion
* @param
* @return Matrix2x2
*/
Matrix2x2 Matrix2x2::inverse() {
	if (values[0] * values[3] == values[1] * values[2]) {
		cout << "Inversion failed" << endl;
		return *this;
	}
	double x = 1 / (values[0] * values[3] - values[1] * values[2]);
	Matrix2x2 ret;
	ret.values[0] = x * values[3];
	ret.values[1] = 0 - (x * values[1]);
	ret.values[2] = 0 - (x * values[2]);
	ret.values[3] = x * values[0];
	return ret;
}

/*
* function Transpose
* @param
* @return Matrix2x2
*/
Matrix2x2 Matrix2x2::transpose()
{
	if ((values[0] * values[3] - values[1] * values[2])<1.0e-6 
		|| (values[0] * values[3] - values[1] * values[2]) == 1.0e-6){
		throw std::overflow_error("Inverse undefined");
	}		
	else{
		Matrix2x2 ret;
		ret.values[0] = values[0];
		ret.values[1] = values[2];
		ret.values[2] = values[1];
		ret.values[3] = values[3];
		return ret;
	}
}

/*
* function determinant
* @param
* @return int
*/
int Matrix2x2::determinant() {
	if ((values[0] * values[3] - values[1] * values[2])<1.0e-6 
		|| (values[0] * values[3] - values[1] * values[2]) == 1.0e-6){
		throw std::overflow_error("Inverse undefined");
	}	
	else{
		return static_cast<int>(values[0] * values[3] - values[1] * values[2]);
	}
}

/*
* calculate the determinant of Matrix
* @param Matrix2x2 m
* @return int
*/
int Matrix2x2::det(Matrix2x2 m){
	return (m.determinant());
}

/*
* function trace
* @param
* @return int
*/
int Matrix2x2::trace() {
	if ((values[0] * values[3] - values[1] * values[2])<1.0e-6 
		|| (values[0] * values[3] - values[1] * values[2]) == 1.0e-6){
		throw std::overflow_error("Inverse undefined");
	}	
	else{
		return static_cast<int>(values[0] + values[3]);
	}
}

/*
* M is symmetric if b = c
* @param
* @return bool
*/
bool Matrix2x2::isSymmetric(){
	if ((values[0] * values[3] - values[1] * values[2])<1.0e-6 
		|| (values[0] * values[3] - values[1] * values[2]) == 1.0e-6){
		throw std::overflow_error("Inverse undefined");
	}
	else{
		return (values[1] == values[2]);
	}
}

/*
* M is similar to M¡ä if det(M) = det(M¡ä) and trace(M) = trace(M¡ä)
* @param
* @return bool
*/
bool Matrix2x2::isSimilar(Matrix2x2 m)
{
	if ((values[0] * values[3] - values[1] * values[2])<1.0e-6 
		|| (values[0] * values[3] - values[1] * values[2]) == 1.0e-6){
		throw std::overflow_error("Inverse undefined");
	}
	else{
		return (det(*this) == det(m) && trace() == m.trace());
	}
}

/*
* print Eigenvalues function
* @param const vector<double>& root
* @param const int i
* @return bool
*/
void Matrix2x2::printEigenvalues(const vector<double>& root, const int i){
	if (i == 1) {	
		cout << "root  " << i << ":  " << static_cast<int>(root[0]) 
		<< " +" << static_cast<int>(root[1]) << "i" << endl;
	}
	else{
		cout << "root  " << i << ":  " << static_cast<int>(root[0]) 
		<< " " << static_cast<int>(root[1]) << "i" << endl;
	}
}

/*
* Overloaded operator () for M() to get the determinant of M
* @param
* @return int: det value of matrix
*/
int Matrix2x2::operator() (){
	return det(*this);
}

/*
* Overloaded operator () for M(int) (in driver for root1 = m1(1) )
* @param const int& a
* @return vector<double>
*/
vector<double> Matrix2x2::operator()(const int& x){
	vector<double > root;
	const auto sq_root = trace()*trace() - 4 * determinant();
	if (x < 1 || x > 3)
		throw invalid_argument("index out of bounds");
	else
	{
		if (x==1)
		{
			if (sq_root >= 0)
			{
				root.reserve(1);
				root.push_back((trace() + sqrt(sq_root)) / 2);
			}
			else
			{
				root.reserve(2);
				root.push_back(static_cast<double>(trace()) / 2);
				root.push_back(sqrt(-sq_root) / 2);
			}
			return root;
		}
		else if (x==2)
		{
			if (sq_root >= 0)
			{
				root.reserve(1);
				root.push_back((trace() - sqrt(sq_root)) / 2);
			}
			else
			{
				root.reserve(2);
				root.push_back(static_cast<double>(trace()) / 2);
				root.push_back(-sqrt(-sq_root) / 2);
			}
			return root;
		}
		else
		{
			return root;
		}			
	}
}

/*
* Overloaded operator[] for M[i]
* @param const int& i
* @return double&
*/
double& Matrix2x2::operator[] (const int& i)
{
	if (i < 0 || i > 4)
		throw invalid_argument("index out of bounds");
	else
		return values[i];
}

/*
* Overloaded operator + for M+M', set as non-member as requested
* @param const Matrix2x2 m1
* @param const Matrix2x2 m2
* @return Matrix2x2
*/
Matrix2x2 operator + (const Matrix2x2 m1, const Matrix2x2 m2){
	Matrix2x2 ret;
	for (int i = 0; i < 4; ++i) {
		ret.values[i] = m1.values[i] + m2.values[i];
	}
	return ret;
}

/*
* Overloaded operator + for M+d, set as non-member as requested
* @param const Matrix2x2 n
* @param const double d
* @return Matrix2x2
*/
Matrix2x2 operator + (const Matrix2x2 n, const double d){
	Matrix2x2 ret;
	for (int i = 0; i < 4; ++i) {
		ret.values[i] = n.values[i] + d;
	}
	return ret;
}

/*
* Overloaded operator + for d+M, set as non-member as requested
* @param const Matrix2x2 n
* @param const double d
* @return Matrix2x2
*/
Matrix2x2 operator + (const double d, const Matrix2x2 n){
	return (n + d);
}

/*
* Overloaded operator + for M-M', set as non-member
* @param const Matrix2x2 m1
* @param const Matrix2x2 rhs
* @return Matrix2x2
*/
Matrix2x2 operator - (const Matrix2x2 m1, const Matrix2x2 m2){
	Matrix2x2 ret;
	for (int i = 0; i < 4; ++i) {
		ret.values[i] = m1.values[i] - m2.values[i];
	}
	return ret;
}

/*
* Overloaded operator + for M-d, set as non-member
* @param const Matrix2x2 n
* @param const double d
* @return Matrix2x2
*/
Matrix2x2 operator - (const Matrix2x2 n, const double d){
	Matrix2x2 ret;
	for (int i = 0; i < 4; ++i) {
		ret.values[i] = n.values[i] - d;
	}
	return ret;
}

/*
* Overloaded operator + for d-M, set as non-member
* @param const Matrix2x2 n
* @param const double d
* @return Matrix2x2
*/
Matrix2x2 operator - (const double d, const Matrix2x2 n){
	Matrix2x2 ret;
	for (int i = 0; i < 4; ++i) {
		ret.values[i] = d - n.values[i];
	}
	return ret;
}

/*
* Overloaded operator * for M*M', set as non-member as requested
* @param const Matrix2x2 m1
* @param const Matrix2x2 m2
* @return Matrix2x2
*/
Matrix2x2 operator * (const Matrix2x2 m1, const Matrix2x2 m2)
{
	Matrix2x2 ret;
	ret.values[0] = m1.values[0] * m2.values[0] + m1.values[1] * m2.values[2];
	ret.values[1] = m1.values[0] * m2.values[1] + m1.values[1] * m2.values[3];
	ret.values[2] = m1.values[2] * m2.values[0] + m1.values[3] * m2.values[2];
	ret.values[3] = m1.values[2] * m2.values[1] + m1.values[3] * m2.values[3];
	return ret;
}

/*
* Overloaded operator * for M*d, set as non-member as requested
* @param const Matrix2x2 n
* @param const double d
* @return Matrix2x2
*/
Matrix2x2 operator * (const Matrix2x2 n, const double d){
	Matrix2x2 ret;
	for (int i = 0; i < 4; ++i) {
		ret.values[i] = n.values[i] * d;
	}
	return ret;
}

/*
* Overloaded operator * for d*M, set as non-member as requested
* @param const Matrix2x2 n
* @param const double d
* @return Matrix2x2
*/
Matrix2x2 operator * (const double d, const Matrix2x2 n){
	return (n * d);
}

/*
* Overloaded operator + for M/M', set as non-member as requested
* @param const Matrix2x2 m1
* @param const Matrix2x2 m2
* @return Matrix2x2
*/
Matrix2x2 operator / (const Matrix2x2 m1, Matrix2x2 m2){
	return m1 * m2.inverse();
}

/*
* Overloaded operator + for M/d, set as non-member as requested
* @param const Matrix2x2 n
* @param const double d
* @return Matrix2x2
*/
Matrix2x2 operator / (const Matrix2x2 n, const double d){
	if (d == 0)
		throw std::overflow_error("Division by zero");
	else{
		Matrix2x2 ret;
		for (int i = 0; i < 4; ++i) {
			ret.values[i] = n.values[i] / d;
		}
		return ret;
	}
}

/*
* Overloaded operator + for d/M, set as non-member as requested
* @param const Matrix2x2 m1
* @param const double d
* @return Matrix2x2
*/
Matrix2x2 operator / (const double d, Matrix2x2 &m){
	Matrix2x2 ret = m.inverse();
	return (ret * d);
}

/*
* Overloaded operator == to tell if 2 matrix have same values
* @param const Matrix2x2 m1
* @param const Matrix2x2 m2
* @return bool
*/
bool operator == (const Matrix2x2 m1, const Matrix2x2 m2)
{
	return 
	abs(m2.values[0] - m1.values[0]) < 1.e-6 
		&& abs(m2.values[1] - m1.values[1]) < 1.e-6
		&& abs(m2.values[2] - m1.values[2]) < 1.e-6 
		&& abs(m2.values[3] - m1.values[3]) < 1.e-6;
}

/*
* Overloaded operator != to tell if 2 matrix have same values
* @param const Matrix2x2 m1
* @param const Matrix2x2 m2
* @return bool
*/
bool operator!=(const Matrix2x2 m1, const Matrix2x2 m2){
	return !(m1 == m2);
}

/*
* print Eigenvalues
* @param const vector<double>&
* @param int
* @return void
*/
void printEigenvalues(const vector<double>& root, int i)
{
	if (i == 1)
		std::cout << "root  " << i << ":  " << static_cast<int>(root[0]) << " +" << static_cast<int>(root[1]) << "i" << endl;
	else
		std::cout << "root  " << i << ":  " << static_cast<int>(root[0]) << " " << static_cast<int>(root[1]) << "i" << endl;
}

/*
* Overloaded input operator for reading Mat2x2 values
* @param Matrix2x2 & n
* @return istream &
*/
istream & operator>>(istream & is, Matrix2x2 & n){
	cout << "To create the following 2x2 matrix :\n";
	cout << "|a  b|\n";
	cout << "|    |\n";
	cout << "|c  d|\n";
	cout << "enter four numbers a, b, c, d, in that order :\n";
	for (auto& value : n.values){
		is >> value;
	}
	return is;
}


/*
* Overloaded output operator for writing Mat2x2 values
* @param const Matrix2x2 &n
* @return out
*/
ostream & operator<< (ostream &os, const Matrix2x2 &n){
	const int a = max(getLength(n.values[0]), getLength(n.values[2]));
	const int b = max(getLength(n.values[1]), getLength(n.values[3]));
	cout << fixed;
	cout << setprecision(2);
	cout << "|" << setw(a) << n.values[0] << " " << setw(b) << n.values[1] << "|" << endl;
	cout << "|" << setw(a) << " " << " " << setw(b) << " " << "|" << endl;
	cout << "|" << setw(a) << n.values[2] << " " << setw(b) << n.values[3] << "|" << endl;
	return os;
}
