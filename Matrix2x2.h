/*
* @ClassName: Matrix2x2.h
* @Function: header file of Matrix2x2.cpp, using methods and overloaded operators to implement the Mat2x2 ADT described in assignment 3
* @author: Lin Li 
* @version: 1.0
* @date: Jul 14, 2018
*/

#pragma once
#include <string>
#include <iostream>
#include <vector>
using namespace std;

class Matrix2x2 {
public:
	//Default constructor: sets all four entries to zero.
	Matrix2x2();
	
	//Normal constructor: accepts initial values for all four entries as parameters.
	Matrix2x2(double a, double b, double c, double d);

	//Default copy constructor
	Matrix2x2(const Matrix2x2& n);
	
	//default assignment operator
	Matrix2x2 & operator =(const Matrix2x2& n);

	//default destructor
	virtual ~Matrix2x2();

	//over loaded operator + and - for positive M / negative M
	Matrix2x2 operator + ();
	Matrix2x2 operator - ();
	

	//Compound Arithmetic Assignment Operator Overloads
	Matrix2x2& operator += (const Matrix2x2& n);	// M+=M'
	Matrix2x2& operator += (const double& d);		// M+=d
	Matrix2x2& operator -= (const Matrix2x2& n);	// M-=M'
	Matrix2x2& operator -= (const double& d);		// M-=d
	Matrix2x2& operator *= (const Matrix2x2& n);	// M*=M'
	Matrix2x2& operator *= (const double& d);		// M*=d
	Matrix2x2& operator /= (Matrix2x2& n);		// M/=M'
	Matrix2x2& operator /= (const double& d);		// M/=d

	//overloaded operator ++ for M++
	Matrix2x2 operator ++ (int);	
	//overloaded operator ++ for ++M
	Matrix2x2& operator ++ ();
	//overloaded operator -- for M--
	Matrix2x2 operator -- (int);
	//overloaded operator -- for --M
	Matrix2x2& operator -- ();

	//function Inversion
	Matrix2x2 inverse();
	//function Transpose
	Matrix2x2 transpose();
	//function Determinant
	int determinant();
	//calculate the determinant of Matrix
	int det(const Matrix2x2 m);
	//function trace
	int trace();
	//M is symmetric if b = c
	bool isSymmetric();
	//M is similar to M′ if det(M) = det(M′) and trace(M) = trace(M′)
	bool isSimilar(Matrix2x2 m);
	//print Eigenvalues function
	static void printEigenvalues(const vector<double>& root, const int i);

	//Overloaded operator () for m1()
	int operator() ();

	//Overloaded operator () for root1 = m1(1)
	vector<double> operator()(const int& i);

	//Overloaded operator[] for m9[0]
	double& operator [] (const int& i);

	//Overloaded operator + for M+M', M+d, d+M, set as non-member
	friend Matrix2x2 operator + (const Matrix2x2 m1, const Matrix2x2 m2);//M+M'
	friend Matrix2x2 operator + (const Matrix2x2 n, const double d);//M+d	
	friend Matrix2x2 operator + (const double d, const Matrix2x2 n);//d+M

	//Overloaded operator + for M-M', M-d, d-M, set as non-member
	friend Matrix2x2 operator - (const Matrix2x2 lhs, const Matrix2x2 rhs);//M-M'
	friend Matrix2x2 operator - (const Matrix2x2 n, const double d);//M-d	
	friend Matrix2x2 operator - (const double d, const Matrix2x2 n);//d-M
	
	//Overloaded operator * for M*M', M*d, d*M, set as non-member
	friend Matrix2x2 operator * (const Matrix2x2 m1, const Matrix2x2 m2);//M*M'
	friend Matrix2x2 operator * (const Matrix2x2 n, const double d);//M*d	
	friend Matrix2x2 operator * (const double d, const Matrix2x2 n);//d*M'

	//Overloaded operator + for M/M', M/d, d/M, set as non-member																//binary Division
	friend Matrix2x2 operator / (const Matrix2x2 m1, Matrix2x2 m2);//M/M'
	friend Matrix2x2 operator / (const Matrix2x2 n, const double d);//M/d	
	friend Matrix2x2 operator / (const double d, Matrix2x2 &n);//d/M

	//Overloaded operator == to tell if 2 matrix have same values
	friend bool operator == (const Matrix2x2 m1, const Matrix2x2 m2);

	//Overloaded operator != to tell if 2 matrix have same values
	friend bool operator != (const Matrix2x2 m1, const Matrix2x2 m2);

	
	//Overloaded input operator for reading Mat2x2 values
	friend std::istream & operator >> (std::istream &is, Matrix2x2 &n);
	
	//Overloaded output operator for writing Mat2x2 values
	friend std::ostream & operator << (std::ostream &os, const Matrix2x2 &n);

private:
	vector <double> values;
};

//print Eigenvalues
void printEigenvalues(const vector<double>& root, int i);
