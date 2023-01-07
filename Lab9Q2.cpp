//============================================================================
// Name        : Lab9Q2.cpp
// Author      :Ahmed Sabzwari
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include <iostream>
#include <cmath>
#include "Complex.h"
#include <cmath>
using namespace std;

Complex::Complex() {
	real = 0.0;
	imaginary = 0.0;
}
Complex::Complex(double real_part) {
	imaginary = real_part;
	real = imaginary;
}
Complex::Complex(double real, double imaginary) {
	this->real = real;
	this->imaginary = imaginary;
}
Complex operator-(const Complex &compa, const Complex &compb) {
	double realDiff = compa.real - compb.real; //the real numbers of the first and second equations are subtracted
	double imaginaryDiff = compa.imaginary - compb.imaginary; //the imaginary numbers of the first and second equations are subtracted
	Complex ahmed(realDiff, imaginaryDiff);
	return ahmed;
}
Complex operator+(const Complex &compa, const Complex &compb) {
	double realDiff = compa.real + compb.real; //the real numbers of the first and second equations are added
	double imaginaryDiff = compa.imaginary + compb.imaginary; //the imagniary numbers of the first and seconc equations are added
	Complex ahmed(realDiff, imaginaryDiff);
	return ahmed;
}
Complex operator*(const Complex &compa, const Complex &compb) {
	double realNum = (compa.real * compb.real)
			- (compa.imaginary * compb.imaginary); // (a+bi)*(c+di)= (a*c -b*d) + (a*d+b*c)*i
	double imaginaryPart1 = (compa.real * compb.imaginary);
	double imaginaryPart2 = (compa.imaginary * compb.real);
	Complex ahmed(realNum, imaginaryPart1 + imaginaryPart2);
	return ahmed;
}

ostream& operator <<(ostream &ost, const Complex &comp) {
	ost << comp.getReal() << " + " << comp.getImaginary() << "i";
	return ost;
}
double Complex::getImaginary() const {
	return imaginary;
}
double Complex::getReal() const {
	return real;
}
int main() {
	Complex x(2, -3);
	Complex y(-7, 5);
	cout << x * y;
}
//output
//1 + 31i
