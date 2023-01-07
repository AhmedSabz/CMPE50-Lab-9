/*
 * complex.h
 *
 *  Created on: Nov 4, 2022
 *      Author: ahmed
 */

#ifndef COMPLEX_H_
#define COMPLEX_H_

#include <iostream>
using namespace std;
class Complex {
private:
	double imaginary;
	double real;
public:
	Complex(double imaginary, double real);
	Complex(double real_part);
	Complex();
	friend Complex operator-(const Complex &compa, const Complex &compb); //subtracts two imaginary equations
	friend Complex operator+(const Complex &compa, const Complex &compb); //adds two imaginary equations
	friend Complex operator*(const Complex &compa, const Complex &compb); //multiplies two imaginary equations
	friend ostream& operator <<(ostream &ost, const Complex &pol); //outputs the sum,difference, or product of two imagniary equations
	double getReal() const; //obtains the real number
	double getImaginary() const; //gets the imaginary number

};

#endif /* COMPLEX_H_ */
