/*
 * Lab9q1.cpp
 *
 *  Created on: Nov 3, 2022
 *      Author: ahmed
 */
#include <iostream>
#include <cmath>
#include "polynomial.h"
#include <cmath>
using namespace std;
Polynomial::Polynomial() {
	degree = 10;
	coeff = new double[degree + 1];
	for (int i = 0; i <= degree; i++) {
		coeff[i] = 0;
	}
}
Polynomial::Polynomial(int degr) {
	degree = degr;
	coeff = new double[degree + 1];
	for (int i = 0; i <= degree; i++) {
		coeff[i] = 0;
	}
}
Polynomial::Polynomial(const Polynomial &poly) {
	degree = poly.get_degree();
	coeff = new double[degree + 1];
	for (int i = 0; i <= degree; i++) {
		coeff[i] = poly.get_coeff(i);
	}
}
Polynomial::Polynomial(double cf[], int deg) {
	degree = deg;
	coeff = new double[degree + 1];
	for (int i = 0; i <= degree; i++) {
		coeff[i] = cf[i];
	}
}
Polynomial::~Polynomial() {
	delete[] coeff;
}
int Polynomial::get_degree() const {
	return degree;
}
double Polynomial::get_coeff(int deg) const {
	if (degree < deg) {
		return 0;
		// The input degree is larger than the polynomial degree
	}
	return coeff[deg];
}
void Polynomial::set_coeff(int degr, double val) {
	if (degree < degr) {
		cout << "Degree exceeded." << endl;
		return;
	}
	coeff[degr] = val;
}
// Evaluate the polynomial
double Polynomial::evaluate(double val) { // plugs in val into the polynomial to get a value in type double
	double answer = 0;
	for (int i = 0; i <= degree; i++) {
		answer = answer + (coeff[i] * pow(val, i));
	}
	return answer;

}
// Assignment operator
void Polynomial::operator =(const Polynomial &poly) {
	if (this == &poly) {
		// Copy to itself. Nothing to be done.
		return;
	}
}
// Overloaded operator +
Polynomial operator+(const Polynomial &pola, const Polynomial &polb) {
	int d = max(pola.degree, polb.degree);
	double *co;
	co = new double[d + 1];
	for (int i = 0; i <= pola.degree; i++) { //dynamic array first gets the coefficients of the first polynomial
		co[i] = pola.get_coeff(i);
	}
	for (int j = 0; j <= polb.degree; j++) { //the coeficcients of the second polynomial are added to the dynamic array
		co[j] += polb.get_coeff(j);
	}

	Polynomial ahmed(co, d);
	return ahmed;

}
Polynomial operator+(int constant, const Polynomial &polb) { //the constant just affects the dynamic array's first index
	int degr = polb.get_degree();
	double *c;
	c = new double[degr + 1];
	for (int i = 0; i <= degr; i++) { //the first index is the constant. the first index is of 0 degree
		if (i == 0) {
			c[i] = polb.get_coeff(i) + constant;
			continue;
		}
		c[i] = polb.get_coeff(i);
	}
	Polynomial a(c, degr);
	return a;
}
Polynomial operator+(const Polynomial &polb, int constant) { //same as the function above
	int degr = polb.get_degree();
	double *c;
	c = new double[degr + 1];
	for (int i = 0; i <= degr; i++) {
		if (i == 0) {
			c[i] = polb.get_coeff(i) + constant;
			continue;
		}
		c[i] = polb.get_coeff(i);
	}
	Polynomial a(c, degr);
	return a;
}
Polynomial operator-(const Polynomial &polb, int constant) {
	int degr = polb.get_degree();
	double *c;
	c = new double[degr + 1];
	for (int i = 0; i <= degr; i++) {
		if (i == 0) { //the coefficient of the zero degree is subtracted by the constant
			c[i] = polb.get_coeff(i) - constant;
			continue;
		}
		c[i] = polb.get_coeff(i);
	}
	Polynomial a(c, degr);
	return a;
}
Polynomial operator-(int constant, const Polynomial &polb) {
	int degr = polb.get_degree();
	double *c;
	c = new double[degr + 1];
	for (int i = 0; i <= degr; i++) {
		if (i == 0) {
			c[i] = constant - polb.get_coeff(i); //the constant subtracts the coefficient of the 0 degree
			continue;
		}
		c[i] = -1 * polb.get_coeff(i); //the coefficients of the polynomial switch sign
	}
	Polynomial a(c, degr);
	return a;
}
Polynomial operator-(const Polynomial &pola, const Polynomial &polb) {
	int d = max(pola.degree, polb.degree);
	double *co;
	co = new double[d + 1];
	for (int i = 0; i <= pola.degree; i++) { //the coefficients of the first polynomial is stored in the dynamic array co
		co[i] = pola.get_coeff(i);
	}
	for (int j = 0; j <= polb.degree; j++) {
		co[j] -= polb.get_coeff(j); //the coefficient values of co is subtracted by the coefficient values in the second polynomial
	}

	Polynomial ahmed(co, d);
	return ahmed;

}
Polynomial operator*(int constant, const Polynomial &polb) {
	int degr = polb.get_degree();
	double *c;
	c = new double[degr + 1];
	for (int i = 0; i <= degr; i++) {
		c[i] = constant * polb.get_coeff(i); //each coeficcient is multiplied by a factor of constant
	}
	Polynomial a(c, degr);
	return a;

}
Polynomial operator*(const Polynomial &polb, int constant) { //same as the *function above
	int degr = polb.get_degree();
	double *c;
	c = new double[degr + 1];
	for (int i = 0; i <= degr; i++) {
		c[i] = constant * polb.get_coeff(i);
	}
	Polynomial a(c, degr);
	return a;

}
Polynomial operator*(const Polynomial &pola, const Polynomial &polb) {
	double *co;
	co = new double[pola.degree + polb.degree];
	for (int i = 0; i <= pola.degree + polb.degree; i++) {
		co[i] = 0; //initialize the values of the dynamic array to 0
	}
	for (int aDegree = 0; aDegree <= pola.degree; aDegree++) {
		for (int bDegree = 0; bDegree <= polb.degree; bDegree++) {
			co[aDegree + bDegree] += pola.get_coeff(aDegree) //each coefficient is multiplied using
			* polb.get_coeff(bDegree);
		}
	}
	Polynomial a(co, pola.degree + polb.degree);
	return a;
}
// Overloaded << operator
ostream& operator <<(ostream &ost, const Polynomial &pol) { //5 + x + 3x^2
	//double a [pol.get_degree()];
	for (int i = 0; i < pol.get_degree(); i++) {
		if (i == 0) {//this coefficient is of the zero degree with no factor of x
			ost << pol.get_coeff(i) << " + ";
			continue;
		}
		if (i == 1) {	//this is the first degree
			ost << pol.get_coeff(i) << "x + ";
			continue;
		}
		ost << pol.get_coeff(i) << "x^" << i << " + ";

	}
	ost << pol.get_coeff(pol.get_degree()) << "x^" << pol.get_degree();
	return ost;

}
int main() {
	double c[3] = { 2, 3, 4 };
	double f[4] = { 1, 2, 3, 4 };
	Polynomial ahmed(c, 2);
	Polynomial ahsan(f, 3);
	cout << ahsan + ahmed;
}
//output
//3 + 5x + 7x^2 + 4x^3
