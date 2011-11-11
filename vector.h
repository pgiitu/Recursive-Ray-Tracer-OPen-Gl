/*
 * vector.h
 *
 *  Created on: 23-Oct-2011
 *      Author: prateek
 */

#ifndef VECTOR_H_
#define VECTOR_H_

class vector {

public:
	double x;
	double y;
	double z;

	vector();

	vector(double a,double b,double c);

	//vector(const vector &a);

	double dot_product(vector b);

	vector cross_product(vector b);

	vector normalize();
	vector add(vector a);
	vector subtract(vector b);
	vector multiply_scalar(double t);
	double magnitude();

	virtual ~vector();
};

#endif /* VECTOR_H_ */
