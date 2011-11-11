/*
 * intersection_point.h
 *
 *  Created on: 23-Oct-2011
 *      Author: prateek
 */

#ifndef INTERSECTION_POINT_H_
#define INTERSECTION_POINT_H_
#include "vector.h"
#include "ray.h"
class intersection_point {
public:
	ray r;
	vector v;
	int object; //0=no intersection,1=triangle;2=sphere
	int index; // to represent the index in the triangles or spheres array
	intersection_point();
	void print_intersection_point();
	virtual ~intersection_point();
};

#endif /* INTERSECTION_POINT_H_ */
