/*
 * ray.h
 *
 *  Created on: 23-Oct-2011
 *      Author: prateek
 */

#ifndef RAY_H_
#define RAY_H_
#include "vector.h"

class ray {
public:
	vector start;
	vector end;

	ray();

	ray(vector a,vector b);
	virtual ~ray();
};

#endif /* RAY_H_ */
