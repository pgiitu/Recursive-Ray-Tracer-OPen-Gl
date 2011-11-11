/*
 * ray.cpp
 *
 *  Created on: 23-Oct-2011
 *      Author: prateek
 */

#include "ray.h"
#include "vector.h"

ray::ray()
{
	start.x=0.0;
	start.y=0.0;
	start.z=0.0;
	end.x=0.0;
	end.y=0.0;
	end.z=0.0;
}
ray::ray(vector a,vector b)
{
	start.x=a.x;
	start.y=a.y;
	start.z=a.z;
	end.x=b.x;
	end.y=b.y;
	end.z=b.z;
}

ray::~ray() {
	// TODO Auto-generated destructor stub
}
