/*
 * vector.cpp
 *
 *  Created on: 23-Oct-2011
 *      Author: prateek
 */

#include "vector.h"
#include <math.h>


vector::vector()
{
	x=0.0;
	y=0.0;
	z=0.0;
}

vector::vector(double a,double b,double c)
{
	x=a;
	y=b;
	z=c;
}
/*
vector::vector(vector a)
{
	x=a.x;
	y=a.y;
	z=a.z;
}
*/

double vector::dot_product(vector a)
{
	double dotp=a.x*x+a.y*y+a.z*z;
	return dotp;
}
vector vector::cross_product(vector b)
{

	double i=y*b.z-z*b.y;
	double j=z*b.x-x*b.z;
	double k=x*b.y-y*b.x;
	vector v(i,j,k);
	return v;

}
vector vector::normalize()
{
	double fac=sqrt(x*x+y*y+z*z);
	vector v(x/fac,y/fac,z/fac);
	return v;
}

vector vector::add(vector a)
{
	vector v(a.x+x,a.y+y,a.z+z);
	return v;
}

vector vector::subtract(vector a)
{
	vector v(x-a.x,y-a.y,z-a.z);
	return v;
}

double vector::magnitude()
{
	double d=sqrt(x*x+y*y+z*z);
	return d;
}

vector vector::multiply_scalar(double t)
{
	vector v(t*x,t*y,t*z);
	return v;
}
vector::~vector() {
	// TODO Auto-generated destructor stub
}
