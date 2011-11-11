/*
 * intersection_point.cpp
 *
 *  Created on: 23-Oct-2011
 *      Author: prateek
 */

#include <iostream>
#include <stdio.h>
#include "intersection_point.h"
using namespace std;

intersection_point::intersection_point() {
	// TODO Auto-generated constructor stub
	r.end=vector();
	r.start=vector();

	v.x=0.0;
	v.y=0.0;
	v.z=0.0;
	object=0;
	index=-1;

}
void intersection_point::print_intersection_point()
{
	cout<<"Found an intersection point \n";
	cout<<"the coordinates are :- \n";
	cout<<"x="<<v.x<<"   y="<<v.y<<"    z="<<v.z<<"\n\n";
}
intersection_point::~intersection_point() {
	// TODO Auto-generated destructor stub
}
