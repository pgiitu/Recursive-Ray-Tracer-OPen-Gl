/*
Author: Frank Pfenning
*/

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <string.h>
#include <math.h>
#include "pic.h"
#include <iostream>
#include <stdio.h>
#include "ray.h"
#include "vector.h"
#include "intersection_point.h"

using namespace std;
#define MAX_TRIANGLES 10000
#define MAX_SPHERES 10000
#define MAX_LIGHTS 10000
#define no_intersection 99999999.23372637
#define sqrt3 1.732050808
#define max_intersection 99999999

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 75.0

#define PLANE_DIST 1.0

unsigned char buffer[HEIGHT][WIDTH][3];

enum shadow { NOT_BLOCKED = 1,
             BLOCKED=0 };

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct ray_color
{
	vector color;
	ray r;
};

struct plane
{
	  // for equation of a plane ax+by+cz=d
	double a;
	double b;
	double c;
	double d;
};

typedef struct _Triangle
{
  struct Vertex v[3];
  // for equation of a plane ax+by+cz=d
  struct plane p;
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;
int Depth=1;

bool flag=true;
vector RT_shading(intersection_point p,ray r,int depth);
vector Recursive_RT(ray r,int depth);



/*
 * function to store the coefficients of a plane equation ax+by+cz=d for a
 * triangle with index i in the array of traingles
 */
void calculate_plane_coefficient(int i)
{
	vector v0(triangles[i].v[0].position[0],triangles[i].v[0].position[1],triangles[i].v[0].position[2]);
	vector v1(triangles[i].v[1].position[0],triangles[i].v[1].position[1],triangles[i].v[1].position[2]);
	vector v2(triangles[i].v[2].position[0],triangles[i].v[2].position[1],triangles[i].v[2].position[2]);
	vector t1=v1.subtract(v0);
	vector t2=v2.subtract(v0);
	vector n=t2.cross_product(t1);
	n=n.normalize();

	triangles[i].p.a=n.x;
	triangles[i].p.b=n.y;
	triangles[i].p.c=n.z;
	triangles[i].p.d=-v0.dot_product(n);
}

/*
 * Given a ray function calculates the intersection with a triangle
 * with index i in the array of triangles
 */

double intersection_with_triangle(ray r,int i)
{
	//first we will find intersection with plane
	vector d=(r.end).subtract(r.start);  //direction vector

	vector v0(triangles[i].v[0].position[0],triangles[i].v[0].position[1],triangles[i].v[0].position[2]);
	vector v1(triangles[i].v[1].position[0],triangles[i].v[1].position[1],triangles[i].v[1].position[2]);
	vector v2(triangles[i].v[2].position[0],triangles[i].v[2].position[1],triangles[i].v[2].position[2]);


	vector n(triangles[i].p.a,triangles[i].p.b,triangles[i].p.c);

	double a=n.dot_product(r.start);
	double b=n.dot_product(d);
	if(b==0.0)
	{
		return no_intersection;
	}
	//here i have missed a condition need to take care of it
	double t=-(triangles[i].p.d+a)/b;  // t gives us the point of intersection

	if(t>0.0)  //if t is greater than one then only we will check whether the point lies in the triangle or not
	{
		vector u=v1.subtract(v0);
		vector v=v2.subtract(v0);
		vector intersection_point(r.start.x+d.x*t,r.start.y+d.y*t,r.start.z+d.z*t);
		//P_I=s*u+t1*v+P0
		vector diff=intersection_point.subtract(v0);
		//diff=s*u+t1*v
/*
		double t1=(diff.y*u.x-diff.x)/(v.y*u.x-v.x);
		double s=(diff.x-t1*v.x)/u.x;
		*/
		double den=v.dot_product(v)*u.dot_product(u)-pow(v.dot_product(u),2);
		double s=(u.dot_product(u)*diff.dot_product(v)-u.dot_product(v)*diff.dot_product(u))/den;
		double t1=(v.dot_product(v)*diff.dot_product(u)-v.dot_product(u)*diff.dot_product(v))/den;

			if(!(s<0 || s>1.0 || t1<0 || t1>1.0 || ((s+t1)>1.0)))
			{
				//the ray intersect the triangle
				return t;
			}
	}
	return no_intersection;
}


/*
 * function to return the interpolated value of the various coefficients and normal for a triangle with index i in the array of triangles
 * index=1  Interpolate Normal
 * index=2  Interpolate diffusion constant
 * index=3  Interpolate Specular constant
 */
vector interpolate_for_triangle(int i,vector v,int index) //v is the point in the triagular plane
{
	vector n0,n1,n2;
	//vector k0,k1,k2;
	switch(index)
	{
		case 1:
		{
			vector k0(triangles[i].v[0].normal[0],triangles[i].v[0].normal[1],triangles[i].v[0].normal[2]);
			vector k1(triangles[i].v[1].normal[0],triangles[i].v[1].normal[1],triangles[i].v[1].normal[2]);
			vector k2(triangles[i].v[2].normal[0],triangles[i].v[2].normal[1],triangles[i].v[2].normal[2]);
			n0=k0;n1=k1;n2=k2;
			break;
		}
		case 2:
		{
			vector k0(triangles[i].v[0].color_diffuse[0],triangles[i].v[0].color_diffuse[1],triangles[i].v[0].color_diffuse[2]);
			vector k1(triangles[i].v[1].color_diffuse[0],triangles[i].v[1].color_diffuse[1],triangles[i].v[1].color_diffuse[2]);
			vector k2(triangles[i].v[2].color_diffuse[0],triangles[i].v[2].color_diffuse[1],triangles[i].v[2].color_diffuse[2]);
			n0=k0;n1=k1;n2=k2;
			break;
		}
		case 3:
		{
			vector k0(triangles[i].v[0].color_specular[0],triangles[i].v[0].color_specular[1],triangles[i].v[0].color_specular[2]);
			vector k1(triangles[i].v[1].color_specular[0],triangles[i].v[1].color_specular[1],triangles[i].v[1].color_specular[2]);
			vector k2(triangles[i].v[2].color_specular[0],triangles[i].v[2].color_specular[1],triangles[i].v[2].color_specular[2]);
			n0=k0;n1=k1;n2=k2;
			break;
		}
	}

	vector v0(triangles[i].v[0].position[0],triangles[i].v[0].position[1],triangles[i].v[0].position[2]);
	vector v1(triangles[i].v[1].position[0],triangles[i].v[1].position[1],triangles[i].v[1].position[2]);
	vector v2(triangles[i].v[2].position[0],triangles[i].v[2].position[1],triangles[i].v[2].position[2]);


	ray r1(v1,v2);  // r1=v1+t1*(v2-v1)
	ray r2(v0,v);	// r2=v0+t2*(v-v0)
	vector d1=v2.subtract(v1);
	vector d2=v.subtract(v0);
/*
 *   vector in=r1.intersect
 *   for intersection r1=r2 ie  v1+t1*(v2-v1)=v0+t2*(v-v0)
 *   v1+t1*d1=v0+t2*d2
 *   cross multiplying by vector d2 on both sides to get t1
 *   cross multiplying by vector d1 on both sides to get t2
 *   the value of t1 must be <1.0
 *	 the value of t2 must be >1.0
 * */
	double t1=(((v0.cross_product(d2)).subtract(v1.cross_product(d2))).magnitude())/((d1.cross_product(d2)).magnitude());
	double t2=(((v1.cross_product(d1)).subtract(v0.cross_product(d1))).magnitude())/((d2.cross_product(d1)).magnitude());
	vector np=n2.multiply_scalar(t1).add(n1.multiply_scalar(1.0-t1));
	vector n=n0.multiply_scalar((t2-1.0)/t2).add(np.multiply_scalar(1.0/t2));
	return n;

}

/*
 * function to interpolate shininess for a triangle in with index i in the array of triangles
 */

double interpolate_shininess(int i,vector v)
{
	vector v0(triangles[i].v[0].position[0],triangles[i].v[0].position[1],triangles[i].v[0].position[2]);
	vector v1(triangles[i].v[1].position[0],triangles[i].v[1].position[1],triangles[i].v[1].position[2]);
	vector v2(triangles[i].v[2].position[0],triangles[i].v[2].position[1],triangles[i].v[2].position[2]);

	double n0=triangles[i].v[0].shininess;
	double n1=triangles[i].v[1].shininess;
	double n2=triangles[i].v[2].shininess;

	ray r1(v1,v2);  // r1=v1+t1*(v2-v1)
	ray r2(v0,v);	// r2=v0+t2*(v-v0)
	vector d1=v2.subtract(v1);
	vector d2=v.subtract(v0);

	double t1=(((v0.cross_product(d2)).subtract(v1.cross_product(d2))).magnitude())/((d1.cross_product(d2)).magnitude());
	double t2=(((v1.cross_product(d1)).subtract(v0.cross_product(d1))).magnitude())/((d2.cross_product(d1)).magnitude());
	//cout<<"t1 is  "<<t1<<"    "<<"t2 is"<<t2<<endl;
	double np=n2*t1 + n1*(1.0-t1);
	double n=n0*(t2-1.0)/t2+np*(1.0/t2);
	return n;

}





/*
 *  this function computes the intersection of a ray with a sphere with index i in the array of spheres
 *
 */
double find_intersection_with_sphere(ray r,int i)
{
	vector d=r.end.subtract(r.start);
//	d=d.normalize();
	double A=d.x*d.x+d.y*d.y+d.z*d.z;
	double B=2*(d.x*(r.start.x-spheres[i].position[0]) + d.y*(r.start.y-spheres[i].position[1]) + d.z*(r.start.z-spheres[i].position[2]));
	
	double C=(r.start.x-spheres[i].position[0])*(r.start.x-spheres[i].position[0]) + (r.start.y-spheres[i].position[1])*(r.start.y-spheres[i].position[1]) + (r.start.z-spheres[i].position[2])*(r.start.z-spheres[i].position[2])-(spheres[i].radius*spheres[i].radius);
	//cout<<"A is  "<<A<<"   B is  "<<B<<"  C is  "<<C<<endl;
	double D=B*B-4*A*C;
	if(D<0)
	{
		return no_intersection;
	}
	else
	{
		if(D==0)
		{
			double k23=-B/(2*A);
			if(k23>1.0)
				return k23;
		}
		else
		{
			double t1=((-B-sqrt(D))/(2.0*A));
//			double t2=((-B+sqrt(D))/(2.0*A));
			if(t1>0.0)
			{
				return t1;
			}
/*
			else if(t2>1.0)
			{
//				return t2;
				return 1.0;
			}
*/
		}
	}
	return no_intersection ;
}


/*
 * function to find the nearest intersection point for a ray. The function iterates through all the
 * objects in the scene ie the array of triangles and spheres
 */

intersection_point find_nearest_intersection(ray r)
{
	int i;
	intersection_point I_P;
	I_P.index=0;
	I_P.r=r;
	double min_i=max_intersection;
	vector v0(r.start.x,r.start.y,r.start.z);
	vector v1(r.end.x,r.end.y,r.end.z);
	vector d=v1.subtract(v0);
//	d=d.normalize();
	double temp;


	for(i=0;i<num_triangles;i++)
	{
		temp=intersection_with_triangle(r,i);
		if(temp<min_i && temp!=no_intersection)
		{
			//cout<<"the intersection parameter is  "<<temp<<"\n";
			min_i=temp;
			vector t(v0.x+temp*d.x,v0.y+temp*d.y,v0.z+temp*d.z);
			I_P.v=t;
			I_P.object=1;
			I_P.index=i;
			//I_P.print_intersection_point();
		}
	}

	for(i=0;i<num_spheres;i++)
	{
		temp=find_intersection_with_sphere(r,i);
//		cout<<"temp is: "<<temp<<endl;
		if(temp<min_i && temp!=no_intersection)
		{
			//cout<<"the intersection parameter is  "<<temp<<"\n";
			min_i=temp;
			vector t(v0.x+temp*d.x,v0.y+temp*d.y,v0.z+temp*d.z);
			I_P.v=t;
			I_P.object=2;
			I_P.index=i;
		}
	}
	return I_P;

}



/*
 * function to calculate the normal vector on a point on a sphere with index i
 */
vector normal_for_sphere(int i,vector v)
{
	vector center(spheres[i].position[0],spheres[i].position[1],spheres[i].position[2]);
	vector n=(v.subtract(center)).multiply_scalar(1.0/spheres[i].radius);
	n=n.normalize();
	return n;
}

/*
 * function to find distance between two points
 */
double find_distance(vector v1,vector v2)
{
	vector diff=v1.subtract(v2);
	return v1.magnitude();
}


/*
 * Function to calculate shadow for a intersection point
 * This function returns whether the ray to light intersects any object or not.
 * a return value of 0 represents blocked
 * a return value of 1 represents not blocked
 */

int find_intersection_light(ray r,intersection_point p)
{
	int i;
	double temp;

	for(i=0;i<num_triangles;i++)
	{
		if(p.object==1 && p.index==i)
		{
			//not to find intersection with the same triangle
		}
		else
		{
			temp=intersection_with_triangle(r,i);
			if(temp!=no_intersection && temp<1.0 && temp>0.0)  //checking for blocking
			{
				return 0;
			}
		}
	}

	for(i=0;i<num_spheres;i++)
	{
		temp=find_intersection_with_sphere(r,i);
		if(temp!=no_intersection && temp<1.0 && temp>0.0)
		{
			return 0;
		}
	}

	return 1;

}

/*
 * function to calculate diffuse color at an intersection point p
 * l:- light vector
 * Id= intensity of light
 * The function calculated on the basis of the intersection point(whether it is a part of a sphere or triangle)
 */
vector diffuse_color(intersection_point p,vector l,vector Id)
{
	if(p.object==1)
	{
		vector diffuse_color;
		vector n=interpolate_for_triangle(p.index,p.v,1);
		n=n.normalize();
		vector Kd=interpolate_for_triangle(p.index,p.v,2);
//		double dis=find_distance(l,p.v);
		double cos_theta=l.dot_product(n);
		double fac=cos_theta;
		if(cos_theta<0.0)
			cos_theta=0.0;
		diffuse_color.x=Id.x*Kd.x*fac;
		diffuse_color.y=Id.y*Kd.y*fac;
		diffuse_color.z=Id.z*Kd.z*fac;
		return diffuse_color;
	}
	else if(p.object==2)
	{

		vector n=normal_for_sphere(p.index,p.v);
		n=n.normalize();

		//for diffusion
		vector diffuse_color;
		double cos_theta=l.dot_product(n);
		if(cos_theta<0.0)
			cos_theta=0.0;
		double fac=cos_theta;
		vector Kd(spheres[p.index].color_diffuse[0],spheres[p.index].color_diffuse[1],spheres[p.index].color_diffuse[2]);
		diffuse_color.x=Id.x*Kd.x*fac;
		diffuse_color.y=Id.y*Kd.y*fac;
		diffuse_color.z=Id.z*Kd.z*fac;
		return diffuse_color;
	}
}

/*
 * function to calculate specular color
 * l:- light vector
 * Id= intensity of light
 * The function calculated on the basis of the intersection point(whether it is a part of a sphere or triangle)
 */
vector specular_color(intersection_point p,vector l,vector Id,ray r)
{
//	struct ray_color rc;
	if(p.object==1)
	{

		//for specular reflection
		vector n=interpolate_for_triangle(p.index,p.v,1);
		vector specular_color;
		vector proj=n.multiply_scalar(n.dot_product(l));
		vector reflection=(proj.multiply_scalar(2.0)).subtract(l);
		reflection=reflection.normalize();


		vector V=r.start.subtract(r.end);
		V=V.normalize();
		vector Ks=interpolate_for_triangle(p.index,p.v,3);

		//double fac1=s*pow(V.dot_product(reflection),shine)/(dis*dis);
		double shine=interpolate_shininess(p.index,p.v);

		double cos_theta=V.dot_product(reflection);
		if(cos_theta<0.0)
			cos_theta=0.0;
		double fac1=pow(cos_theta,shine);
		specular_color.x=Id.x*Ks.x*fac1;
		specular_color.y=Id.y*Ks.y*fac1;
		specular_color.z=Id.z*Ks.z*fac1;
		return specular_color;
	}
	else if(p.object==2)
	{
		//for specular reflection
		vector n=normal_for_sphere(p.index,p.v);
		n=n.normalize();
		vector specular_color;
		vector proj=n.multiply_scalar(n.dot_product(l));   //proj=ncos(0)
		vector reflection=(proj.multiply_scalar(2.0)).subtract(l);
		reflection=reflection.normalize();
//			vector V=(r.end).subtract(p.v);


		vector V=r.start.subtract(r.end);
		V=V.normalize();
		vector Ks(spheres[p.index].color_specular[0],spheres[p.index].color_specular[1],spheres[p.index].color_specular[2]);
//			double fac1=s*pow(V.dot_product(reflection),spheres[p.index].shininess)/(dis*dis);
		double cos_theta=V.dot_product(reflection);
		double fac1=0.0;
		if(cos_theta<0.0)
		{
			cos_theta=0.0;
		}
		else
		{
			fac1=pow(cos_theta,spheres[p.index].shininess);
		}
		specular_color.x=Id.x*Ks.x*fac1;
		specular_color.y=Id.y*Ks.y*fac1;
		specular_color.z=Id.z*Ks.z*fac1;
		return specular_color;
	}
}

/*
 * function to calculate the ambient color at an intersection point
 */

vector ambient_color(intersection_point p)
{
	vector I_a(ambient_light[0],ambient_light[1],ambient_light[2]);

	if(p.object==1)  //for triangle
	  {
			vector K_a=interpolate_for_triangle(p.index,p.v,2);
			vector ambient_color(I_a.x*K_a.x,I_a.y*K_a.y,I_a.z*K_a.z);
			return ambient_color;
	  }
	  else if(p.object==2)  //for sphere
	  {
			//for ambient light
		  vector K_a(spheres[p.index].color_diffuse[0],spheres[p.index].color_diffuse[1],spheres[p.index].color_diffuse[2]);
		  vector ambient_color(I_a.x*K_a.x,I_a.y*K_a.y,I_a.z*K_a.z);
		  return ambient_color;
	  }

}

/*
 * function to calculate the refleted ray for an intersection point
 */
ray find_reflective_ray(intersection_point p)
{
	vector l=p.r.start.subtract(p.r.end);
	l=l.normalize();
	vector n;
	if(p.object==1)  //for triangle
	{
		n=interpolate_for_triangle(p.index,p.v,1);
		n=n.normalize();

	}
	else if(p.object==2)  //for sphere
	{

		n=normal_for_sphere(p.index,p.v);
		n=n.normalize();
	}
	vector proj=n.multiply_scalar(n.dot_product(l));
	vector reflection=(proj.multiply_scalar(2.0)).subtract(l);
	ray r(p.v,p.v.add(reflection));
	return r;
}

/*
 * Function to calculate shininess at a intersection point p
 */
double find_shininess(intersection_point p)
{
	double shine=0.0;
	if(p.object==1)  //for triangle
	{
		shine=interpolate_shininess(p.index,p.v);
	}
	else if(p.object==2)  //for sphere
	{
		shine=spheres[p.index].shininess;
	}
	return shine;
}

/*
 * function to calculate the specular constant at an intersection point as requied for recursive ray tracing
 */
vector find_specular_coeffcient(intersection_point p)
{
	if(p.object==1)  //for triangle
	{
		vector ks=interpolate_for_triangle(p.index,p.v,3);
		return ks;
	}
	else if(p.object==2)  //for sphere
	{
		vector v(spheres[p.index].color_specular[0],spheres[p.index].color_specular[1],spheres[p.index].color_specular[2]);
		return v;
	}

}

/*
 * function which return the color for a ray and the depth of the ray. this function uses all the above mentioned functions and gives
 * us the total color at  the point of intersection
 */
vector Recursive_RT(ray r,int depth)
{
	vector color;
	intersection_point p=find_nearest_intersection(r);
	if(p.object>0)
	{
		  	  color=RT_shading(p,r,depth);
	}
	else
	{
			color.x=0.0;
			color.y=0.0;
			color.z=0.0;
	}

	  return color;
}

/*
 * Function to calculate the shading at an intersection point p , for a given ray.
 * This function calls recursively RT_shading function
 */
vector RT_shading(intersection_point p,ray r,int depth)
{

	vector color=ambient_color(p);
	//ray reflected_ray;

	for(int i=0;i<num_lights;i++)
	{
		vector l1(lights[i].position[0],lights[i].position[1],lights[i].position[2]);
		ray l_ray(p.v,l1);

		int s=find_intersection_light(l_ray,p);  //for shadow
		vector Id(lights[i].color[0],lights[i].color[1],lights[i].color[2]);
		vector l=l1.subtract(p.v);
		l=l.normalize();

		vector d_c= diffuse_color(p,l,Id);
		vector s_c= specular_color(p,l,Id,r);
		color=color.add((d_c.add(s_c)).multiply_scalar(s));

	}

	ray reflected_ray=find_reflective_ray(p);

	vector ks=find_specular_coeffcient(p);
	//double shine=find_shininess(p);
	vector c(0.0,0.0,0.0);

	if(depth<Depth && ks.x>0.0 && ks.y>0.0 && ks.z>0.0)
	{
		c=Recursive_RT(reflected_ray,depth+1);
		c.x=c.x*pow(ks.x,depth);
		c.y=c.x*pow(ks.y,depth);
		c.z=c.x*pow(ks.z,depth);
		//color=color.add(c.multiply_scalar(shine));
		color.x=color.x*(1-ks.x);
		color.y=color.y*(1-ks.y);
		color.z=color.z*(1-ks.z);

	}

	color=color.add(c);

	return color;

}


void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

/*
 * function which initially calculates the dimensions of the screen given the plane distance and then generate rays
 * the camera/eye point and the color to be shown at that pixel(as returned by other functions) is stored in the buffer
 * and then lated displayed in one go.
 *
 */
void fill_buffer()
{
	double height=PLANE_DIST*tan((fov/2)*M_PI/180);
	double width=PLANE_DIST*tan((fov/2)*M_PI/180);
	double increment_height=2*height/HEIGHT;
	double increment_width=2*width/WIDTH;
	vector V0,V1;
	V0.x=V0.y=V0.z=0.0;
	V1.z=-PLANE_DIST;

	unsigned int x,y;
	  //simple output
    for(x=0;x < WIDTH;x++)
	  {
        glBegin(GL_POINTS);
	      V1.x = -width+x*increment_width;
	    for(y=0;y < HEIGHT;y++)
	    {

	          vector color; //it denotes the total color
		      V1.y = -height+y*increment_height;
		      ray r(V0,V1);

		      color=Recursive_RT(r,1);
		      if(color.x==0.0 && color.y==0.0 && color.z==0.0)
		      {
		    	  color.x=1.0;
		    	  color.y=1.0;
		    	  color.z=1.0;
		      }
			  if(color.x>1.0)
				  color.x=1.0;
			  if(color.y>1.0)
				  color.y=1.0;
			  if(color.z>1.0)
				  color.z=1.0;

			  buffer[y][x][0]=color.x*255;
			  buffer[y][x][1]=color.y*255;
			  buffer[y][x][2]=color.z*255;

		      plot_pixel_display(x,y,buffer[y][x][0]%256,buffer[y][x][1]%256,buffer[y][x][2]%256);
	    }
	    glEnd();
	    glFlush();
	  }
    printf("Done!\n");
//    fflush(stdout);

}


//MODIFY THIS FUNCTION
void draw_scene()
{

	fill_buffer();
/*
  unsigned int x,y;
  //simple output
  for(x=0;x < WIDTH;x++)
  {
	//glPointSize(5.0);
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
      //plot_pixel(x,y,x%256,y%256,(x+y)%256);
        plot_pixel_display(x,y,buffer[y][x][0]%256,buffer[y][x][1]%256,buffer[y][x][2]%256);
    }
    glEnd();
    glFlush();

  }
*/
//  printf("Done!\n");
 // fflush(stdout);
//  printf("Done!\n");

}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glPointSize(1.0);
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{// this function computes the intersection of a ray with a triangle with index i in the array of triangles

  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;
/*
  in = pic_alloc(640, 480, 3, NULL);
*/
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
/*
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);
*/
}


void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	  calculate_plane_coefficient(num_triangles-1);
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  //glOrtho(-1*WIDTH,WIDTH,-HEIGHT,HEIGHT,1,-1);
  glOrtho(0, WIDTH, 0, HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
//	  fill_buffer();
/*
      if(mode == MODE_JPEG)
	save_jpg();
*/
    }
  once=1;
}


int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> <Recursive(0)/NonRecursive(1) Ray Tracing> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      if(atoi(argv[2])==0)
      {
    	  Depth=1;
      }
      else if(atoi(argv[2])==1)
      {
    	  Depth=3;
      }
      else
      {
    	  Depth=1;
      }

      filename = argv[3];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);
//  display();


  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
