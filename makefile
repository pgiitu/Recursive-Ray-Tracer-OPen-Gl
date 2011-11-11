all: vector.o ray.o intersection_point.o assign3.o
	g++ -lglut -lGLU -lGL -ltiff -ljpeg vector.o ray.o intersection_point.o assign3.o -o trace_ray

vector: vector.cpp
	g++ -c vector.cpp
ray: ray.cpp
	g++ -c ray.cpp
intersection_point: intersection_point.cpp
	g++ -c intersection_point.cpp
assign3: assign3.cpp
	g++ -c assign3.cpp
clean:
	rm -rf *.o trace_ray
