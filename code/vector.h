/**
 * 3D Vector helper class. 
 * Vectors are allocated with new, but because operator overloading uses a return statement, 
 * freeing memory can't be done using the class destructor. Because of this, you 
 * **must use the free_memory function**
 * failing to do so will result on memory leakage, because it's the only way this class 
 * implements it. If you have a solution for this please refer to
 * www.github.com/jpvanegasc/UsefulPrograms , this project main repository, and check if 
 * it has been solved already, or make a pull request.
 * 
 * Author : Jose Daniel Munoz et al.
 * Modified: Juan Pablo Vanegas
 */
#include <iostream>
#include <cmath>

class Vector3D{
    double *v = NULL;
    //double v[3];
    public:
        Vector3D(double =0.0, double =0.0, double =0.0);
        void free_memory(void);
        void load(double x0, double y0, double z0);
        void show(void);
        /* @return x component*/
        double x(void){return v[0];};
        /* @return y component*/
        double y(void){return v[1];};
        /* @return z component*/
        double z(void){return v[2];};
        /* @return nth component*/
        double & operator[](int n){return v[n];};
        Vector3D operator= (Vector3D v2);
        Vector3D operator+ (Vector3D v2);
        Vector3D operator+=(Vector3D v2);
        Vector3D operator- (Vector3D v2);
        Vector3D operator-=(Vector3D v2);
        Vector3D operator* (double a);
        Vector3D operator*=(double a);
        Vector3D operator/ (double a);
        double operator* (Vector3D v2);
        Vector3D operator^ (Vector3D v2);
        friend Vector3D operator* (double a, Vector3D v1);
        friend double norm2(Vector3D v1);
        friend double norm(Vector3D v1);
};

/* Initialize vector. Defaults to zero */
Vector3D::Vector3D(double x, double y, double z){
    v = new double[3];
    load(x,y,z);
}
void Vector3D::free_memory(void){
    delete[] v;
}
/* Load vector values */
void Vector3D::load(double x0, double y0, double z0){
    v[0] = x0; v[1] = y0; v[2] = z0;
}
void Vector3D::show(void){
    std::cout << "(" << v[0]<< "," << v[1]<< "," << v[2]<< ")" << std::endl;
}
Vector3D Vector3D::operator=(Vector3D v2){
    v[0] = v2.v[0];
    v[1] = v2.v[1];
    v[2] = v2.v[2];
    
    return *this;
}
Vector3D Vector3D::operator+(Vector3D v2){
    v[0] += v2.v[0];
    v[1] += v2.v[1];
    v[2] += v2.v[2];
    
    return *this;
}
Vector3D Vector3D::operator+=(Vector3D v2){
    *this = *this + v2;
    return *this;
}
Vector3D Vector3D::operator-(Vector3D v2){
    v[0] -= v2.v[0];
    v[1] -= v2.v[1];
    v[2] -= v2.v[2];
    
    return *this;
}
Vector3D Vector3D::operator-=(Vector3D v2){
    *this = *this - v2;
    return *this;
}
/* Vector times scalar */
Vector3D Vector3D::operator*(double a){
    v[0] *= a;
    v[1] *= a;
    v[2] *= a;
    
    return *this;
}
/* Vector times scalar */
Vector3D Vector3D::operator*=(double a){
    *this = (*this)*a;
    return *this;
}
/* Vector divided by scalar */
Vector3D Vector3D::operator/(double a){
    double inverse = 1.0/a;
    v[0] *= inverse;
    v[1] *= inverse;
    v[2] *= inverse;
    
    return *this;
}
/* Dot product */
double Vector3D::operator*(Vector3D v2){
    return v[0]*v2.v[0] + v[1]*v2.v[1] +v [2]*v2.v[2];;
}
/* Cross product */
Vector3D Vector3D::operator^(Vector3D v2){
    Vector3D result;
    result.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
    result.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
    result.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];
    *this = result;
    result.free_memory;
    return *this;
}
/* Dont use this. Memory leakage */
Vector3D operator*(double a, Vector3D v1){
    Vector3D result;
    result = v1*a;
    return result;
}
/* @return vector norm squared */
double norm2(Vector3D v1){
    return v1.v[0]*v1.v[0] + v1.v[1]*v1.v[1] + v1.v[2]*v1.v[2];
}
/* @return vector norm */
double norm(Vector3D v1){
    return std::sqrt(norm2(v1));
}
