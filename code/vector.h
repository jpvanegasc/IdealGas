/**
 * 3D Vector helper class
 * Author : Jose Daniel Munoz et al.
 */
#include <iostream>
#include <cmath>

class Vector3D{
    double *v = NULL;
    public:
        Vector3D(double =0.0, double =0.0, double =0.0);
        ~Vector3D();
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
Vector3D::~Vector3D(){
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
    for(int i=0;i<3;i++)
        v[i] = v2.v[i];
    return *this;
}
Vector3D Vector3D::operator+(Vector3D v2){
    Vector3D result;
    for(int i=0;i<3;i++)
        result.v[i] = v[i] + v2.v[i];
    return result;
}
Vector3D Vector3D::operator+=(Vector3D v2){
    *this = *this + v2;
    return *this;
}
/* Vector times scalar */
Vector3D Vector3D::operator*(double a){
    Vector3D result;
    for(int i=0;i<3;i++)
            result.v[i] = a*v[i];
    return result;
}
/* Vector times scalar */
Vector3D Vector3D::operator*=(double a){
    *this = (*this)*a;
    return *this;
}
/* Vector divided by scalar */
Vector3D Vector3D::operator/(double a){
    double inver = 1.0/a;
    Vector3D result;
    for(int i=0;i<3;i++)
            result.v[i] = inver*v[i];
    return result;
}
Vector3D Vector3D::operator-(Vector3D v2){
    return *this + v2*(-1); 
}
Vector3D Vector3D::operator-=(Vector3D v2){
    *this = *this - v2;
    return *this;
}
/* Dot product */
double Vector3D::operator*(Vector3D v2){
    double dot = 0;
    for(int i=0;i<3;i++)
            dot += v[i]*v2.v[i];
    return dot;
}
/* Cross product */
Vector3D Vector3D::operator^(Vector3D v2){
    Vector3D result;
    result.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
    result.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
    result.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];
    return result;
}
Vector3D operator*(double a, Vector3D v1){
    Vector3D result;
    result = v1*a;	
    return result;
}
/* @return vector norm squared */
double norm2(Vector3D v1){
    double norm = 0;
    for(int i=0;i<3;i++)
            norm += v1.v[i]*v1.v[i];
    return norm;
}
/* @return vector norm */
double norm(Vector3D v1){
    return std::sqrt(norm2(v1));
}
