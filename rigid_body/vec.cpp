#include "stdafx.h"
#include "rigid.h"

bool operator==(vec &a,vec &b)
{
	return (a.x==b.x) && (a.y==b.y) && (a.z==b.z);
}
vec operator*(number_t a,vec &v)
{
	return vec(a*v.x,a*v.y,a*v.z);
}
vec operator*(vec &v,number_t a)
{
	return vec(a*v.x,a*v.y,a*v.z);
}
vec operator/(vec &v,number_t a)
{
	return vec(v.x/a,v.y/a,v.z/a);
}
vec cross(vec &a,vec &b)
{
	return vec(a.y*b.z-a.z*b.y,-a.x*b.z+a.z*b.x,a.x*b.y-a.y*b.x);
}
quaternion quaternion::normalize()
{
	number_t n=len();
	return quaternion(s/n,vx/n,vy/n,vz/n);
}
quaternion operator*(quaternion &a,quaternion &b)
{
	return quaternion(a.s*b.s-a.vx*b.vx-a.vy*b.vy-a.vz*b.vz,
					  a.s*b.vx+a.vx*b.s+a.vy*b.vz-a.vz*b.vy,
					  a.s*b.vy-a.vx*b.vz+a.vy*b.s+a.vz*b.vx,
					  a.s*b.vz+a.vx*b.vy-a.vy*b.vx+a.vz*b.s);
}
quaternion operator*(quaternion &a,number_t b)
{ 
	return quaternion(a.s*b,a.vx*b,a.vy*b,a.vz*b); 
}
quaternion operator*(number_t a,quaternion &b)
{ 
	return quaternion(b.s*a,b.vx*a,b.vy*a,b.vz*a); 
}
matrix3::matrix3(quaternion q)
{
	number_t s=q.s,vx=q.vx,vy=q.vy,vz=q.vz;
	m[0][0]=1-2*vy*vy-2*vz*vz;
	m[0][1]=2*vx*vy-2*s*vz;
	m[0][2]=2*vx*vz+2*s*vy;

	m[1][0]=2*vx*vy+2*s*vz;
	m[1][1]=1-2*vx*vx-2*vz*vz;
	m[1][2]=2*vy*vz-2*s*vx;

	m[2][0]=2*vx*vz-2*s*vy;
	m[2][1]=2*vy*vz+2*s*vx;
	m[2][2]=1-2*vx*vx-2*vy*vy;
}
number_t matrix3::det()
{
	return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
		  -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
		  +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
}
matrix3::matrix3(const matrix3 &a)
{
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			m[i][j]=a.m[i][j];
}
matrix3& matrix3::operator=(const matrix3 &a)
{
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			m[i][j]=a.m[i][j];
	return *this;
}
matrix3 matrix3::operator/(number_t a)
{
	matrix3 ma;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			ma.m[i][j]=m[i][j]/a;
	return ma;
}
matrix3 matrix3::inv()
{
	matrix3 m_inv;
	m_inv.m[0][0]=m[1][1]*m[2][2]-m[1][2]*m[2][1];
	m_inv.m[0][1]=m[0][2]*m[2][1]-m[0][1]*m[2][2];
	m_inv.m[0][2]=m[0][1]*m[1][2]-m[0][2]*m[1][1];

	m_inv.m[1][0]=m[1][2]*m[2][0]-m[1][0]*m[2][2];
	m_inv.m[1][1]=m[0][0]*m[2][2]-m[0][2]*m[2][0];
	m_inv.m[1][2]=m[0][2]*m[1][0]-m[0][0]*m[1][2];

	m_inv.m[2][0]=m[1][0]*m[2][1]-m[1][1]*m[2][0];
	m_inv.m[2][1]=m[0][1]*m[2][0]-m[0][0]*m[2][1];
	m_inv.m[2][2]=m[0][0]*m[1][1]-m[0][1]*m[1][0];

	return m_inv/this->det();
}
matrix3 matrix3::operator*(const matrix3 &a)
{
	matrix3 product;
	number_t sum=0;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			sum=0;
			for(int k=0;k<3;k++)
				sum+=m[i][k]*a.m[k][j];
			product.m[i][j]=sum;
		}
	}
	return product;
}
matrix3 matrix3::transpose()
{
	matrix3 trans;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			trans.m[i][j]=m[j][i];
	return trans;
}
vec matrix3::operator*(const vec &a)
{
	return vec(a.x*m[0][0]+a.y*m[0][1]+a.z*m[0][2],
			   a.x*m[1][0]+a.y*m[1][1]+a.z*m[1][2],
			   a.x*m[2][0]+a.y*m[2][1]+a.z*m[2][2]);

}