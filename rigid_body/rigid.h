#include "stdafx.h"
#include <gl/glut.h>
#include <math.h>
#include <vector>
#include <queue>
#include <list>
#include <iostream>
using namespace std;

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))
typedef float number_t;
struct vec{
	number_t x,y,z;
	vec(){}
	vec(number_t x_,number_t y_,number_t z_)
		:x(x_),y(y_),z(z_){}
	vec(const vec &v)
	{ x=v.x;y=v.y;z=v.z;}
	vec& operator=(const vec &v)
	{ x=v.x;y=v.y;z=v.z;return *this;}
	number_t operator*(const vec &b)
	{ return x*b.x+y*b.y+z*b.z; }
	number_t len()
	{ return sqrt(x*x+y*y+z*z); }
	number_t len2()
	{ return x*x+y*y+z*z; }
	vec normalize()
	{
		number_t n=len();
		return vec(x/n,y/n,z/n);
	}
	vec operator+(vec &a)
	{ return vec(x+a.x,y+a.y,z+a.z); }
	vec operator-(vec &a)
	{ return vec(x-a.x,y-a.y,z-a.z); }
};
bool operator==(vec &a,vec &b);
vec operator*(number_t a,vec &v);
vec operator*(vec &v,number_t a);
vec operator/(vec &v,number_t a);
vec cross(vec &a,vec &b);

struct quaternion{
	number_t s,vx,vy,vz;
	quaternion(){}
	quaternion(number_t s_,number_t vx_,number_t vy_,number_t vz_)
		:s(s_),vx(vx_),vy(vy_),vz(vz_){ }
	quaternion(number_t vx_,number_t vy_,number_t vz_)
		:s(0),vx(vx_),vy(vy_),vz(vz_){ }
	quaternion(const vec &a)
		:s(0),vx(a.x),vy(a.y),vz(a.z){ }
	bool operator==(quaternion &b)
	{ return (s==b.s) && (vx==b.vx) && (vy==b.vy) && (vz==b.vz); }

	quaternion operator+(quaternion &b)
	{ return quaternion(s+b.s,vx+b.vx,vy+b.vy,vz+b.vz); }
	quaternion operator-(quaternion &b)
	{ return quaternion(s-b.s,vx-b.vx,vy-b.vy,vz-b.vz); }
	quaternion operator/(number_t b)
	{ return quaternion(s/b,vx/b,vy/b,vz/b); }
	number_t len()
	{ return sqrt(s*s+vx*vx+vy*vy+vz*vz); }
	number_t len2()
	{ return s*s+vx*vx+vy*vy+vz*vz; }
	quaternion normalize();
	vec xyz()
	{ return vec(vx,vy,vz); }
};
quaternion operator*(quaternion &a,quaternion &b);
quaternion operator*(quaternion &a,number_t b);
quaternion operator*(number_t a,quaternion &b);

struct matrix3{
	number_t m[3][3];
	matrix3(){}
	matrix3(const matrix3 &a);
	matrix3(quaternion q);
	number_t det();
	matrix3 operator/(number_t a);
	matrix3& operator=(const matrix3 &a);
	matrix3 inv();
	matrix3 operator*(const matrix3 &a);
	vec     operator*(const vec &a);
	matrix3 transpose();
};



class rigid_sys{//3D

public:
	
	struct body_info{
		vector<vec> pos;//relative position
		vector<number_t> mass;
		number_t total_mass;
		matrix3 Ibody_Inv;
		void generate_Ibody_inv();
	};
	struct dynamic_info{
		vec x,P,L;
		quaternion q;
		
		dynamic_info operator/(number_t a)
		{ dynamic_info di;di.x=x/a;di.P=P/a;di.L=L/a;di.q=q/a;return di;}
		dynamic_info operator+(dynamic_info &a)
		{ dynamic_info di;di.x=x+a.x;di.P=P+a.P;di.L=L+a.L;di.q=q+a.q;return di;}
		dynamic_info operator-(dynamic_info &a)
		{ dynamic_info di;di.x=x-a.x;di.P=P-a.P;di.L=L-a.L;di.q=q-a.q;return di;}
	};
	friend dynamic_info operator*(dynamic_info &b,number_t a)
	{ dynamic_info di;di.x=b.x*a;di.P=b.P*a;di.L=b.L*a;di.q=b.q*a;return di;}
	friend dynamic_info operator*(number_t a,dynamic_info &b)
	{ dynamic_info di;di.x=b.x*a;di.P=b.P*a;di.L=b.L*a;di.q=b.q*a;return di;}

	vector<body_info> bodies;     //one element one body
	vector<dynamic_info> dynamics;//one element one body
	
	virtual bool collision_dectect(vector<dynamic_info> status)=0;
	virtual void collision_handler(vector<dynamic_info> &status)=0;
	virtual vec extern_force(number_t t, int body_no,int subbody_no){ return vec(0,0,0); }
	
	number_t step_max,step_min,err_max,err_min;
	number_t precision;
	
	int iteration_conllision_rollback;

	void update(number_t t0,number_t tn);
	vector<body_info> current_body(vector<dynamic_info> &status);
	vector<body_info> current_body();

	rigid_sys(vector<body_info> b,vector<dynamic_info> dynamic_info0);
	virtual ~rigid_sys();
	

private:
	
	struct RK{
		vector<dynamic_info> RK4;
		vector<dynamic_info> RK5;
	};
	vector<dynamic_info> get_dynamic(number_t t,vector<dynamic_info> start);
	RK delta(vector<dynamic_info> &start,number_t t0,number_t h);

};
class squares:public rigid_sys{

public:
	//for simplicity: pos=position of the vertices,counter-/clockwise
	squares(vector<body_info> b,vector<dynamic_info> dynamic_info0)
		:rigid_sys(b,dynamic_info0)	
	{
		position=b;
		restitution_coef=0.5;
	}
	bool collision_dectect(vector<dynamic_info> status);
	void collision_handler(vector<dynamic_info> &status);
	vec extern_force(number_t t, int body_no,int subbody_no);
	void square_flyin(body_info new_body,vec velo,vec center,number_t sight_size);
	~squares()
	{
		collision_pair.clear();
		position.clear();
	}

private:

	number_t restitution_coef;
	struct pair{
		int a,b;
		pair(int a_,int b_):a(a_),b(b_){}
		void change(){int tmp=a;a=b;b=tmp;}
	};
	vector<pair> collision_pair;
	struct points{
		number_t posi_x;
		number_t posi_y;
		int square_no;
		bool end;
	};

	struct points_compare_x{
		bool operator()(points &a, points &b)
		{ return (a.posi_x>b.posi_x); }
	};
	struct points_compare_y{
		bool operator()(points &a,points &b)
		{ return (a.posi_y>b.posi_y); }
	};
	vector<pair> bound_box_detect();
	bool apply_impulse(vector<dynamic_info> &status_temp,
					   vector<dynamic_info> &status,pair curr,vec p,vec n);

	vector<body_info> position;
	bool squares_collision(int a_no,int b_no,vec center_a,vec center_b);
};
//bound_box
//fly-away
