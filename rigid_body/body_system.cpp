#include "stdafx.h"
#include "rigid.h"
rigid_sys::rigid_sys(vector<body_info> b,vector<dynamic_info> dynamic_info0)
{
	bodies=b;
	dynamics=dynamic_info0;
	step_max=2e-1;
	step_min=1e-2;
	err_max=1e-2;
	err_min=1e-4;
	precision=1e-4;
	iteration_conllision_rollback=2;
}

rigid_sys::~rigid_sys()
{
	bodies.clear();
	dynamics.clear();
}
rigid_sys::RK rigid_sys::delta(vector<dynamic_info> &start,number_t t0,number_t h)
{
	int i,N=start.size();

	vector<dynamic_info> temp=start;

	vector<dynamic_info> k1=get_dynamic(t0,temp);
	for(i=0;i<N;i++)
		k1[i]=k1[i]*h;

	
	for(i=0;i<N;i++)
		temp[i]=start[i]+k1[i]/4;
	vector<dynamic_info> k2=get_dynamic(t0+h/4,temp);
	for(i=0;i<N;i++)
		k2[i]=k2[i]*h;


	for(i=0;i<N;i++)
		temp[i]=start[i]+3.0/32*k1[i]+9.0/32*k2[i];
	vector<dynamic_info> k3=get_dynamic(t0+3.0/8*h,temp);
	for(i=0;i<N;i++)
		k3[i]=k3[i]*h;



	for(i=0;i<N;i++)
		temp[i]=start[i]+1932.0/2197*k1[i]-7200.0/2197*k2[i]+7296.0/2197*k3[i];
	vector<dynamic_info> k4=get_dynamic(t0+12.0/13*h,temp);
	for(i=0;i<N;i++)
		k4[i]=k4[i]*h;


	for(i=0;i<N;i++)
		temp[i]=start[i]+439.0/216*k1[i]-8*k2[i]+3680.0/513*k3[i]-845.0/4104*k4[i];
	vector<dynamic_info> k5=get_dynamic(t0+h,temp);
	for(i=0;i<N;i++)
		k5[i]=k5[i]*h;


	for(i=0;i<N;i++)
		temp[i]=start[i]-8.0/27*k1[i]+2*k2[i]-3544.0/2565*k3[i]+1859.0/4104*k4[i]-11.0/40*k5[i];
	vector<dynamic_info> k6=get_dynamic(t0+h/2,temp);
	for(i=0;i<N;i++)
		k6[i]=k6[i]*h;
	

	
	RK rk;
	for(i=0;i<N;i++)
	{
		temp[i]=25.0/216*k1.at(i)+1408.0/2565*k3.at(i)+2197.0/4104*k4.at(i)-1.0/5*k5.at(i);
		rk.RK4.push_back(temp[i]);
		
		temp[i]=16.0/135*k1.at(i)+6656.0/12825*k3.at(i)+28561.0/56430*k4.at(i)-9.0/50*k5.at(i)+2.0/55*k6.at(i);
		rk.RK5.push_back(temp[i]);
	}

	return rk;
}
void rigid_sys::update(number_t t0,number_t tn)
{

	int i;
	int N=dynamics.size();
	number_t h=step_max;
	number_t err=0;
	
	vector<dynamic_info> y=dynamics,z=dynamics,y_temp=dynamics,z_temp=dynamics,y_back;
	RK temp;
	//adaptive
	while(t0<=tn)
	{
		while(1)
		{
			err=0;
			temp=delta(y,t0,h);
			for(i=0;i<N;i++)
			{
				y_temp[i]=y[i]+temp.RK4[i];
				z_temp[i]=z[i]+temp.RK5[i];
				
				dynamic_info diff=y_temp[i]-z_temp[i];
				err+=diff.L.len2()+diff.P.len2()+diff.q.len2()+diff.x.len2();
			}
			err=sqrt(err);
			//y_temp.at(N)=y.at(N);
			if(err<=err_max && err>=err_min)
			{
				y_back=y;
				y=y_temp;
				z=z_temp;
				break;
			}

			if(err>=err_max)
			{
				if(h==step_min)
				{
					cout<<"not accurate at t="<<t0+h<<endl;
					y_back=y;
					y=y_temp;
					z=z_temp;
					break;
				}
				h=max(step_min,h/2.0);
				continue;
			}

			if(err<err_min)
			{
				h=min(step_max,h*2);
				y_back=y;
				y=y_temp;
				z=z_temp;
				break;
			}
		}
		if(collision_dectect(y))
		{
			//bisection to find the time
			//int itr=2;
			number_t dtl=0,dtr=h;
			vector<dynamic_info> rk4;
			for(int i=0;i<iteration_conllision_rollback;i++)
			{
				//y_back=delta(y_back,t0,(dtl+dtr)/2).RK4;
				rk4=delta(y_back,t0,(dtl+dtr)/2).RK4;
				for(int j=0;j<rk4.size();j++)
					rk4[j]=rk4[j]+y_back[j];

				if(collision_dectect(rk4))
					dtr=(dtl+dtr)/2;
				else
					dtl=(dtl+dtr)/2;

			}
			rk4=delta(y_back,t0,dtr).RK4;
			for(int j=0;j<rk4.size();j++)
			{
				y[j]=rk4[j]+y_back[j];
			}
			t0+=dtr;
			
			collision_handler(y);
			
			//start over
			z=y;
			y_temp=y;
			z_temp=y;
		}
		else
		{
			t0+=h;
		}
		for(i=0;i<y.size();i++)
		{
			if(y[i].L.len()<precision)
				y[i].L=vec(0,0,0);
			if(y[i].P.len()<precision)
				y[i].P=vec(0,0,0);
		}
	}
	dynamics=y;
		

}
vector<rigid_sys::body_info> rigid_sys::current_body(vector<dynamic_info> &status)
{
	vector<body_info> cb;
	for(int i=0;i<bodies.size();i++)
	{
		cb.push_back(bodies[i]);
		for(int j=0;j<bodies[i].pos.size();j++)
			cb[i].pos[j]=matrix3(status[i].q)*bodies[i].pos[j];
	}
	return cb;
}
vector<rigid_sys::body_info> rigid_sys::current_body()
{
	return current_body(dynamics);
}
vector<rigid_sys::dynamic_info> rigid_sys::get_dynamic(number_t t,vector<dynamic_info> start)
{
	vector<rigid_sys::dynamic_info> dynamic_dot;
	for(int i=0;i<start.size();i++)
	{
		rigid_sys::dynamic_info change;
		rigid_sys::dynamic_info origin;

		origin=start[i];

		//dx
		change.x=origin.P/bodies[i].total_mass;
	
		//dq
		matrix3 R(origin.q);
		matrix3 I_inv=R*bodies[i].Ibody_Inv*R.transpose();
		vec omega=I_inv*origin.L;

		change.q=quaternion(omega)*origin.q;

		//dP
		vec total_force(0,0,0);
		for(int j=0;j<bodies[i].mass.size();j++)
			total_force=total_force+extern_force(t,i,j);
		change.P=total_force;

		//dL
		vec total_torque(0,0,0);
		for(int j=0;j<bodies[i].mass.size();j++)
			total_torque=total_torque+cross(matrix3(origin.q)*bodies[i].pos[j]-origin.x,extern_force(t,i,j));
		change.L=total_torque;
		
		dynamic_dot.push_back(change);
	}
	return dynamic_dot;
}
void rigid_sys::body_info::generate_Ibody_inv()
{
	int i;
	number_t sum=0;
	for(i=0;i<pos.size();i++)
		sum+=mass[i]*(pos[i].y*pos[i].y+pos[i].z*pos[i].z);
	Ibody_Inv.m[0][0]=sum;
	
	sum=0;
	for(i=0;i<pos.size();i++)
		sum+=-mass[i]*pos[i].x*pos[i].y;
	Ibody_Inv.m[0][1]=sum;
	Ibody_Inv.m[1][0]=sum;

	sum=0;
	for(i=0;i<pos.size();i++)
		sum+=-mass[i]*pos[i].x*pos[i].z;
	Ibody_Inv.m[0][2]=sum;
	Ibody_Inv.m[2][0]=sum;

	sum=0;
	for(i=0;i<pos.size();i++)
		sum+=mass[i]*(pos[i].x*pos[i].x+pos[i].z*pos[i].z);
	Ibody_Inv.m[1][1]=sum;

	sum=0;
	for(i=0;i<pos.size();i++)
		sum+=-mass[i]*pos[i].y*pos[i].z;
	Ibody_Inv.m[1][2]=sum;
	Ibody_Inv.m[2][1]=sum;

	sum=0;
	for(i=0;i<pos.size();i++)
		sum+=mass[i]*(pos[i].x*pos[i].x+pos[i].y*pos[i].y);
	Ibody_Inv.m[2][2]=sum;

	Ibody_Inv=Ibody_Inv.inv();
}

number_t line_equation(vec a1,vec a2,vec x)
{
	return (x.x-a2.x)/(a1.x-a2.x)-(x.y-a2.y)/(a1.y-a2.y);
}
bool squares::squares_collision(int a_no,int b_no,vec center_a,vec center_b)
{
	int i,j;
	for(i=0;i<4;i++)
	{
		bool ref=(line_equation(position[a_no].pos[i]+center_a,
								position[a_no].pos[(i+1)%4]+center_a,position[a_no].pos[(i+2)%4]+center_a)>=0);
		for(j=0;j<4;j++)
		{
			if(ref==(line_equation(position[a_no].pos[i]+center_a,
								position[a_no].pos[(i+1)%4]+center_a,position[b_no].pos[j]+center_b)>=0))
				break;
		}
		if(j==4)
			return false;
	}
	for(i=0;i<4;i++)
	{
		bool ref=(line_equation(position[b_no].pos[i]+center_b,
								position[b_no].pos[(i+1)%4]+center_b,position[b_no].pos[(i+2)%4]+center_b)>=0);
		for(j=0;j<4;j++)
		{
			if(ref==(line_equation(position[b_no].pos[i]+center_b,
								position[b_no].pos[(i+1)%4]+center_b,position[a_no].pos[j]+center_a)>=0))
				break;
		}
		if(j==4)
			return false;
	}
	return true;
}
number_t max4(number_t a1,number_t a2,number_t a3,number_t a4)
{
	return max(a1,max(a2,max(a3,a4)));
}
number_t min4(number_t a1,number_t a2,number_t a3,number_t a4)
{
	return min(a1,min(a2,min(a3,a4)));
}
vector<squares::pair> squares::bound_box_detect()
{
	vector<pair> collision_susp;
	priority_queue<points,vector<points>,points_compare_x> posi_x;
	vector<points> boundbox_info;
	for(int i=0;i<position.size();i++)
	{
		points sq;
		sq.posi_x=max4(position[i].pos[0].x,position[i].pos[1].x,position[i].pos[2].x,
				  position[i].pos[3].x)+dynamics[i].x.x;
		sq.posi_y=max4(position[i].pos[0].y,position[i].pos[1].y,position[i].pos[2].y,
					   position[i].pos[3].y)+dynamics[i].x.y;
		sq.square_no=i;
		sq.end=true;
		posi_x.push(sq);
		boundbox_info.push_back(sq);

		sq.posi_x=min4(position[i].pos[0].x,position[i].pos[1].x,position[i].pos[2].x,
					   position[i].pos[3].x)+dynamics[i].x.x;
		sq.posi_y=min4(position[i].pos[0].y,position[i].pos[1].y,position[i].pos[2].y,
					   position[i].pos[3].y)+dynamics[i].x.y;
		sq.square_no=i;
		sq.end=false;
		posi_x.push(sq);
		boundbox_info.push_back(sq);
		
	}

	list<int> conflict_list;
	list<int>::iterator conflict_list_itr=conflict_list.begin();
	priority_queue<points,vector<points>,points_compare_x> posi_y;

	while(posi_x.size())
	{
		if(posi_x.top().end==false)
			conflict_list.push_back(posi_x.top().square_no);

		if(posi_x.top().end==true)
		{
			conflict_list_itr=conflict_list.begin();
			for(int i=0;i<conflict_list.size();i++)
			{
				posi_y.push(boundbox_info[*conflict_list_itr*2]);
				posi_y.push(boundbox_info[*conflict_list_itr*2+1]);
				conflict_list_itr++;
			}
			list<int> conflict_list_y;
			list<int>::iterator conflict_list_y_itr=conflict_list_y.begin();
			pair pair_temp(posi_x.top().square_no,0);
			while(posi_y.size())
			{
				if(posi_y.top().end==false)
					conflict_list_y.push_back(posi_y.top().square_no);

				if(posi_y.top().end==true)
				{
					conflict_list_y.push_back(posi_y.top().square_no);
					while(conflict_list_y_itr!=conflict_list_y.end())
					{
						pair_temp.b=*conflict_list_y_itr;
						collision_susp.push_back(pair_temp);
						conflict_list_y_itr++;
					}
				}

				posi_y.pop();
			}
			conflict_list_itr=conflict_list.begin();
			conflict_list.remove(posi_x.top().square_no);

		}

		posi_x.pop();
	}

	return collision_susp;
}
bool squares::collision_dectect(vector<dynamic_info> status)
{
	position=current_body(status);

	vector<pair> vp;//=bound_box_detect();
	for(int i=0;i<(int)vp.size();i++)
		cout<<vp[i].a<<"   "<<vp[i].b<<endl;

	
	bool collision=false;
	collision_pair.clear();
	for(int i=0;i<(int)status.size();i++)
	{
		for(int j=i+1;j<(int)status.size();j++)
		{
			if(squares_collision(i,j,status[i].x,status[j].x))
			{
				collision_pair.push_back(pair(i,j));
				collision=true;
			}
		}
	}
	return collision;
}
vec squares::extern_force(number_t t, int body_no,int subbody_no)
{
	matrix3 R(dynamics[body_no].q);
	matrix3 I_inv=R*position[body_no].Ibody_Inv*R.transpose();
	vec omega=I_inv*dynamics[body_no].L;
	vec p_dot=dynamics[body_no].P/position[body_no].total_mass+cross(omega,position[body_no].pos[subbody_no]);
	if(p_dot.len()>1e-6)
		return (-0.0001)*p_dot.normalize();
	return vec(0,0,0);
}
bool squares::apply_impulse(vector<dynamic_info> &status_temp,
					   vector<dynamic_info> &status,pair curr,vec p,vec n)
{
	if(n==vec(0,0,0))
		return false;
	//set impulse
	n=n.normalize();
	vec ra=p-status[curr.a].x;
	vec rb=p-status[curr.b].x;

	matrix3 R(status[curr.a].q);
	matrix3 I_inv_a=R*position[curr.a].Ibody_Inv*R.transpose();
	vec omega_a=I_inv_a*status[curr.a].L;
	vec p_dot_a=status[curr.a].P/position[curr.a].total_mass+cross(omega_a,ra);
		
	R=matrix3(status[curr.b].q);
	matrix3 I_inv_b=R*position[curr.b].Ibody_Inv*R.transpose();
	vec omega_b=I_inv_b*status[curr.b].L;
	vec p_dot_b=status[curr.b].P/position[curr.b].total_mass+cross(omega_b,rb);

	number_t v_rel=n*(p_dot_a-p_dot_b);
	if(v_rel>-1e-4)
		return false;

	number_t j=-(1+restitution_coef)*v_rel/(1/position[curr.a].total_mass+1/position[curr.b].total_mass
										+n*cross((I_inv_a*cross(ra,n)),ra)
										+n*cross((I_inv_b*cross(rb,n)),rb));
	cout<<j<<endl;
	vec P=j*n;
	status_temp[curr.a].P=status_temp[curr.a].P+P;
	status_temp[curr.b].P=status_temp[curr.b].P-P;

	status_temp[curr.a].L=status_temp[curr.a].L+cross(ra,P);
	status_temp[curr.b].L=status_temp[curr.b].L-cross(rb,P);
	return true;

}
void squares::collision_handler(vector<dynamic_info> &status)
{
	int j,k;
	vector<dynamic_info> status_temp=status;
	for(int i=0;i<collision_pair.size();i++)
	{
		//get n and p: further collision detection
		vec n(0,0,0),p,p_prev;
		pair curr=collision_pair[i];
		bool has_added=false;
		bool finish=false;
		for(j=0;j<4;j++)
		{//vertex a edge b
			vec point=status[curr.a].x+position[curr.a].pos[j];
			for(k=0;k<4;k++)
			{
				//vec edge=position[curr.b].posi[(k+1)%4]-position[curr.b].posi[k];
				vec diff1=point-position[curr.b].pos[k]-status[curr.b].x;
				vec diff2=point-position[curr.b].pos[(k+1)%4]-status[curr.b].x;
				number_t temp=abs(diff1*diff2+diff1.len()*diff2.len());
				if(abs(diff1*diff2+diff1.len()*diff2.len())<1e-4)
				{

					n=position[curr.b].pos[k]-position[curr.b].pos[(k+3)%4];
					p=point;
					has_added=apply_impulse(status_temp,status,curr,p,n);
					finish=has_added;
					
					if(finish)
					{
						cout<<"impulse"<<endl;
						break;
						
					}
				}
			}
			if(finish)
				break;
		}

		p_prev=p;
		finish=false;
		curr.change();
		for(j=0;j<4;j++)
		{//vertex b edge a
			vec point=status[curr.a].x+position[curr.a].pos[j];
			for(k=0;k<4;k++)
			{
				//vec edge=position[curr.a].posi[(k+1)%4]-position[curr.a].posi[k];
				vec diff1=point-position[curr.b].pos[k]-status[curr.b].x;
				vec diff2=point-position[curr.b].pos[(k+1)%4]-status[curr.b].x;
				if(abs(diff1*diff2+diff1.len()*diff2.len())<1e-4)
				{
					n=position[curr.b].pos[k]-position[curr.b].pos[(k+3)%4];
					p=point;
					
					if((p-p_prev).len()>1e-4)
					{
						
						finish=apply_impulse(status_temp,status,curr,p,n);


						if(has_added)
						{
							status_temp[curr.a].P=(status_temp[curr.a].P+status[curr.a].P)/2;
							status_temp[curr.b].P=(status_temp[curr.b].P+status[curr.b].P)/2;
							status_temp[curr.a].L=(status_temp[curr.a].L+status[curr.a].L)/2;
							status_temp[curr.b].L=(status_temp[curr.b].L+status[curr.b].L)/2;
						}
						
						//curr.change();
						//status=status_temp;
						if(finish)
						{
							cout<<"impulse"<<endl;
							break;
						}
					}
				}
			}
			if(finish)
				break;
		}


	}
	status=status_temp;//optimizable
}
void squares::square_flyin(body_info new_body,vec velo,vec center,number_t sight_size=1.0)
{
	bodies.push_back(new_body);
	dynamic_info new_sq;
	new_sq.P=velo*new_body.total_mass;
	int sep0=1;
	new_sq.x=center-velo.normalize()*sight_size*sep0;
	bool have_colision=true;
	int i;
	while(have_colision)
	{

		for(i=0;i<dynamics.size();i++)
		{
			if((new_sq.x-dynamics[i].x).len()<(new_body.pos[0]-new_body.pos[1]).len()
				+(position[i].pos[0]-position[i].pos[1]).len())
			{
				new_sq.x=new_sq.x-velo.normalize()*sight_size*sep0;
				break;
			}
		}
		if(i==dynamics.size())
			have_colision=false;
		
	}
	new_sq.q=quaternion(1,0,0,0);
	new_sq.L=vec(0,0,0);
	dynamics.push_back(new_sq);
	position.push_back(new_body);
}