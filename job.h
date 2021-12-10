//#define _CRT_SECURE_NO_WARNINGS
#ifndef __JOB__
#define __JOB__
//#define _CRT_NONSTDC_NO_WARNINGS
#include<iostream>
#include<vector>
#include<list>
#include"constant.h"

using namespace std;


extern int Job_Num;
extern int Iterater_NUM;
extern int Instances_Num;
extern int Times;
extern int N;
extern int M_s1;
extern int M_s2;
extern int M_s3;
extern int M_s4;
extern int M_s5;
extern int L_1;
extern int L_2;
extern int L_3;
extern int L_4;
extern int L_5;
extern int M_1;
extern int M_2;
extern int M_3;
extern int M_4;
extern int M_5;
extern double v_1;
extern double v_2;
extern double v_3;
extern double v_4;
extern double v_5;

template<typename T1, typename T2>
class two_point
{
public:
	T1 x;
	T2 y;
	two_point(const two_point& p)
	{
		x = p.x;
		y = p.y;
	}
	two_point(T1 x = 0.0, T2 y = 0.0) :x(x), y(y) {}

	~two_point() = default;
};
class two
{
public:
	int A;
	int B;
	two(int a = 0, int b = 0) :A(a), B(b)
	{}
};
class NO_I//num-obj(int)
{
public:
	NO_I(int n = 0, int o = 0) :num(n), O(o)
	{

	}
	int num;
	int O;
};
class NO_B//num-obj(bool)
{
public:
	NO_B(int n = 0, bool o = false) :num(n), O(o)
	{

	}
	int num;
	bool O;
};
class NO//num-obj(double)
{
public:
	NO(int n = 0, double o = 0.0) :num(n), O(o)
	{

	}
	int num;
	double O;
};
class NOO//num-obj-obj
{
public:
	NOO(int n=0, double a=0.0, double b=0.0) :num(n),A(a),B(b)
	{

	}
	int num;
	double A;
	double B;
};
class Two
{
public:
	Two() :Max_Lmax(0.0), Max_TC(0.0), Lmax(0.0), TC(0.0), T_Lmax(0.0), T_TC(0.0) {}
	double Max_Lmax;
	double Max_TC;
	double Lmax;
	double TC;
	double T_Lmax;
	double T_TC;
};
class job
{
private:
	int J_Inum;
	int J_Ip;
	int J_Is;
	int J_Ir;
	int J_Id;
	double T_j;
	double C_j;
	double rand_key;
public:
	double selection_feature;
	bool flag;
	job(int I_num = 0, int I_p = 0, int I_s = 0, int I_r = 0, int I_d = 0, double T = 0.0, double C = 0.0, double key = 0.0, double feature = 0.0, bool f = false) :
		J_Inum(I_num), J_Ip(I_p), J_Is(I_s), J_Ir(I_r), J_Id(I_d),
		rand_key(key), T_j(T), C_j(C), selection_feature(feature), flag(f)
	{

	}

	~job() = default;
	void set_num(const int& I_num) { J_Inum = I_num; }
	int get_num() const { return J_Inum; }			
	void set_p(const int& I_p) { J_Ip = I_p; }
	int get_p() const { return J_Ip; };		
	void set_s(const int& I_s) { J_Is = I_s; }
	int get_s() const { return J_Is; };			
	void set_r(const int& I_r) { J_Ir = I_r; }
	int get_r() const { return J_Ir; };		                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ;
	void set_d(const int& I_d) { J_Id = I_d; }
	int get_d() const { return J_Id; };			
	void set_C_j(const double& D_C_j) { C_j = D_C_j; }
	double get_C_j() const { return C_j; }			
	void set_T_j(const double& D_T_j) { T_j = D_T_j; }
	double get_T_j() const { return T_j; }			
	void set_rand_key(const double& D_rand_key) { rand_key = D_rand_key; }
	double get_rand_key() const { return rand_key; }
};
class batch
{

private:
	int B_Inum;
	int B_Is;
	double B_DPb;
	double B_DRb;
	double B_DCb;
	double B_DDb;
	double B_DSb;
	double B_DTb;
	int B_IMs;
	int B_IMnum;
	int B_Irs;
	int B_Jcount; 
public:

	vector<job> v_Bjobs;
	void update_batch(const job& j, const double& M_v);
	void cal_B_Pb(const double& M_v);
	void cal_B_Rb();
	void cal_B_Db();
	//void Init(int B_id);
	void Init(int B_id, int M_id);
	batch() :B_Inum(0), B_Is(0), B_DPb(0.0),  B_DRb(0.0), B_DCb(0.0), B_DDb(INF), B_DSb(0.0), B_DTb(INF), B_IMs(0), B_IMnum(0), B_Irs(0), B_Jcount(0)
	{

	};
	batch(int Inum) :B_Inum(Inum), B_Is(0), B_DPb(0.0),  B_DRb(0.0), B_DCb(0.0), B_DDb(INF), B_DSb(0.0), B_DTb(INF), B_IMs(0), B_IMnum(0), B_Irs(0), B_Jcount(0)
	{

	};
	~batch() = default;
	void set_num(const int& I_num) { B_Inum = I_num; }
	int get_num() const { return B_Inum; }			
	void set_s(const int& I_s) { B_Is = I_s; }	
	int get_s() const { return B_Is; }				
	void set_Rb(const double& D_Rb) { B_DRb = D_Rb; }	
	double get_Rb() const { return B_DRb; }				
	//void set_Pb(const double& I_Pb) { B_DPb = D_Pb; }
	double get_Pb() const { return B_DPb; }				
	void set_Cb(const double& D_Cb) { B_DCb = D_Cb; }	
	double get_Cb() const { return B_DCb; }
	void set_Db(const double& D_Db) { B_DDb = D_Db; }	
	double get_Db() const { return B_DDb; }			
	void set_Sb(const double& D_Sb) { B_DSb = D_Sb; }	
	double get_Sb() const { return B_DSb; }			
	void set_Tb(const double& D_Tb) { B_DTb = D_Tb; }	
	double get_Tb() const { return B_DTb; }			
	void set_Ms(const int& I_Ms) { B_IMs = I_Ms; }	
	int get_Ms() const { return B_IMs; }				
	void set_Mnum(const int& I_Mnum) { B_IMnum = I_Mnum; }
	int get_Mnum() const { return B_IMnum; }			
	void set_rs(const int& I_rs) { B_Irs = I_rs; }
	int get_rs() const { return B_Irs; }		 
	void set_J_num(const int& I_Jnum) { B_Jcount = I_Jnum; }	
	int get_J_num() const { return B_Jcount; }			
};

class machine
{
private:
	int M_Inum;
	int M_IL;
	double M_Dv;
	int M_Is;
	double M_D_Lmax;
	double M_D_TC;
	int M_Bcount;
public:
	vector<batch> v_Mbatches;
	void cal_M();
	machine() :M_Inum(0), M_IL(0), M_Dv(0.0), M_Is(0), M_D_Lmax(static_cast<double>(-MAX)), M_D_TC(0.0), M_Bcount(0)
	{

	};
	machine(int Inum) :M_Inum(Inum), M_IL(0), M_Dv(0.0), M_D_Lmax(static_cast<double>(-MAX)), M_D_TC(0.0), M_Bcount(0)
	{
		if (Inum > 0 && Inum <= M_1)
		{
			M_IL = L_1;
			M_Dv = v_1;
			M_Is = M_s1;
		}
		else if (Inum > M_1 && Inum <= (M_1 + M_2))
		{
			M_IL = L_2;
			M_Dv = v_2;
			M_Is = M_s2;
		}
		else if (Inum > (M_1 + M_2) && Inum <= (M_1 + M_2 + M_3))
		{
			M_IL = L_3;
			M_Dv = v_3;
			M_Is = M_s3;
		}
		else if (Inum > (M_1 + M_2 + M_3) && Inum <= (M_1 + M_2 + M_3 + M_4))
		{
			M_IL = L_4;
			M_Dv = v_4;
			M_Is = M_s4;
		}
		else if (Inum > (M_1 + M_2 + M_3 + M_4) && Inum <= (M_1 + M_2 + M_3 + M_4 + M_5))
		{
			M_IL = L_5;
			M_Dv = v_5;
			M_Is = M_s5;
		}
	};
	~machine() = default;
	void set_num(const int& I_num) { M_Inum = I_num; }	
	int get_num() const { return M_Inum; }				
	void set_L(const int& I_L) { M_IL = I_L; }		
	int get_L() const { return M_IL; }					
	void set_v(const double& D_v) { M_Dv = D_v; }	
	double get_v() const { return M_Dv; }					
	int get_s() const { return M_Is; }				
	void set_Lmax(const double& D_Lmax) { M_D_Lmax = D_Lmax; }	
	double get_Lmax() const { return M_D_Lmax; }				
	void set_TC(const double& D_TC) { M_D_TC = D_TC; }	
	double get_TC() const { return M_D_TC; }			

	void set_B_num(const int& I_Bnum) { M_Bcount = I_Bnum; }	
	int get_B_num() const { return M_Bcount; }		

};

class scheme
{
private:
	double Lmax;
	double TC;
public:
	vector<machine> v_Smachines;
	void cal_Fitness();
	double get_Lmax()const { return Lmax; };
	double get_TC()const { return TC; };

	void set_Lmax(const double& D_Lmax) { Lmax = D_Lmax; }
	void set_TC(const double& D_TC) { TC = D_TC; }
	scheme() :Lmax(0), TC(0) {};
	~scheme() = default;

};
class chrom
{
private:
	int number;
	int ChromLength;
public:
	int unNumber;
	int ass_reg;
	double ass_angle1;
	int No_Rank;
	double ass_angle2;
	double ass_angle3;
	double first_index;
	double distance2o;
	vector<int> Set_p;
	int N_p;
	bool flag;
	scheme chro_sche(vector<job>& jobs);
	machine S_G(vector<job> jobs, const int& mids, int bids);
	Two ChromFitness;
	vector<job> v_chromejobs;
	chrom(int num=0,int unnum=0,int ass_reg=0, double ass_angle1 = 0.0, int No_Rank = 0, double ass_angle2 = 0.0, double ass_angle3=0.0, double first_index=0.0, double distance2o=0.0, int N_p=0, bool flag=false) :number(num), ChromLength(Job_Num), unNumber(unnum), ass_reg(ass_reg), ass_angle1(ass_angle1), No_Rank(No_Rank), ass_angle2(ass_angle2), ass_angle3(ass_angle3), first_index(first_index), distance2o(distance2o), N_p(N_p), flag(flag)
	{
		v_chromejobs.reserve(Job_Num);
	};
	~chrom() = default;
	void set_number(const int& num) { number = num; }
	int get_number() const { return number; }
	void set_Chro_Length(const int& Length) { ChromLength = Length; }
	int get_Chro_Length()const { return ChromLength; }
	bool operator==(chrom& v1)
	{
		if (this->ChromFitness.Lmax == v1.ChromFitness.Lmax && this->ChromFitness.TC == v1.ChromFitness.TC)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator!=(chrom& v1)
	{
		if (this->ChromFitness.Lmax != v1.ChromFitness.Lmax || this->ChromFitness.TC != v1.ChromFitness.TC)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

#endif // !__JOB__

