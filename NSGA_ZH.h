#ifndef __NSGA_ZH__
#define __NSGA_ZH__
#include"job.h"
//#include"constant.h"
#include<algorithm>
#include<random>
#include<cmath>
#include <fstream> 
//#include"AGNES.h"
namespace aa
{
	void test_selection_indidvial();
}
class NSGA
{
public:
	vector<vector<chrom>> archive;
	vector<chrom>Population;
	Two z;
	void cal_z(vector<chrom>& T);
	void Transform1(vector<chrom>& T);
	void Transform2(vector<chrom>& T, const double& a, const double& b);
	vector<int> Find_Extreme_Point(vector<chrom>& T);
	void Normalize(vector<chrom>& T);

	double cal_distance(const chrom& c,const NOO& n);

	void Fast_Non_Dominated_Sort(vector<chrom> &P, int& R_num);

	vector<chrom> Init(vector<job>& J);//Initialization
	chrom generateSolution(vector<job>& J,const int& num);
	vector<chrom> Init_Pure(const vector<job>& J);//Initialization_Pure

	chrom& tournament(chrom &ind1,chrom &ind2);
	//vector<chrom> recombination(vector<chrom>& T, const int Mate_Num);//Recombination
	void selction_crossover(vector<chrom>& T1, vector<chrom>& T2);//SBX
	void crossover(chrom& c1, chrom& c2, chrom& c3, chrom& c4);//SBX
	void mutation(vector<chrom>& T);//polynomial mutation
	vector<chrom> merge(vector<chrom>& P1, vector<chrom>& P2);

	vector<chrom> E_select(vector<chrom>& merged_population_T, int iter, int max_iter);//Environmental Selection
	void Encode(vector<job>& Jcode);
	void Decode(vector<chrom>& T);
	void NSGA_schedule(vector<job> ZH_instacne, int max_iterator, int _instance, vector<chrom>& Parent, const int& Count);
	chrom& find_p(const int& p, vector<chrom> &_P);
	//NO& find_rho_j(const int &i,vector<NO>& N);
	vector<chrom> find_Q(const int& i, vector<chrom>& _P);
	void Print_Pop(const vector<chrom>& T);


	int Roulette_Wheel_Selection(vector<NOO>& I);
	int Machine_Selection(vector<machine>& _M,vector<NO_B>& machineList);
	int Machine_Selection(vector<machine>& _M, int& M_num_x);
	vector<job> Cons_Cand_List(const vector<job>& J,const int& S);
	bool check_state(const vector<job>& J);
	job choose_job(vector<job>& _J);
	NSGA() {
	};
	~NSGA() = default;
	/*~NSGA() {};*/
};

class ANGLE
{
public:
	vector<vector<chrom>> A_Set;
	vector<chrom> origin_solution;
	ANGLE() = default;;
	ANGLE(const int& k, const vector<chrom>& S);
	~ANGLE() = default;;
	vector<vector<chrom>> cluster(const int& k);
	vector<chrom> selection_indidvial(vector<vector<chrom>>& div_S, const int& F_num,const int& iter, const int& max_iter,Two z);
};

#endif // !__NSGA_ZH__


