#include"job.h"
#include<algorithm>
#include <cassert>
#define SET_NAME_BASE(text1,text2) text1##text2
#define SET_NAME(text1,text2) SET_NAME_BASE(text1,text2)


void sort_Rand_key(vector<job>& jobs)
{
	sort(jobs.begin(), jobs.end(), [](const job& j1, const job& j2)
		{
			return j1.get_rand_key() < j2.get_rand_key();
		});
}
/************************************************************************************************/

void batch::update_batch(const job& j, const double& M_v)
{
	B_Irs -= j.get_s();
	B_Is += j.get_s();
	B_Jcount++;
	cal_B_Pb(M_v);
	cal_B_Rb();
	cal_B_Db();
	/*if (B_Inum != 1)
	{
		B_ICb = (this - 1)->B_ISb + B_IPb;
	}
	else
	{
		B_ICb = B_IPb;
	}*/

	if (B_Inum != 1)
	{
		if ((this - 1)->B_DCb > B_DRb)
		{	
			B_DSb = (this - 1)->B_DCb;
			
		}
		else
		{
			B_DSb = B_DRb;
		}	
	}
	else
	{
		B_DSb = B_DRb;
	}
	B_DCb = this->B_DSb + B_DPb;
	B_DTb = B_DCb - B_DDb;
}


void batch::cal_B_Pb(const double& M_v)
{
	int temp = 0;
	vector<job>::const_iterator iter_job;
	for (iter_job = v_Bjobs.cbegin(); iter_job != v_Bjobs.cend(); ++iter_job)
	{
		if (iter_job->get_p() > temp)
		{
			temp = iter_job->get_p();
		}
	}
	B_DPb = temp/ M_v;
}

void batch::cal_B_Rb()
{
	double temp = 0.0;
	vector<job>::const_iterator iter_job;
	for (iter_job = v_Bjobs.cbegin(); iter_job != v_Bjobs.cend(); ++iter_job)
	{
		if (iter_job->get_r() > temp)
		{
			temp = iter_job->get_r();
		}
	}
	B_DRb = temp;
}

void batch::cal_B_Db()
{
	double temp = INF;
	vector<job>::const_iterator iter_job;
	for (iter_job = v_Bjobs.cbegin(); iter_job != v_Bjobs.cend(); ++iter_job)
	{
		if (iter_job->get_d() < temp)
		{
			temp = iter_job->get_d();
		}
	}
	B_DDb = temp;
}



void batch::Init(int B_id, int M_id)
{
	B_IMnum = M_id;
	B_Inum = B_id;
	if (M_id > 0 && M_id <= M_1)
	{
		B_IMs = M_s1;
	}
	else if (M_id > M_1 && M_id <= (M_1 + M_2))
	{
		B_IMs = M_s2;
	}
	else if (M_id > (M_1 + M_2) && M_id <= (M_1 + M_2 + M_3))
	{
		B_IMs = M_s3;
	}
	else if (M_id > (M_1 + M_2 + M_3) && M_id <= (M_1 + M_2 + M_3 + M_4))
	{
		B_IMs = M_s4;
	}
	else if (M_id > (M_1 + M_2 + M_3 + M_4) && M_id <= (M_1 + M_2 + M_3 + M_4 + M_5))
	{
		B_IMs = M_s5;
	}
	B_Is = 0;
	B_DCb = 0;
	B_DSb = 0;
	B_DDb = INF;
	B_DTb = INF;
	B_DPb = 0;
	B_Irs = B_IMs;
	B_Jcount = 0;
	v_Bjobs.clear();
}

/************************************************************************************/

void machine::cal_M()
{

	for (auto i = v_Mbatches.begin(); i != v_Mbatches.end(); ++i)
	{
		i->set_J_num(static_cast<int>(i->v_Bjobs.size()));
		i->cal_B_Pb(M_Dv);
		i->cal_B_Rb();
	}
	for (auto i = v_Mbatches.begin(); i != v_Mbatches.end(); ++i)
	{
		if (i == v_Mbatches.begin())
		{
			i->set_Sb(i->get_Rb());

		}
		else
		{
			i->set_Sb(((i - 1)->get_Cb() > i->get_Rb()) ? (i - 1)->get_Cb() : i->get_Rb());
			//i->set_Cb(i->get_Sb() + i->get_Pb());
		}
		i->set_Cb(i->get_Pb() + i->get_Sb());
		for (auto j = i->v_Bjobs.begin(); j != i->v_Bjobs.end(); ++j)
		{
			j->set_C_j(i->get_Cb());
		}
		for (auto j = i->v_Bjobs.begin(); j != i->v_Bjobs.end(); ++j)
		{
			if (j->get_d() < i->get_Db())
			{
				i->set_Db(j->get_d());
			}
		}

		i->set_Tb(i->get_Cb() - i->get_Db());
		for (auto j = i->v_Bjobs.begin(); j != i->v_Bjobs.end(); ++j)
		{
			j->set_T_j(i->get_Tb());
		}
	}
	for (auto i = v_Mbatches.begin(); i != v_Mbatches.end(); ++i)
	{
		M_D_TC += i->get_J_num() * i->get_Pb() * M_IL;
		if (i->get_Tb() > M_D_Lmax)
		{
			M_D_Lmax = i->get_Tb();
		}

	}
	M_Bcount = static_cast<int>(v_Mbatches.size());
}




void scheme::cal_Fitness()//计算总误工时间,总误工数量,能耗
{
	for (auto i = v_Smachines.begin(); i != v_Smachines.end(); ++i)
	{
		//i->cal_M();
		//Lmax = k->get_Lmax() > Lmax?k->get_Lmax():Lmax;
		if (i->get_Lmax() > Lmax)
		{
			Lmax = i->get_Lmax();
		}
		TC += i->get_TC();
	}
	if (TC == 0 && Lmax == 0)
	{
		for (auto i = v_Smachines.cbegin(); i != v_Smachines.cend(); ++i)
		{
			cout << i->get_num() << " " << i->get_Lmax() << " " << i->get_TC() << endl;
		}
		cout << "cal_Fitness错误" << endl;
		system("pause");
	}
}

scheme chrom::chro_sche(vector<job>& jobs)
{
	scheme TempScheme;
	//scheme TempScheme2; scheme TempScheme3; scheme TempScheme4; scheme TempScheme5;
	using V_job = vector<vector<job>>;
	V_job Dg_job_name;
	vector<job>::iterator iter_Sjob;
	vector<machine>::iterator iter_Sm;
	for (int i = 1; i <= M; ++i)
	{
		vector<job> SET_NAME(job, i);
		Dg_job_name.push_back(jobi);//set the number of job set, eg.job1,job2
	}
	for (iter_Sjob = jobs.begin(); iter_Sjob != jobs.end(); iter_Sjob++)
	{
		//V_job::iterator iter_Dg = Dg_job_name.begin();
		if (iter_Sjob->get_rand_key() > int(M + 1))
		{
			cout << "chro_sche____" << iter_Sjob->get_num();
		}
		
		Dg_job_name.at(int(iter_Sjob->get_rand_key()) - 1).push_back(*iter_Sjob);
		
		
	}
	
	for (V_job::iterator iter_Dg = Dg_job_name.begin(); iter_Dg != Dg_job_name.end(); ++iter_Dg)
	{
		

		if ((*iter_Dg).size() == 0)
		{
			continue;

		}
		else
		{
			TempScheme.v_Smachines.push_back(S_G(*iter_Dg, int((*iter_Dg).begin()->get_rand_key()), 1));//有多少台机器，就运行多少次
			
		}

	}
	
	TempScheme.cal_Fitness();

	return TempScheme;
}

machine chrom::S_G(vector<job> jobs, const int& mids, int bids)
{
	
	int batchcapacity = 0;
	batch tempbatch;
	
	tempbatch.Init(bids, mids);
	vector<job> ::iterator iter_job1;
	machine tempMachine(mids);
	sort_Rand_key(jobs);

	for (iter_job1 = jobs.begin(); iter_job1 != jobs.end(); ++iter_job1)
	{
		
		if (tempbatch.get_rs() - iter_job1->get_s() >= 0)
		{
			tempbatch.v_Bjobs.push_back(*iter_job1);
			tempbatch.set_rs(tempbatch.get_rs() - iter_job1->get_s());
			tempbatch.set_s(tempbatch.get_s() + iter_job1->get_s());

		}
		else
		{
			--iter_job1;
			tempbatch.set_J_num(static_cast<int>(tempbatch.v_Bjobs.size()));
			tempMachine.v_Mbatches.push_back(tempbatch);
			tempbatch.Init(++bids, mids);
			
		}
		if (iter_job1 == jobs.end() - 1)
		{
			tempbatch.set_J_num(static_cast<int>(tempbatch.v_Bjobs.size()));
			tempMachine.v_Mbatches.push_back(tempbatch);
		}
	}
	
	tempMachine.set_B_num(static_cast<int>(tempMachine.v_Mbatches.size()));
	tempMachine.cal_M();
	return tempMachine;

}