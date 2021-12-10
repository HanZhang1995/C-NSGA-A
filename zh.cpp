//#pragma warning(disable:4996)

#include"job.h"
#include"NSGA_ZH.h"
#include"constant.h"
#include<algorithm>
//#include "gnuplot_i.hpp"
#include <time.h> 
#include<string>
#include <fstream>
#include <sstream>

int Job_Num;
int Iterater_NUM;
int Instances_Num;
int Times;
int N;
int M_s1;
int M_s2;
int M_s3;
int M_s4;
int M_s5;
int L_1;
int L_2;
int L_3;
int L_4;
int L_5;
int M_1;
int M_2;
int M_3;
int M_4;
int M_5;
double v_1;
double v_2;
double v_3;
double v_4;
double v_5;

bool fun(Two A, Two B)
{

	return (A.TC < B.TC);
}
bool fun1(Two A, Two B)
{
	return (A.TC == B.TC && A.Lmax == B.Lmax);
}



int main()
{
	
//	aa::test_selection_indidvial();
	ifstream in("\\machine_information.txt");
	string line;
	vector<string> machineInfo;
	if (in) 
	{
		while (getline(in, line)) 
		{
			vector<string> res;
			string result;
			stringstream input(line);
			while (input >> result)
			{
				res.push_back(result);
			}
			for (int i = 0; i < res.size(); i++)
			{
				int mas = 0;
				string machineSize = res[0];
				if (machineSize == "Job_Num") mas = -4;
				if (machineSize == "Iterater_NUM") mas = -3;
				if (machineSize == "Instances_Num") mas = -2;
				if (machineSize == "Times") mas = -1;
				if (machineSize == "N") mas = 0;
				if (machineSize == "M_s1") mas = 1;
				if (machineSize == "M_s2") mas = 2;
				if (machineSize == "M_s3") mas = 3;
				if (machineSize == "M_s4") mas = 4;
				if (machineSize == "M_s5") mas = 5;
				if (machineSize == "M_1") mas = 6;
				if (machineSize == "M_2") mas = 7;
				if (machineSize == "M_3") mas = 8;
				if (machineSize == "M_4") mas = 9;
				if (machineSize == "M_5") mas = 10;
				if (machineSize == "L_1") mas = 11;
				if (machineSize == "L_2") mas = 12;
				if (machineSize == "L_3") mas = 13;
				if (machineSize == "L_4") mas = 14;
				if (machineSize == "L_5") mas = 15;
				if (machineSize == "v_1") mas = 16;
				if (machineSize == "v_2") mas = 17;
				if (machineSize == "v_3") mas = 18;
				if (machineSize == "v_4") mas = 19;
				if (machineSize == "v_5") mas = 20;
				switch (mas)
				{
				case -4:
					Job_Num = stoi(res[1], nullptr);
					//cout << Job_Num << endl;
					break;
				case -3:
					Iterater_NUM = stoi(res[1], nullptr);
					//cout << Iterater_NUM << endl;
					break;
				case -2:
					Instances_Num = stoi(res[1], nullptr);
					//cout << Instances_Num << endl;
					break;
				case -1:
					Times = stoi(res[1], nullptr);
					//cout << Times << endl;
					break;
				case 0:
					N = stoi(res[1], nullptr);
					//cout << N << endl;
					break;
				case 1:
					M_s1 = stoi(res[1], nullptr);
					//cout << m_s1 << endl;
					break;
				case 2:
					M_s2 = stoi(res[1], nullptr);
					//cout << m_s2 << endl;
				case 3:
					M_s3 = stoi(res[1], nullptr);
					//cout << m_s3 << endl;
					break;
				case 4:
					M_s4 = stoi(res[1], nullptr);
					//cout << m_s4 << endl;
					break;
				case 5:
					M_s5 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 6:
					M_1 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 7:
					M_2 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 8:
					M_3 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 9:
					M_4 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 10:
					M_5 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 11:
					L_1 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 12:
					L_2 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 13:
					L_3 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 14:
					L_4 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 15:
					L_5 = stoi(res[1], nullptr);
					//cout << m_s5 << endl;
					break;
				case 16:
					v_1 = stod(res[1], nullptr);
					//cout << v_1 << endl;
					break;
				case 17:
					v_2 = stod(res[1], nullptr);
					//cout << v_2 << endl;
					break;
				case 18:
					v_3 = stod(res[1], nullptr);
					//cout << v_3 << endl;
					break;
				case 19:
					v_4 = stod(res[1], nullptr);
					//cout << v_4 << endl;
					break;
				case 20:
					v_5 = stod(res[1], nullptr);
					//cout << v_5 << endl;
					break;
				}
			}
		}

	}
	else
	{
		cout << "no such file" << endl;
	}
	for (int instances = 1; instances <= Instances_Num; instances++)
	{
		double start = clock();
		ifstream ifile;
		char filename[100];
		sprintf_s(filename, "\\Experiment_data\\%d-%d.txt", Job_Num, instances);
		ifile.open(filename, ifstream::in);//打开文件
		if (!ifile)
		{
			cout << "open file" << Job_Num << "-" << instances << "fail" << endl;
			exit(1);
		}
		vector<job> vtemp_job;
		int Inum, Ip, Is, Ir, Id;
		//int i = 0;
		while (!ifile.eof())
		{
			ifile >> Inum >> Is >> Ip >> Ir >> Id;
			job tempjob;
			tempjob.set_num(Inum); tempjob.set_p(Ip); tempjob.set_s(Is); tempjob.set_r(Ir); tempjob.set_d(Id);
			vtemp_job.push_back(tempjob);
		}
		ifile.close();
		double TIME{ 0 };
		char str[100];
		char str1[100];
		char str2[100];
		sprintf_s(str, "NSGA_A_%d_%d.xls", Job_Num, instances);
		sprintf_s(str1, "NSGA_A_%d_%d_non.xls", Job_Num, instances);
		sprintf_s(str2, "\\Solution_data\\NSGA_A_%d_%d.txt", Job_Num, instances);
		cout << Job_Num << "-" << instances << " in NSGA_A doing..." << endl;
		ofstream out;
		ofstream out1;
		ofstream out2;
		out.open(str, fstream::app);
		out1.open(str1, fstream::app);
		out2.open(str2, fstream::out);
		vector<Two> avgSolution;
		for (int count = 1; count <= Times; ++count)
		{
			
			NSGA NSGA_ZH;
			NSGA_ZH.NSGA_schedule(vtemp_job, Iterater_NUM, instances, NSGA_ZH.Population, count);


			vector<Two> Save_solution;
			for (auto i = NSGA_ZH.Population.begin(); i != NSGA_ZH.Population.end(); ++i)
			{
				Save_solution.push_back(i->ChromFitness);

			}

			


			out1 << "jobnum\t" << "times\n" << Job_Num << '\t' << TIME << '\t' << "TC" << '\t';
			for (auto i = Save_solution.begin(); i != Save_solution.end(); ++i)
			{
				out1 << i->TC << "\t";
			}
			out1 << '\n';
			out1 << '\t' << '\t' << "TT" << '\t';
			for (auto i = Save_solution.begin(); i != Save_solution.end(); ++i)
			{
				out1 << i->Lmax << "\t";
			}
			out1 << '\n';
			sort(Save_solution.begin(), Save_solution.end(), fun);
			Save_solution.erase(unique(Save_solution.begin(), Save_solution.end(), fun1), Save_solution.end());
			for (auto o = Save_solution.begin(); o != Save_solution.end(); ++o)
			{
				double temp = o->Lmax, temp1 = o->TC;
				for (auto p = o + 1; p != Save_solution.end(); ++p)
				{
					if (p->Lmax > temp && p->TC > temp1)
					{
						p->TC = MAX;
						p->Lmax = MAX;
					}
					if (p->Lmax == temp)
					{
						if(p->TC > temp1)
						{
							p->TC = MAX;
							p->Lmax = MAX;
						}else
						{
							o->TC = MAX;
							o->Lmax = MAX;
							break;
						}
					}
					
					if ( p->TC == temp1)
					{

						if (p->Lmax > temp)
						{
							p->TC = MAX;
							p->Lmax = MAX;
						}
						else
						{
							o->TC = MAX;
							o->Lmax = MAX;
							break;
						}
					}
				}
			}
			for (auto it = Save_solution.begin(); it != Save_solution.end();)
			{
				if (it->TC == MAX && it->Lmax == MAX)
					it = Save_solution.erase(it);    
				else
					++it;    
			}

			double stop = clock();
			++TIME;
			for (const auto& s : Save_solution)
			{
				avgSolution.push_back(s);
			}
			out << "jobnum\t" << "times\n" << Job_Num << '\t' << TIME << '\t' << "TC" << '\t';
			for (auto i = Save_solution.begin(); i != Save_solution.end(); ++i)
			{
				out << i->TC << "\t";
			}
			out << '\n';
			out << '\t' << '\t' << "TT" << '\t';
			for (auto i = Save_solution.begin(); i != Save_solution.end(); ++i)
			{
				out << i->Lmax << "\t"; 
			}
			out << '\n';
			out << '\t' << (stop - start) / CLOCKS_PER_SEC << endl;


		}
		
		sort(avgSolution.begin(), avgSolution.end(), fun);

		for (auto o = avgSolution.begin(); o != avgSolution.end(); ++o)
		{
			double temp = o->Lmax, temp1 = o->TC;
			for (auto p = o + 1; p != avgSolution.end(); ++p)
			{
				if (p->Lmax > temp && p->TC > temp1)
				{
					p->TC = MAX;
					p->Lmax = MAX;
				}
				if (p->Lmax == temp)
				{
					if(p->TC > temp1)
					{
						p->TC = MAX;
						p->Lmax = MAX;
					}else
					{
						o->TC = MAX;
						o->Lmax = MAX;
						break;
					}
				}
				
				if ( p->TC == temp1)
				{

					if (p->Lmax > temp)
					{
						p->TC = MAX;
						p->Lmax = MAX;
					}
					else
					{
						o->TC = MAX;
						o->Lmax = MAX;
						break;
					}
				}
			}
		}
		for (auto it = avgSolution.begin(); it != avgSolution.end();)
		{
			if (it->TC == MAX && it->Lmax == MAX)
				it = avgSolution.erase(it);    
			else
				++it;    
		}
		for (auto i = avgSolution.begin(); i != avgSolution.end(); ++i)
		{
			out2 << i->TC << "\t" << i->Lmax << "\n";
		}
		out.close();
		out1.close();
		out2.close();
		cout << endl;
	}

	return 0;
}