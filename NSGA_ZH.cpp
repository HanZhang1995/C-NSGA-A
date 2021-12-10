#include "NSGA_ZH.h"
#include <cassert>
//#include "gnuplot_i.hpp"
static random_device rd;   // non-deterministic generator  
mt19937_64 gen(rd());  // to seed mersenne twister.  
static uniform_real_distribution<> dist1(1, M + 1); // distribute results between 1 and M inclusive.  
static uniform_real_distribution<> dist2(0, N); // distribute results between 0 and N inclusive.

static uniform_real_distribution<> dist3(0, 1); // distribute results between 0 and 1 inclusive.

static double random_x(const int& i, const int& j)
{
	uniform_real_distribution<> distX(i, static_cast<int>(j + 1));
	return distX(gen);
}

static double random_y(const int& i, const int& j)
{
	uniform_real_distribution<> distY(i, j);
	return distY(gen);
}

void sort_Lmax(vector<chrom>& chs)
{
	sort(chs.begin(), chs.end(), [](const chrom& c1, const chrom& c2)
		{
			return c1.ChromFitness.Lmax < c2.ChromFitness.Lmax;
		});
}


void sort_TC(vector<chrom>& chs)
{
	sort(chs.begin(), chs.end(), [](const chrom& c1, const chrom& c2)
		{
			return c1.ChromFitness.TC < c2.ChromFitness.TC;
		});
}

void sort_angle(vector<chrom>& chs)
{
	sort(chs.begin(), chs.end(), [](const chrom& c1, const chrom& c2)
		{
			return c1.ass_angle1 < c2.ass_angle1;
		});
}
double max(double x, double y, double z)
{
	if (x > y)
		y = x;
	if (y > z)
		y = z;
	return z;
}


job& find_j(const int& j, vector<job>& _J)
{
	for (job& m : _J)
	{
		if (m.get_num() == j)
		{
			return m;
		}
	}
}


int Dominated_Compared(chrom P_1, chrom P_2)
{
	if (((P_1.ChromFitness.Lmax < P_2.ChromFitness.Lmax) && (P_1.ChromFitness.TC <= P_2.ChromFitness.TC)) ||
		((P_1.ChromFitness.Lmax <= P_2.ChromFitness.Lmax) && (P_1.ChromFitness.TC < P_2.ChromFitness.TC)))
	{
		return 1;
	}
	else if (((P_1.ChromFitness.Lmax > P_2.ChromFitness.Lmax) && (P_1.ChromFitness.TC >= P_2.ChromFitness.TC)) ||
		((P_1.ChromFitness.Lmax >= P_2.ChromFitness.Lmax) && (P_1.ChromFitness.TC > P_2.ChromFitness.TC)))
	{
		return -1;
	}
	else
	{
		return 0;
	}
}


ANGLE::ANGLE(const int& k, const vector<chrom>& S)
{
	A_Set.resize(k);
	origin_solution = S;
}

vector<vector<chrom>> ANGLE::cluster(const int& k)
{
	for(auto s=origin_solution.begin();s!=origin_solution.end();++s)
	{
		two_point<double, double> O{};
		two_point<double, double> P2(s->ChromFitness.T_Lmax, s->ChromFitness.T_TC);
		two_point<double, double> P3(0.0, 1.0);
		double A = P2.y - O.y;
		double B = O.x - P2.x;
		double A_plus = P3.y - O.y;
		double B_plus = O.x - P3.x;
		double value = abs((A * A_plus + B * B_plus) / (pow(pow(A, 2) + pow(B, 2), 0.5) * pow(pow(A_plus, 2) + pow(B_plus, 2), 0.5)));
		//value = -1;
		s->ass_angle1 = ceil(acos(value) * 180.0 / M_PI * 100.0) / 100.0;
		if(s->ass_angle1==90)
		{
			s->ass_reg = k;
		}else
		{
			s->ass_reg = static_cast<int>(floor(s->ass_angle1 / (90.0 / k) + 1));
		}
		

	}
	sort_angle(origin_solution);
	for(int i=1;i<=k;++i)
	{
		for (const auto& s : origin_solution)
		{
			if(s.ass_reg==i)
			{
				A_Set.at(i-1).push_back(s);
			}
			
		}
		
	}
	
	
	return A_Set;
}


void aa::test_selection_indidvial()
{
	cout << "--- test begin ----"<<endl;
	vector<chrom> ch;
	vector <vector<chrom>> div;
	vector<chrom> temp1;
	chrom chr_1(1, 0, 1, 20, 0, 0.0, 0.0, 0.0, 0);
	chr_1.ChromFitness.T_Lmax = 0.1;
	chr_1.ChromFitness.T_TC = 0.9;
	chrom chr_2(2, 0, 1, 15, 0, 0.0, 0.0, 0.0, 0);
	chr_2.ChromFitness.T_Lmax = 0.3;
	chr_2.ChromFitness.T_TC = 0.7;
	temp1.push_back(chr_1);
	temp1.push_back(chr_2);
	div.push_back(temp1);
	vector<chrom> temp2;
	chrom chr_3(3, 0, 2, 35, 0, 0.0, 0.0, 0.0, 0);
	chr_3.ChromFitness.T_Lmax = 0.4;
	chr_3.ChromFitness.T_TC = 0.6;
	chrom chr_4(4, 0, 2, 55, 0, 0.0, 0.0, 0.0, 0);
	chr_4.ChromFitness.T_Lmax = 0.6;
	chr_4.ChromFitness.T_TC = 0.4;
	temp2.push_back(chr_3);
	temp2.push_back(chr_4);
	div.push_back(temp2);
	vector<chrom> temp3;
	chrom chr_5(5, 0, 3, 68, 0, 0.0, 0.0, 0.0, 0);
	chr_5.ChromFitness.T_Lmax = 0.8;
	chr_5.ChromFitness.T_TC = 0.1;
	chrom chr_6(6, 0, 3, 75, 0, 0.0, 0.0, 0.0, 0);
	chr_6.ChromFitness.T_Lmax = 0.7;
	chr_6.ChromFitness.T_TC = 0.3;
	
	temp3.push_back(chr_5);
	temp3.push_back(chr_6);
	div.push_back(temp3);

	ANGLE A;
	A.origin_solution.push_back(chr_1);
	A.origin_solution.push_back(chr_2);
	A.origin_solution.push_back(chr_3);
	A.origin_solution.push_back(chr_4);
	
	//ch = A.selection_indidvial(div);
	cout << "--- test end ----" << endl;
}
vector<chrom> ANGLE::selection_indidvial(vector<vector<chrom>>& div_S,const int& F_num, const int& iter, const int& max_iter, Two z)
{
	vector<chrom> choose_solution;
	choose_solution.resize(div_S.size());
	for(auto d=div_S.begin();d!=div_S.end();++d)
	{
		
		if(!d->empty())
		{
			
			for(auto s=d->begin();s!=d->end();++s)
			{
				
				double value = (180.0 * (d - div_S.begin()) + 90.0) / (2 * div_S.size());
				s->ass_angle2 = abs(s->ass_angle1 - value);
				
				s->first_index = s->ass_angle2;
				
				
			}
			sort(d->begin(), d->end(), [](const chrom& ch1, const chrom& ch2)
				{
					return ch1.first_index < ch2.first_index;
				});
			vector<chrom> vec_chr;
			for(auto temp=d->begin();temp!=d->end();++temp)
			{
				if(temp->first_index == d->begin()->first_index)
				{
					vec_chr.push_back(*temp);
				}
			}
			chrom selected_ind;
			if(vec_chr.size()==1)
			{
				selected_ind = vec_chr.at(0);
			}
			else
			{

				for(auto &v:vec_chr)
				{
					v.distance2o = pow(pow(v.ChromFitness.T_Lmax, 2) + pow(v.ChromFitness.T_TC, 2), 0.5);//* (cos(v.ass_angle2) + sin(v.ass_angle2))
				}
				sort(vec_chr.begin(), vec_chr.end(), [](const chrom& ch1, const chrom& ch2)
					{
						return ch1.distance2o < ch2.distance2o;
					});
				vector<chrom> vec_chr_dis;
				
				for (auto temp = vec_chr.begin(); temp != vec_chr.end(); ++temp)
				{
					if (temp->distance2o == vec_chr.begin()->distance2o)
					{
						vec_chr_dis.push_back(*temp);
					}
				}
				if (vec_chr_dis.size() == 1)
				{
					selected_ind = vec_chr_dis.at(0);
				}
				else
				{
					sort(vec_chr_dis.begin(), vec_chr_dis.end(), [](const chrom& ch1, const chrom& ch2)
						{
							return ch1.ChromFitness.T_Lmax < ch2.ChromFitness.T_Lmax;
						});
					vector<chrom> vec_val;
					for(auto v= vec_chr_dis.begin();v!=vec_chr_dis.end();++v)
					{
						if (v->ChromFitness.T_Lmax== vec_chr_dis.begin()->ChromFitness.T_Lmax&& v->ChromFitness.T_TC == vec_chr_dis.begin()->ChromFitness.T_TC)
						{
							vec_val.push_back(*v);
						}
					}

					assert(vec_val.size() == 1);
					if (vec_val.size() > 1)
					{
						vec_chr_dis.erase(unique(vec_chr_dis.begin(), vec_chr_dis.end(), [](const chrom& ch1, const chrom& ch2)
							{
								return ch1.ChromFitness.T_Lmax == ch2.ChromFitness.T_Lmax && ch1.ChromFitness.T_TC == ch2.ChromFitness.T_TC;
							}), vec_chr_dis.end());
						vec_chr_dis.push_back(vec_val.at(static_cast<int>(random_y(0, static_cast<int>(vec_val.size())))));
					}
					if (vec_chr_dis.size() == 1)
					{
						selected_ind = vec_val.at(0);
					}
					else
					{
						two_point<double, double> near1{0.0,1.0};
						two_point<double, double> near2{1.0,0.0};
						if(d!= div_S.begin())
						{
							for (auto i = d - div_S.begin() - 1; i >= 0; --i)
							{
								if (!div_S.at(i).empty())
								{
									sort(div_S.at(i).begin(), div_S.at(i).end(), [](const chrom& ch1, const chrom& ch2)
										{
											return ch1.ass_angle1 < ch2.ass_angle1;
										});
									near1.x = div_S.at(i).at(div_S.at(i).size() - 1).ChromFitness.T_Lmax;
									near1.y = div_S.at(i).at(div_S.at(i).size() - 1).ChromFitness.T_TC;
								}
							}
						}
						if (d != div_S.end())
						{
							for (auto i = d - div_S.begin() + 1; i < static_cast<int>(div_S.size()); ++i)
							{
								if (!div_S.at(i).empty())
								{
									sort(div_S.at(i).begin(), div_S.at(i).end(), [](const chrom& ch1, const chrom& ch2)
										{
											return ch1.ass_angle1 < ch2.ass_angle1;
										});
									near2.x = div_S.at(i).at(0).ChromFitness.T_Lmax;
									near2.y = div_S.at(i).at(0).ChromFitness.T_TC;
								}
							}
						}

							
						two_point<double, double> O{};//原点
							
						two_point<double, double> P3(0.0, 1.0);

							
						double A = near1.y - O.y;
						double B = O.x - near1.x;
						double A_plus = near2.y - O.y;
						double B_plus = O.x - near2.x;
							
						double A_plus_plus = P3.y - O.y;
						double B_plus_plus = O.x - P3.x;
						double near1_value = abs((A * A_plus_plus + B * B_plus_plus) / (pow(pow(A, 2) + pow(B, 2), 0.5) * pow(pow(A_plus_plus, 2) + pow(B_plus_plus, 2), 0.5)));
						double near2_value = abs((A_plus * A_plus_plus + B_plus * B_plus_plus) / (pow(pow(A_plus, 2) + pow(B_plus, 2), 0.5) * pow(pow(A_plus_plus, 2) + pow(B_plus_plus, 2), 0.5)));							
						//value = -1;
						double dif_angle= (acos(near1_value) + acos(near2_value)) / 2 * 180.0 / M_PI;

							
						for (auto s = vec_chr_dis.begin(); s != vec_chr_dis.end(); ++s)
						{
							two_point<double, double> P2(s->ChromFitness.T_Lmax, s->ChromFitness.T_TC);
							double _A = P2.y - O.y;
							double _B = O.x - P2.x;
							double s_value = abs((_A * A_plus_plus + _B * B_plus_plus) / (pow(pow(_A, 2) + pow(_B, 2), 0.5) * pow(pow(A_plus_plus, 2) + pow(B_plus_plus, 2), 0.5)));
							double s_angle = acos(s_value) * 180.0 / M_PI;
							s->ass_angle3 = ceil(abs(dif_angle - s_angle)* 100000)/ 100000;

						}
						sort(vec_chr_dis.begin(), vec_chr_dis.end(), [](const chrom& ch1, const chrom& ch2)
							{
								return ch1.ass_angle3 < ch2.ass_angle3;
							});
						vector<chrom> vec_ang_dis;
						for (auto v = vec_chr_dis.begin(); v != vec_chr_dis.end(); ++v)
						{

							if(v->ass_angle3 == vec_chr_dis.begin()->ass_angle3)
							{
								vec_ang_dis.push_back(*v);
							}
						}
						assert(!vec_ang_dis.empty());
						if(vec_ang_dis.size()==1)
						{
							selected_ind = vec_ang_dis.at(0);
						}
						if (vec_ang_dis.size() > 1)
						{
							selected_ind = vec_ang_dis.at(static_cast<int>(random_y(0, static_cast<int>(vec_ang_dis.size()))));
						}
					}
				}
			}

			for(auto& o:origin_solution)
			{
				if(o.get_number()== selected_ind.get_number())
				{
					o.flag = true;
				}
			}
			choose_solution.at(d - div_S.begin())= selected_ind;
		}
		
	}
	for (auto d = div_S.begin(); d != div_S.end(); ++d)
	{
		for(auto& _d:*d)
		{
			_d.ass_angle2 = 0.0;
			_d.first_index = 0.0;
			_d.ass_angle3 = 0.0;
		}
		if (d->empty())
		{

			chrom selected_s;

			for (auto o = origin_solution.begin(); o != origin_solution.end(); ++o)
			{
				
				double value = (180.0 * (d - div_S.begin()) + 90.0) / (2 * div_S.size());
				o->ass_angle2 = abs(o->ass_angle1 - value);

				o->first_index = o->ass_angle2;
				
			}
			sort(origin_solution.begin(), origin_solution.end(), [](const chrom& ch1, const chrom& ch2)
				{
					return ch1.first_index < ch2.first_index;
				});
			chrom temp1;
			for(auto o = origin_solution.begin(); o != origin_solution.end(); ++o)
			{
				if(o->flag==false)
				{
					temp1 = *o;
					break;
				}

			}
			vector<chrom> duplicate_chroms, duplicate_angle_chroms;
			for (auto o = origin_solution.begin(); o != origin_solution.end(); ++o)
			{
				if(o->flag == false)
				{
					if (temp1.ChromFitness.T_Lmax == o->ChromFitness.T_Lmax && temp1.ChromFitness.T_TC == o->ChromFitness.T_TC)
					{
						duplicate_chroms.push_back(*o);
					}
					
				}
				
			}
			assert(duplicate_chroms.size() == 1);
			if (duplicate_chroms.size() != 1)
			{
				sort(origin_solution.begin(), origin_solution.end(), [](const chrom& ch1, const chrom& ch2)
					{
						return ch1.ChromFitness.T_Lmax < ch2.ChromFitness.T_Lmax;
					});

				origin_solution.erase(unique(origin_solution.begin(), origin_solution.end(), [](const chrom& ch1, const chrom& ch2)
					{
						return ch1.ChromFitness.T_Lmax == ch2.ChromFitness.T_Lmax && ch1.ChromFitness.T_TC == ch2.ChromFitness.T_TC;
					}), origin_solution.end());
				origin_solution.push_back(duplicate_angle_chroms.at(static_cast<int>(random_y(0, static_cast<int>(duplicate_angle_chroms.size())))));
			}
			temp1 = duplicate_chroms.at(0);
			for (auto o = origin_solution.begin(); o != origin_solution.end(); ++o)
			{
				if (o->flag == false)
				{
					if (temp1.first_index == o->first_index && o->flag==false)
					{
						duplicate_angle_chroms.push_back(*o);
					}
				}

			}
			if(duplicate_angle_chroms.size()==1)
			{

				selected_s = duplicate_angle_chroms.at(0);
			}
			else
			{

				for (auto& v : duplicate_angle_chroms)
				{
					v.distance2o = pow(pow(v.ChromFitness.T_Lmax, 2) + pow(v.ChromFitness.T_TC, 2), 0.5) ;//* (cos(v.ass_angle2)+ sin(v.ass_angle2))
				}
				sort(duplicate_angle_chroms.begin(), duplicate_angle_chroms.end(), [](const chrom& ch1, const chrom& ch2)
					{
						return ch1.distance2o < ch2.distance2o;
					});
				vector<chrom> vec_chr_dis;

				for (auto temp = duplicate_angle_chroms.begin(); temp != duplicate_angle_chroms.end(); ++temp)
				{
					if (temp->distance2o == duplicate_angle_chroms.begin()->distance2o)
					{
						vec_chr_dis.push_back(*temp);
					}
				}
				assert(!vec_chr_dis.empty());
				if (vec_chr_dis.size() == 1)
				{

					selected_s = vec_chr_dis.at(0);
				}
				else
				{					

						two_point<double, double> near1{ 0.0,1.0 };
						two_point<double, double> near2{ 1.0,0.0 };
						int d_i = 1;
						int d_j = 1;
						while(choose_solution.at(d - div_S.begin() - d_i).get_number()==0)
						{
							++d_i;
						}
						while (choose_solution.at(d - div_S.begin() + d_j).get_number() == 0)
						{
							++d_j;
						}
						near1.x= choose_solution.at(d - div_S.begin() - d_i).ChromFitness.T_Lmax;
						near1.y = choose_solution.at(d - div_S.begin() - d_i).ChromFitness.T_TC;
						near2.x = choose_solution.at(d - div_S.begin() + d_j).ChromFitness.T_Lmax;
						near2.y = choose_solution.at(d - div_S.begin() + d_j).ChromFitness.T_TC;

						two_point<double, double> O{};
							
						two_point<double, double> P3(0.0, 1.0);



						double A = near1.y - O.y;
						double B = O.x - near1.x;
						double A_plus = near2.y - O.y;
						double B_plus = O.x - near2.x;
							
						double A_plus_plus = P3.y - O.y;
						double B_plus_plus = O.x - P3.x;
						double near1_value = abs((A * A_plus_plus + B * B_plus_plus) / (pow(pow(A, 2) + pow(B, 2), 0.5) * pow(pow(A_plus_plus, 2) + pow(B_plus_plus, 2), 0.5)));
						double near2_value = abs((A_plus * A_plus_plus + B_plus * B_plus_plus) / (pow(pow(A_plus, 2) + pow(B_plus, 2), 0.5) * pow(pow(A_plus_plus, 2) + pow(B_plus_plus, 2), 0.5)));
							
						double dif_angle = (acos(near1_value) + acos(near2_value)) / 2 * 180.0 / M_PI ;

							
						for (auto s = vec_chr_dis.begin(); s != vec_chr_dis.end(); ++s)
						{
							two_point<double, double> P2(s->ChromFitness.T_Lmax, s->ChromFitness.T_TC);
							double _A = P2.y - O.y;
							double _B = O.x - P2.x;
							double s_value = abs((_A * A_plus_plus + _B * B_plus_plus) / (pow(pow(_A, 2) + pow(_B, 2), 0.5) * pow(pow(A_plus_plus, 2) + pow(B_plus_plus, 2), 0.5)));
							//value = -1;
							double s_angle = acos(s_value) * 180.0 / M_PI;
							s->ass_angle3 = ceil(abs(dif_angle - s_angle)*100000)/100000;
							
						}
						sort(vec_chr_dis.begin(), vec_chr_dis.end(), [](const chrom& ch1, const chrom& ch2)
							{
								return ch1.ass_angle3 < ch2.ass_angle3;
							});
						vector<chrom> vec_ang_dis;
						for (auto v = vec_chr_dis.begin(); v != vec_chr_dis.end(); ++v)
						{
							if (v->ass_angle3 == vec_chr_dis.begin()->ass_angle3)
							{
								vec_ang_dis.push_back(*v);
							}
						}
						assert(!vec_ang_dis.empty());
						if (vec_ang_dis.size() == 1)
						{
							selected_s = vec_ang_dis.at(0);
						}
						if (vec_ang_dis.size() > 1)
						{
							int x = static_cast<int>(random_y(0, static_cast<int>(vec_ang_dis.size())));
							selected_s = vec_ang_dis.at(x);
						}
					//}
				}
			}
			for(auto &j:origin_solution)
			{
				if(j.get_number()== selected_s.get_number())
				{
					j.flag = true;
				}
			}
			choose_solution.at(d - div_S.begin()) = selected_s;
		}
	}
	return choose_solution;
}

void NSGA::cal_z(vector<chrom>& T)
{
	double a_TC = INF, b_Lmax = INF, a_M_TC = 0, b_M_Lmax = 0;
	z.TC = INF; z.Lmax = INF;
	for (auto iter = T.begin(); iter != T.end(); ++iter)
	{
		if (iter->ChromFitness.TC < a_TC)
		{
			a_TC = iter->ChromFitness.TC;
		}
		if (iter->ChromFitness.Lmax < b_Lmax)
		{
			b_Lmax = iter->ChromFitness.Lmax;
		}
		if (iter->ChromFitness.TC > a_M_TC)
		{
			a_M_TC = iter->ChromFitness.TC;
		}
		if (iter->ChromFitness.Lmax > b_M_Lmax)
		{
			b_M_Lmax = iter->ChromFitness.Lmax;
		}
	}
	if (z.TC > a_TC)
	{
		z.TC = a_TC;
	}
	if (z.Lmax > b_Lmax)
	{
		z.Lmax = b_Lmax;
	}
	if (z.Max_TC < a_M_TC)
	{
		z.Max_TC = a_M_TC;
	}
	if (z.Max_Lmax < b_M_Lmax)
	{
		z.Max_Lmax = b_M_Lmax;
	}

}
//
void NSGA::Transform1(vector<chrom>& T)
{
	for (auto iter = T.begin(); iter != T.end(); ++iter)
	{
		iter->ChromFitness.T_TC = (iter->ChromFitness.TC - z.TC);
		iter->ChromFitness.T_Lmax = (iter->ChromFitness.Lmax - z.Lmax);
	}
}
void NSGA::Transform2(vector<chrom>& T, const double& a, const double& b)
{
	for (auto iter = T.begin(); iter != T.end(); ++iter)
	{
		iter->ChromFitness.T_Lmax /= a;
		iter->ChromFitness.T_TC /= b;
	}
}

vector<int> NSGA::Find_Extreme_Point(vector<chrom>& T)
{
	vector<NO> Ex_T1, Ex_T2;
	vector<int> F_Ind;
	for (auto i = T.begin(); i != T.end(); ++i)
	{
		Ex_T1.push_back({ i->get_number(),max(i->ChromFitness.T_Lmax / 1.0,i->ChromFitness.T_TC / 1.0e-6) });
		Ex_T2.push_back({ i->get_number(),max(i->ChromFitness.T_Lmax / 1.0e-6,i->ChromFitness.T_TC / 1.0) });
	}
	auto ind1 = *min_element(Ex_T1.begin(), Ex_T1.end(), [](const NO& n1, const NO& n2)
		{
			return n1.O > n2.O;
		});
	auto ind2 = *min_element(Ex_T2.begin(), Ex_T2.end(), [](const NO& n1, const NO& n2)
		{
			return n1.O > n2.O;
		});
	F_Ind.push_back({ ind1.num });
	F_Ind.push_back({ ind2.num });
	return F_Ind;
}
void NSGA::Normalize(vector<chrom>& T)
{
	cal_z(T);
	Transform1(T);
	
}


double NSGA::cal_distance(const chrom& c, const NOO& n)
{
	double temp = abs(n.A * c.ChromFitness.T_TC - n.B * c.ChromFitness.T_Lmax) * 1.0 / pow((pow(n.A, 2) + pow(n.B, 2)), 0.5);
	return temp;
}


void NSGA::Encode(vector<job>& Jcode)//-----------------------------------------------------------给每个工件编码
{
	for (auto iter_jobs = Jcode.begin(); iter_jobs != Jcode.end(); ++iter_jobs)
	{
		double randomKey = 0.0;
		if (iter_jobs->get_s() <= M_s1)
		{
			randomKey = random_y(1, M_1 + M_2 + M_3 + M_4 + M_5 + 1);
		}
		else if (iter_jobs->get_s() <= M_s2)
		{
			randomKey = random_y(M_1 + 1, M_1 + M_2 + M_3 + M_4 + M_5 + 1);
		}
		else if (iter_jobs->get_s() <= M_s3)
		{
			randomKey = random_y(M_1 + M_2 + 1, M_1 + M_2 + M_3 + M_4 + M_5 + 1);
		}
		else if (iter_jobs->get_s() <= M_s4)
		{
			randomKey = random_y(M_1 + M_2 + M_3 + 1, M_1 + M_2 + M_3 + M_4 + M_5 + 1);
		}
		else
		{
			randomKey = random_y(M_1 + M_2 + M_3 + M_4 + 1, M_1 + M_2 + M_3 + M_4 + M_5 + 1);
		}
		iter_jobs->set_rand_key(randomKey);

	}
}

void NSGA::Decode(vector<chrom>& T)
{
	for (auto& i : T)
	{
		auto Temp = i.chro_sche(i.v_chromejobs);
		i.ChromFitness.TC = Temp.get_TC();
		i.ChromFitness.Lmax = Temp.get_Lmax();
	}
	
}


void NSGA::NSGA_schedule(vector<job> ZH_instacne, int max_iterator, int _instance, vector<chrom>& Parent, const int& Count)//--------------------------------迭代
{


	for (int iter_times = 1; iter_times <= max_iterator; ++iter_times)
	{
		if (iter_times == 1)
		{
			Parent = Init(ZH_instacne);
		}

		vector<chrom> Children = Init_Pure(ZH_instacne);

		selction_crossover(Parent, Children);
		mutation(Children);
		Decode(Children);
		auto Mix_pop = merge(Parent, Children);
		Parent = E_select(Mix_pop,iter_times,max_iterator);
	}
}



chrom& NSGA::find_p(const int& p, vector<chrom>& _P)
{
	for (chrom& m : _P)
	{
		if (m.get_number() == p)
		{
			return m;
		}
	}
}


vector<chrom> NSGA::find_Q(const int& i, vector<chrom>& _P)
{
	vector<chrom> temp_Q;
	for (auto& j : _P)
	{
		if (j.No_Rank == i)
		{
			temp_Q.push_back(j);
		}
	}
	return temp_Q;
}

void NSGA::Print_Pop(const vector<chrom>& T)
{
	int j = 1;
	for (auto& i : T)
	{
		cout << j << "-" << i.get_number() << " " << i.ChromFitness.Lmax << " " << i.ChromFitness.TC << " " << endl;
		++j;
	}
}

int NSGA::Roulette_Wheel_Selection(vector<NOO>& I)
{
	double s_val{0.0};
	for (auto& i : I)
	{
		if (i.A <= 0)
		{

		}
		s_val += i.A;
	}
	double temp{ 0.0 };
	for (auto& i : I)
	{
		i.B=temp+i.A/s_val;
		temp = i.B;
	}
	double fSlice = random_y(0,1);
	for (auto i = I.begin();i!=I.end();++i)
	{
		if (fSlice < i->B)
		{
			return i->num;
		}
	}
}

int NSGA::Machine_Selection(vector<machine>& _M, vector<NO_B>& machineList)
{
	int i = 1;
	vector<NOO> I;
	

	double bigTC = 0.0, bigLmax = 0.0;

	for (auto m = _M.begin(); m != _M.end(); ++m)
	{
		//m->cal_M();

		double val1 = 0.000001, val2 = 0.000001;

		if (m->v_Mbatches.size() != 0)
		{
			if (z.Lmax != 0.0 && z.TC != 0.0)
			{
				val2 = 1.0 / (pow(abs(m->v_Mbatches.rbegin()->get_Cb() - m->v_Mbatches.rbegin()->get_Db()) /z.Lmax, 0) +pow(m->get_L() * m->v_Mbatches.rbegin()->get_Cb() / z.TC,1));
			}
			else
			{
				val2 = 1.0 / m->get_L();
			}
		}
		else
		{
			val2 = 1.0 / m->get_L();
		}
		I.push_back({ m->get_num(),val2 ,0.0 });

	}
	int num = Roulette_Wheel_Selection(I);
	machineList.at(num - 1).O = true;
	return num;
}

int NSGA::Machine_Selection(vector<machine>& _M, int& M_num_x)
{
	int i = 1;
	vector<NOO> I;
	

	double bigTC = 0.0, bigLmax = 0.0;

	for (auto m = _M.begin(); m != _M.end(); ++m)
	{
			double val1 = 0.000001, val2 = 0.000001;
			val2 = 1.0 / m->get_L();
			I.push_back({ m->get_num(),val2 ,0.0 });

	}
	const int num = Roulette_Wheel_Selection(I);
	M_num_x = num;
	return num;
}

vector<job> NSGA::Cons_Cand_List(const vector<job>& J, const int& S)
{
	vector<job> temp_j;
	for (auto& j : J)
	{
		if (j.get_s() <= S&&j.flag==false)
		{
			temp_j.push_back(j);
		}
	}
	return temp_j;
}

void NSGA::Fast_Non_Dominated_Sort(vector<chrom>& P, int& R_num)
{
	vector <vector<chrom>> F(N - 1);
	for (auto iter_P = P.begin(); iter_P != P.end(); ++iter_P)
	{
		//iter_P->Set_p.clear();
		//iter_P->N_p;
		for (auto iter_Q = P.begin(); iter_Q != P.end(); ++iter_Q)
		{
			if (*iter_Q != *iter_P)
			{
				if (Dominated_Compared(*iter_P, *iter_Q) == 1)
				{
					iter_P->Set_p.push_back(iter_Q->get_number());
				}
				else if (Dominated_Compared(*iter_Q, *iter_P) == 1)
				{
					++(iter_P->N_p);
				}
			}
		}
		if (iter_P->N_p == 0)
		{
			iter_P->No_Rank = 1;
			F.at(0).push_back(*iter_P);
		}
	}

	int i = 1;
	while (F.at(i - 1).size() != 0)
	{
		vector<chrom> Q;
		for (auto& _iter_P : F.at(i - 1))
		{
			for (auto& _iter_q : _iter_P.Set_p)
			{
				find_p(_iter_q, P).N_p--;
				if (find_p(_iter_q, P).N_p == 0)
				{
					find_p(_iter_q, P).No_Rank = i + 1;
					Q.push_back(find_p(_iter_q, P));
				}
			}
		}
		++i;
		F.at(i - 1) = Q;
	}
	R_num = i - 1;
	
}


bool NSGA::check_state(const vector<job>& J)
{
	int total = 0;
	for (auto& j : J)
	{
		if (j.flag == true)
		{
			++total;
		}
	}
	if (total == Job_Num)
	{
		return false;
	}
	else
	{
		return true;
	}
}

job NSGA::choose_job(vector<job>& _J)
{
	double val{ 0.0 };
	vector<NOO> temp_N;
	int i = 1;
	for (auto& j : _J)
	{
		val += j.get_s();
		temp_N.push_back({i++,val,0.0});
	}

	return Roulette_Wheel_Selection(temp_N);
}


vector<chrom> NSGA::Init(vector<job>& J)
{
	vector<chrom> temp_P;
	for (int i = 1; i <= N; ++i)
	{
		auto TempChrom = generateSolution(J, i);
		auto Temp = TempChrom.chro_sche(TempChrom.v_chromejobs);
		TempChrom.ChromFitness.TC = Temp.get_TC();
		TempChrom.ChromFitness.Lmax = Temp.get_Lmax();
		temp_P.push_back(TempChrom);
		cal_z(temp_P);
	}
	return temp_P;

}

chrom NSGA::generateSolution(vector<job>& J, const int& num)
{
	for (auto& j : J)
	{
		j.flag = false;
	}

	scheme S;
	for (int k = 1; k <= M; ++k)
	{
		S.v_Smachines.push_back({ machine(k) });
	}

	vector<NO_B> mList;
	for (int i = 1; i <= M; ++i)
	{
		mList.push_back(NO_B(i));
	}
	int x = -1;
	while(check_state(J))
	{
		vector<machine> MS,MS_de;
		for(auto j=S.v_Smachines.begin();j!=S.v_Smachines.end();++j)
		{
			if(Cons_Cand_List(J, j->get_s()).size() != 0)
			{
				MS.push_back(*j);
			}
		}
		for(const auto& i:MS)
		{
			MS_de.push_back(i);
		}
		int M_num = Machine_Selection(MS_de, x);
		
		batch temp;
		int v_B_size = static_cast<int>(S.v_Smachines.at(M_num - 1).v_Mbatches.size());
		temp.Init(v_B_size + 1, M_num);
		S.v_Smachines.at(M_num - 1).v_Mbatches.push_back(temp);

		while (Cons_Cand_List(J, S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).get_rs()).size() != 0)
		{

			vector<job> temp_L = Cons_Cand_List(J, S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).get_rs());

			if (S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).get_rs() == S.v_Smachines.at(M_num - 1).get_s())
			{
				for (auto j = temp_L.begin(); j != temp_L.end(); ++j)
				{
					j->selection_feature = (static_cast<double>(j->get_r()) + j->get_p())* abs(static_cast<double>(j->get_r()) + j->get_p() - j->get_d());
				}
				sort(temp_L.begin(), temp_L.end(), [](const job& job1, const job& job2)
					{
						return job1.selection_feature < job2.selection_feature;
					});
				vector<job> select_first_jobs;
				for (auto j = temp_L.begin(); j != temp_L.end(); ++j)
				{
					if (j->selection_feature == temp_L.begin()->selection_feature)
					{
						select_first_jobs.push_back(*j);
					}
				}
				job a = select_first_jobs.at(min_element(select_first_jobs.begin(), select_first_jobs.end(), [S, M_num, v_B_size](const job& job1, const job& job2)
					{
						return job1.get_d() < job2.get_d();
					}) - select_first_jobs.begin());
				job& temp_j = find_j(a.get_num(), J);
				temp_j.flag = true;
				if (temp_j.get_num() < 0 || temp_j.get_num() > 550)
				{
					cout << endl;
				}
				S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).v_Bjobs.push_back(a);
				S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).update_batch(temp_j, S.v_Smachines.at(M_num - 1).get_v());
				
			}
			else
			{

				vector<job> temp_L1;
				for (const auto& j : temp_L)
				{

					if ((S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).get_Cb() -
						max(static_cast<double>(j.get_r()), S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).get_Sb())) - j.get_p() >= 0)
					{
						temp_L1.push_back(j);
					}
				}
				//begin
				vector<NOO> sel_jobL;
				double val{ 0.0 };
				///////////////启发式信息//////////////////
				if (temp_L1.size() == 0)
				{
					temp_L1 = temp_L;
				}
				int max_size = 0;
				for(const auto& j : temp_L1)
				{
					if(j.get_s()>max_size)
					{
						max_size = j.get_s();
					}
				}
				for (const auto& j : temp_L1)
				{
					double val_t = 1.0 / (abs(S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).get_Pb() - j.get_p())/ S.v_Smachines.at(M_num - 1).get_v() * (j.get_s()) + 1);
					
					sel_jobL.push_back({ j.get_num(),val_t,0.0 });
				}
				////////////////////////////////////////
				int j_num = Roulette_Wheel_Selection(sel_jobL);
				//end
				if (j_num > Job_Num)
				{
					cout << j_num << endl;
				}
				job& temp_j = find_j(j_num, J);
				temp_j.flag = true;
				if (temp_j.get_num() < 0 || temp_j.get_num() > 550)
				{
					cout << endl;
				}
				S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).v_Bjobs.push_back(temp_j);
				S.v_Smachines.at(M_num - 1).v_Mbatches.at(v_B_size).update_batch(temp_j, S.v_Smachines.at(M_num - 1).get_v());
			}



		}

	}
	for (auto s = S.v_Smachines.begin(); s != S.v_Smachines.end(); ++s)
	{
	
		s->cal_M();
	}
	S.cal_Fitness();


	chrom TempChrom;
	TempChrom.set_number(num);
	for (auto s = S.v_Smachines.begin(); s != S.v_Smachines.end(); ++s)
	{
		vector<job> temp_jobs;
		for (auto k = s->v_Mbatches.begin(); k != s->v_Mbatches.end(); ++k)
		{
			for (auto j = k->v_Bjobs.begin(); j != k->v_Bjobs.end(); ++j)
			{
				if (j->get_s() > s->get_s())
				{
					cout << endl;
				}
				temp_jobs.push_back(*j);
			}
		}
		vector<double> val_set;
		for (int i = 0; i < temp_jobs.size(); ++i)
		{
			val_set.push_back({ random_y(0, 1) });
		}
		sort(val_set.begin(), val_set.end());
	
		vector<job>::iterator iter_job;
		vector<double>::iterator iter_d;
		for (iter_job = temp_jobs.begin(), iter_d = val_set.begin();
			iter_job != temp_jobs.end() && iter_d != val_set.end(); ++iter_job, ++iter_d)
		{
			iter_job->set_rand_key(s->get_num()+*iter_d);
			TempChrom.v_chromejobs.push_back(*iter_job);
		}
	
	}
	return TempChrom;
}

vector<chrom> NSGA::Init_Pure(const vector<job>& J)
{
	vector<chrom> temp_P;
	for (int i = 1; i <= N; ++i)
	{
		chrom TempChrom;
		TempChrom.set_number(i);

		TempChrom.v_chromejobs = J;

		temp_P.push_back(TempChrom);
	}
	return temp_P;
}

chrom& NSGA::tournament(chrom& ind1, chrom& ind2)
{
	int B = Dominated_Compared(ind1, ind2);
	if (B == 1)
	{
		return ind1;
	}
	if (B == -1)
	{
		return ind2;
	}

	if (dist3(gen) <= 0.5)
	{
		return ind1;
	}
	else
	{
		return ind2;
	}
}



void NSGA::selction_crossover(vector<chrom>& T1, vector<chrom>& T2)//SBX
{
	int rand, temp;
	vector<int> t1, t2;
	for (int i = 1; i <= N; ++i)
	{
		t1.push_back(i);
		t2.push_back(i);
	}
	for (auto i = 0; i < N; ++i)
	{
		rand = static_cast<int>(dist2(gen));
		temp = t1[rand];
		t1[rand] = t1[i];
		t1[i] = temp;
		rand = static_cast<int>(dist2(gen));
		temp = t2[rand];
		t2[rand] = t2[i];
		t2[i] = temp;
	}
	double r{ 0.0 };
	chrom parent1, parent2;
	for (int i1 = 0; i1 < N; i1 += 4)
	{

		parent1 = tournament(find_p(t1[i1], T1), find_p(t1[i1 + 1], T1));
		parent2 = tournament(find_p(t1[i1 + 2], T1), find_p(t1[i1 + 3], T1));
		crossover(parent1, parent2, T2[i1], T2[i1 + 1]);
		parent1 = tournament(find_p(t2[i1], T1), find_p(t2[i1 + 1], T1));
		parent2 = tournament(find_p(t2[i1 + 2], T1), find_p(t2[i1 + 3], T1));
		crossover(parent1, parent2, T2[i1 + 2], T2[i1 + 3]);
		
	}
}

void NSGA::crossover(chrom& c1, chrom& c2, chrom& c3, chrom& c4)
{
	double r{ 0.0 };
	double u = dist3(gen);

	if (u <= 0.5)
	{
		r = pow(2 * u, (1.0 / SBX_Index));
	}
	else
	{
		r = pow(2 * (1 - u), (1.0 / SBX_Index));
	}

	for (int J_i = 0; J_i < Job_Num; ++J_i)
	{
		c3.v_chromejobs[J_i].set_rand_key(0.5 * (c1.v_chromejobs[J_i].get_rand_key() * (1 + r) +
			c2.v_chromejobs[J_i].get_rand_key() * (1 - r)));
		c4.v_chromejobs[J_i].set_rand_key(0.5 * (c1.v_chromejobs[J_i].get_rand_key() * (1 + r) +
			c2.v_chromejobs[J_i].get_rand_key() * (1 - r)));

		if (c3.v_chromejobs[J_i].get_s() <= M_s1)
		{
			while (c3.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 +M_4 +M_5 + 1))
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() - (M_1 + M_2 + M_3 + M_4 + M_5));
				
			}
			while (c3.v_chromejobs[J_i].get_rand_key() < 1.0)
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}
		else if (c3.v_chromejobs[J_i].get_s() <= M_s2)
		{
			while (c3.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() - (M_2 + M_3 + M_4 + M_5));
				
			}
			while (c3.v_chromejobs[J_i].get_rand_key() < M_1 + 1.0)
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}
		else if (c3.v_chromejobs[J_i].get_s() <= M_s3)
		{
			while (c3.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() - (M_3 + M_4 + M_5));

			}
			while (c3.v_chromejobs[J_i].get_rand_key() < M_1 + M_2 + 1.0)
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}
		else if (c3.v_chromejobs[J_i].get_s() <= M_s4)
		{
			while (c3.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() - (M_4 + M_5));

			}
			while (c3.v_chromejobs[J_i].get_rand_key() < M_1 + M_2 + M_3 + 1.0)
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}
		else
		{
			while (c3.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() - M_5);
				
			}
			while (c3.v_chromejobs[J_i].get_rand_key() < M_1 + M_2 + M_3 + M_4 + 1.0)
			{
				c3.v_chromejobs[J_i].set_rand_key(c3.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}
		/******************************************************************************************/
		if (c4.v_chromejobs[J_i].get_s() <= M_s1)
		{
			while (c4.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() - (M_1 + M_2 + M_3 + M_4 + M_5));
				
			}
			while (c4.v_chromejobs[J_i].get_rand_key() < 1.0)
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}
		else if (c4.v_chromejobs[J_i].get_s() <= M_s2)
		{
			while (c4.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() - (M_2 + M_3 + M_4 + M_5));
				
			}
			while (c4.v_chromejobs[J_i].get_rand_key() < M_1 + 1.0)
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}
		else if (c4.v_chromejobs[J_i].get_s() <= M_s3)
		{
			while (c4.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() - (M_3 + M_4 + M_5));

			}
			while (c4.v_chromejobs[J_i].get_rand_key() < M_1 + M_2 + 1.0)
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}
		else if (c4.v_chromejobs[J_i].get_s() <= M_s4)
		{
			while (c4.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() - (M_4 + M_5));

			}
			while (c4.v_chromejobs[J_i].get_rand_key() < M_1 + M_2 + M_3 + 1.0)
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}
		else
		{
			while (c4.v_chromejobs[J_i].get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() - M_5);
				
			}
			while (c4.v_chromejobs[J_i].get_rand_key() < M_1 + M_2 + M_3 + M_4 + 1.0)
			{
				c4.v_chromejobs[J_i].set_rand_key(c4.v_chromejobs[J_i].get_rand_key() + 1);
			}
		}

	}
}

void NSGA::mutation(vector<chrom>& T)//--------------------------------------polynomial mutation
{
	double B, u;
	for (auto i = T.begin(); i != T.end(); ++i)
	{
		u = dist3(gen);
		if (u < 0.5)
		{
			B = pow(2 * u, (1.0 / (P_mutation_Index + 1)) - 1);
		}
		else
		{
			B = 1 - pow(2 * (1 - u), (1.0 / (P_mutation_Index + 1)));
		}
		for (auto& j : i->v_chromejobs)
		{
			j.set_rand_key(j.get_rand_key() + B);
			if (j.get_s() <= M_s1)
			{
				while (j.get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
				{
					j.set_rand_key(j.get_rand_key() - (M_1 + M_2 + M_3 + M_4 + M_5));
					while (j.get_rand_key() < 1.0)
					{
						j.set_rand_key(j.get_rand_key() + 1);
					}
				}
			}
			else if (j.get_s() <= M_s2)
			{
				while (j.get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
				{
					j.set_rand_key(j.get_rand_key() - (M_2 + M_3 + M_4 + M_5));
					while (j.get_rand_key() < M_1 + 1.0)
					{
						j.set_rand_key(j.get_rand_key() + 1);
					}
				}
			}
			else if (j.get_s() <= M_s3)
			{
				while (j.get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
				{
					j.set_rand_key(j.get_rand_key() - (M_3 + M_4 + M_5));
					while (j.get_rand_key() < M_1 + M_2 + 1.0)
					{
						j.set_rand_key(j.get_rand_key() + 1);
					}
				}
			}
			else if (j.get_s() <= M_s4)
			{
				while (j.get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
				{
					j.set_rand_key(j.get_rand_key() - (M_4 + M_5));
					while (j.get_rand_key() < M_1 + M_2 + M_3 + 1.0)
					{
						j.set_rand_key(j.get_rand_key() + 1);
					}
				}
			}
			else
			{
				while (j.get_rand_key() > (M_1 + M_2 + M_3 + M_4 + M_5 + 1))
				{
					j.set_rand_key(j.get_rand_key() - M_5);
					while (j.get_rand_key() < M_1 + M_2 + M_3 + M_4 + 1.0)
					{
						j.set_rand_key(j.get_rand_key() + 1);
					}
				}
			}
		}
	}
	
}

vector<chrom> NSGA::merge(vector<chrom>& P1, vector<chrom>& P2)//P1父代，P2子代
{
	vector<chrom> temp_P(P1);
	for (auto iter_i = P2.cbegin(); iter_i != P2.cend(); ++iter_i)
	{
		temp_P.push_back(*iter_i);
	}
	for (int i = 0; i < temp_P.size(); ++i)
	{
		temp_P.at(i).set_number(i + 1);
	}
	return temp_P;
}


vector<chrom> NSGA::E_select(vector<chrom>& merged_population_T,int iter, int max_iter)//-------------------------------------------------------------Environmental Selection
{

	vector<chrom> copy_merged_population_T = merged_population_T;
	sort(copy_merged_population_T.begin(), copy_merged_population_T.end(), [](chrom& A, chrom& B)
		{
		
			return (A.ChromFitness.Lmax < B.ChromFitness.Lmax||(A.ChromFitness.Lmax == B.ChromFitness.Lmax&& A.ChromFitness.TC < B.ChromFitness.TC));
		});
	copy_merged_population_T.erase(unique(copy_merged_population_T.begin(), copy_merged_population_T.end(), [](chrom& A, chrom& B)
		{
			return (A.ChromFitness.TC == B.ChromFitness.TC && A.ChromFitness.Lmax == B.ChromFitness.Lmax);
		}), copy_merged_population_T.end());

	merged_population_T = copy_merged_population_T;

	
	vector<vector<chrom>> Q;

	int R_Num{ 0 };//分层总数
	//auto temp_T = T;
	vector<chrom> TT;
	for (auto& o : merged_population_T)
	{
		o.ass_angle1 = 0;
		o.ass_angle2 = 0;
		o.No_Rank = 0;
		o.N_p = 0;
		o.Set_p.clear();
	}
	Fast_Non_Dominated_Sort(merged_population_T, R_Num);
	int C_Num1 = 0, C_Num2 = 0;
	for (int i = 1; i <= R_Num; ++i)
	{
		vector<chrom> temp_T = find_Q(i, merged_population_T);
		vector<chrom> selectionT= temp_T;
		C_Num1 += static_cast<int>(temp_T.size());
		if (C_Num1 > N)
		{

			Normalize(temp_T);
			vector<chrom> temp_tran_T = temp_T;
			int k = N - C_Num2;

			selectionT.clear();
			
			ANGLE angle_cluster(k, temp_tran_T);
			
			vector<vector<chrom>> div_S=angle_cluster.cluster(k);
			
			
			selectionT = angle_cluster.selection_indidvial(div_S, static_cast<int>(temp_tran_T.size()),iter,max_iter,z);
			
		}
		Q.push_back(selectionT);
		C_Num2 += static_cast<int>(selectionT.size());
		if (C_Num2 == N)
		{
			break;
		}

	}
	int snum = 1;
	for (auto& k : Q)
	{
		for (auto& l : k)
		{
			l.set_number(snum);
			TT.push_back(l);
			++snum;
		}
	}
	return TT;
}
