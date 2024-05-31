#include "Defines.h"
#include <algorithm>
#include <thread>
#include <stdlib.h>
#include <vector>
#include <mutex>
/*全局变量*/
int thread_num;//线程数
DataType ach_re0 = 0.00, ach_re1 = 0.00, ach_re2 = 0.00, ach_im0 = 0.00, ach_im1 = 0.00, ach_im2 = 0.00;//从循环里面提出来
std::mutex the_mutex;//互斥量
inline void correntess(ComplexType result1, ComplexType result2, ComplexType result3)
{
	double re_diff, im_diff;

	re_diff = fabs(result1.real() - -264241151.454552);
	im_diff = fabs(result1.imag() - 1321205770.975190);
	re_diff += fabs(result2.real() - -137405397.758745);
	im_diff += fabs(result2.imag() - 961837795.884157);
	re_diff += fabs(result3.real() - -83783779.241634);
	im_diff += fabs(result3.imag() - 754054017.424472);
	printf("%f,%f\n",re_diff,im_diff);
	printf("%f, %f\n", fabs(-264241151.454552 + -137405397.758745 + -83783779.241634) * 1e-6, fabs(1321205770.975190 + 961837795.884157 + 754054017.424472) * 1e-6);

	if (re_diff < fabs(result1.real() + result2.real() + result3.real()) * 1e-6 && im_diff < fabs(result1.imag() + result2.imag() + result3.imag()) * 1e-6)
	{
		printf("\n!!!! SUCCESS - !!!! Correctness test passed :-D :-D\n\n");
	}
	else
	{
		printf("\n!!!! FAILURE - Correctness test failed :-( :-(  \n");
	}
}
int main(int argc, char **argv) 
{
	int number_bands = 0, nvband = 0, ncouls = 0, nodes_per_group = 0;
	int npes = 1;
	if (argc == 2) {
		number_bands = 512;
		nvband = 2;
		ncouls = 32768;
		nodes_per_group = 20;
		thread_num = atoi(argv[1]);
		if (thread_num <= 0 || thread_num >100) 
		{
			std::cout << "wrong thread_num!\n";
			exit(0);
		}
	} else if (argc == 6) 
	{
		thread_num = atoi(argv[1]);
		if (thread_num <= 0 || thread_num >100) 
		{
			std::cout << "wrong thread_num!\n";
			exit(0);
		}
		number_bands = atoi(argv[2]);
		nvband = atoi(argv[3]);
		ncouls = atoi(argv[4]);
		nodes_per_group = atoi(argv[5]);
	} 
	else 
	{
		std::cout << "The correct form of input is : " << endl;
		std::cout << "./main.exe <thread_num>"<< endl;//稍微把输入错误时候的提示改了改
		std::cout << " ./main.exe <thread_num> <number_bands> <number_valence_bands> "
			"<number_plane_waves> <nodes_per_mpi_group> "
			<< endl;
		exit(0);
	}
	int ngpown = ncouls / (nodes_per_group * npes);

	// Constants that will be used later
	const DataType e_lk = 10;
	const DataType dw = 1;
	const DataType to1 = 1e-6;
	const DataType limittwo = pow(0.5, 2);
	const DataType e_n1kq = 6.0;

	// Using time point and system_clock
	time_point<system_clock> start, end, k_start, k_end;
	start = system_clock::now();
	double elapsedKernelTimer;

	// Printing out the params passed.
	std::cout << "Sizeof(ComplexType = "
		<< sizeof(ComplexType) << " bytes" << std::endl;
	std::cout << "number_bands = " << number_bands << "\t nvband = " << nvband
		<< "\t ncouls = " << ncouls
		<< "\t nodes_per_group  = " << nodes_per_group
		<< "\t ngpown = " << ngpown << "\t nend = " << nend
		<< "\t nstart = " << nstart << endl;

	size_t memFootPrint = 0.00;

	// ALLOCATE statements .
	ARRAY1D achtemp(nend - nstart);
	memFootPrint += (nend - nstart) * sizeof(ComplexType);

	ARRAY2D aqsmtemp(number_bands, ncouls);
	ARRAY2D aqsntemp(number_bands, ncouls);
	memFootPrint += 2 * (number_bands * ncouls) * sizeof(ComplexType);

	ARRAY2D I_eps_array(ngpown, ncouls);
	ARRAY2D wtilde_array(ngpown, ncouls);
	memFootPrint += 2 * (ngpown * ncouls) * sizeof(ComplexType);

	ARRAY1D_DataType vcoul(ncouls);
	memFootPrint += ncouls * sizeof(DataType);

	ARRAY1D_int inv_igp_index(ngpown);
	ARRAY1D_int indinv(ncouls + 1);
	memFootPrint += ngpown * sizeof(int);
	memFootPrint += (ncouls + 1) * sizeof(int);

	ARRAY1D_DataType wx_array(nend - nstart);
	memFootPrint += 3 * (nend - nstart) * sizeof(DataType);

	// Print Memory Foot print
	cout << "Memory Foot Print = " << memFootPrint / pow(1024, 3) << " GBs"
		<< endl;

	ComplexType expr(.5, .5);
	for (int i = 0; i < number_bands; i++)
		for (int j = 0; j < ncouls; j++) {
			aqsmtemp(i, j) = expr;
			aqsntemp(i, j) = expr;
		}

	for (int i = 0; i < ngpown; i++)
		for (int j = 0; j < ncouls; j++) {
			I_eps_array(i, j) = expr;
			wtilde_array(i, j) = expr;
		}

	for (int i = 0; i < ncouls; i++)
		vcoul(i) = 1.0;

	for (int ig = 0; ig < ngpown; ++ig)
		inv_igp_index(ig) = (ig + 1) * ncouls / ngpown;

	for (int ig = 0; ig < ncouls; ++ig)
		indinv(ig) = ig;
	indinv(ncouls) = ncouls - 1;

	for (int iw = nstart; iw < nend; ++iw) {
		wx_array(iw) = e_lk - e_n1kq + dw * ((iw + 1) - 2);
		if (wx_array(iw) < to1)
			wx_array(iw) = to1;
	}

	k_start = system_clock::now();
	noflagOCC_solver(number_bands, ngpown, ncouls, inv_igp_index, indinv,
			wx_array, wtilde_array, aqsmtemp, aqsntemp, I_eps_array,
			vcoul, achtemp);

	k_end = system_clock::now();
	duration<double> elapsed = k_end - k_start;
	elapsedKernelTimer = elapsed.count();

	// Check for correctness
	correntess(achtemp(0),achtemp(1),achtemp(2));	

	printf("\n Final achtemp\n");
	ComplexType_print(achtemp(0));
	ComplexType_print(achtemp(1));
	ComplexType_print(achtemp(2));

	end = system_clock::now();
	elapsed = end - start;

	cout << "********** Kernel Time Taken **********= " << elapsedKernelTimer
		<< " secs" << endl;
	cout << "********** Total Time Taken **********= " << elapsed.count()
		<< " secs" << endl;

	return 0;
}
//核心函数，每个线程都调用这个函数
void thread_function(size_t ngpown_start, size_t ngpown_end,size_t number_bands, size_t ncouls,
		ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
		ARRAY1D_DataType &wx_array, ARRAY2D &wtilde_array,
		ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
		ARRAY2D &I_eps_array, ARRAY1D_DataType &vcoul,
		ARRAY1D &achtemp)
{
	register  DataType thread_ach_re[3], thread_ach_im[3];//尽量不要用通过内存寻址访问，而是通过寄存器
	thread_ach_im[0] = thread_ach_im[1] = thread_ach_im[2] = thread_ach_re[0] = thread_ach_re[1] = thread_ach_re[2] = 0.00;
	bool can_accelerate=true;
	if(ngpown_end-ngpown_start>=256)
		can_accelerate=false;
	if(can_accelerate)
	{
		ComplexType conj_array[number_bands][ngpown_end-ngpown_start];
		int igp_array[ngpown_end-ngpown_start];
		for(int i=0;i<number_bands;++i)
		{
			for(int j=0;j<ngpown_end-ngpown_start;j++)
			{
				igp_array[j] = indinv(inv_igp_index(j+ngpown_start));
				conj_array[i][j]=ComplexType_conj(aqsmtemp(i,igp_array[j]))*aqsntemp(i,igp_array[j]);
			}
		}
		for (int thread_igp = ngpown_start; thread_igp < ngpown_end; ++thread_igp) 
		{
			for (int n1 = 0; n1 < number_bands; ++n1) 
			{
				for (int ig = 0; ig < ncouls; ++ig) 
				{
					DataType achtemp_re_loc[nend - nstart]={0.00}, achtemp_im_loc[nend - nstart]={0.00};
					ComplexType sch_store1 =conj_array[n1][thread_igp-ngpown_start] * wtilde_array(thread_igp,igp_array[thread_igp-ngpown_start])*I_eps_array(thread_igp, ig);
					for (int iw = nstart; iw < nend; ++iw)
					{
						ComplexType wdiff =	wx_array(iw) - wtilde_array(thread_igp, ig);
						DataType wdiff_real=wdiff.real();
						DataType wdiff_imag = wdiff.imag();
						DataType asas =1.0/ (wdiff_real*wdiff_real + wdiff_imag* wdiff_imag);
						ComplexType delw =ComplexType_conj(wdiff) ;
						ComplexType sch_array =	delw * sch_store1;
						achtemp_re_loc[iw] += (sch_array).real()* 0.5* asas *	vcoul(igp_array[thread_igp-ngpown_start]);
						achtemp_im_loc[iw] += (sch_array).imag()* 0.5* asas *	vcoul(igp_array[thread_igp-ngpown_start]);
					}
					thread_ach_re[0] += achtemp_re_loc[0];
					thread_ach_re[1] += achtemp_re_loc[1];
					thread_ach_re[2] += achtemp_re_loc[2];
					thread_ach_im[0] += achtemp_im_loc[0];
					thread_ach_im[1] += achtemp_im_loc[1];
					thread_ach_im[2] += achtemp_im_loc[2];
				}
			} 
		}
	}
	else
	{
		for (int thread_igp = ngpown_start; thread_igp < ngpown_end; ++thread_igp) 
		{
			for (int n1 = 0; n1 < number_bands; ++n1) 
			{
				for (int ig = 0; ig < ncouls; ++ig) 
				{
					DataType achtemp_re_loc[nend - nstart]={0.00}, achtemp_im_loc[nend - nstart]={0.00};
					int indigp = inv_igp_index(thread_igp);
					int igp = indinv(indigp);
					ComplexType sch_store1 = ComplexType_conj(aqsmtemp(n1, igp)) * aqsntemp(n1, igp) * wtilde_array(thread_igp,igp)*I_eps_array(thread_igp, ig);
					for (int iw = nstart; iw < nend; ++iw)
					{
						ComplexType wdiff =	wx_array(iw) - wtilde_array(thread_igp, ig);
						DataType wdiff_real=wdiff.real();
						DataType wdiff_imag = wdiff.imag();
						DataType asas =1.0/ (wdiff_real*wdiff_real + wdiff_imag* wdiff_imag);
						ComplexType delw =ComplexType_conj(wdiff) ;
						ComplexType sch_array =	delw * sch_store1;
						achtemp_re_loc[iw] += (sch_array).real()* 0.5* asas *	vcoul(igp);
						achtemp_im_loc[iw] += (sch_array).imag()* 0.5* asas *	vcoul(igp);
					}
					thread_ach_re[0] += achtemp_re_loc[0];
					thread_ach_re[1] += achtemp_re_loc[1];
					thread_ach_re[2] += achtemp_re_loc[2];
					thread_ach_im[0] += achtemp_im_loc[0];
					thread_ach_im[1] += achtemp_im_loc[1];
					thread_ach_im[2] += achtemp_im_loc[2];
				}
			} 
		}
	}
	//进入临界区
	the_mutex.lock();
	ach_re0 += thread_ach_re[0];
	ach_re1 += thread_ach_re[1];  
	ach_re2 += thread_ach_re[2];
	ach_im0 += thread_ach_im[0];
	ach_im1 += thread_ach_im[1];
	ach_im2 += thread_ach_im[2];
	the_mutex.unlock();
	//离开临界区
}
void noflagOCC_solver(size_t number_bands, size_t ngpown, size_t ncouls,
		ARRAY1D_int &inv_igp_index, ARRAY1D_int &indinv,
		ARRAY1D_DataType &wx_array, ARRAY2D &wtilde_array,
		ARRAY2D &aqsmtemp, ARRAY2D &aqsntemp,
		ARRAY2D &I_eps_array, ARRAY1D_DataType &vcoul,
		ARRAY1D &achtemp)
{
	time_point<system_clock> start, end;
	start = system_clock::now();
	std::cout << " ----------duoxiancheng begin----------\n";
	std::cout << "----------the number of xiancheng is " << thread_num << "----------\n";
	if (thread_num == 1)//单线程
	{
		std::cout<<"WARNING: danxiancheng running!\n";
		ach_re0 = 0.00, ach_re1 = 0.00, ach_re2 = 0.00, ach_im0 = 0.00, ach_im1 = 0.00, ach_im2 = 0.00;
		thread_function(0, ngpown,number_bands, ncouls, 
				std::ref(inv_igp_index ), std::ref(indinv), std::ref(wx_array), std::ref(wtilde_array), 
				std::ref(aqsmtemp), std::ref(aqsntemp), std::ref(I_eps_array), std::ref(vcoul),  std::ref(achtemp));
	}
	else//多线程了
	{
		std::vector<std::thread> thread_vector;
		time_point<system_clock>  time_start, time_end; 
		duration<double> time_duration ;
		time_start = system_clock::now();
		size_t ngpown_size = (ngpown / thread_num) + 1, ngpown_start = 0, ngpown_end = ngpown_size;
		for (int i = 0; i < thread_num; ++i) 
		{
			thread_vector.push_back(std::thread(thread_function, ngpown_start, ngpown_end, number_bands,  ncouls,std::ref(inv_igp_index ), 
				std::ref(indinv), std::ref(wx_array), std::ref(wtilde_array),std::ref(aqsmtemp), std::ref(aqsntemp),
				 std::ref(I_eps_array), std::ref(vcoul),  std::ref(achtemp)));
			ngpown_start += ngpown_size, ngpown_end += ngpown_size;
			if (ngpown_end > ngpown) ngpown_end = ngpown;//如果不能整除分给各个线程，最后一个end要缩小到ngpown才不会溢出
			if (ngpown_start >= ngpown_end) break;
		}
		ach_re0 = ach_re1 = ach_re2 = ach_im0 = ach_im1 = ach_im2 = 0.00;
		for (auto &each_thread : thread_vector)
			each_thread.join();
		thread_vector.clear();
		time_end = system_clock::now();
		time_duration = time_end - time_start;
		std::cout << "##########duoxiancheng time:" << time_duration.count() << "\n";
	}
	std::cout << "-----------duoxiancheng end----------\n";
	achtemp(0) = ComplexType(ach_re0, ach_im0);
	achtemp(1) = ComplexType(ach_re1, ach_im1);
	achtemp(2) = ComplexType(ach_re2, ach_im2);
}