#include "Defines.h"
#include <omp.h>

inline void correntess(ComplexType result1, ComplexType result2, ComplexType result3)
{
    double re_diff, im_diff;
    int numThreads;
#pragma omp parallel
    {
        int ttid = omp_get_thread_num();
        if (ttid == 0)
            numThreads = omp_get_num_threads();
    }
    printf("here are %d threads \n", numThreads);
    if (numThreads <= 64)
    {
        re_diff = fabs(result1.real() - -264241151.454552);
        im_diff = fabs(result1.imag() - 1321205770.975190);
        re_diff += fabs(result2.real() - -137405397.758745);
        im_diff += fabs(result2.imag() - 961837795.884157);
        re_diff += fabs(result3.real() - -83783779.241634);
        im_diff += fabs(result3.imag() - 754054017.424472);
        printf("%f,%f\n", re_diff, im_diff);
    }
    else
    {
        re_diff = fabs(result1.real() - -264241151.200123);
        im_diff = fabs(result1.imag() - 1321205763.246570);
        re_diff += fabs(result2.real() - -137405398.773852);
        im_diff += fabs(result2.imag() - 961837794.726070);
        re_diff += fabs(result3.real() - -83783779.939936);
        im_diff += fabs(result3.imag() - 754054018.099450);
    }
    if (re_diff < 10 && im_diff < 10)
        printf("\n!!!! SUCCESS - !!!! Correctness test passed :-D :-D\n\n");
    else
        printf("\n!!!! FAILURE - Correctness test failed :-( :-(  \n");
}

// 主函数
int main(int argc, char **argv)
{

    int number_bands = 0, nvband = 0, ncouls = 0, nodes_per_group = 0;
    int npes = 1;
    if (argc == 1)
    {
        number_bands = 512;
        nvband = 2;
        ncouls = 32768;
        nodes_per_group = 20;
    }
    else if (argc == 5)
    {
        number_bands = atoi(argv[1]);
        nvband = atoi(argv[2]);
        ncouls = atoi(argv[3]);
        nodes_per_group = atoi(argv[4]);
    }
    else
    {
        std::cout << "The correct form of input is : " << endl;
        std::cout << " ./main.exe <number_bands> <number_valence_bands> "
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

    // ================Total Time Start===============================================

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
    int size = sizeof(ComplexType);
    int stepLength = nend - nstart;
    
    ComplexType achtemp[stepLength];
    memFootPrint += (stepLength) * size;

    int tmpLength1 = number_bands * ncouls;
    ComplexType *aqsmtemp = new ComplexType[tmpLength1];
    ComplexType *aqsntemp = new ComplexType[tmpLength1];
    memFootPrint += 2 * tmpLength1 * size;

    int tmpLength2 = ngpown * ncouls;
    ComplexType* I_eps_array = new ComplexType[tmpLength2];
    ComplexType* wtilde_array = new ComplexType[tmpLength2];
    memFootPrint += 2 * tmpLength2 * size;

    
    DataType vcoul[ncouls];
    memFootPrint += ncouls * sizeof(DataType);

    int inv_igp_index[ngpown];
    int indinv[ncouls + 1];
    memFootPrint += ngpown * sizeof(int);
    memFootPrint += (ncouls + 1) * sizeof(int);

    DataType wx_array[stepLength];
    memFootPrint += 3 * (stepLength) * sizeof(DataType);

    // Print Memory Foot print
    cout << "Memory Foot Print = " << memFootPrint / pow(1024, 3) << " GBs"
         << endl;

    ComplexType expr(.5, .5);
    // 初始化赋值
    for (int i = 0; i < number_bands; i++)
        for (int j = 0; j < ncouls; j++)
        {
            aqsmtemp[i * ncouls + j] = expr;
            aqsntemp[i * ncouls + j] = expr;
        }

    // 初始化赋值
    for (int i = 0; i < ngpown; i++)
        for (int j = 0; j < ncouls; j++)
        {
            I_eps_array[i * ncouls + j] = expr;
            wtilde_array[i * ncouls + j] = expr;
        }

    for (int i = 0; i < ncouls; i++)
        vcoul[i] = 1.0;

    for (int ig = 0; ig < ngpown; ++ig)
        inv_igp_index[ig] = (ig + 1) * ncouls / ngpown;

    for (int ig = 0; ig < ncouls; ++ig)
        indinv[ig] = ig;
    indinv[ncouls] = ncouls - 1;

    for (int iw = nstart; iw < nend; ++iw)
    {
        wx_array[iw] = e_lk - e_n1kq + dw * ((iw + 1) - 2);
        if (wx_array[iw] < to1)
            wx_array[iw] = to1;
    }

    // ================Kernel Time Start=============
    k_start = system_clock::now();
    noflagOCC_solver(number_bands, ngpown, ncouls, inv_igp_index, indinv,
                     wx_array, wtilde_array, aqsmtemp, aqsntemp, I_eps_array,
                     vcoul, achtemp);

    // ================Kernel Time End=============
    k_end = system_clock::now();
    duration<double> elapsed = k_end - k_start;
    elapsedKernelTimer = elapsed.count();

    // Check for correctness
    correntess(achtemp[0], achtemp[1], achtemp[2]);
    printf("\n Final achtemp\n");
    ComplexType_print(achtemp[0]);
    ComplexType_print(achtemp[1]);
    ComplexType_print(achtemp[2]);

    // ================Total Time End=============
    end = system_clock::now();
    elapsed = end - start;

    cout << "********** Kernel Time Taken **********= " << elapsedKernelTimer
         << " secs" << endl;
    cout << "********** Total Time Taken **********= " << elapsed.count()
         << " secs" << endl;
    delete[] aqsmtemp;
    delete[] aqsntemp;
    delete[] I_eps_array;
    delete[] wtilde_array;
    return 0;
}

void noflagOCC_solver(size_t number_bands, size_t ngpown, size_t ncouls,
                      int *inv_igp_index, int *indinv,
                      DataType *wx_array, ComplexType *wtilde_array,
                      ComplexType *aqsmtemp, ComplexType *aqsntemp,
                      ComplexType *I_eps_array, DataType *vcoul,
                      ComplexType *achtemp)
{
    // time_point<system_clock> start, end;
    // start = system_clock::now();
    // Vars to use for reduction
    DataType ach_re0 = 0.00, ach_re1 = 0.00, ach_re2 = 0.00, ach_im0 = 0.00,
             ach_im1 = 0.00, ach_im2 = 0.00;

#pragma omp parallel for reduction(+ \
                                   : ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2)
    for (int n1 = 0; n1 < number_bands; ++n1)
    {
        for (int ig = 0; ig < ncouls; ++ig)
        {
            for (int my_igp = 0; my_igp < ngpown; ++my_igp)
            {
                int indigp = inv_igp_index[my_igp];
                int igp = indinv[indigp];
                DataType achtemp_re_loc[nend - nstart], achtemp_im_loc[nend - nstart];
                for (int iw = nstart; iw < nend; ++iw)
                {
                    achtemp_re_loc[iw] = 0.00;
                    achtemp_im_loc[iw] = 0.00;
                }
                ComplexType sch_store1 =
                    ComplexType_conj(aqsmtemp[n1 * ncouls + igp]) * aqsntemp[n1 * ncouls + igp] * 0.5 *
                    vcoul[igp] * wtilde_array[my_igp * ncouls + igp];

                for (int iw = nstart; iw < nend; ++iw)
                {
                    ComplexType wdiff =
                        wx_array[iw] - wtilde_array[my_igp * ncouls + ig];
                    ComplexType delw =
                        ComplexType_conj(wdiff) *
                        (1 / (wdiff * ComplexType_conj(wdiff)).real());
                    ComplexType sch_array =
                        delw * I_eps_array[my_igp * ncouls + ig] * sch_store1;

                    achtemp_re_loc[iw] += (sch_array).real();
                    achtemp_im_loc[iw] += (sch_array).imag();
                }
                ach_re0 += achtemp_re_loc[0];
                ach_re1 += achtemp_re_loc[1];
                ach_re2 += achtemp_re_loc[2];
                ach_im0 += achtemp_im_loc[0];
                ach_im1 += achtemp_im_loc[1];
                ach_im2 += achtemp_im_loc[2];
            }
        }
    }

    achtemp[0] = ComplexType(ach_re0, ach_im0);
    achtemp[1] = ComplexType(ach_re1, ach_im1);
    achtemp[2] = ComplexType(ach_re2, ach_im2);
}
