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
    memFootPrint += (stepLength)*size;

    int tmpLength1 = number_bands * ncouls;
    // ComplexType *aqsmtemp = new ComplexType[tmpLength1];
    // 2^6 = 64
    // DataType *aqsmtemp_a = new DataType[tmpLength1];
    // DataType *aqsmtemp_b = new DataType[tmpLength1];
    DataType *aqsmtemp_a = (DataType*)aligned_alloc(32, sizeof(DataType) * tmpLength1 );
    DataType *aqsmtemp_b = (DataType*)aligned_alloc(32, sizeof(DataType) * tmpLength1 );
    // ComplexType *aqsntemp = new ComplexType[tmpLength1];
    // DataType *aqsntemp_a = new DataType[tmpLength1];
    // DataType *aqsntemp_b = new DataType[tmpLength1];
    DataType *aqsntemp_a = (DataType*)aligned_alloc(32, sizeof(DataType) * tmpLength1 );
    DataType *aqsntemp_b = (DataType*)aligned_alloc(32, sizeof(DataType) * tmpLength1 );
    memFootPrint += 2 * tmpLength1 * size;

    int tmpLength2 = ngpown * ncouls;
    // ComplexType* I_eps_array = new ComplexType[tmpLength2];
    // DataType *I_eps_array_a = new DataType[tmpLength2];
    // DataType *I_eps_array_b = new DataType[tmpLength2];
    DataType* I_eps_array_a = (DataType*)aligned_alloc(32, sizeof(DataType) * tmpLength2 );
    DataType* I_eps_array_b = (DataType*)aligned_alloc(32, sizeof(DataType) * tmpLength2 );
    // ComplexType* wtilde_array = new ComplexType[tmpLength2];
    // DataType *wtilde_array_a = new DataType[tmpLength2];
    // DataType *wtilde_array_b = new DataType[tmpLength2];
    DataType* wtilde_array_a = (DataType*)aligned_alloc(32, sizeof(DataType) * tmpLength2 );
    DataType* wtilde_array_b = (DataType*)aligned_alloc(32, sizeof(DataType) * tmpLength2 );
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
            int pos = i * ncouls + j;
            aqsmtemp_a[pos] = .5;
            aqsmtemp_b[pos] = .5;
            // aqsmtemp[i * ncouls + j] = expr;
            aqsntemp_a[pos] = .5;
            aqsntemp_b[pos] = .5;
        }

    // 初始化赋值
    double t = ncouls / ngpown;
    for (int i = 0; i < ngpown; i++)
    {
        inv_igp_index[i] = (i + 1) * t;
        for (int j = 0; j < ncouls; j++)
        {
            int pos = i * ncouls + j;
            I_eps_array_a[pos] = .5;
            I_eps_array_b[pos] = .5;
            wtilde_array_a[pos] = .5;
            wtilde_array_b[pos] = .5;
        }
    }

    for (int i = 0; i < ncouls; i++)
        vcoul[i] = 1.0;

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
                     wx_array, wtilde_array_a, wtilde_array_b,
                     aqsmtemp_a, aqsmtemp_b,
                     aqsntemp_a, aqsntemp_b,
                     I_eps_array_a, I_eps_array_b,
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
    free(aqsmtemp_a);
    free(aqsmtemp_b);
    free(aqsntemp_a);
    free(aqsntemp_b);
    free(I_eps_array_a);
    free(I_eps_array_b);
    free(wtilde_array_a);
    free(wtilde_array_b);
    return 0;
}

void noflagOCC_solver(size_t number_bands, size_t ngpown, size_t ncouls,
                      int *inv_igp_index, int *indinv,
                      DataType *wx_array,
                      DataType *wtilde_array_a, DataType *wtilde_array_b,
                      DataType *aqsmtemp_a, DataType *aqsmtemp_b,
                      DataType *aqsntemp_a, DataType *aqsntemp_b,
                      DataType *I_eps_array_a, DataType *I_eps_array_b,
                      DataType *vcoul,
                      ComplexType *achtemp)
{
    // time_point<system_clock> start, end;
    // start = system_clock::now();
    // Vars to use for reduction
    DataType ach_re0 = 0.00, ach_re1 = 0.00, ach_re2 = 0.00, ach_im0 = 0.00,
             ach_im1 = 0.00, ach_im2 = 0.00;

#pragma omp parallel for reduction(+ \
                                   : ach_re0, ach_re1, ach_re2, ach_im0, ach_im1, ach_im2)
    for (int my_igp = 0; my_igp < ngpown; ++my_igp)
    {
        int indigp = inv_igp_index[my_igp];
        // int indigp = (my_igp + 1) * ncouls / ngpown
        int igp = indinv[indigp];
        // int igp = indigp & ~ncouls ? indigp : ncouls -1;
        int pos_base = my_igp * ncouls;
        int pos2 = pos_base + igp;

        for (int n1 = 0; n1 < number_bands; ++n1)
        {
            int pos1 = n1 * ncouls + igp;
            for (int ig = 0; ig < ncouls; ++ig)
            {
                int pos3 = pos_base + ig;

                // printf("---igp: %d, pos1: %d, pos2: %d, pos3: %d\n", igp, pos1, pos2, pos3);
                DataType r0 = 0.5 * vcoul[igp] * wtilde_array_a[pos2];
                DataType i0 = 0.5 * vcoul[igp] * wtilde_array_b[pos2];
                DataType r1 = aqsntemp_a[pos1] * r0 - aqsntemp_b[pos1] * i0;
                DataType i1 = aqsntemp_a[pos1] * i0 + aqsntemp_b[pos1] * r0;
                DataType r2 = aqsmtemp_a[pos1] * r1 + aqsmtemp_b[pos1] * i1;
                DataType i2 = aqsmtemp_a[pos1] * i1 - aqsmtemp_b[pos1] * r1;

                DataType achtemp_re_loc[nend - nstart], achtemp_im_loc[nend - nstart];
                // for (int iw = nstart; iw < nend; ++iw)
                // {
                // achtemp_re_loc[iw] = 0.00;
                // achtemp_im_loc[iw] = 0.00;
                // }

                DataType b2 = wtilde_array_b[pos3] * wtilde_array_b[pos3];
#pragma unroll(3)
                for (int iw = nstart; iw < nend; ++iw)
                {
                    DataType wdiff_a = wx_array[iw] - wtilde_array_a[pos3];
                    DataType fm = wdiff_a * wdiff_a + b2;
                    DataType delw_a = wdiff_a / fm;
                    DataType delw_b = wtilde_array_b[pos3] / fm;

                    achtemp_re_loc[iw] = (delw_a * I_eps_array_a[pos3] - delw_b * I_eps_array_b[pos3]) * r2 - (delw_a * I_eps_array_b[pos3] + delw_b * I_eps_array_a[pos3]) * i2;
                    achtemp_im_loc[iw] = (delw_a * I_eps_array_a[pos3] - delw_b * I_eps_array_b[pos3]) * i2 + (delw_a * I_eps_array_b[pos3] + delw_b * I_eps_array_a[pos3]) * r2;
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
