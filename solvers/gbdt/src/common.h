#pragma GCC diagnostic ignored "-Wunused-result"

#ifndef _COMMON_H_
#define _COMMON_H_

#define flag { printf("\nLINE: %d\n", __LINE__); fflush(stdout); }

#include <cstdio>
#include <string>
#include <vector>
#include <cmath>

#include <pmmintrin.h>

struct Node
{
    Node() : i(0), v(0) {}
    Node(uint32_t const i, float const v) : i(i), v(v) {}
    uint32_t i;
    float v;
};

struct Problem
{
    //X: 保存按各个field由小至大排序
    //Z: 原始的各个instance的结果，未排序
    Problem() : nr_instance(0), nr_field(0), nr_sparse_field(0) {}
    Problem(uint32_t const nr_instance, uint32_t const nr_field) 
        : nr_instance(nr_instance), nr_field(nr_field), nr_sparse_field(0),
          X(nr_field, std::vector<Node>(nr_instance)), 
          Z(nr_field, std::vector<Node>(nr_instance)), 
          Y(nr_instance) {}
    uint32_t const nr_instance, nr_field;
    uint32_t nr_sparse_field;
    std::vector<std::vector<Node>> X, Z;

    //SJ保存稀疏特征
    //SJ：所有instance为1的位置如    [2,4,5,  1,5, 2]
    //SJP: 各个instance的边界，对以上[0,      3,   5,6]
    std::vector<uint32_t> SI, SJ;
    std::vector<uint64_t> SIP, SJP;
    std::vector<float> Y;
};

inline std::vector<float> 
construct_instance(Problem const &prob, uint32_t const i)
{
    uint32_t const nr_field = prob.nr_field; 
    uint32_t const nr_sparse_field = prob.nr_sparse_field;
    std::vector<uint32_t> const &SJ = prob.SJ;
    std::vector<uint64_t> const &SJP = prob.SJP;

    std::vector<float> x(nr_field+nr_sparse_field, 0);
    for(uint32_t j = 0; j < prob.nr_field; ++j)
        x[j] = prob.Z[j][i].v;
    for(uint64_t p = SJP[i]; p < SJP[i+1]; ++p)
        x[SJ[p]+nr_field] = 1;

    return x;
}

Problem read_data(std::string const &dense_path, 
    std::string const &sparse_path);

FILE *open_c_file(std::string const &path, std::string const &mode);

std::vector<std::string> 
argv_to_args(int const argc, char const * const * const argv);

#endif // _COMMON_H_
