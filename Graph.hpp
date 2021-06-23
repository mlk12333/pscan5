//
//  Graph.hpp
//  pscan-index5
//
//  Created by 孟令凯 on 2021/6/22.
//

#ifndef Graph_hpp
#define Graph_hpp

#include <iostream>
#include <stdio.h>
#include "string"
#include <ext/hash_map>
#include <vector>
#include <list>
#include <fstream>
#include <queue>
#include <set>

using namespace std;
using namespace __gnu_cxx;


class Graph {
private:
    string dir; //输入的图的地址
    
    int n, in_m, out_m,m,degree_m,edge_m;; //顶点数，入边数，出边数，总边数，最大度

    int eps_a1, eps_b1,eps_a2, eps_b2, miu; // eps_a2/eps_b2 = eps^2
    
    unsigned int *pstart,*out_pstart,*in_pstart;// 节点相邻点偏移量
    unsigned int *edges;// 相邻边
    unsigned int *reverse; // 反向边在边中的位置
    int *min_cn1, *min_cn2; //minimum common neighbor: -2 means not similar; -1 means similar; 0 means not sure; > 0 means the minimum common neighbor
    unsigned int *degree,*in_degree,*out_degree;//总度，入度，出度
    
    int *similar_degree, *effective_degree;//sd和ed
    
    int *pa,*rank;//为每个点找到所属的集群
    
    int *cid;
    
    unsigned int *out_edges,*in_edges;
    
    vector<pair<int,int>> noncore_cluster;//非核心点

public:
    Graph(const string _dir) ;
    ~Graph();
    void read_graph() ;//读图
    void pSCAN(const char *eps_s1,const char *eps_s2, int miu) ;//算法主流程
    void getDegree();//获取每个点的度数
    void get_eps(const char *eps_s1,const char *eps_s2);//获取eps
    void prune_and_cross_link(int eps_a1, int eps_b1,int eps_a2, int eps_b2, int miu, int *cores, int &cores_e) ;//预处理
//    int compute_common_neighbor_lowerbound(int u,int v,int eps_a2,int eps_b2, int k);
    unsigned int binary_search(unsigned int *array, int b, int e, int val);//array数组中是否存在val
    int find_root(int u);//找到u的根结点
    int check_common_neighbor(int u, int v, int c1,int c2);//获取公共邻居数

    int similar_check_OP(int u, int idx, int eps_a1, int eps_b1, int eps_a2, int eps_b2);//检查是否相似
    vector<int> get_f_and_c_num(int a,int b, int c);//找到abc之间的关系
    unsigned int BinarySearch(unsigned int *array, int val, int b, int e);//找到value位置
    void my_union(int u,int v);//将u和v归结到一个集群
    
    void cluster_noncore_vertices(int eps_a1, int eps_b1,int eps_a2, int eps_b2, int mu);//对非核心点聚类
    
    void output();//输出
};


#endif /* Graph_hpp */


