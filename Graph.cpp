//
//  Graph.cpp
//  pscan-index5
//
//  Created by 孟令凯 on 2021/6/22.
//

#include "Graph.hpp"
Graph::Graph(const string _dir) {//初始化一些数据
    dir = _dir;

    n = in_m = out_m =m=degree_m=edge_m= 0;

    eps_a1 = eps_b1 = eps_a2 = eps_b2 = miu = 0;
    pstart = NULL;
    edges = NULL;
    in_edges = NULL;
    out_edges = NULL;
    reverse = NULL;
    min_cn1 = NULL;
    min_cn2 = NULL;
    cid = NULL;
    
    out_pstart=NULL;
    in_pstart=NULL;

    degree = NULL;
    in_degree = NULL;
    out_degree = NULL;
    effective_degree = NULL;
    similar_degree = NULL;

    pa = NULL;
    rank = NULL;
    
}

Graph::~Graph() {
    if(pstart != NULL) {
        delete[] pstart;
        pstart = NULL;
    }
    if(in_pstart != NULL) {
        delete[] in_pstart;
        in_pstart = NULL;
    }
    if(out_pstart != NULL) {
        delete[] out_pstart;
        out_pstart = NULL;
    }
    if(in_edges != NULL) {
        delete[] in_edges;
        in_edges = NULL;
    }
    if(out_edges != NULL) {
        delete[] out_edges;
        out_edges = NULL;
    }
    if(edges != NULL) {
        delete[] edges;
        edges = NULL;
    }
    if(reverse != NULL) {
        delete[] reverse;
        reverse = NULL;
    }
    if(min_cn1 != NULL) {
        delete[] min_cn1;
        min_cn1 = NULL;
    }
    if(min_cn2 != NULL) {
        delete[] min_cn2;
        min_cn2 = NULL;
    }
    if(cid != NULL) {
        delete[] cid;
        cid = NULL;
    }
    if(degree != NULL) {
        delete[] degree;
        degree = NULL;
    }
    if(in_degree != NULL) {
        delete[] in_degree;
        in_degree = NULL;
    }
    if(out_degree != NULL) {
        delete[] out_degree;
        out_degree = NULL;
    }
    if(effective_degree != NULL) {
        delete[] effective_degree;
        effective_degree = NULL;
    }
    if(similar_degree != NULL) {
        delete[] similar_degree;
        similar_degree = NULL;
    }
    if(pa != NULL) {
        delete[] pa;
        pa = NULL;
    }
    if(rank != NULL) {
        delete[] rank;
        rank = NULL;
    }
}

void Graph::read_graph(){
    ifstream infile;   //输入流
    
    infile.open(dir+"/degree.txt", ios::in);
    if (!infile.is_open()){
        cout<<"Open degree file failure"<<endl;
        exit(0);
    }
    infile>>n>>m;
    
    if(out_pstart == NULL) out_pstart = new unsigned int[n+1];
    if(in_pstart == NULL) in_pstart = new unsigned int[n+1];
    if(out_edges == NULL) out_edges = new unsigned int[m];
    if(in_edges == NULL) in_edges = new unsigned int[m];
    
     
    out_pstart[0] = 0;
    in_pstart[0] = 0;
    
    int outD, inD, starti = 0;
    
    while (!infile.eof())            // 若未到文件结束一直循环
    {
        infile >> outD >> inD;
        out_pstart[starti+1] = outD + out_pstart[starti];
        in_pstart[starti+1] = inD + in_pstart[starti];
        starti++;
    }
    infile.close();
    
    infile.open(dir+"/out_edges.txt", ios::in);
    if (!infile.is_open()){
        cout<<"Open out_edges file failure"<<endl;
        exit(0);
    }
    int mm_;
    infile >> mm_;
    int outi = 0;
    while (!infile.eof())            // 若未到文件结束一直循环
    {
        infile >> out_edges[outi];
        outi++;
    }
    infile.close();
    
    infile.open(dir+"/in_edges.txt", ios::in);
    
    if (!infile.is_open()){
        cout<<"Open in_edges file failure"<<endl;
        exit(0);
    }
    int ini = 0;
    infile >> mm_;
    while (!infile.eof())            // 若未到文件结束一直循环
    {
        infile >> in_edges[ini];
        ini++;
    }
    infile.close();
    
    cout<<"read!"<<endl;
    getDegree();
    
    cout<<"getDegree!"<<endl;
    
}

void Graph::getDegree(){
    if(pstart == NULL) pstart = new unsigned int[n+1];
    if(degree == NULL) degree = new unsigned int[n];
    pstart[0] = 0;
    if(edges == NULL)  edges = new unsigned int[(long)2*m];
//    edges = new unsigned int[(long)2*m];
    long edgesi = 0;
    for(int i = 0;i<n;i++){
        int l1 = out_pstart[i+1] - out_pstart[i];
        int l2 = in_pstart[i+1] - in_pstart[i];
        int i1 = 0, i2 = 0, num = 0;
        while(i1<l1 && i2<l2){
            
            if(out_edges[out_pstart[i] + i1] == in_edges[in_pstart[i] + i2]){
                edges[edgesi] = out_edges[out_pstart[i] + i1];
                edgesi++;
                i1++;
                i2++;
                num++;
            }else if(out_edges[out_pstart[i] + i1] < in_edges[in_pstart[i] + i2]){
                edges[edgesi] = out_edges[out_pstart[i] + i1];
                edgesi++;
                i1++;
                num++;
            }else{
                edges[edgesi] = in_edges[in_pstart[i] + i2];
                edgesi++;
                i2++;
                num++;
            }
        }
        while(i1<l1){
            edges[edgesi] = out_edges[out_pstart[i] + i1];
            edgesi++;
            i1++;
            num++;
        }
        while(i2<l2){
            edges[edgesi] = in_edges[in_pstart[i] + i2];
            edgesi++;
            i2++;
            num++;
        }
        
        degree[i] = num;
        edge_m = edge_m + num;
        pstart[i+1] = pstart[i] + num;
    }
    
    reverse = new unsigned int[(long)edge_m];
    if(min_cn1 == NULL) min_cn1 = new int[(long)edge_m];
    memset(min_cn1, 0, sizeof(int)*(long)edge_m);
    if(min_cn2 == NULL) min_cn2 = new int[(long)edge_m];
    memset(min_cn2, 0, sizeof(int)*(long)edge_m);
    
    if(cid == NULL) cid = new int[n];//用cid分离出有几个群体
    for(unsigned int i = 0;i < n;i ++) cid[i] = n;
}

void  Graph::pSCAN(const char *eps_s1,const char *eps_s2, int _miu){
    get_eps(eps_s1,eps_s2);
    miu = _miu;
    
    if(similar_degree == NULL) similar_degree= new int[n];
    memset(similar_degree, 0, sizeof(int)*n);
    if(effective_degree == NULL) effective_degree= new int[n];
    for(int i = 0;i < n;i ++){
//        similar_degree[i] = 0;
        effective_degree[i] = degree[i]-1;
    }
    if(pa == NULL) pa = new int[n];//pa和rank用来找core
    if(rank == NULL) rank = new int[n];
    for(int i = 0;i < n;i ++) {
        pa[i] = i;
        rank[i] = 0;
    }
    
    unsigned int *edge_buf = new unsigned int[n];
    int *cores = new int[n];
    int cores_n = 0;
    
    
    clock_t startTime,endTime;
    
    startTime = clock();//计时开始
   
    prune_and_cross_link(eps_a1, eps_b1, eps_a2, eps_b2, miu, cores, cores_n);//预处理，并得到达到相似要达到的要求
    printf("\t*** Finished prune and cross link!\n");
    
    endTime = clock();//计时结束
    cout << "The prune_and_cross_link time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    
    int *bin_head = new int[n];//记录第一个ed值为n的点
    int *bin_next = new int[n];//bin_next的每一个值指向下一个ed值相同的点，遇到-1结束
    for(int i = 0;i < n;i ++) bin_head[i] = -1;//初始值为-1
    int max_ed = 0;
    for(int i = 0;i < n;i++) if(effective_degree[i] >= miu) {//根据每个点的ed大小进行排序，从大到小开始探索
        int ed = effective_degree[i];
        if(ed > max_ed) max_ed = ed;//出现新的更大的ed值，标记一下
        bin_next[i] = bin_head[ed];
        bin_head[ed] = i;
    }
    while(true) {
        int u = -1;
        if(cores_n) u = cores[-- cores_n];//先找已经确认的核心点
        else {//根据之前的排序找到ed最大的点
            while(max_ed >= miu&&u == -1) {
                for(int x = bin_head[max_ed];x != -1;) {
                    int tmp = bin_next[x];
                    int ed = effective_degree[x];
                    if(ed == max_ed) {
                        u = x;
                        bin_head[max_ed] = bin_next[x];
                        break;
                    }
                    else if(ed >= miu) {
                        bin_next[x] = bin_head[ed];//超了，强制变回max_ed
                        bin_head[ed] = x;
                    }
                    x = tmp;
                }
                if(u == -1) {//度最大的找完了
                    bin_head[max_ed] = -1;
                    -- max_ed;
                }
            }
        }
        if(u == -1) break;//最大的度已经不满足要求了
        int edge_buf_n = 0;
        for(int j = pstart[u];j < pstart[u+1];j ++) {//查看当前点的所有邻居
            if(min_cn1[j] == -2) continue;
//            if(find_root(u) != find_root(edges[j])) edge_buf[edge_buf_n ++] = j;
            if(similar_degree[u] < miu||find_root(u) != find_root(edges[j]))edge_buf[edge_buf_n++] = j;            //edge_buf存储还没有确定关系的点以及他的邻居，edge_buf_n是个数
        }

        int i = 0;
        while(similar_degree[u] < miu&&effective_degree[u] >= miu&&i < edge_buf_n) {//还有机会相似的点
            unsigned int idx = edge_buf[i];//u的邻居的一部分
            if(min_cn1[idx] != -1) {//-1已经合格 不许再谈
#ifdef _DEBUG_
                if(min_cn1[idx] == 0) printf("WA min_cn!\n");
#endif
                int v = edges[idx];//具体的点
                min_cn1[idx] = min_cn1[reverse[idx]] = min_cn2[idx] = min_cn2[reverse[idx]] = similar_check_OP(u, idx, eps_a1, eps_b1, eps_a2, eps_b2);
                if(min_cn1[idx] == -1) ++ similar_degree[u];
                else -- effective_degree[u];

                if(effective_degree[v] >= 0) {
                    if(min_cn1[idx] == -1) {
                        ++ similar_degree[v];
                        if(similar_degree[v] == miu) cores[cores_n ++] = v;
                    }
                    else -- effective_degree[v];
                }
            }
            ++ i;
        }
        effective_degree[u] = -1;//标记为已检查过的点
        if(similar_degree[u] < miu) continue;//没希望成为核心了
        for(int j = 0;j < edge_buf_n;j++) {
            unsigned int idx = edge_buf[j];
            if(min_cn1[idx] == -1&&similar_degree[edges[idx]] >= miu)
                my_union(u, edges[idx]);//核心抱团
        }
        while(i < edge_buf_n) {//检查剩余的点
            unsigned int idx = edge_buf[i];
            int v = edges[idx];
            if(min_cn1[idx] < 0||similar_degree[v] < miu||find_root(u) == find_root(v)) {//没有确定相似，但已经是核心并且没抱在一起才能跳过continue，是核心，但是还不确定想不相思
                ++ i;
                continue;
            }
            min_cn1[idx] = min_cn1[reverse[idx]] = min_cn2[idx] = min_cn2[reverse[idx]] = similar_check_OP(u, idx, eps_a1, eps_b1, eps_a2, eps_b2);
            if(effective_degree[v] >= 0) {//说明还没有检查过
                if(min_cn1[idx] == -1) {
                    ++ similar_degree[v];
                    if(similar_degree[v] == miu) cores[cores_n ++] = v;
                }
                else -- effective_degree[v];
            }
            if(min_cn1[idx] == -1) my_union(u,v);//已经是核心并且相似 那就抱一起吧
            ++ i;
        }
        //printf(")\n");
    }
    delete[] edge_buf; edge_buf = NULL;
    delete[] cores; cores = NULL;
    delete[] bin_head; bin_head = NULL;
    delete[] bin_next; bin_next = NULL;
    cluster_noncore_vertices(eps_a1,eps_b1,eps_a2,eps_b2,miu);
}


void Graph::get_eps(const char *eps_s1,const char *eps_s2) {
    int i = 0, eps_a = 0, eps_b = 1;
    while(eps_s1[i] != '\0'&&eps_s1[i] != '.') {
        eps_a = eps_a*10 + (eps_s1[i]-'0');
        ++ i;
    }

    if(eps_s1[i] == '.') {
        ++ i;
        while(eps_s1[i] != '\0') {
            eps_a = eps_a*10 + (eps_s1[i]-'0');
            eps_b *= 10;
            ++ i;
        }
    }

    if(eps_a > eps_b||eps_b > 100||eps_a <= 0) {
        printf("??? Wrong eps format: %d/%d\n", eps_a, eps_b);
        exit(1);
    }

    eps_a1 = eps_a * eps_a;
    eps_b1 = eps_b * eps_b;
//    cout<<eps_a1<<" "<<eps_b1<<endl;
    
    i = 0;
    eps_a = 0;
    eps_b = 1;
   while(eps_s2[i] != '\0'&&eps_s2[i] != '.') {
       eps_a = eps_a*10 + (eps_s2[i]-'0');
       ++ i;
   }

   if(eps_s2[i] == '.') {
       ++ i;
       while(eps_s2[i] != '\0') {
           eps_a = eps_a*10 + (eps_s2[i]-'0');
           eps_b *= 10;
           ++ i;
       }
   }

   if(eps_a > eps_b||eps_b > 100||eps_a <= 0) {
       printf("??? Wrong eps format: %d/%d\n", eps_a, eps_b);
       exit(1);
   }

   eps_a2 = eps_a * eps_a;
   eps_b2 = eps_b * eps_b;
//    cout<<eps_a2<<" "<<eps_b2<<endl;

}

void Graph::prune_and_cross_link(int eps_a1, int eps_b1,int eps_a2, int eps_b2, int miu, int *cores, int &cores_e){
    
    for(int i = 0;i < n;i ++) { //must be iterating from 0 to n-1
        for(int j = pstart[i];j < pstart[i+1];j ++) {//i点所有邻居
            if(edges[j] < i) {//避免重复检查，邻居节点顺序是排好的
                if(min_cn1[j] == 0) min_cn1[j] = -2;//此时还没计算说明这条边已经不相似
                if(min_cn2[j] == 0) min_cn2[j] = -2;

                continue; //节点已检查 跳过
            }
            
            int v = edges[j];//v是i的邻居
            int a = degree[i], b = degree[v];//此时度是加一的
            if(a > b) swap(a, b);
            
            if((((long long)a)*eps_b2 < ((long long)b)*eps_a2) || (((long long)a)*eps_b1 < ((long long)b)*eps_a1)) {//近似计算，判断是否没有相似的可能
                min_cn1[j] = -2;
                min_cn2[j] = -2;
                -- effective_degree[i];
                -- effective_degree[v];
            }
            else {
//                int c1=9,c2=9;
                int c1 = sqrtl(((long double)a*(long double)b*(long double)eps_a1)/(long double)eps_b1) * (long double)2;//计算要达到相似所需要的公共邻居数
                int c2 = sqrtl(((long double)a*(long double)b*(long double)eps_a2)/(long double)eps_b2) * (long double)6;
                if(c1 <= 4 && c2<=12) {//直接认定相似
                    min_cn1[j] = -1;
                    min_cn2[j] = -1;
                    ++ similar_degree[i];//sd对应增加
                    ++ similar_degree[v];
                    if(similar_degree[i] == miu)
                        cores[cores_e ++] = i;
                    if(similar_degree[v] == miu)
                        cores[cores_e ++] = v;
                }else{
                    min_cn1[j] = c1;
                    min_cn2[j] = c2;
                }
            }
            if(min_cn1[j] != -2) {//若还有机会相似，则将该边的反向边也赋值
                int r_id = binary_search(edges, pstart[v], pstart[v+1], i);
                reverse[j] = r_id;//记录反向边的位置
                reverse[r_id] = j;

                min_cn1[r_id] = min_cn1[j];//最低要求是一样的
                min_cn2[r_id] = min_cn2[j];//最低要求是一样的
            }
        }
    }
}

unsigned int Graph::binary_search(unsigned int *array, int b, int e, int val) {
#ifdef _DEBUG_
    if(e < b) printf("??? WA1 in binary_search\n");
#endif
    -- e;
    if(array[e] < val) return e+1;
    while(b < e) {
        int mid = b + (e-b)/2;
        if(array[mid] >= val) e = mid;
        else b = mid+1;
    }
#ifdef _DEBUG_
    if(array[e] < val) printf("??? WA2 in binary_search\n");
#endif
    return e;
}

int Graph::find_root(int u) {
    int x = u;
    while(pa[x] != x) x = pa[x];
//    cout<<endl;
//    for(int iu = 0;iu<n;iu++) cout<<pa[iu]<<" ";
//    cout<<endl;
//    for(int iu = 0;iu<n;iu++) cout<<rank[iu]<<" ";
//    cout<<endl;

    while(pa[u] != x) {
        int tmp = pa[u];
        pa[u] = x;
        u = tmp;
    }

    return x;
}


int Graph::check_common_neighbor(int u, int v, int c1,int c2) {
    
    int cn1 = 4, cn2 = 12;
    

    if(degree[u] > degree[v]) swap(u,v);

    int du = degree[u]+1, dv = degree[v]+1;
    int i = pstart[u], j = pstart[v];
    while(i < pstart[u+1]&&j < pstart[v+1]&&(cn1 < c1 || cn2 < c2)&&(2*du>=c1&&2*dv>=c1)&&(6*du>=c2&&6*dv>=c2)) {
        if(edges[i] < edges[j]) {
            -- du;
            ++ i;
        }
        else if(edges[i] > edges[j]) {
            -- dv;
            ++ j;
        }
        else {
            
            vector<int> cn11 = get_f_and_c_num(u,v,edges[j]);

            cn1 = cn1 + cn11[0];
            cn2 = cn2 + cn11[1];
            ++ i;
            ++ j;
        }
    }

    if(cn1 >= c1&&cn2 >= c2) return -1;
    return -2;
}



int Graph::similar_check_OP(int u, int idx, int eps_a1, int eps_b1, int eps_a2, int eps_b2) {
    int v = edges[idx];

#ifdef _DEBUG_
    if(min_cn1[idx] == -1||min_cn1[idx] == -2) printf("??? WA in similar_check\n");
#endif

    if(min_cn1[idx] == 0) {//肯定还没计算最低标准
        int du = degree[u], dv = degree[v];
        min_cn1[idx] = min_cn1[reverse[idx]] = ceil(sqrtl(((long double)du*(long double)dv*(long double)eps_a1)/(long double)eps_b1) * (long double)2);
        min_cn2[idx] = min_cn2[reverse[idx]] = ceil(sqrtl(((long double)du*(long double)dv*(long double)eps_a2)/(long double)eps_b2) * (long double)6);
    

#ifdef _DEBUG_
        if(du < c||dv < c) return -2;
#endif

        if(min_cn1[idx] <= 4 && min_cn2[idx]<=12) return -1;
        
    }
    
    return check_common_neighbor(u, v, min_cn1[idx], min_cn2[idx]);
}

vector<int> Graph::get_f_and_c_num(int a,int b, int c){
//    int a_to_b = 1;
//    int b_to_a = 0;
//    int a_to_c = 1;
//    int c_to_a = 0;
//    int c_to_b = 1;
//    int b_to_c = 0;

    int a_to_b = BinarySearch(out_edges,b,out_pstart[a],out_pstart[a+1]);
    int b_to_a = BinarySearch(out_edges,a,out_pstart[b],out_pstart[b+1]);
    int a_to_c = BinarySearch(out_edges,c,out_pstart[a],out_pstart[a+1]);
    int c_to_a = BinarySearch(out_edges,a,out_pstart[c],out_pstart[c+1]);
    int c_to_b = BinarySearch(out_edges,b,out_pstart[c],out_pstart[c+1]);
    int b_to_c = BinarySearch(out_edges,c,out_pstart[b],out_pstart[b+1]);

    int all =a_to_b+b_to_a+a_to_c+c_to_a+c_to_b+b_to_c;
    
    int cNum = 0;
    int fNum = 0;
    
    if(all == 6){
        cNum = cNum + 2;
        fNum = fNum + 6;    }
    else if(all == 5){
        cNum = cNum + 1;
        fNum = fNum + 3;
    }else{
        if(a_to_b && b_to_c && c_to_a) cNum++;
        if(a_to_c && c_to_b && b_to_a) cNum++;
        if(a_to_b && c_to_b && c_to_a) fNum++;
        if(a_to_b && a_to_c && c_to_b) fNum++;
        if(a_to_b && b_to_c && a_to_c) fNum++;
        if(a_to_c && b_to_c && b_to_a) fNum++;
        if(b_to_a && b_to_c && c_to_a) fNum++;
        if(b_to_a && c_to_b && c_to_a) fNum++;
    }
    vector<int> cn;
    cn.push_back(cNum);
    cn.push_back(fNum);
    return cn;
}

unsigned int Graph::BinarySearch(unsigned int *array, int val ,int b, int e)
{
    int low, high, mid;
    low = b;
    high = e-1;
    while(low<=high)
    {
        mid = (low+high)/2;
        if(array[mid]==val)
            return 1;
        if(array[mid]>val)
            high = mid-1;
        if(array[mid]<val)
            low = mid+1;
    }
    return 0;
}

void Graph::my_union(int u, int v) {
    int ru = find_root(u);
    int rv = find_root(v);

    if(ru == rv) return ;

    if(rank[ru] < rank[rv]) pa[ru] = rv;
    else if(rank[ru] > rank[rv]) pa[rv] = ru;
    else {
        pa[rv] = ru;
        ++ rank[ru];
    }
}

void Graph::cluster_noncore_vertices(int eps_a1, int eps_b1,int eps_a2, int eps_b2, int mu){
    
    for(int i = 0;i < n;i ++) {//先遍历一次，核心的相似的邻居归到同一个聚类里
        if(similar_degree[i] >= mu) {
            int x = find_root(i);//寻找根节点
            if(i < cid[x]) cid[x] = i;
        }
    }
    noncore_cluster.clear();
    noncore_cluster.reserve(n);
    for(int i = 0;i < n;i ++)
    if(similar_degree[i] >= mu) {//只有核心才有资格聚类自己的邻居
        for(int j = pstart[i];j < pstart[i+1];j ++)
        if(similar_degree[edges[j]] < mu) {//只需要找非核心的
            if(min_cn1[j] >= 0) {//还没有相似性检查
                min_cn1[j] = min_cn2[j] = similar_check_OP(i, j, eps_a1, eps_b1, eps_a2, eps_b2);//进行相似性检查
//                if(reverse[reverse[j]] != j) printf("WA cluster_noncore\n");
                min_cn1[reverse[j]] = min_cn1[j];
                min_cn2[reverse[j]] = min_cn2[j];
                if(min_cn1[j] == -1) {
                    ++ similar_degree[i];
                    ++ similar_degree[edges[j]];
                }
            }
            if(min_cn1[j] == -1) noncore_cluster.push_back(make_pair(cid[pa[i]], edges[j]));//把游离的点记录下来
        }
    }
}

void Graph::output(){
    printf("\t*** Start write result into disk!\n");
    printf("c/n vertex_id cluster_id\n");

    int mu = miu;
//    for(int i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
//        printf("c %d %d\n", i, cid[pa[i]]);//先输出核心
//    }

    sort(noncore_cluster.begin(), noncore_cluster.end());
    noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
//    for(int i = 0;i < noncore_cluster.size();i ++) {
//        printf("n %d %d\n", noncore_cluster[i].second, noncore_cluster[i].first);
//    }

}

