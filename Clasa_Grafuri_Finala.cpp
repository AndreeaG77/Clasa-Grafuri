#include <iostream>
#include <fstream>
#include <bits/stdc++.h>

using namespace std;

ofstream g("rezultate.out");

#define nmax 505
#define nmax_dfs 10005
#define mmax 200010
#define cmax 110005
const int inf = 0x3f3f3f3f;


class Graf{
private:
    int nr_noduri, nr_muchii;
    vector <int> la[nmax];
    vector <pair<int,int>> lac[nmax];

    //Functii ajutatoare
    void dfs(const int nod, int viz[nmax_dfs], vector<int> la[nmax_dfs]);
    int Tarjan(int nod, int nr, int id, int indx[nmax], int viz[nmax], int dmin[nmax], vector<vector<int>> &conex, stack<int> &s);
    int reprez(int u, int tata[nmax]);
    void reuneste(const int u, const int v, int tata[nmax], int h[nmax]);
    vector<int> bfs_arb(vector<int> la[nmax], const int s);
public:
    Graf();
    // Constructor cu paramentrii pentru problemele cu lista de adiacenta
    Graf(const int nr_noduri, const int nr_muchii, const vector <int> la[nmax]);
    // Constructor cu paramentrii pentru probleme cu lista de adiacenta si costuri pe muchii
    Graf(const int nr_noduri, const int nr_muchii, const vector <pair<int,int>> lac[nmax]);
    ~Graf();
    vector<int> BFS(const int ns);
    int DFS();
    int HHakimi(int n, vector<int> v);
    void CTC();
    vector<int> Sortare_topologica();
    void APM();
    vector<int> Dijkstra();
    void BellmanFord();
    void Disjoint();
    int Flux_Max(int c[nmax][nmax]);
    void FloydWarshall(const int n, int m[nmax][nmax]);
    int Dm_Arb();

};

Graf::Graf(){
        this->nr_noduri=0;
        this->nr_muchii=0;
    }

Graf::Graf(const int nr_noduri, const int nr_muchii, const vector <int> la[nmax]){
        this->nr_noduri = nr_noduri;
        this->nr_muchii = nr_muchii;
        for(int i=1; i<=nr_noduri; i++){
            this->la[i].clear();
        }
        for(int i=1; i<=nr_noduri; i++){
            for(int j=0; j<la[i].size(); j++){
                this->la[i].push_back(la[i][j]);
            }
        }
    }

Graf::Graf(const int nr_noduri, const int nr_muchii, const vector <pair<int,int>> lac[nmax]){
        this->nr_noduri = nr_noduri;
        this->nr_muchii = nr_muchii;
        for(int i=1; i<=nr_noduri; i++){
            this->lac[i].clear();
        }
        for(int i=1; i<=nr_noduri; i++){
            for(int j=0; j<lac[i].size(); j++){
                this->lac[i].push_back(lac[i][j]);
            }
        }
}

Graf::~Graf(){
        for(int i=1; i<=nr_noduri; i++){
            this->la[i].clear();
            this->lac[i].clear();
        }

}

vector<int> Graf::BFS(const int ns){

    int coada[nmax];
    vector<int> cost(nmax+1, -1);

    // Introduc nodul de start in coada
    coada[1] = ns;
    cost[ns] = 0;
    int index_coada = 1;

    // Se parcurg nodurile din coada
    for(int i=1; i<=index_coada; i++){
    // Pentru fiecare nod din coada se verifica vecinii
        for(int j=0; j<la[coada[i]].size(); j++){
    // Daca vecinul nu a fost vizitat inca, se adauga in coada si i se seteaza costul
            if(cost[la[coada[i]][j]] == -1){
                cost[la[coada[i]][j]] = cost[coada[i]] + 1;
                index_coada++;
                coada[index_coada] = la[coada[i]][j];
            }
        }
    }

    //cout<<"\nBFS: ";
    //for(int j=1; j<=nr_noduri; j++)
   // cout<<cost[j]<<" ";
   return cost;

}

void Graf::dfs(const int nod, int viz[nmax_dfs], vector<int> la[nmax_dfs]){
    // Se marcheaza nodul curent ca fiind vizitat
    viz[nod]=1;
    // Se viziteaza toate nodurile la care poate ajunge nodul curent,
    // apelandu-se recursiv functia DFS
    for(int j=0; j<la[nod].size(); j++){
        if(viz[la[nod][j]]==0)
            dfs(la[nod][j],viz,la);
    }
}

int Graf::DFS(){

        int viz[nmax_dfs]={0}, conex=0;

        for(int j=1; j<=nr_noduri; j++){
        // Pentru fiecare nod nevizitat se marcheaza toate nodurile
        // componentei conexe din care face parte apelandu-se functia DFS
            if(viz[j]==0){
                conex++;
                dfs(j,viz,la);
            }
        }

        return conex;

}

int Graf::HHakimi(int n, vector<int> v){
    int i,x,ok=0;

    while(ok==0){
        sort(v.begin(), v.end(), greater<int>());
        x = v.front();
        v.erase(v.begin());
        n--;

        if(x==0) ok=2;

        for(i=0;i<n && x>0;i++){
            if(v[i]!=0) {
                v[i]--;
                x--;
            }
        }

        if(x!=0) ok=1;
    }

    if(ok==2) return 1;
    if(ok==1) return 0;
}

int Graf::Tarjan(int nod, int nr, int id, int indx[nmax], int viz[nmax], int dmin[nmax], vector<vector<int>> &conex, stack<int> &s){

    indx[nod] = id;
    dmin[nod] = id;
    viz[nod] = 1;
    id++;
    s.push(nod);

    // Pentru fiecare vecin nevizitat al nodului curent se apeleaza recursiv Tarjan
    for(int j=0; j<la[nod].size(); j++){
        if(indx[la[nod][j]]==-1){
            nr=Tarjan(la[nod][j],nr,id,indx,viz,dmin,conex,s);
            dmin[nod]=min(dmin[nod], dmin[la[nod][j]]);
        }
    // Daca vecinul a fost deja vizitat, actualizam distanta minima (exista drum de intoarcere in arborele DFS)
        else if(viz[la[nod][j]]==1){
            dmin[nod]=min(dmin[nod], indx[la[nod][j]]);
        }
    }

    if(indx[nod]==dmin[nod]){
        int n;
        nr++;
        vector<int> aux;
        do{
            n=s.top();
            s.pop();
            viz[n]=0;
            aux.push_back(n);
            //conex[nr].push_back(n);

        }while(n!=nod);
        conex.push_back(aux);
    }
    return nr;
}

void Graf::CTC(){
        vector <vector<int>> conex;
        int indx[nmax], viz[nmax]={0}, dmin[nmax], id=0, nr=0;
        stack<int> s;

        // Vectorul indx tine ordinea in care sunt vizitate nodurile intr-o parcurgere DFS
        memset(indx, -1, sizeof(indx));
        // Vectorul dmin tine distanta minima fata de radacina arborelui DFS
        memset(dmin, -1, sizeof(dmin));

        // Pentru fiecare nod nevizitat se aplica algoritmul Tarjan
        for(int j=1; j<=nr_noduri; j++){
            if(indx[j]==-1){
                nr=Tarjan(j,nr,id,indx,viz,dmin,conex,s);
            }
        }

        g<<"\nComponente conexe: ";
        g<<nr<<"\n";
        for(int k=0; k<nr; k++){
            for(int l=0; l<conex[k].size(); l++){
                g<<conex[k][l]<<" ";
            }
            g<<"\n";
        }

    }

vector<int> Graf::Sortare_topologica(){

        vector<int> sortare_top;
        int grd_int[nmax]={0};
        queue<int> nstart;

        for(int i=1; i<=nr_noduri; i++){
            for(int j=0; j<la[i].size(); j++){
                // In vectorul grd_int se introduce gradul interior pentru fiecare nod
                grd_int[la[i][j]]++;
            }
        }

        // In vectorul nstart se pun nodurile cu gradul interior egal cu 0
        for(int j=1; j<=nr_noduri; j++){
            if(grd_int[j]==0){
                nstart.push(j);
            }
        }

        int nod;
        while (!nstart.empty()){
            // Se ia din coada primul nod si se adauga la solutie
            sortare_top.push_back(nstart.front());
            nod = nstart.front();
            nstart.pop();

            // Se elimina toate muchiile care pornesc din nodul curent
            for(int j=0; j<la[nod].size(); j++){
            // Se actualizeaza gradul interior pentru fiecare nod unde s-a sters muchia (nod curent, nod)
                grd_int[la[nod][j]]--;
            // Daca gradul interior actualizat este egal cu 0 atunci se adauga in coada
                if(grd_int[la[nod][j]]==0){
                    nstart.push(la[nod][j]);
                }
            }
            la[nod].clear();
        }

        //g<<"\nSortare topologica: ";
        //for(int k=0; k<nr_noduri; k++){
        //    g<<sortare_top[k]<<" ";
       // }
    return sortare_top;
}

void Graf::APM(){
    //vector <pair<int,int>> la[nmax];
    int viz[nmax]={0};
    vector <pair<int,int>> apm;

    int s=0,nrmuchii=0, idx=0;


    // Implementare priority queue cu struct: https://www.geeksforgeeks.org/stl-priority-queue-for-structure-or-class/

    // Am creat o structura care tine nodurile unei muchii si costul ei
     struct muchie {
        int x,y,c;
        muchie(int x, int y, int c): x(x), y(y), c(c) {}
    };

    // Am defenit o functie care compara costurile a doua muchii
    struct CompareC {
    bool operator()(muchie const& m1, muchie const& m2)
    {
        return m1.c > m2.c;
    }
    };

    priority_queue<muchie, vector<muchie>, CompareC> pq;
    int v=1;
    viz[1]=1;
    // Introduc in pq muchiile care pornesc din nodul 1
    for(int j=0; j<lac[v].size(); j++){
        pq.push(muchie(v,lac[v][j].first, lac[v][j].second));
    }

    // Cat timp coada nu este goala verific daca mai pot adauga muchii la apm
    while(!pq.empty()){
        // Se extrage muchia cu costul minim
        muchie top = pq.top();
        pq.pop();
        // Daca nodul in care ajunge muchia nu a fost inca vizitat, se adauga la apm
        if(viz[top.y]==0){
            viz[top.y]=1;
            s=s+ top.c;
            nrmuchii++;
            // Se adauga muchia la solutie
            apm.push_back(make_pair(top.x,top.y));

            // Se adauga in pq muchiile nodului curent
            for(int k=0; k<lac[top.y].size(); k++){
                pq.push(muchie(top.y,lac[top.y][k].first, lac[top.y][k].second));
            }
        }
    }

    g<<"\nAPM:\n";
    g<<s<<"\n"<<nrmuchii<<"\n";

    for(int l=0; l<apm.size(); l++){
        g<<apm[l].first<<" "<<apm[l].second<<"\n";
    }
}

vector<int> Graf::Dijkstra(){
    vector<int> d(nmax+1, inf);
    set <pair<int,int>> s;
    int idx=0;

    d[1]=0;

    // Se adauga nodul de start in set
    s.insert(make_pair(d[1],1));

    while(!s.empty()){
        int nod = s.begin()->second;
        int dmin = s.begin()->first;
        s.erase(s.begin());
        // Pentru fiecare vecin al nodului curent se actualizeaza distanta minima
        for(int k=0; k<lac[nod].size(); k++){
            int nod2 = lac[nod][k].first;
            int cost = lac[nod][k].second;
            if(dmin+cost< d[nod2]){
                if(d[nod2] != inf){
                    s.erase(s.find(make_pair(d[nod2], nod2)));
                }
                d[nod2] = dmin+cost;
                s.insert(make_pair(d[nod2], nod2));
            }
        }
    }


    return d;
}

void Graf::BellmanFord(){

    int d[nmax], frecv[nmax]={0};
    queue <int> q;
    bool inq[nmax];
    int cn=0;
    int idx=0;

    // Se seteaza vectorul de distante minime cu inf
    memset(d, inf, sizeof(d));
    d[1]=0;

    // Se adauga nodul de start in coada
    q.push(1);

    while(!q.empty() && cn==0){
        int nod = q.front();
        q.pop();
        inq[nod]=false;
        // Pentru fiecare vecin al nodului curent se actualizeaza distanta minima
        for(int k=0; k<lac[nod].size(); k++){
            int nod2 = lac[nod][k].first;
            int cost = lac[nod][k].second;
            if(d[nod]+cost< d[nod2]){
                d[nod2] = d[nod]+cost;
                if(!inq[nod2]){
                    if(frecv[nod2]>nr_noduri) cn=1;
                    else {
                        inq[nod2]=true;
                        q.push(nod2);
                        frecv[nod2]++;
                    }
                }
            }
        }
    }

    // Afisare
    g<<"\nBellman-Ford: ";
    if(cn==1) g<<"Ciclu negativ!";
    else
        for(int l=2; l<=nr_noduri; l++){
            if(d[l]==inf) g<<0<<" ";
            else g<<d[l]<<" ";
        }
}

int Graf::reprez(int u, int tata[nmax]) {
    while (tata[u] != 0)
        u = tata[u];
    return u;
}

void Graf::reuneste(const int u, const int v, int tata[nmax], int h[nmax]) {
    // Gasesc radacinile arborilor din care fac parte cele doua valori
    int ru=reprez(u,tata), rv=reprez(v,tata);
    // Se leaga printr-o muchie cele doua radacini
    if (h[ru] > h[rv])
        tata[rv] = ru;
    else {
        tata[ru] = rv;
        // In caz de egalitate se incrementeaza inaltimea unui arbore cu 1
        if (h[ru] == h[rv])
            h[rv]++;
    }
}

void Graf::Disjoint(){
    ifstream f("disjoint.in");
    //Memorez multimile ca arbori
    // Folosesc un vector de tati si un vector pentru inaltimile arborilor
    int tata[nmax]={0}, h[nmax]={0};

    int n,m,op,x,y;
    f>>n>>m;

    g<<"\nDisjoint: ";
    for(int j=0; j<m; j++){
        f>>op>>x>>y;
        if(op==1){
            reuneste(x,y,tata,h);
        }
        else {
            if(reprez(x,tata)==reprez(y,tata)) g<<"DA\n";
            else g<<"NU\n";
        }
    }
}

int Graf::Flux_Max(int c[nmax][nmax]){

    int i,j=0,flow=0;
    int n = nr_noduri;

    int nc;
    int ok=0;

    while(ok==0){
        int p[n+1]={0};
        queue <int> q;
        q.push(1);

        while(!q.empty()){
            nc = q.front();
            q.pop();
            for(i=0; i<la[nc].size(); i++){
                int naux = la[nc][i];
                if(p[naux]==0 && c[nc][naux]>0){
                    p[naux]=nc;
                    q.push(naux);
                }
            }

        }

        if(p[n]!=0){
            int fm=cmax;
            int t;
            nc = n;
            while(nc!=1){
                t=p[nc];
                if(c[t][nc]<fm) fm = c[t][nc];
                nc=t;
            }
            nc = n;
            while(nc!=1){
                t=p[nc];
                c[t][nc] -= fm;
                nc=t;
            }
            flow += fm;
        }
    if(p[n]==0) ok=1;
    }

    return flow;
}

void Graf::FloydWarshall(const int n, int m[nmax][nmax]){
    int i,j,k;

    for(k=1; k<=n; k++){
        for(i=1;i<=n;i++){
            for(j=1;j<=n;j++){
               if(m[i][j]>m[i][k]+m[k][j])
                    m[i][j] = m[i][k]+m[k][j];
            }
        }
    }

    g<<"\nFloyd Warshall:\n";
    for(i=1;i<=n;i++){
        for(j=1;j<=n;j++){
            g<<m[i][j]<<" ";
        }
        g<<"\n";
    }
}

vector<int> Graf::bfs_arb(vector<int> la[nmax], const int s){
    queue <int> q;
    int v[nmax]={0},last,nc,dist[nmax]={0},dm=0;

    v[s]=1;
    q.push(s);
    dist[s]=1;

    while(!q.empty()){

        nc = q.front();
        last = q.front();
        q.pop();
        v[nc] = 1;

        for(int i=0; i<la[nc].size(); i++){
            if(v[la[nc][i]]==0){
                v[la[nc][i]]=1;
                q.push(la[nc][i]);
                dist[la[nc][i]] = dist[nc] + 1;
                dm = dist[la[nc][i]];
            }
        }
    }
    vector<int> sol;
    sol.push_back(last);
    sol.push_back(dm);
    return sol;
}

int Graf::Dm_Arb(){

    vector<int> sol_aux = bfs_arb(la,1);
    vector<int> sol = bfs_arb(la,sol_aux[0]);
    return sol[1];
}


void init_BFS(){

    ifstream f("bfs.in");

    vector <int> la[nmax];
    int n,m,s,x,y;
    f>>n>>m>>s;
    for(int i=0; i<m; i++){
        f>>x>>y;
        la[x].push_back(y);
    }

    Graf gr(n,m,la);
    vector<int> c=gr.BFS(s);

    g<<"\nBFS: ";
    for(int j=1; j<=n; j++)
    g<<c[j]<<" ";

}

void init_DFS(){
    ifstream f("dfs.in");

    vector <int> la[nmax_dfs];
    int n,m,x,y;
    f>>n>>m;
    for(int i=0; i<m; i++){
        f>>x>>y;
        la[x].push_back(y);
        la[y].push_back(x);
    }

    Graf gr(n,m,la);
    int conex = gr.DFS();

    g<<"\nDFS";
    g<<"\nNumarul de elemente conexe: "<<conex;
}

void init_HHakimi(){

    ifstream f("hhakimi.in");

    vector <int> v;
    int n,x;
    f>>n;
    for(int i=0; i<n; i++){
        f>>x;
        v.push_back(x);
    }

    Graf gr;
    int r = gr.HHakimi(n,v);

    g<<"\nHavel Hakimi: ";
    if(r) g<<"Da";
    else g<<"Nu";
}

void init_CTC(){
    ifstream f("ctc.in");

    vector <int> la[nmax];
    int n,m,x,y;
    f>>n>>m;
    for(int i=0; i<m; i++){
        f>>x>>y;
        la[x].push_back(y);
    }

    Graf gr(n,m,la);
    gr.CTC();
}

void init_Sortare_Top(){

    ifstream f("sortaret.in");

    vector <int> la[nmax];
    int n,m,x,y;
    f>>n>>m;
    for(int i=0; i<m; i++){
        f>>x>>y;
        la[x].push_back(y);
    }

    Graf gr(n,m,la);
    vector<int> st=gr.Sortare_topologica();

    g<<"\nSortare topologica: ";
    for(int k=0; k<n; k++){
        g<<st[k]<<" ";
    }
}

void init_APM(){
    ifstream f("apm.in");

    vector <pair<int,int>> lac[nmax];
    int n,m,c,x,y;
    f>>n>>m;
    for(int i=0; i<m; i++){
        f>>x>>y>>c;
        lac[x].push_back(make_pair(y,c));
        lac[y].push_back(make_pair(x,c));
    }

    Graf gr(n,m,lac);
    gr.APM();

}

void init_Dijkstra(){
    ifstream f("dijkstra.in");

    vector <pair<int,int>> lac[nmax];
    int n,m,c,x,y;
    f>>n>>m;
    for(int i=0; i<m; i++){
        f>>x>>y>>c;
        lac[x].push_back(make_pair(y,c));
    }

    Graf gr(n,m,lac);
    vector<int> d=gr.Dijkstra();

    g<<"\nDijkstra: ";
    for(int i=2; i<=n; i++){
       if(d[i]==inf) g<<0<<" ";
       else g<<d[i]<<" "; }
}

void init_BellmanF(){
    ifstream f("bellmanford.in");

    vector <pair<int,int>> lac[nmax];
    int n,m,c,x,y;
    f>>n>>m;
    for(int i=0; i<m; i++){
        f>>x>>y>>c;
        lac[x].push_back(make_pair(y,c));
    }

    Graf gr(n,m,lac);
    gr.BellmanFord();

}

void init_Flux_Max(){
    ifstream f("maxflow.in");

    vector <int> la[nmax];
    int c[nmax][nmax];
    int n,m,cst,x,y;
    f>>n>>m;
    for(int i=0; i<m; i++){
        f>>x>>y>>cst;
        la[x].push_back(y);
        c[x][y]=cst;
    }

    Graf gr(n,m,la);
    int fl=gr.Flux_Max(c);

    g<<"\nFlux maxim: ";
    g<<fl;
}

void init_Floyd_Warshall(){
    ifstream f("royfloyd.in");

    int n,i,j,x,k;
    int m[nmax][nmax];

    f>>n;
    for(i=1;i<=n;i++){
        for(j=1;j<=n;j++){
            f>>x;
            if(x==0 && i!=j)
                m[i][j]=inf;
            else
                m[i][j]=x;
        }
    }

    Graf gr;
    gr.FloydWarshall(n,m);

}

void init_Dm_Arb(){
    ifstream f("darb.in");

    vector <int> la[nmax];
    int n,i,x,y;

    f>>n;
    for(i=1; i<n; i++){
        f>>x>>y;
        la[x].push_back(y);
        la[y].push_back(x);
    }

    Graf gr(n,n-1,la);
    int dm=gr.Dm_Arb();

    g<<"\nDiamentrul unui arbore: ";
    g<<dm;
}

int main()
{

    //init_BFS();
    //init_DFS();
    //init_HHakimi();
    //init_CTC();
    //init_Sortare_Top();
    //init_APM();
    //init_Dijkstra();
    //init_BellmanF();
    //Graf g;
    //g.Disjoint();
    //init_Flux_Max();
    //init_Floyd_Warshall();
    //init_Dm_Arb();


    return 0;
}
