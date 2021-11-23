#include <iostream>
#include <fstream>
#include <bits/stdc++.h>

using namespace std;

#define nmax 100010
#define mmax 200010
const int inf = 0x3f3f3f3f;

class Graf{
private:
    int nr_noduri, nr_muchii;
    int* m;
    int* cst;
public:
    Graf();
    Graf(int nr_noduri, int nr_muchii, int* m, int* cst);
    ~Graf();
    void BFS();
    void DFS();
    void dfs(int nod, int viz[nmax], vector<int> la[nmax]);
    void CTC();
    int Tarjan(int nod, int nr, int id, int indx[nmax], int viz[nmax], int dmin[nmax], vector<int> conex[nmax], vector<int> la[nmax], stack<int> s);
    void Sortare_topologica();
    void APM();
    void Dijkstra();
    void Bellmanford();
    void Disjoint();
    int reprez(int u, int tata[nmax]);
    void reuneste(int u, int v, int tata[nmax], int h[nmax]);

    // Supraincarcarea operatorului >>
    friend istream& operator>>(istream& in, Graf& g){
        cout<<"\nNumar noduri: "; in>>g.nr_noduri;
        cout<<"\nNumar muchii: "; in>>g.nr_muchii;
        cout<<"\nMuchii: ";
        if(g.m!=NULL)
            delete[] g.m;
        if(g.cst!=NULL)
            delete[] g.cst;
        if (g.nr_muchii==0){
            g.m=NULL;
            g.cst=NULL;
        }
        else {
            int j=0;
            g.m=new int[2*g.nr_muchii];
            g.cst=new int[g.nr_muchii];
            for(int i=0; i<2*g.nr_muchii; i=i+2){
                in>>g.m[i]>>g.m[i+1]>>g.cst[j++];
            }
        }
        return in;
    }

    // Supraincarcarea operatorului <<
    friend ostream& operator<<(ostream& out, const Graf& g){
        out<<"\nNumar noduri: "<<g.nr_noduri;
        out<<"\nNumar muchii: "<<g.nr_muchii;
        out<<"\nMuchii:\n";
        int j=0;
        for(int i=0; i<2*g.nr_muchii; i=i+2){
            out<<g.m[i]<<" "<<g.m[i+1]<<" "<<g.cst[j++]<<"\n";
        }
        return out;
    }

};

Graf::Graf(){
        this->nr_noduri=0;
        this->nr_muchii=0;
        this->m=NULL;
        this->cst=NULL;
    }

Graf::Graf(int nr_noduri, int nr_muchii, int* m, int* cst){
        this->nr_noduri = nr_noduri;
        this->nr_muchii = nr_muchii;
        if(this->m!=NULL)
            delete[] this->m;
        this->m=new int[2*nr_muchii];
        for(int i=0; i<2*nr_muchii; i++){
            this->m[i]=m[i];
        }
        if(this->cst!=NULL)
            delete[] this->cst;
        this->cst=new int[nr_muchii];
        for(int i=0; i<nr_muchii; i++){
            this->cst[i]=cst[i];
        }
    }

Graf::~Graf(){
        if(this->m!=NULL)
            delete[] this->m;
        if(this->cst!=NULL)
            delete[] this->cst;
}

void Graf::BFS(){
    int nod_start = 1;
    vector <int> la[nmax];
    // Se introduc muchiile in liste de adiacenta
    for(int i=0; i<2*nr_muchii; i=i+2){
        la[m[i]].push_back(m[i+1]);
    }
    int coada[nmax], cost[nmax];
    memset(cost, -1, sizeof(cost));

    // Introduc nodul de start in coada
    coada[1] = nod_start;
    cost[nod_start] = 0;
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
    cout<<"\nBFS: ";
    for(int j=1; j<=nr_noduri; j++)
    cout<<cost[j]<<" ";
}

void Graf::dfs(int nod, int viz[nmax], vector<int> la[nmax]){
    // Se marcheaza nodul curent ca fiind vizitat
    viz[nod]=1;
    // Se viziteaza toate nodurile la care poate ajunge nodul curent,
    // apelandu-se recursiv functia DFS
    for(int j=0; j<la[nod].size(); j++){
        if(viz[la[nod][j]]==0)
            dfs(la[nod][j],viz,la);
    }
}

void Graf::DFS(){
        vector <int> la[nmax];
        // Se introduc muchiile in liste de adiacenta
        for(int i=0; i<2*nr_muchii; i=i+2){
            la[m[i]].push_back(m[i+1]);
            la[m[i+1]].push_back(m[i]);
        }
        int viz[nmax]={0}, conex=0;

        for(int j=1; j<=nr_noduri; j++){
        // Pentru fiecare nod nevizitat se marcheaza toate nodurile
        // componentei conexe din care face parte apelandu-se functia DFS
            if(viz[j]==0){
                conex++;
                dfs(j,viz,la);
            }
        }

        cout<<"Numar elemente conexe: "<<conex;

}

int Graf::Tarjan(int nod, int nr, int id, int indx[nmax], int viz[nmax], int dmin[nmax], vector<int> conex[nmax], vector<int> la[nmax], stack<int> s){

    indx[nod] = id;
    dmin[nod] = id;
    viz[nod] = 1;
    id++;
    s.push(nod);

    // Pentru fiecare vecin nevizitat al nodului curent se apeleaza recursiv Tarjan
    for(int j=0; j<la[nod].size(); j++){
        if(indx[la[nod][j]]==-1){
            nr=Tarjan(la[nod][j],nr,id,indx,viz,dmin,conex,la,s);
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
        do{
            n=s.top();
            s.pop();
            viz[n]=0;
            conex[nr].push_back(n);

        }while(n!=nod);
    }
    return nr;
}

void Graf::CTC(){
        vector <int> la[nmax], conex[nmax];
        int indx[nmax], viz[nmax]={0}, dmin[nmax], id=0, nr=0;
        stack<int> s;

        // Se introduc muchiile in liste de adiacenta
        for(int i=0; i<2*nr_muchii; i=i+2){
            la[m[i]].push_back(m[i+1]);
        }

        // Vectorul indx tine ordinea in care sunt vizitate nodurile intr-o parcurgere DFS
        memset(indx, -1, sizeof(indx));
        // Vectorul dmin tine distanta minima fata de radacina arborelui DFS
        memset(dmin, -1, sizeof(dmin));

        // Pentru fiecare nod nevizitat se aplica algoritmul Tarjan
        for(int j=1; j<=nr_noduri; j++){
            if(indx[j]==-1){
                nr=Tarjan(j,nr,id,indx,viz,dmin,conex,la,s);
            }
        }

        cout<<"\nComponente conexe: ";
        cout<<nr<<"\n";
        for(int k=1; k<=nr; k++){
            for(int l=0; l<conex[k].size(); l++){
                cout<<conex[k][l]<<" ";
            }
            cout<<"\n";
        }

    }

void Graf::Sortare_topologica(){
        vector <int> la[nmax];
        int sortare_top[nmax], grd_int[nmax]={0};
        queue<int> nstart;

        // Se introduc muchiile in liste de adiacenta
        for(int i=0; i<2*nr_muchii; i=i+2){
            la[m[i]].push_back(m[i+1]);
            // In vectorul grd_int se introduce gradul interior pentru fiecare nod
            grd_int[m[i+1]]++;
        }

        // In vectorul nstart se pun nodurile cu gradul interior egal cu 0
        for(int j=1; j<=nr_noduri; j++){
            if(grd_int[j]==0){
                nstart.push(j);
            }
        }

        int l = 0, nod;
        while (!nstart.empty()){
            // Se ia din coada primul nod si se adauga la solutie
            sortare_top[l] = nstart.front();
            nod = nstart.front();
            nstart.pop();
            l++;
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

        cout<<"\nSortare topologica: ";
        for(int k=0; k<nr_noduri; k++){
            cout<<sortare_top[k]<<" ";
        }

}

void Graf::APM(){
    vector <pair<int,int>> la[nmax];
    int viz[nmax]={0};
    vector <pair<int,int>> apm;

    int s=0,nrmuchii=0, idx=0;

    // Se introduc perechile (nod, cost) in liste de adiacenta
    for(int i=0; i<2*nr_muchii; i=i+2){
        la[m[i]].push_back(make_pair(m[i+1],cst[idx]));
        la[m[i+1]].push_back(make_pair(m[i],cst[idx++]));
    }

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
    for(int j=0; j<la[v].size(); j++){
        pq.push(muchie(v,la[v][j].first, la[v][j].second));
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
            for(int k=0; k<la[top.y].size(); k++){
                pq.push(muchie(top.y,la[top.y][k].first, la[top.y][k].second));
            }
        }
    }

    cout<<"\nAPM:\n";
    cout<<s<<"\n"<<nrmuchii<<"\n";

    for(int l=0; l<apm.size(); l++){
        cout<<apm[l].first<<" "<<apm[l].second<<"\n";
    }
}

void Graf::Dijkstra(){
    vector <pair<int,int>> la[nmax];
    int d[nmax];
    set <pair<int,int>> s;
    int idx=0;

    // Se introduc perechile (nod, cost) in liste de adiacenta
    for(int i=0; i<2*nr_muchii; i=i+2){
        la[m[i]].push_back(make_pair(m[i+1],cst[idx++]));
    }

    // Se seteaza vectorul de distante minime cu inf
    memset(d, inf, sizeof(d));
    d[1]=0;

    // Se adauga nodul de start in set
    s.insert(make_pair(d[1],1));

    while(!s.empty()){
        int nod = s.begin()->second;
        int dmin = s.begin()->first;
        s.erase(s.begin());
        // Pentru fiecare vecin al nodului curent se actualizeaza distanta minima
        for(int k=0; k<la[nod].size(); k++){
            int nod2 = la[nod][k].first;
            int cost = la[nod][k].second;
            if(dmin+cost< d[nod2]){
                if(d[nod2] != inf){
                    s.erase(s.find(make_pair(d[nod2], nod2)));
                }
                d[nod2] = dmin+cost;
                s.insert(make_pair(d[nod2], nod2));
            }
        }
    }

    // Afisare
    for(int l=2; l<=nr_noduri; l++){
        if(d[l]==inf) cout<<0<<" ";
        else cout<<d[l]<<" ";
    }
}

void Graf::Bellmanford(){
    vector <pair<int,int>> la[nmax];
    int d[nmax], frecv[nmax]={0};
    queue <int> q;
    bool inq[nmax];
    int cn=0;
    int idx=0;

    // Se introduc perechile (nod, cost) in liste de adiacenta
    for(int i=0; i<2*nr_muchii; i=i+2){
        la[m[i]].push_back(make_pair(m[i+1],cst[idx++]));
    }

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
        for(int k=0; k<la[nod].size(); k++){
            int nod2 = la[nod][k].first;
            int cost = la[nod][k].second;
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
    if(cn==1) cout<<"Ciclu negativ!";
    else
        for(int l=2; l<=nr_noduri; l++){
            if(d[l]==inf) cout<<0<<" ";
            else cout<<d[l]<<" ";
        }
}

int Graf::reprez(int u, int tata[nmax]) {
    while (tata[u] != 0)
        u = tata[u];
    return u;
}

void Graf::reuneste(int u, int v, int tata[nmax], int h[nmax]) {
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
    //Memorez multimile ca arbori
    // Folosesc un vector de tati si un vector pentru inaltimile arborilor
    int tata[nmax]={0}, h[nmax]={0};

    int n,m,op,x,y;
    cin>>n>>m;

    for(int j=0; j<m; j++){
        cin>>op>>x>>y;
        if(op==1){
            reuneste(x,y,tata,h);
        }
        else {
            if(reprez(x,tata)==reprez(y,tata)) cout<<"DA\n";
            else cout<<"NU\n";
        }
    }
}

int main()
{
   /* Graf g;
    cin>>g;
    cout<<g;
    g.Disjoint(); */
    return 0;
}
