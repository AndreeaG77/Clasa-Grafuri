#include <iostream>
#include <fstream>
#include <bits/stdc++.h>

using namespace std;

#define nmax 100010
#define mmax 200010

void dfs(int nod, int viz[nmax], vector<int> la[nmax]){

    // Se marcheaza nodul curent ca fiind vizitat
    viz[nod]=1;
    // Se viziteaza toate nodurile la care poate ajunge nodul curent,
    // apelandu-se recursiv functia DFS
    for(int j=0; j<la[nod].size(); j++){
        if(viz[la[nod][j]]==0)
            dfs(la[nod][j],viz,la);
    }
}


int Tarjan(int nod, int nr, int id, int indx[nmax], int viz[nmax], int dmin[nmax], vector<int> conex[nmax], vector<int> la[nmax], stack<int> s){

    indx[nod] = id;
    dmin[nod] = id;
    viz[nod] = 1;
    id++;
    s.push(nod);

    // Pentru fiecare vecin nevizitat al nodului curent se apeleaza recursiv Tarjan
    for(int j=0; j<la[nod].size(); j++){
        if(indx[la[nod][j]]==-1){
            Tarjan(la[nod][j],nr,id,indx,viz,dmin,conex,la,s);
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




class Graf{
private:
    int nr_noduri, nr_muchii;
    int* m;
public:
    // Constructorul fara parametrii
    Graf(){
        this->nr_noduri=0;
        this->nr_muchii=0;
        this->m=NULL;
    }
    // Constructorul cu toti parametrii
    Graf(int nr_noduri, int nr_muchii, int* m){
        this->nr_noduri = nr_noduri;
        this->nr_muchii = nr_muchii;
        if(this->m!=NULL)
            delete[] this->m;
        this->m=new int[2*nr_muchii];
        for(int i=0; i<2*nr_muchii; i++){
            this->m[i]=m[i];
        }
    }
    // Supraincarcarea operatorului >>
    friend istream& operator>>(istream& in, Graf& g){
        cout<<"\nNumar noduri: "; in>>g.nr_noduri;
        cout<<"\nNumar muchii: "; in>>g.nr_muchii;
        cout<<"\nMuchii: ";
        if(g.m!=NULL)
            delete[] g.m;
        if (g.nr_muchii==0)
            g.m=NULL;
        else {
            g.m=new int[2*g.nr_muchii];
            for(int i=0; i<2*g.nr_muchii; i++){
                in>>g.m[i];
            }
        }
        return in;
    }
    // Supraincarcarea operatorului <<
    friend ostream& operator<<(ostream& out, const Graf& g){
        out<<"\nNumar noduri: "<<g.nr_noduri;
        out<<"\nNumar muchii: "<<g.nr_muchii;
        out<<"\nMuchii:\n";
        for(int i=0; i<2*g.nr_muchii; i=i+2){
            out<<g.m[i]<<" "<<g.m[i+1]<<"\n";
        }
        return out;
    }
    // Destructorul
    ~Graf(){
        if(this->m!=NULL)
            delete[] this->m;
    }

    void BFS(){
        cout<<"ok1";
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

    void DFS(){
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

    void CTC(){
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

    void Sortare_topologica(){
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


};



int main()
{
  /*  Graf g;
    cin>>g;
    cout<<g; */
    return 0;
}
