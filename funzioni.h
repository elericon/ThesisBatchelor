#ifndef funzioni_h
#define funzioni_h

using namespace std;
using namespace arma;

extern int N;//number of spins
extern double magn;//value of magnation
extern int dimensione;//cluster dimension for wolff
extern int pv;//number of first nearest neighbours


double pow1(double base, int esp) {//same as pow but faster
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
double mod(rowvec &r){//calculates the module of the relative distance 
    double mod = sqrt(pow1(r(0),2)+pow1(r(1),2));
    return mod;
}
double mod(cube &r, double riga, double colonna){
//calculates the module of the relative distance
    double mod = sqrt(pow1(r(riga,colonna,0),2)+pow1(r(riga,colonna,1),2));
    return mod;
}
void plot_reticolo() {//executes the plot of the lattice
    string comando;
    comando = "gnuplot";
    comando += " plot_reticolo.plt";
    cout << comando << endl;
    system(comando.c_str());
    return;
}
void blocking_plot(){//executes the plot of the blocking
    //ora faccio il plot
    string comando;

    comando = "gnuplot";
    comando += " plot_blocking.plt";
    cout << comando << endl;
    system(comando.c_str());
}
void plot_osservabili_T(int c_Tm){
//executes the plot of the observables as functions of T
    if(c_Tm!=-1){
        string comando;
        comando = "gnuplot";
        comando += " plot_osservabili.plt";
        cout << comando << endl;
        system(comando.c_str());
        return;
    }
    else{
        string comando;
        comando = "gnuplot";
        comando += " plot_osservabili_T.plt";
        cout << comando << endl;
        system(comando.c_str());
        return;
    }
}
void stampa_reticolo(mat &r){//prints the lattice on a file
    ofstream reticolo;
    reticolo.open("out/reticolo.txt");
    for (int i = 0; i < N; ++i){
        reticolo << r(i, 0) << "\t" << r(i,1) << "\t" << r(i,2) << endl;
        //0 posizione x, 1 posizione y, 2 spin
    }
    reticolo.close();
    // plot_reticolo();
    return;
}
void crea_reticolo_quadrato(mat &r, double L){// creates the square lattice    
    int n = sqrt(N);
    double L_cella = L / n;
    int cont = 0;
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            r.row(cont) = {(double)i * L_cella, (double)j * L_cella, 1};
            cont++;
        }
    }
    stampa_reticolo(r);
    return;
}
double aggiunta_y(int cont_y, double H){
//needed for the creation of the exagonal lattice
    if(cont_y==0){
        return -H;
    }
    else if (cont_y==1 || cont_y==3){
        return 0;
    }
    else{
        return H;
    }
}
double aggiunta_x(int cont_x){//needed for the creation of the exagonal lattice
    if(cont_x==0){
        return 0.5;
    }
    else{
        return 1;
    }
}
void crea_reticolo_esagonale(mat &r, int L){// creates the exagonal lattice
    double H = sqrt(3.)/2.;//"height" of half a hexagon
    int cont_x = 0;
    int cont_y = 0;
    int cont_part = 1;

    double y = H, x = 0;
    do{
        r.row(cont_part-1) = {x, y, 1};
        x += aggiunta_x(cont_x);
        y += aggiunta_y(cont_y, H);
        if(cont_y==3){
            cont_y=0;
        }
        else{
            cont_y++;
        }
        if(cont_x==1){
            cont_x--;
        }
        else{ 
            cont_x++;
        }
        if(cont_part%(L)==0){//at the start of the chain
            x = 0;
            y += 2 * H;
            cont_x = 0;
            cont_y = 0;
        }
        cont_part++;
    }
    while(cont_part!=L*L+1);
    stampa_reticolo(r);
    return;
}
double potenziale_ising(mat &r, double J, int i, mat &primivicini){
//calculates the actual potential 
    double V = 0;

    rowvec p_v = primivicini.row(i);

    for (int j = 0; j < pv; ++j){
        V += r((int)p_v(j),2) * r(i,2) * J;
    }
    return -V;
}
double modifica_potenziale(mat &r, double J, int i, mat &primivicini){
//calculates the modified potential 
    double V = 0;

    rowvec p_v = primivicini.row(i);

    for (int j = 0; j < pv; ++j){
        V += r((int)p_v(j),2);
    }
    V = J * 2 * r(i,2) * V;
    return V;
}
void crea_primi_vicini(double L, string p_v_path){
//creates a file with the first nearest neighbours for each spin
    double L_y = sqrt(3)*(L), L_x = L*3/4;
    
    mat prim_vic(N,pv), r(N,3);
    rowvec dr(2);

    ifstream reticolo;
    reticolo.open("out/reticolo.txt");

    for (int i = 0; i < N; ++i){
        reticolo >> r(i,0) >> r(i,1) >> r(i,2);
    }
    reticolo.close();

    ofstream primivicini;
    primivicini.open(p_v_path);

    for (int i = 0; i < N; ++i){
    //creates the matrix for the first nearest neighbours
        int cont = 0;
        
        for (int j = 0; j < N; ++j){
            if(i!=j){
                for (int k = 0; k < 2; ++k){
                    dr(k) = r(j,k) - r(i,k);
                    if(pv==4){
                        dr(k) -= L * rint(dr(k)/L);//moves in [-L/2,+L/2]
                    }
                }
                if(pv==3){
                    dr(0) -= L_x * rint(dr(0)/L_x);//moves in [-L/2,+L/2]
                    dr(1) -= L_y * rint(dr(1)/L_y);//moves in [-L/2,+L/2]
                }
                
                if(cont==pv){
                    break;
                }
                
                if(mod(dr) <= 1.1){
                //if they are close as much as 1 (1.1 to be sure even with
                //approximation errors) they are first nearest neighbours
                    prim_vic(i,cont) = j;
                    cont++;
                }
            }
        }
        for (int j = 0; j < pv; ++j){
            if(j!=pv-1){
                primivicini << prim_vic(i,j) << "\t";
            }
            else{
                primivicini << prim_vic(i,j) << endl;
            }
        }
    }
    primivicini.close();
}
void blocking(int N_t, int num_oss){//blocking method
    mat O(num_oss, N_t);

    ifstream dati_blocking;
    dati_blocking.open("out/dati_blocking.txt");
    for (int i = 0; i < N_t; ++i){
        for (int j = 0; j < num_oss; ++j){
            dati_blocking >> O(j, i);
        }
    }
    dati_blocking.close();

    ofstream blocking;
    blocking.open("out/blocking.txt");

    int N_B_prev=0;
    for (int B = 10; B < N_t/4; B+=10){

        int N_B = floor(N_t / B);
        
        if(N_B!=N_B_prev){
            
            mat O_mB(num_oss, N_B, fill::value(0));
            //mean observable of the block
            rowvec var_O(num_oss, fill::value(0));
            //variance of each observable
            rowvec O_mean(num_oss,fill::value(0));
            //complessive mean for each observable
            
            for (int i = 0; i < N_B; ++i){//cycle over the blocks
                for (int j = 0; j < B; ++j){
                    for (int k = 0; k < num_oss; ++k){
                        O_mB(k, i) += O(k, (i * B + j)) / B;
                    }
                }
                for (int k = 0; k < num_oss; ++k){
                    O_mean(k) += O_mB(k,i) / N_B;
                    //calculates the mean over all the blocks
                }
            }
            for (int k = 0; k < num_oss; ++k){
                for (int i = 0; i < N_B; ++i){
                    var_O(k) += pow1(O_mB(k,i) - O_mean(k),2) / N_B;
                }
            }
            blocking << B << "\t";
            for (int i = 0; i < num_oss; ++i){
                if(i==num_oss-1){
                    blocking << sqrt(var_O(i)/N_B) << endl;//error of a mean
                }
                else{
                    blocking << sqrt(var_O(i)/N_B) << "\t";
                }
            }
            N_B_prev=N_B;
        }
    }
    blocking.close();
    // blocking_plot();
}
void jackknife(int N_t, int num_oss){//jackknife method
    mat O(num_oss, N_t);

    ifstream dati_blocking;
    ofstream jackknife;
    
    jackknife.open("out/jackknife.txt");
    dati_blocking.open("out/dati_blocking.txt");
    
    for (int i = 0; i < N_t; ++i){
        for (int j = 0; j < num_oss; ++j){
            dati_blocking >> O(j, i);
        }
    }
    dati_blocking.close();

    int N_B_prev=1;//scarto il blocco singolo

    for (int B = (int)(N_t / 1e4) + 1; B < N_t/2; B+=10){
    //takes away blocks which are too small 
        int N_B = floor(N_t / B);
        
        if(N_B!=N_B_prev){
            mat O_m(num_oss, N_B, fill::value(0));
            //mean of the observables of each block
            mat O_m_jack(num_oss, N_B, fill::value(0));
            //means found by excluding a block

            rowvec var_O(num_oss, fill::value(0));
            //variance of each observable
            rowvec O_media(num_oss,fill::value(0));
            //complessive mean for each observable

            //calculate the means of each block
            for (int i = 0; i < N_B; ++i){//cycle over the blocks
                for (int j = 0; j < B; ++j){
                //cycle over the points of the block
                    for (int k = 0; k < num_oss; ++k){
                        O_m(k, i) += O(k, (i * B + j)) / B;
                    }
                }
            }

            for (int i = 0; i < N_B; ++i){
                for (int j = 0; j < N_B; ++j){
                    if(j!=i){//exclude one different block each cycle
                        for (int k = 0; k < num_oss; ++k){
                            O_m_jack(k, i) += O_m(k, j) / (N_B - 1.0);
                        }
                    }
                }
                for (int k = 0; k < num_oss; ++k){
                    O_media(k) += O_m_jack(k, i) / N_B;
                }
            }

            for (int i = 0; i < N_B; ++i){
                for (int k = 0; k < num_oss; ++k){
                    var_O(k) += pow1(O_m_jack(k,i) - O_media(k), 2);
                }
            }
            jackknife << B << "\t";
            for (int i = 0; i < num_oss; ++i){
                if (i!=num_oss-1){
                    jackknife << sqrt((N_B - 1) * var_O(i) / N_B) << "\t";
                }
                else{
                    jackknife << sqrt((N_B - 1) * var_O(i) / N_B) << endl;
                }
            }
            N_B_prev=N_B;
        }
    }

    jackknife.close();
}
void MRT2(mat &r, double *V, double T, double J, mat &primivicini){
//random metropolis algorithm 
    double s= rand() / (RAND_MAX + 1.0);
    int m = (int)rint((N-1) * s);//find the spin for MRT2
  
    double V_mod = modifica_potenziale(r,J,m,primivicini);
    //changes the contribution to the potential for the modified spin
    
    if(V_mod>0){
        if(exp(-V_mod/T)>(rand()/((double)RAND_MAX+1.0))){
            r(m,2) *= -1;
            magn +=  2 * r(m, 2) / N;
            *V += V_mod;
        }
    }
    else{//if the change in energy is less than zero then always accept
        r(m,2) *= -1;
        magn +=  2 * r(m, 2) / N;
        *V += V_mod;
    }
    return; 
}
void MRT2(mat &r, double *V, double T, double J, mat &primivicini, int i){
//sequential metropolis algorithm 
    int m = i%N;//find the spin for MRT2
    double V_mod = modifica_potenziale(r,J,m,primivicini);
    //changes the contribution to the potential for the modified spin
    
    if(V_mod>0){
        if(exp(-V_mod/T)>(rand()/((double)RAND_MAX+1.0))){
            r(m,2) *= -1;
            magn +=  2 * r(m, 2) / N;
            *V += V_mod;
        }
    }
    else{//if the change in energy is less than zero then always accept
        r(m,2) *= -1;
        magn +=  2 * r(m, 2) / N;
        *V += V_mod;
    }
    return; 
}
void spin_inverter(mat &r, double sp_c, mat &p_v, double prob, rowvec &spin_m){
//selects the spin for the cluster update
    spin_m(sp_c) = 1;

    for (int i = 0; i < pv; ++i){//cycle over the first nearest neighburs
        int n = p_v(sp_c, i);//spin to add to the cluster?
        if (r(sp_c,2) == r(n, 2) && spin_m(n) == -1){
            if(prob > (rand()/(RAND_MAX + 1.))){
                spin_inverter(r, n, p_v, prob, spin_m);
            }
        }
    }
}
void Wolff(mat &r, double *V, double T, double J, mat &primivicini){
// wolff algorithm
    double s= rand() / (RAND_MAX + 1.0);
    int m = (int)rint((N-1) * s);//find the first spin for the cluster
    
    rowvec spin_da_mod(N, fill::value(-1));

    double prob = 1 - exp(- 2 * J / T);
    dimensione = 0;
    
    spin_inverter(r, m, primivicini, prob, spin_da_mod);//selects the cluster

    for (int i = 0; i < N; ++i){
        if(spin_da_mod(i) == 1){//invert the cluster
            dimensione++;
            r(i, 2) *= -1;
            magn += 2 * r(i,2) / N;
        }
    }
    *V = 0;
    for (int i = 0; i < N; ++i){//recalculates the potential
        *V += potenziale_ising(r, J, i, primivicini);
    }
    *V/=2;
}
void err_max_blocking(int N_t, rowvec &sigma, int num_oss){
//blocking method to find the max error to use on the errorbars for the program
    mat O(num_oss, N_t);

    ifstream dati_blocking;
    dati_blocking.open("out/dati_blocking.txt");
    for (int i = 0; i < N_t; ++i){
        for (int j = 0; j < num_oss; ++j){
            dati_blocking >> O(j, i);
        }
    }
    dati_blocking.close();

    int N_B_prev=0;
    sigma.fill(0);
    for (int B = 1; B < N_t/4; B+=10){

        int N_B = floor(N_t / B);
        
        if(N_B!=N_B_prev){
            
            mat O_mB(num_oss, N_B, fill::value(0));
            //mean observable of the block
            rowvec var_O(num_oss, fill::value(0));
            //variance of each observable
            rowvec O_mean(num_oss,fill::value(0));
            //complessive mean for each observable
            
            for (int i = 0; i < N_B; ++i){//cycle over the blocks
                for (int j = 0; j < B; ++j){
                    for (int k = 0; k < num_oss; ++k){
                        O_mB(k, i) += O(k, (i * B + j)) / B;
                    }
                }
                for (int k = 0; k < num_oss; ++k){
                    O_mean(k) += O_mB(k,i) / N_B;
                    //calculates the mean over all the blocks
                }
            }
            for (int k = 0; k < num_oss; ++k){
                for (int i = 0; i < N_B; ++i){
                    var_O(k) += pow1(O_mB(k,i) - O_mean(k),2) / N_B;
                }
            }
            for (int i = 0; i < num_oss; ++i){
                if (sigma(i)<sqrt(var_O(i)/N_B)){
                    sigma(i) = sqrt(var_O(i)/N_B);
                }
            }
            
            N_B_prev=N_B;
        }
    }
}
void plot_config() {//executes the plot of the configurations
    string comando;
    comando = "gnuplot";
    comando += " out/config.plt";
    cout << comando << endl;
    system(comando.c_str());
    return;
}
void configurations(mat &r, int iter, double L){
//configurations in iterations at the same temperature
    ofstream config;
    config.open("out/config.plt");
    config<<"set term png"<<endl;
    config<< "set output \"out/config/ising_2d_iter_"<<iter+1<<".png\""<<endl;
    config<<"set xrange [ 0 :   "<<L<<" ]"<<endl;
    config<<"set yrange [ 0 :   "<<L<<" ]"<<endl;
    config<<"set nokey"<<endl;
    config<<"set title \"Configuration at "<<iter+1<<" iterations for N="<<L;
    config<<"x"<<L<<"\""<<endl;
    config<<"unset tics"<<endl;
    config<<"set size ratio    1.00000"<<endl;
    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j){
            if(r(i*L+j,2)==1){
                config<<"set object rectangle from "<<i<<", "<< j <<" to " <<i+1;
                config<<", "<< j+1 << "fc rgb \"white\""<<endl;
            }
            else{
                config<<"set object rectangle from "<<i<<", "<< j <<" to " <<i+1;
                config<<", "<< j+1 << "fc rgb \"black\""<<endl;
            }
        }
    }
    config<<"plot 1"<<endl;
    config<<"quit"<<endl;
    config.close();
    plot_config();
}
void configurations(mat &r, double T, double L){//configurations in temperature
    ofstream config;
    config.open("out/config.plt");
    config<<"set term png"<<endl;
    config<< "set output \"out/config/ising_2d_T_"<<T<<".png\""<<endl;
    config<<"set xrange [ 0 :   "<<L<<" ]"<<endl;
    config<<"set yrange [ 0 :   "<<L<<" ]"<<endl;
    config<<"set nokey"<<endl;
    config<<"set title \"Configuration at temperature "<<T<<" for N="<<L;
    config<<"x"<<L<<"\""<<endl;
    config<<"unset tics"<<endl;
    config<<"set size ratio    1.00000"<<endl;
    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j){
            if(r(i*L+j,2)==1){
                config<<"set object rectangle from "<<i<<", "<< j <<" to " <<i+1;
                config<<", "<< j+1 << "fc rgb \"white\""<<endl;
            }
            else{
                config<<"set object rectangle from "<<i<<", "<< j <<" to " <<i+1;
                config<<", "<< j+1 << "fc rgb \"black\""<<endl;
            }
        }
    }
    config<<"plot 1"<<endl;
    config<<"quit"<<endl;
    config.close();
    plot_config();
}
#endif 
