/**In this header file, the SPH class is going to be defined
 * with the appropriate constructors and destructors
 * and the member functions. The SPH class is going to admit a matrix
 * as an input, as well as the number of particles.
 * Details on the functions, the arrays and the overloadings
 * can be found in the file : class.cpp
 **/

using namespace std;


class SPH{

private:

    unsigned int n;                         //size of the Matrix

    int t;                                  //time

    double dt;                              //timestep

    double h;                               //Radius of influence

    int rank;                               //Rank of processor

    int np;                                 //Total number of processors

    double ** vMatrix = new double*[4];     //AlloCating space memory for the Matrix

    //Constants of the problem
    double k =2000.0;
    double roo = 1000.0; 
    double m = 1.0;
    double miou = 1.0;
    double e = 0.5;
    const double g = 9.81;
    const int root = 0;

public:
 
    double * x;
    double * y ;
    double * vx;
    double * vy ;
    double * r ;
    double * x_new;
    double * y_new ;
    double * vx_new;
    double * vy_new ;
    int * npRank;
    double * ro;
    double * p ;
    double * vi;
    double * q;
    double * Freduce = new double[4];
    double * Ftake = new double[4];
    double * xyv ;
    double E_k,E_p;                                 //Energies
    double Fip,Fiv,Fig,Fipx,Fivx,Fipy,Fivy,Figy;    //Forces
    double Figx=0.0;
    int i,j;


    /******** CONSTRUCTORS/DESTRUCTOR********/

    SPH() = default;                                           //Default constructor

    ~SPH();                                                    //Destructor
    
    SPH(const unsigned n_new);                                 //User defined constructor for allocating the dimensions of the Matrix


    /**********OVERLOADINGS**********/
    
    double & operator () (unsigned row, unsigned col);

    int & operator > (unsigned t);                            

    double & operator >> (double dt);  

    double & operator < (double h_new);

    int & operator | (int rank_new);  

    int & operator ^ (int np_new);  

    /**********MEMBER-FUNCTIONS*********/
    
    void x0();
    
    void y0();
    
    void vx0();
    
    void vy0();

    void rVec();

    void den();

    void pressure();

    double FiP(double * x_y);

    double FiV(double * v);

    double FiG();

    double vInit(double * v,double & Fip, double & Fiv, double & Fig);

    double v_x_y(double * v, double & Fip, double & Fiv, double & Fig);

    void spatial();

    double rety(int l);

    double retx(int l);

    double Ek();

    double Ep();

    void mass();

    void nprank();

    void getdata();

};