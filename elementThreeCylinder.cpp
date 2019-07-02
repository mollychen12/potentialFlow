#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<cmath>

//#include"elementPotentialFlow.h"

using namespace std;

const int in=1;
const int out=2;
const int up=3;
const int down=4;
const int cylinder1=5;
const int cylinder2=6;
const int cylinder3=7;

int boundaryNum=7;
int lboundary=2;
int numOfSeeds=9;
const int lineElementN=63;
const int n_max=145;
const int elem_maxt=294;
const int elem_max=231;
double x_coordinate[n_max];
double y_coordinate[n_max];
int n[3][elem_max];
int nr[lineElementN];
double ur[lineElementN];
double A_t[n_max][n_max];
double f[n_max];

//double A[3][3];
double kesai[n_max];
double vlctXOfElement[elem_max];
double vlctYOfElement[elem_max];
double vx[n_max],vy[n_max],p[n_max];

double eps=10e-6;//eos shi error tolerance



void readAnyArray(int arrayN){
    int i=0;
ifstream myfile("mesh1.msh");
//ifstream ftmp;
ofstream fout("mesh1.txt");
if(!myfile.is_open()){
    cout<<"Unable to open my file";
    system("pause");
    exit(1);
}
vector<string> vec;
string temp;
while(getline(myfile,temp)){//use getline to read every line
    vec.push_back(temp);
}
vector<double> needed;
for(auto it=vec.begin();it !=vec.end();it++){
    istringstream is(*it);//use every line as a string stream
    string s;
    int pam=0;
    while(is>>s){//以空格为界，把istringstream中数据格式取出放入到依次s中
            if(pam==(arrayN-1))//get the  array data
            {
                double r=atof(s.c_str());//datatype casting ,convert the string to double
                needed.push_back(r);
            }
        pam++;
    }
}
for(auto it=needed.begin();it!=needed.end();it++){
  // cout<<"the "<<arrayN<<" Array is "<<*it<<endl;
    fout<<*it<<endl;
}
//system("pause");
}


void read_coordinate(){
    ifstream ftmp1;
    readAnyArray(2);
    ftmp1.open("mesh1.txt");
vector<double> xcoordinate;
while(ftmp1){
   double idata;
   ftmp1>>idata;
   xcoordinate.push_back(idata);
   }
for(int i=0;i<n_max;i++){
    x_coordinate[i]=xcoordinate[i];
    cout<<"x["<<i<<"]"<<x_coordinate[i]<<endl;
   }


    ifstream ftmp2;
    readAnyArray(3);
    ftmp2.open("mesh1.txt");
vector<double> ycoordinate;
while(ftmp2){
   double idata;
   ftmp2>>idata;
   ycoordinate.push_back(idata);
   }
for(int i=0;i<n_max;i++){
    y_coordinate[i]=ycoordinate[i];
    cout<<"y["<<i<<"]"<<y_coordinate[i]<<endl;
   }
}

void read_element_node(){
    ifstream ftmpn1;
    ifstream ftmpn2;
    ifstream ftmpn3;

    readAnyArray(6);
    ftmpn1.open("mesh1.txt");
vector<int> node1;
while(ftmpn1){
   double idata1;
   ftmpn1>>idata1;
   node1.push_back(idata1);
   }
for(int i=0;i<elem_max;i++){
    n[0][i]=node1[i+lineElementN];
   cout<<"n[0]["<<i<<"]"<<n[0][i]<<endl;
   }
ftmpn1.close();
    readAnyArray(7);
    ftmpn2.open("mesh1.txt");
vector<int> node2;
while(ftmpn2){
   double idata2;
   ftmpn2>>idata2;
   node2.push_back(idata2);
   }
for(int i=0;i<elem_max;i++){
    n[1][i]=node2[i+lineElementN];
    cout<<"n[1]["<<i<<"]"<< n[1][i]<<endl;
   }
ftmpn2.close();

   readAnyArray(8);
    ftmpn3.open("mesh1.txt");
vector<int> node3;
while(ftmpn3){
   double idata3;
   ftmpn3>>idata3;
   node3.push_back(idata3);
   }
for(int i=0;i<elem_max;i++){
    n[2][i]=node3[i];
    cout<<"n[2]["<<i<<"]"<<  n[2][i]<<endl;
   }
ftmpn3.close();
}



void read_firstBc(){
    ifstream ftmpr1;
    ifstream ftmpr2;
readAnyArray(6);
    ftmpr1.open("mesh1.txt");
vector<int> noder1;
while(ftmpr1){
   double idatar1;
   ftmpr1>>idatar1;
   noder1.push_back(idatar1);
   }
   //should pick out the lboundary node keep the rboundary node
   //should give the value to each rboundary node,
   //in: ur=y_coordinate,out:ur=x_coordinate,up&down :ur=3.5,cylinders:ur=0

for(int i=0;i<lineElementN;i++){
    nr[i]=noder1[i];

   if(i<out*numOfSeeds) ur[i]=y_coordinate[nr[i]];
else if(i<down*numOfSeeds) ur[i]=3.5;
       else if(i<lineElementN) ur[i]=0;
       cout<<"nr["<<i<<"]"<<nr[i]<<endl;
       cout<<"ur["<<i<<"]"<<ur[i]<<endl;
   }

ftmpr1.close();

}



void elmentMatrix(int elem_num,double A_e[][3]){
    int el_n=elem_num;
    double b[3],c[3];
    double x1,x2,x3,y1,y2,y3;
    double A[elem_max];
   // double A_e[3][3];


    x1=x_coordinate[n[0][el_n]];
    x2=x_coordinate[n[1][el_n]];
    x3=x_coordinate[n[2][el_n]];

    y1=y_coordinate[n[0][el_n]];
    y2=y_coordinate[n[1][el_n]];
    y3=y_coordinate[n[2][el_n]];

    A[el_n]=0.5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));

    b[0]=(y2-y3)/(2*A[el_n]);
    b[1]=(y3-y1)/(2*A[el_n]);
    b[2]=(y1-y2)/(2*A[el_n]);

    c[0]=(x3-x2)/(2*A[el_n]);
    c[1]=(x1-x3)/(2*A[el_n]);
    c[2]=(x2-x1)/(2*A[el_n]);

    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            A_e[i][j]=A[el_n]*(b[i]*b[j]+c[i]*c[j]);
            }}

}

//made up of the monothilic matrix
//the monothilic matrix responding to the nodeNumber
void monothilicMatrix(double A_t[][n_max]){//n_max is the total number of the node
int a_i,a_j;
double A[3][3];

//double A_t[n_max][n_max];
for(int k=1;k<elem_max;k++){
        //elementMatrix(k,A[3][3]);
        elmentMatrix(k,A);
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            a_i=n[i][k];
            a_j=n[j][k];
            A_t[a_i][a_j]+=A[i][j];
        }
      }
}
}


void eliminateLineForBc(double **a,double *b){//double A_t[n_max][n_max];
    //double b[n_max];
    //int nr[nr_max];
    //first BC//double ur[nr_max];
    //the value of the first BC
int nr_max=lineElementN;
   for(int i=0;i<nr_max;i++){
        for(int j=0;j<nr_max;j++){b[j]=b[j]-ur[i]*a[j][nr[i]];}
        //f[nr[i]]=f[nr[i]]-ur[i]*A_t[i][nr[i];
       for(int i=0;i<n_max;i++){a[i][nr[i]]=0;}
       for(int i=0;i<n_max;i++){a[nr[i]][i]=0;}
       a[nr[i]][nr[i]]=1;

}
}

void eliminateEnlargeForBc(double a[][n_max],double *b,int nr_max){
    for(int i=0;i<nr_max;i++){
    b[nr[i]]=a[nr[i]][nr[i]]*ur[i]*pow(10.0,20.0);
    a[nr[i]][nr[i]]=a[nr[i]][nr[i]]*pow(10.0,20.0);
    }
}


void Gauss_seidel(double* x, double a[][n_max],double* b,int n){    //x need give the preliminary cvalue
    int flag;
    int iteration=0;
    double xold,sum,error;

    do{                                                   //can use  a loop to initialz the x
        flag=0;
        iteration++;
        for(int i=0;i<n;i++){
            xold=x[i];
            sum=0.0;
            for(int j=0;j<n;j++){
                if(j!=i){sum+=a[i][j]*x[j];}
            }
            x[i]=(b[i]-sum)/a[i][i];
            error=fabs(xold-x[i])/x[i];
            if(error>=eps){flag=1;}
         }
   }while(flag==1);
for(int i=0;i<n_max;i++)
kesai[i]=x[i];s
}

void velocityOfElement(int elemN){//use the linear interpolation
int j,k;
double xOfE[3],yOfE[3],uOfE[3];
double area;//the area of the elem,use the determinants
double x_interpolation[3],y_interpolation[3];
for(int i=0;i<3;i++){
    k=n[i][elemN];
    xOfE[i]=x_coordinate[k];
    yOfE[i]=y_coordinate[k];
    uOfE[i]=f[k];
}
area=xOfE[1]*yOfE[2]-xOfE[2]*yOfE[1]+xOfE[2]*yOfE[0]-xOfE[0]*yOfE[2]+xOfE[0]*yOfE[1]-xOfE[1]*yOfE[0];
for(int i=0;i<3;i++){
   j=i+1;
     if(j>=3) j=j-3;
   k=j+1;
     if(k>=3) k=k-3;
    y_interpolation[i]=(yOfE[j]-yOfE[k])/area;
    x_interpolation[i]=(xOfE[k]-xOfE[j])/area;
}
vlctXOfElement[elemN]=uOfE[0]*x_interpolation[0]+uOfE[1]*x_interpolation[1]+uOfE[2]*x_interpolation[2];
vlctYOfElement[elemN]=uOfE[0]*y_interpolation[0]+uOfE[1]*y_interpolation[1]+uOfE[2]*y_interpolation[2];
}



void vlctyAndPressure(int node){
int elemN;
int j=0;
vx[node]=0;
  for(elemN=0;elemN<elem_max;elemN++){
    for(int i=0;i<3;i++){
        if(node==n[i][elemN]){
            vx[node]=vx[node]+vlctXOfElement[elemN];
            vy[node]=vy[node]+vlctYOfElement[elemN];
            j++;
        }
    }
  }
vx[node]=vx[node]/j;
vy[node]=vy[node]/j;
p[node]=0.5*(1-vx[node]*vx[node]-vy[node]*vy[node]);
}




void putOutData(){
    ofstream fout;
    fout.open("output_mesh1.plt");
    fout<<" x "<<" y "<<" kesai "<<endl;
    cout<<" x "<<" y "<<" kesai "<<endl;
  //  fout<<"zoneNumber="<<n_max<<"elemenNumber"<<elem_max;
   // cout<<"zoneNumber="<<n_max<<"elemenNumber"<<elem_max;
    for(int i=0;i<n_max;i++){
        fout<<x_coordinate[i]<<" "<< y_coordinate[i]<<" ";
        cout<<x_coordinate[i]<<" "<< y_coordinate[i]<<" ";
        fout<<kesai[i]<<endl;
        cout<<kesai[i]<<endl;
    }
    for(int i=0;i<elem_max;i++){
        for(int j=0;j<3;j++){
            fout<<n[j][i]<<" ";
            cout<<n[j][i]<<" ";
        }
    }fout<<endl;
    cout<<endl;


}


//for test
int main(){
    double initx[n_max];
    double At[145][145];
 read_coordinate();
 read_element_node();
 monothilicMatrix(At);
 eliminateEnlargeForBc(At,f,lineElementN);
 for(int i=0;i<n_max;i++){
    initx[i]=0;
 }
 Gauss_seidel(initx, At,f,n_max);

 for(int i=0;i<elem_max;i++){velocityOfElement(i);}
 for(int i=0;i<n_max;i++){vlctyAndPressure(i);}

 putOutData();




}
