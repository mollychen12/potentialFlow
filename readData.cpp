#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<sstream>

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
double x_coordinate[n_max+1];
double y_coordinate[n_max+1];
int n[3][elem_max];
int nr[lineElementN];
double ur[lineElementN];



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




//for test
int main(){
//read_coordinate();
//read_element_node();
//read_ycoordinate();
//read_element_node();
read_firstBc();
}

