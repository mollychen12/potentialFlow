#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#define n_max 145

using namespace std;
double x_coordinate[145];
void readAnyArray(int arrayN){
    int i=0;
ifstream myfile("mesh1.msh");
ifstream ftmp;
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
vector<double> radius;
cout<<"the data read"<<endl;
for(auto it=vec.begin();it !=vec.end();it++){
    cout<<*it<<endl;
    istringstream is(*it);//use every line as a string stream
    string s;
    int pam=0;

    while(is>>s){//以空格为界，把istringstream中数据格式取出放入到依次s中
            if(pam==(arrayN-1))//get the sixth array data
            {
                double r=atof(s.c_str());//datatype casting ,convert the string to double
                radius.push_back(r);
            }
        pam++;
    }
}
//cout<<"the sixth array of the readDate is"<<endl;
for(auto it=radius.begin();it!=radius.end();it++){
   // cout<<*it<<endl;
    fout<<*it<<endl;

}
system("pause");
ftmp.open("mesh1.txt");
vector<double> x;
while(ftmp){
   double idata;
   ftmp>>idata;
   x.push_back(idata);
   }
for(int i=0;i<n_max;i++){
    x_coordinate[i]=x[i];
    cout<<"x["<<i<<"]"<<x_coordinate[i]<<endl;
   }

}



int main(){
readAnyArray(2);
return 1;
}
