#include"elementPotentialFlow.h"

double** elementPotentialFlowe::ementMatrix(int elem_number){
    int elem_max;
    double b[3],c[3];
    double x1,x2,y1,y2;
    double A[elem_max];
    double A_e[3][3];


    x1=x_coordinate[n1[elem_number]];
    x2=x_coordinate[n2[elem_number]];
    x3=x_coordinate[n3[elem_number]];

    y1=y_coordinate[n2[elem_number]];
    y2=y_coordinate[n2[elem_number]];
    y3=y_coordinate[n3[elem_number]];

    A[elem_number]=((x2-x1)*(y3-y1)-(y2-y1)(x3-x1))/2.0;

    b[0]=(y2-y3)/(2*A[elem_number]);
    b[1]=(y3-y1)/(2*A[elem_number]);
    b[2]=(y1-y2)/(2*A[elem_number]);

    c[0]=(x3-x2)/(2*A[elem_number]);
    c[1]=(x1-x3)/(2*A[elem_number]);
    c[2]=(x2-x1)/(2*A[elem_number]);

    for(int i=0;i<3;i++){
        for(int j=0;j<3:j++){
            A_e[i][j]=A[elem_number]*(b[i]*b[j]+c[i]*c[j]);
            }}
return A_e;
}

//made up of the monothilic matrix
//the monothilic matrix responding to the nodeNumber
double** elementPotentialFlow::monothilicMatrix(int n_max){//n_max is the total number of the node
int a_i,a_j;
double **A;
double A_t[n_max][n_max]={0};
for(int k=1;k<elem_max;k++){
    A[3][3]=elementMatrix(k);
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            a_i=n[i][k];
            a_j=n[j][k];
            A_t[a_i][a_j]+=A[i][j];
        }
      }
      }
      return A_t;
    }


void elementPotentialFlow::eliminateLineForBc(double **a,double *b){
double A_t[n_max][n_max];
double b[n_max];
int nr[nr_max];//first BC
double ur[nr_max];//the value of the first BC
   for(int i=0;i<nr_max;i++){
        for(int j=0;j<nr_max;j++){f[j]=f[j]-ur[i]*A_t[j][nr[i];}
        //f[nr[i]]=f[nr[i]]-ur[i]*A_t[i][nr[i];
       for(int i=0;i<n_max;i++){A[i][nr[i]]=0;}
       for(int i=0;i<n_max;i++){A[nr[i]][i]=0;}
       A_t[nr[i][nr[i]]=1;

}
}

void elementPotentialFlow::eliminateEnlargeForBc(double **a,double *b,nr_max){
double A_t[n_max][n_max];
double b[n_max];
int nr[nr_max];
double ur[nr_max];
for(int i=0;i<nr_max;i++){
    f[nr[i]]=*ur[i]*pow(10.0,20.0);
    A_t[nr[i]][nr[i]]=A_t[nr[i]][nr[i]]*pow(10.0,20.0);
}
}


void elementPotentialFlow::Gauss_seidel(double* x, double** a,double* b,int n){//x need give the preliminary cvalue
    do{                                                                        //can use  a loop to initialz the x
        flag=0;
        iteration++;
        for(int i=0;i<n;i++){
            xold=x[i];
            sum=0.0;
            for(j=0;j<n;j++){
                if(j!=i){sum+=a[i][j]*x[j];}
            }
            x[i]=(b[i]-sum)/a[i][i];
            error=fabs(xold-x[i])/x[i];
            if(error>=eps){flag=1;}
         }
   }while(flag==1);

}


void  elementPotentialFlow::guass(double a[n+1][n+1],double b[n+1],int n){//the parameters are the matrix and the order of the matrix
        //after execution the b is the flow function value matrix
   int i,j,k,i0,ip1,kp1;
	double C,T;

   //double **a=new double* [n+1];
 //  for(i=0;i<n+1;i++) a[i]=new double [n+1];
   // new double [n+1] b;

	for(k=1;k<=(n-1);k++){
		C=0;
		for(i=k;i<=n;i++){
			if(abs(a[i][k]>abs(C))){
				C=a[i][k];
				i0=i;}}
		if(i0!=k){
			for(j=k;j<=n;j++){
				T=a[k][j];
				a[k][j]=a[i0][j];
				a[i0][j]=T;}
			T=b[k];
			b[k]=b[i0];
			b[i0]=T;}
		  kp1=k+1;
		  C=1/C;
		  b[k]=b[k]*C;
		for(j=kp1;j<=n;j++){
			a[k][j]=a[k][j]*C;
			for(i=kp1;i<=n;i++)
				a[i][j]=a[i][j]-a[i][k]*a[k][j];
			b[j]=b[j]-a[j][k]*b[k];}
    }
	b[n]=b[n]/a[n][n];
	for(k=1;k<=(n-1);k++){
		i=n-k;
		C=0;
		ip1=i+1;
		for(j=ip1;j<=n;j++)
			C=C+a[i][j]*b[j];
		b[i]=b[i]-C;
	}
   //for(i=0;i<n+1;i++){delete [] a[i];}
  // delete [] a;
   // delete [] b;
   }







