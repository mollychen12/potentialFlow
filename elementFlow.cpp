#include "stdio.h"
#include "math.h"
#include<cstdio>
#include<cmath>
FILE *fp1;
FILE *fp2;
FILE *fp3;
FILE *fp4;
int TNP;
int TNE;
int TNR;
int TNS;
int N[4][10000];
double X[8000],Y[8000],Z[8000];
int NR[8000],NS[8000];
double UR[5000],GS[5000];
int ND[10000];
double AIJ[4][4];
double FE[4];
int NDE;
double ANM[8000][8000];
double FU[8000];
double UE[4];
double VXE[10000],VYE[10000];
double VXP[8000],VYP[8000];
double P[8000];

void readdata()
{
	int e;
	int i;
	fp1=fopen("显示输入数据.txt","w");
	fp2=fopen("原始数据.txt","r");
	fp3=fopen("计算结果.txt","w");
	fp4=fopen("作图数据.txt","w");

	fscanf(fp2,"%d",&TNP);
	fprintf(fp1,"TNP(节点总数)=%d\n",TNP);

	fscanf(fp2,"%d",&TNE);
	fprintf(fp1,"TNE(单元总数)=%d\n",TNE);

	fscanf(fp2,"%d",&TNR);
	fprintf(fp1,"TNR(本质边界节点总数)=%d\n",TNR);

	fscanf(fp2,"%d",&TNS);
	fprintf(fp1,"TNS(自然边界节点总数)=%d\n",TNS);

	for(e=1;e<=TNE;e++)
		for(i=1;i<=3;i++)
			fscanf(fp2,"%d",&N[i][e]);
	for(e=1;e<=TNE;e++)
		for(i=1;i<=3;i++)
			fprintf(fp1,"第%d个单元第%d个节点的总体节点号：%d\n",e,i,N[i][e]);

	for(i=1;i<=TNP;i++)
		fscanf(fp2,"%lf%lf",&X[i],&Y[i]);
	for(i=1;i<=TNP;i++)
		fprintf(fp1,"第%d个节点的坐标(m)：%f %f\n",i,X[i],Y[i]);

	for(i=1;i<=TNR;i++)
		fscanf(fp2,"%d",&NR[i]);
	fprintf(fp1,"本质边界的节点号为：");
	for(i=1;i<=TNR;i++)
		fprintf(fp1,"%d ",NR[i]);
	fprintf(fp1,"\n");

	for(i=1;i<=TNS;i++)
		fscanf(fp2,"%d",&NS[i]);
	fprintf(fp1,"自然边界的节点号为：");
	for(i=1;i<=TNS;i++)
		fprintf(fp1,"%d ",NS[i]);
	fprintf(fp1,"\n");

	for(i=1;i<=TNR;i++)
		fscanf(fp2,"%lf",&UR[i]);
	fprintf(fp1,"本质边界各节点的值为：");
	for(i=1;i<=TNR;i++)
		fprintf(fp1,"%f ",UR[i]);
	fprintf(fp1,"\n");

	for(i=1;i<=TNS;i++)
		fscanf(fp2,"%lf",&GS[i]);
	fprintf(fp1,"自然边界各节点的值为：");
	for(i=1;i<=TNS;i++)
		fprintf(fp1,"%f ",GS[i]);
	fprintf(fp1,"\n");

}



void elmt(int NE)
{
	int i;
	int j;
	int k;
	int KI;
	int m;
		int ii;
	int K2;
	int K3;
	double XE[4],YE[4];
	double B[4],C[4];
	int KG[4];
	double D;
	double G2,G3;
	double SL;
	for(i=1;i<=3;i++)
	{
		k=N[i][NE];
		XE[i]=X[k];
		YE[i]=Y[k];
		FE[i]=0;
	}
	D=XE[2]*YE[3]-XE[3]*YE[2]+XE[3]*YE[1]-XE[1]*YE[3]+XE[1]*YE[2]-XE[2]*YE[1];
	for(i=1;i<=3;i++)
	{
		j=i+1;
		if(j>3)
			j=j-3;
		k=j+1;
		if(k>3)
			k=k-3;
		B[i]=(YE[j]-YE[k])/D;
		C[i]=(XE[k]-XE[j])/D;
	}
	for(i=1;i<=3;i++)
		for(j=1;j<=3;j++)
			AIJ[i][j]=0.5*(B[i]*B[j]+C[i]*C[j])*D;
	if(NDE==1)
	{
		for(i=1;i<=3;i++)
		{
			KI=N[i][NE];
			for(m=1;m<=TNS;m++)
			{
				if(KI==NS[m])
				{
					KG[i]=m;
					j=1;
					break;
				}
				else
					j=2;
			}
			if(j==1)
				continue;
			else if (j==2)
				ii=i;
		}
		j=ii+1;
		if(j>3)
			j=j-3;
		k=j+1;
		if(k>3)
			k=k-3;
		K2=KG[j];
		K3=KG[k];
		G2=GS[K2];
		G3=GS[K3];
		SL=sqrt((XE[j]-XE[k])*(XE[j]-XE[k])+(YE[j]-YE[k])*(YE[j]-YE[k]));
		FE[j]=SL*(2*G2+G3)/6;
		FE[k]=SL*(G2+2*G3)/6;
	}
}



void gauss(double A[8000][8000],double B[8000],int n)
{
	int nm1;
	int i;
	int j;
	int k;
	int i0;
	int ip1;
	int kp1;
	double C;
	double T;
	nm1=n-1;
	for(k=1;k<=nm1;k++)
	{
		C=0;
		for(i=k;i<=n;i++)
		{
			if(abs(A[i][k]>abs(C)))
			{
				C=A[i][k];
				i0=i;
			}
		}
		if(i0!=k)
		{
			for(j=k;j<=n;j++)
			{
				T=A[k][j];
				A[k][j]=A[i0][j];
				A[i0][j]=T;
			}
			T=B[k];
			B[k]=B[i0];
			B[i0]=T;
		}
		kp1=k+1;
		C=1/C;
		B[k]=B[k]*C;
		for(j=kp1;j<=n;j++)
		{
			A[k][j]=A[k][j]*C;
			for(i=kp1;i<=n;i++)
				A[i][j]=A[i][j]-A[i][k]*A[k][j];
			B[j]=B[j]-A[j][k]*B[k];
		}
	}
	B[n]=B[n]/A[n][n];
	for(k=1;k<=nm1;k++)
	{
		i=n-k;
		C=0;
		ip1=i+1;
		for(j=ip1;j<=n;j++)
			C=C+A[i][j]*B[j];
		B[i]=B[i]-C;
	}
}

void vesln(int NE)
{
	int i;
	int j;
	int k;
	double XE[4],YE[4];
	double D;
	double B[4],C[4];
	for(i=1;i<=3;i++)
	{
		k=N[i][NE];
		XE[i]=X[k];
		YE[i]=Y[k];
		UE[i]=FU[k];
	}
	D=XE[2]*YE[3]-XE[3]*YE[2]+XE[3]*YE[1]-XE[1]*YE[3]+XE[1]*YE[2]-XE[2]*YE[1];
	for(i=1;i<=3;i++)
	{
		j=i+1;
		if(j>3)
			j=j-3;
		k=j+1;
		if(k>3)
			k=k-3;
		B[i]=(YE[j]-YE[k])/D;
		C[i]=(XE[k]-XE[j])/D;
	}
	VYE[NE]=UE[1]*C[1]+UE[2]*C[2]+UE[3]*C[3];
	VXE[NE]=UE[1]*B[1]+UE[2]*B[2]+UE[3]*B[3];
}

void vpsln(int n)
{
	int i;
	int j;
	int NE;
	j=0;
	VXP[n]=0;
	for(NE=1;NE<=TNE;NE++)
	{
		for(i=1;i<=3;i++)
		{
			if(n==N[i][NE])
			{
				VXP[n]=VXP[n]+VXE[NE];
				VYP[n]=VYP[n]+VYE[NE];
				j=j+1;
			}
		}
	}
	VXP[n]=VXP[n]/j;
	VYP[n]=VYP[n]/j;
	P[n]=0.5*(1-VXP[n]*VXP[n]-VYP[n]*VYP[n]);
}

void main()
{
	int i;
	int j;
	int k;
	int e;
	int c;
	int ii;
	int jj;
	int NE;
	readdata();
	for (NE=1;NE<=TNE;NE++)
	{
		c=0;
		for(i=1;i<=3;i++)
		{
			for(j=1;j<=TNS;j++)
			{
				if(N[i][NE]==NS[j])
				{
					c++;
					break;
				}
			}
		}
		if(c==2)
			ND[NE]=1;
		else
			ND[NE]=0;
	}
	for(i=1;i<=TNP;i++)
	{
		FU[i]=0;
		for(j=1;j<=TNP;j++)
			ANM[i][j]=0;
	}
	for(NE=1;NE<=TNE;NE++)
	{
		NDE=ND[NE];
		elmt(NE);
		for(i=1;i<=3;i++)
		{
			ii=N[i][NE];
			FU[ii]=FU[ii]+FE[i];
			for(j=1;j<=3;j++)
			{
				jj=N[j][NE];
				ANM[ii][jj]=ANM[ii][jj]+AIJ[i][j];
			}
		}
	}
	for(i=1;i<=TNR;i++)
	{
		k=NR[i];
		ANM[k][k]=1e20*ANM[k][k];
		FU[k]=ANM[k][k]*UR[i];
	}
	gauss(ANM,FU,TNP);
	for(i=1;i<=TNP;i++)
		fprintf(fp3,"第%d个节点的流函数值为：%f\n",i,FU[i]);
	for(NE=1;NE<=TNE;NE++)
	{
		vesln(NE);
		fprintf(fp3,"第%d个单元内，x方向速度值：%f，y方向速度值：%f\n",NE,VXE[NE],VYE[NE]);
	}
	for(i=1;i<=TNP;i++)
	{
		vpsln(i);
		fprintf(fp3,"第%d个节点，x方向速度值：%f，y方向速度值：%f，压力值：%f\n",i,VXP[i],VYP[i],P[i]);
		fprintf(fp4,"%f %f %f %f %f %f\n",X[i],Y[i],VXP[i],VYP[i],P[i],FU[i]);
	}
	for(e=1;e<=TNE;e++)
	{
		for(i=1;i<=3;i++)
			fprintf(fp4,"%d ",N[i][e]);
		fprintf(fp4,"\n");
	}
}





