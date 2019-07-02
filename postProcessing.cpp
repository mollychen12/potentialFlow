#include"elemntPotentailFlow.h"
//to get the velocity of the element
void vesln(int ne)
{
	int i;
	int j;
	int k;
	double xe[4],ye[4];
	double d;
	double x_grid_d[4],y_grid_d[4];
	for(i=1;i<=3;i++)
	{
		k=N[i][ne];
		xe[i]=X[k];
		ye[i]=Y[k];
		ue[i]=fu[k];
	}
	d=xe[2]*ye[3]-xe[3]*ye[2]+xe[3]*ye[1]-xe[1]*ye[3]+xe[1]*ye[2]-xe[2]*ye[1];
	for(i=1;i<=3;i++)
	{
		j=i+1;
		if(j>3)
			j=j-3;
		k=j+1;
		if(k>3)
			k=k-3;
		x_grid_d[i]=(ye[j]-ye[k])/d;
		y_grid_d[i]=(xe[k]-xe[j])/d;
	}
	VYE[NE]=UE[1]*C[1]+UE[2]*C[2]+UE[3]*C[3];
	VXE[NE]=UE[1]*B[1]+UE[2]*B[2]+UE[3]*B[3];
}


//to get the pressure value of the node
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
