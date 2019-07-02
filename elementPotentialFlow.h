#ifndef ELEMENTPOTENTIALFLOW_H_INCLUDED
#define ELEMENTPOTENTIALFLOW_H_INCLUDED
class elementPotentialFlow{
private:
    int elem_number,n;//n is the rank number of the node
    int elem_max,n_max;//total elem and node number




//   int n[3][elem_max];

   //struct nodeTable{
  // int n[3][elem_number];
   //int n2[elem_number];
  // int n3[elem_number];
   //};

   struct nodeCoordinate{
//   double x_coordinate[n_max];
 //  double y_coordinate[n_max];
   };

   struct firstBCnode{


   };
public:
    elementPotentialFlow(){}
    ~elementPotentialFlow(){}
    void readData(int &nodeTable,double &nodeCoordinate,double &firstBCnode,double );
    double** elementMatrix(int elem_number);
    double** monothilicMatrix(int n);
//    void  guass(double a[n+1][n+1],double b[n+1],int n)




};


#endif // ELEMENTPOTENTIALFLOW_H_INCLUDED
