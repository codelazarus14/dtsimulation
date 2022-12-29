public class Embed{

private int n,total,dplus,dplusplus;
public double[] x,y;
public int[] mark;
    
Embed(){
// simplex_number=pointer_number after tidy()
n=DT.node_number;

    //System.out.println("dplus "+dplus);
x = new double[n];
y = new double[n];
}
    public void ComputeEmbedding(){
        int coord;
        formmatrix();
        coord=0;
        CG_solver(coord);
        coord=1;
        CG_solver(coord);
    }

public void sparse_mult(double[] v, double[] p){
for(int i=0;i<n;i++){
if(mark[i]==0){
p[i]=0.0;
for(int j=DT.nstart[i];j<DT.nstart[i+1];j++){
    if(mark[DT.ncol[j]]==0){
        p[i]=p[i]+DT.nlap[j]*v[DT.ncol[j]];
        }
}
}
}
}

	
    public void formmatrix(){
        int i;
        mark=new int[n];
        
        //System.out.println("Boundary in formmatrix "+ DT.boundary_length);
        for(i=0;i<n;i++){
            mark[i]=0;}
        
        mark[0]=DT.boundary_length+1;
        for(i=0;i<DT.boundary_length;i++){
            mark[DT.boundary[i]]=i+1;
          //  System.out.println("boundary mark "+DT.boundary[i]);
        }
    }
    
    public double dot(double[] a, double[] b){
        double d=0.0;
        for(int i=0;i<n;i++){
        if(mark[i]==0){
                d+=a[i]*b[i];}}
        return(d);
        }

    // remove one vertex from sphere. This makes disk.
    // Pin neighbour triangles around circle and solve for
    // interior points placing each at center of mass of its
    // neighbours
    


public void CG_solver(int coord){
double[] R_0,R_1,B,P_0,P_1;
double[] SS,SOL_1,SOL;
double alpha,beta,rrtmp,rrtmp2,resid,psdot;
int count=0,i,j,itmp,done;
double CG_RESIDUAL= 0.0000001;
    
   // System.out.println("In CG_solver");
    
    R_0=new double[n];
    R_1=new double[n];
    P_0=new double[n];
    P_1=new double[n];
    SS=new double[n];
    SOL=new double[n];
    SOL_1=new double[n];
    B=new double[n];
    
// pin neighbours of simplex 0 to unit circle. Yields RHS/source for
// laplacian problem
        
    for(i=0;i<n;i++){
        B[i]=0.0;
    }
    
    
    for(i=0;i<n;i++){
        if(mark[i]==0){
   // bulk pt
    for(j=DT.nstart[i];j<DT.nstart[i+1];j++){
        itmp=DT.ncol[j];
       //System.out.println("itmp is "+itmp);
        
        if(mark[itmp]!=0){
        // connects to boundary pt
            //System.out.println("here ...");
        if(coord==0){
        B[i]=B[i]+Math.cos((2*Math.PI/DT.boundary_length)*mark[itmp]);}
        if(coord==1){
        B[i]=B[i]+Math.sin((2*Math.PI/DT.boundary_length)*mark[itmp]);}
                }
            }
        }
    }
    
for(i=0;i<n;i++){
SOL[i]=0.0;
P_0[i]=B[i];
    //System.out.println("B is "+B[i]);
R_0[i]=B[i];
}

do{

sparse_mult(P_0,SS);
    //System.out.println("doing sparse mult");

rrtmp=dot(R_0,R_0);
psdot=dot(SS,P_0);
 //   System.out.println("dot= "+psdot);

alpha=rrtmp/psdot;

for(i=1;i<n;i++){
R_1[i]=R_0[i]-alpha*SS[i];}

for(i=1;i<n;i++){
SOL_1[i]=SOL[i]+alpha*P_0[i];}

rrtmp2=dot(R_1,R_1);

beta=rrtmp2/rrtmp;

for(i=1;i<n;i++){
P_1[i]=R_1[i]+beta*P_0[i];}

resid=Math.sqrt(rrtmp2/n);
    //System.out.println("resid is "+resid);

for(i=1;i<n;i++){
R_0[i]=R_1[i];
P_0[i]=P_1[i];
SOL[i]=SOL_1[i];}

count++;
   
}
while((resid>CG_RESIDUAL)&&(count<1000));
    //System.out.println("residual is "+resid);
    
    for(i=1;i<n;i++){
        if(coord==0){x[i]=SOL[i];}
        if(coord==1){y[i]=SOL[i];}
    }
    
    for(j=0;j<DT.boundary_length;j++){
        if(coord==0){
    x[DT.boundary[j]]=Math.cos((2*Math.PI/DT.boundary_length)*(j+1));}
        if(coord==1){
    y[DT.boundary[j]]=Math.sin((2*Math.PI/DT.boundary_length)*(j+1));}
        }
    
    //System.out.println("Output soln:");
   // for(i=0;i<n;i++){
     //   if(coord==0){
       //     System.out.print(x[i]+" ");}
    //    if(coord==1){
        //    System.out.print(y[i]+" ");}
      //  }
    
    
    
    return;
}

}
