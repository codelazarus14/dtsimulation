public class SIMPLEX {

	public int[] vertex;
	public int  sum;
	public int  label;
	public DT.LOGIC flag;
    public SIMPLEX[] neighbour;

    SIMPLEX(int a[], int N){
        int  i,count;
        int DPLUS=N+1;
        vertex = new int[DPLUS];
        neighbour = new SIMPLEX[DPLUS];
        count=0;
        
        for(i=0;i<DPLUS;i++)
        {
            //System.out.println("a[i] is "+a[i]);
            vertex[i]=a[i];
            neighbour[i]=null;
            count+=a[i];
        }
        sum=count;
        //System.out.println("count in SIMPLEX is"+count);
        
        label=DT.pointer_number;
       // System.out.println("pointer_number in SIMPLEX is "+label);
        flag=DT.LOGIC.NO;
    }
}




