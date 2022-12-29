import java.io.File;  // Import the File class
import java.io.FileNotFoundException;  // Import this class to handle errors
import java.util.Scanner; // Import the Scanner class to read text files

public class DT{

	/* some global constants for run */

	public enum LOGIC {NO, YES};
	public double kappa_0=0.0,BETA=3.0,kappa_d=0.0, kappa_0b=0.0,ALPHA=0.25;
	public int  DISKVOL=472;
	public int DISKMINVOL=DISKVOL-30;
	public int DISKMAXVOL=DISKVOL+30;
	public double FRAC=0.75;
	public double TC=7.0;
	public int MARKEDQ=282; //FRAC*(DISKVOL+2)/(2.0-FRAC)

	public static int D=2;
	public int DPLUS=D+1;
	public int DPLUSPLUS=D+2;
	public int VOL=DISKVOL+MARKEDQ;
	public int BIGVOL=4*VOL;
	public int MAXVOL=VOL+30;
	public int MINVOL=DPLUSPLUS;
	public int NUMROW=VOL/2;
	public int NONZEROS=NUMROW*VOL;

	public int THERMALISE=500;
	public int SWEEPS=20000;
	public int TUNE_COUPLING=100;
	public int SEED=1;
	public int GAP=10;
	public int DV=2;

	/* simple global pointers and counters */

	public NODE stack_head = null;
	public int simplex_number=0;

	public static int pointer_number=0;

	public static int node_number=0;
	public static int boundary_length;
	public static int[] boundary;

	public static int[] nstart, ncol;
	public static double[] nlap;

	/* data structures */

	public static SIMPLEX[] simplex_point;
	public int stack_count=0;

	/* measurements etc */

	public int number_measure=0, max_point;
	public int[] legal_subsimplex, try_subsimplex, manifold_subsimplex, go_subsimplex;
	public int vol_monitor=0, num_monitor=0, b_monitor=0, growing_vol;
	public double manifold_check;

	public LOGIC grow;

	/* simple observables -- number of nodes, simplices and size */
	/* and variable versions of couplings */

	public double real_simplex=0.0, real_node=0.0;

	/* routines checks manifold constraint */
	/* this entails examining neighbour simplices for an occurrence */
	/* of the 'opposing' vertex in a simplex which also contains the other*/
	/* new common subsimplex vertices */
	/* this is equivalent to requiring that the new common subsimplex */
	/* is not present anywhere else in the triangulation */
	/* examine all neighbours gotten by moving out on faces */
	/* which contain the other new common subsimplex vertices */

void PrintConfig(){
        int i,j;
        System.out.println("Configuration");
    System.out.println("Number of simplices and nodes "+simplex_number+" "+node_number);
        for(i=0;i<pointer_number;i++){
            System.out.println("Simplex "+i);
            if(simplex_point[i]!=null){
            for(j=0;j<DPLUS;j++){
                System.out.println("vertex "+simplex_point[i].vertex[j]);
                System.out.println("neighbour "+simplex_point[i].neighbour[j].label);}
            }
    }
    }

	public LOGIC allowed_move(SIMPLEX p, int sub, int[] a)
	{

		SIMPLEX[] array1, array2, dum, examine;
		int i,j,number1,number2,search;
		int[] b;
		LOGIC good;
		examine=new SIMPLEX[VOL];
		array1=new SIMPLEX[VOL];
		array2=new SIMPLEX[VOL];
		dum=new SIMPLEX[DPLUS];
		b=new int[DPLUS];

		/* need to teminate when either 1. see opposing vertex a[DPLUS] in some */
		/* simplex or 2. have examined all simplices containing the new common */
		/* vertices derived from original simplex */
		/* so need to flag simplices that have been examined */

		good=LOGIC.YES;

		if(subsimplex==D)
			return(good);

		if (subsimplex==0)
			return(good);

		j=0;
		for(i=0;i<DPLUS;i++)
			if(i>subsimplex)
			{b[j]=a[i];j++;}

		array1[0]=p;
		number1=1;
		examine[0]=p;
		search=1;
		p.flag=LOGIC.YES;

		/* loop while new neighbour simplices to examine */

		while(number1>0)
		{

			number2=0;

			for(i=0;i<number1;i++)
			{

				/* examine to see if contains 'opposing vertex */

				for(j=0;j<DPLUS;j++)
					if(array1[i].vertex[j]==a[DPLUS]){good=LOGIC.NO; break; }

				if(good==LOGIC.NO) break;

				/* find simplices which also in common with this subsimplex */
				/* and neighbour to simplex array1[i] */

				common_simplex(array1[i],b,D-subsimplex,dum);

				/* check to see whether any of these have been seen before */
				/* if not then add to array2 and flag seen. Also make note in examine */

				for(j=0;j<subsimplex+1;j++)
					if((dum[j].flag)==LOGIC.NO)
					{array2[number2]=dum[j];
					number2++;
					dum[j].flag=LOGIC.YES;

					examine[search]=dum[j];
					search++;
					}

			}

			if(good==LOGIC.NO) break;

			for(i=0;i<number2;i++)
				array1[i]=array2[i];

			number1=number2;

		}

		/* depending on value of good we have either found the opposing vertex or */
		/* examined all the associated simplices */

		/* first set all local flags back to zero */

		for(i=0;i<search;i++)
			examine[i].flag=LOGIC.NO;

			manifold_check+=(double)search/VOL;

			return(good);
	}

	/* routine takes n element vector a[] and returns n-1 element in b[] */
	/* thus n combinations possible selected according to count */
	/* count is the index 'left out' in getting n-1 from n */
	/* vertex left out of final vector b[] is returned in b[n-1] */

	public void combo(int[] a, int[] b, int  n, int  count)
	{
		int  i,add;

		add=0;
		for(i=0;i<n;i++)
		{
			if(count!=i)
			{
				b[add]=a[i];
				add++;
			}
			else b[n-1]=a[i];
		}
		return;
	}

	/* routine takes a named subsimplex a[0]...a[n-1] and searches for it */
	/* in a simplex pointed at by p. */
	/* returns a vector containing the local indices of pointers to */
	/* neighbouring simplices */
	/* which share a face also encompassing this subsimplex */

	public void common_simplex(SIMPLEX p, int[] a, int n, SIMPLEX[] face)
	{
        int i, j;
		int[] b, mask = new int[DPLUS];
		LOGIC[] found = new LOGIC[DPLUS];
		b=new int[DPLUS];

		for(i=0;i<n;i++)
			found[i]=LOGIC.NO;

		/* find positions/local indices of subsimplex in this simplex */

		for(i=0;i<n;i++)
			for(j=0;j<DPLUS;j++)
				if(p.vertex[j]==a[i]){b[i]=j;found[i]=LOGIC.YES;break;}

		for(i=0;i<n;i++)
			if(found[i]==LOGIC.NO){System.out.println("Error in common_simplex");System.exit(1);}

		for(j=0;j<DPLUS;j++)
			mask[j]=0;

		for(j=0;j<n;j++)
			mask[b[j]]=1;

		j=0;
		for(i=0;i<DPLUS;i++)
			if(mask[i] == 0){face[j]=p.neighbour[i];j++;}

		return;
	}

	/* routine takes pointer to simplex and a vertex and returns local */
	/* index to conjugate face */

	public int  find_face(SIMPLEX p, int  a)
	{
		int i;
//        System.out.println("looking for face opposite to "+a);
		for(i=0;i<DPLUS;i++)
			if(p.vertex[i]==a) return(i);

		System.out.println("Error in find_face ");
		System.exit(1);
		/* return dummy if get to here error */

		return(VOL);
	}

    // wrap num in single element array to pass by ref !

	public int find_order(int a, SIMPLEX pp){
		int[] dum,num;
		SIMPLEX[] list;
		list=new SIMPLEX[MAXVOL];
		dum=new int[DPLUS];
        num=new int[1];
		/* finds num of d-simps sharing vertex a */

		dum[0]=a;
        int dummy=1;
		find_simplices(pp,dum,dummy,list,num);
       //System.out.println("order of vertex "+a+" is "+num[0]);
		return(num[0]);
	}
	/* finds addresses of all simplices which share a given subsimplex */

	public void find_simplices(SIMPLEX p, int[] a, int sub,
			SIMPLEX[] s_near, int[] num) {

		SIMPLEX[] array1,array2,near;
		int i,j,k,num1,num2;
		array1=new SIMPLEX[MAXVOL];
		array2=new SIMPLEX[MAXVOL];
		near=new SIMPLEX[DPLUS];
        int number;

		array1[0]=p;
		num1=1;
		s_near[0]=p;
		number=1;
		//num = number;
		s_near[0].flag=LOGIC.YES;

		while(num1>0){
			num2=0;

			for(i=0;i<num1;i++){
				common_simplex(array1[i],a,sub,near);

				for(j=0;j<DPLUS-sub;j++)
					if(near[j].flag==LOGIC.NO){
						s_near[number]=near[j];
						array2[num2]=near[j];
						/*printf("simplex %d found with vertices:\n",near[j].label);
						 * for(k=0;k<DPLUS;k++)
						 * printf("%d ",near[j].vertex[k]);
						 * printf("\n");
						 * fflush(stdout);*/

						near[j].flag=LOGIC.YES;
						(number)++;
						num2++;
					}
			}

			for(i=0;i<num2;i++)
				array1[i]=array2[i];

			num1=num2;

		}

		for(i=0;i<(number);i++)
			s_near[i].flag=LOGIC.NO;

        num[0]=number;
        return;
	}
	/* takes simplex and order of subsimplex returns logical flag if legal move */
	/* legal move is equivalent to being exactly d+1-i associated simplices. */
	/* these subsimplex vertices occupy first i+1 entries in a[] */
	/* other d-i 'external' vertices of original simplex occupy end of a[] */
	/* also returns pointers to these d+1-i simplices */
	/* and opposing vertex is placed at end of a[] */

	public LOGIC good_subsimplex(SIMPLEX p, int sub, int[] a,
			SIMPLEX[] isimplex)
	{
		int i,add,temp, seen_already, opposing;
		int[] aind;
		aind=new int[DPLUS];

		/* test whether subsimplex is simplex itself i.e node insertion move */


		if(sub==D)
		{

			for(i=0;i<DPLUS;i++)
				a[i]=p.vertex[i];

				/* opposing vertex in this case is 'new' one obtained off stack */

				if(stack_head==null)
					a[DPLUS]=node_number;
				else
					a[DPLUS]=stack_head.name;

					isimplex[DPLUS]=p;
					return(LOGIC.YES);
		}

		/* otherwise generate subsimplex at random placing its indices in aind */

		add=0;
		while (add<subsimplex+1)
		{
			temp=(int)(Math.random() * DPLUS);
			/* generate random index */

			/* scan existing ones to see if already produced */

			seen_already=0;
			for(i=0;i<add;i++)
				if(temp==aind[i]) {seen_already=1;break;}

			if (seen_already == 0)
			{
				aind[add]=temp;
				add++;
			}
		}


		/* now create array of indices to d-i remaining vertices */

		temp=add;
		for(i=0;i<DPLUS;i++)
		{
			seen_already = 0;
			for(add=0;add<subsimplex+1;add++)
				if(i==aind[add]) {seen_already=1;break;}

			if(seen_already == 0)
			{
				aind[temp]=i;
				temp++;
			}

		}

		if(temp!=(DPLUS)){System.out.println("Error in good_subsimplex");System.exit(1);}

		/* now loop over all possible faces constucted to include this subsimplex */
		/* by selecting d-i-1 out of the d-i remaining indices */

		for(i=0;i<DPLUS;i++)
			a[i]=p.vertex[aind[i]];

        isimplex[subsimplex+1]=(p.neighbour[aind[subsimplex+1]]);
			opposing=(isimplex[subsimplex+1].sum)-sum_face(p,a[subsimplex+1]);

			/* protect following loop if sub=D-1 */

			if(sub<(D-1))
				for(i=sub+2;i<DPLUS;i++)
				{
					isimplex[i]=(p.neighbour[aind[i]]);
					if(((isimplex[i].sum)-sum_face(p,a[i]))!=opposing)
						return(LOGIC.NO);
				}

			a[DPLUS]=opposing;
			isimplex[DPLUS]=p;
			return(LOGIC.YES);

	}


	public void ReadFile() {
			try {
								File myObj = new File("config.txt");
								int  count,s_number,i,j,k,dummy,l_number,c,e1,e2,simp;
								int[] dum;
								int[][] dum2;
								dum=new int[DPLUS];

								double k0,b;
								double temp;

								Scanner myReader = new Scanner(myObj);

								s_number= myReader.nextInt();
								dum2=new int[s_number][DPLUS];
								node_number=myReader.nextInt();
								count=myReader.nextInt();
								temp=myReader.nextFloat();
								temp=myReader.nextFloat();
								temp=myReader.nextFloat();

								System.out.println(s_number+" "+node_number+" "+count);



								simplex_number=0;
								stack_count=0;
								pointer_number=0;
								System.out.println("Reading in existing configuration\n");


								System.out.println("\nConfiguration has volume=" +s_number);
								System.out.println("\nNode number="+node_number);
								System.out.println("\n  k0b coupling, k2 coupling and curvaturesq coupling ="+kappa_0b+" "+kappa_d +" "+BETA);


								for(i=0;i<count;i++)
								{
									dummy=myReader.nextInt();
									push(dummy);
								}


								for(i=0;i<s_number;i++)
								{
										for(j=0;j<DPLUS;j++){
											dum[j]=myReader.nextInt();
											dum2[i][j]=myReader.nextInt();
										}

										simplex_point[i]=new SIMPLEX(dum,D);
										simplex_number++;
										pointer_number++;

								}

								for(i=0;i<s_number;i++)
								{
										for(j=0;j<DPLUS;j++){
											simplex_point[i].neighbour[j]=simplex_point[dum2[i][j]];
										}

								}


								for(i=0;i<s_number;i++)
								{
									//System.out.println(simplex_point[i].label+ ":\t" );
									for(j=0;j<DPLUS;j++){
										System.out.println(simplex_point[i].vertex[j]+ "\t" +simplex_point[i].neighbour[j].label+ "\t" );
									}
										System.out.println("\n");
								}



								//for(i=0;i<VOL;i++)
								//fscanf(fp1,"%d",&localvol[i]);

								growing_vol=simplex_number;


								System.out.printf("\nHave read data successfully\n ");



								myReader.close();

			} catch (FileNotFoundException e) {
				System.out.println("An error occurred.");
				e.printStackTrace();
			}
			return;
	}



    void thermalize(){
    int therm;

		System.out.println("Let's run reading files");
		simplex_point=new SIMPLEX[BIGVOL];

		ReadFile();


    nlap=new double[NONZEROS];
    ncol=new int[NONZEROS];
    nstart=new int[VOL];
    boundary=new int[VOL];

    //System.out.println("About to call initial_config");
		//initial_config();
        //PrintConfig();

		legal_subsimplex=new int[DPLUS];
		try_subsimplex=new int[DPLUS];
		manifold_subsimplex=new int[DPLUS];
		go_subsimplex=new int[DPLUS];

		/* build lattice */
    //    System.out.println("About to grow lattice");

		grow=LOGIC.NO;
		if(grow==LOGIC.YES){
				while(growing_vol<VOL) {
		//			System.out.println(growing_vol);
					trial_change();
					growing_vol+=D;
		            //PrintConfig();
				}
		}
        tidy();

		grow=LOGIC.NO;

		/* thermalise and output info on run */

		header();

		System.out.println("Thermalizing lattice");

		for(therm=1;therm<=THERMALISE;therm++)
		{
			for(int iter=0;iter<VOL;iter++){
				trial_change();
                //tidy();
			}
            tidy();
            //PrintConfig();
						// find pointer to node 0
					 int done=0,v;
					 //System.out.println("simplex number= "+simplex_number);
					 SIMPLEX p=null;
					 int[] nn=null;
					 int[] num=null;
					 num=new int[1];
					 nn=new int[VOL];

					 for(int i=0;i<simplex_number;i++){
					 if(done==1) break;
					 for(int j=0;j<(DT.D+1);j++){
					 if(simplex_point[i].vertex[j]==0){p=simplex_point[i];done=1;break;}
					 }}
					 v=0;

					 getnn(p,v,nn,num);
					 //boundary_length=num[0];


			num_monitor++;
			vol_monitor+=simplex_number;
            b_monitor+= num[0];

			if(therm%TUNE_COUPLING==0){
					System.out.println("V: " + (simplex_number-num[0])  );
					System.out.println("L: " + num[0]);
					System.out.println("N0B: "+(node_number-boundary_length-1));

				shift_coupling();
			}


		}

		init();

		return;
}


// hardwired for D=2 right now ..

public void getnn(SIMPLEX p, int v, int[] nn, int[] vnum){
    int[] dum,seen,num,v1,v2;
    int i,j,k,currentpt,nextpt,working,index;
    SIMPLEX[] list;

    list=new SIMPLEX[BIGVOL];
    dum=new int[DPLUS];
    seen=new int[VOL];
    num=new int[1];
    v1=new int[VOL];
    v2=new int[VOL];

    //System.out.println("in getnn");
    //tidy();
    //relabelnodes();


    dum[0]=v;
    int dummy=1;
    //System.out.println("simps around node");
    find_simplices(p,dum,dummy,list,num);
   // for(i=0;i<num[0];i++){
    //System.out.println("tri "
    //    +list[i].vertex[0]+" "
     //   +list[i].vertex[1]+" "
      //  +list[i].vertex[2]);
   // }

    for(i=0;i<VOL;i++){
    seen[i]=0;}

    k=0;
    for(i=0;i<num[0];i++){
    index=0;

    for(j=0;j<DPLUS;j++){
    if(list[i].vertex[j]==v){
    index=j;}
    }

    v1[k]=list[i].vertex[(index+1)%DPLUS];
    v2[k]=list[i].vertex[(index+2)%DPLUS];
    k++;
    }

    //for(i=0;i<num[0];i++){
    //System.out.println("v1,v2 "+v1[i]+" "+v2[i]);
    //}

    nn[0]=v1[0];
    //boundary[1]=v2[0];
    seen[0]=1;
    k=1;

    // look for v1[0] in rest of v1/v2 arrays
    currentpt=v1[0];
    do{
    nextpt=0;

    for(i=0;i<num[0];i++){
    if(seen[i]==0){
    //System.out.println("examining next triangle with i= "+i);
    if(v1[i]==currentpt){
    nn[k]=v2[i];
    nextpt=v2[i];k++;seen[i]=1;}
    if(v2[i]==currentpt){
    nn[k]=v1[i];
    nextpt=v1[i];k++;seen[i]=1;}
    }

    }
    //System.out.println("nextpt is "+nextpt);
    currentpt=nextpt;
    }while(k<num[0]);


    vnum[0]=num[0];

    //System.out.println("ordered verts");
    //for(i=0;i<vnum[0];i++){
    //System.out.println("neighbor vertices "+nn[i]);
    //}
return;
}

public void laplacian(){
int i,j,v,k,done;
int[] num_in_row, col, start, not_seen,nn,num,q,on_boundary;
double[] lap;
SIMPLEX p=null;
num=new int[1];
nn=new int[VOL];
not_seen=new int[VOL];

//System.out.println("In laplacian");

 // find pointer to node 0

done=0;
for(i=0;i<simplex_number;i++){
if(done==1) break;
for(j=0;j<(DT.D+1);j++){
if(simplex_point[i].vertex[j]==0){p=simplex_point[i];done=1;break;}
}}
v=0;

getnn(p,v,nn,num);
boundary_length=num[0];
System.out.println("Boundary length "+boundary_length);
System.out.println("Total bulk node "+(node_number-boundary_length-1));

//System.out.println("Boundary vertices: ");
for(i=0;i<boundary_length;i++){
boundary[i]=nn[i];
//System.out.println(boundary[i]);
}

num_in_row=new int[VOL];
col=new int[NONZEROS];
start=new int[VOL];
not_seen=new int[VOL];
q=new int[VOL];
nn=new int[VOL];
lap=new double[NONZEROS];
on_boundary=new int[VOL];
num=new int[1];

for(i=0;i<VOL;i++){
not_seen[i]=1;}

for(i=0;i<VOL;i++){
on_boundary[i]=0;}

not_seen[0]=0;
//on_boundary[0]=1;
//for(i=0;i<boundary_length;i++){
//not_seen[boundary[i]]=0;
//on_boundary[boundary[i]]=1;
//System.out.println("boundary pt "+boundary[i]);
//}

    for(k=0;k<VOL;k++){
    start[k]=k*NUMROW;
    num_in_row[k]=0;
    q[k]=0;
    }
    for(k=0;k<NONZEROS;k++){
        col[k]=-1;
    }

int maxv=0;int minv=100000;
int vcount=0;
for(i=0;i<simplex_number;i++){
for(j=0;j<DPLUS;j++){
v=simplex_point[i].vertex[j];

if(v>maxv){maxv=v;}
if(v<minv){minv=v;}

//printf("vertex %d\n",v);
// leave out marked node

if(not_seen[v]==1){
vcount++;
// new vertex
//printf("new vertex\n");
    //System.out.println("vertex "+v);
getnn(simplex_point[i],v,nn,num);
    for(k=0;k<num[0];k++){
       // System.out.println("neighbor is "+nn[k]);
        if(nn[k]!=0){
            lap[start[v]+num_in_row[v]]=-1.0;
            col[start[v]+num_in_row[v]]=nn[k];
            num_in_row[v]++;
        }
    }

    // diagonal piece
    lap[start[v]+num_in_row[v]]=num_in_row[v];
    col[start[v]+num_in_row[v]]=v;
    num_in_row[v]++;
    if(num_in_row[v]>NUMROW){System.out.println("oops - dimension sparse row too large at v="+v);}

    not_seen[v]=0;
  }
}
}
//printf("node_number is %d\n",node_number);
//printf("vcount is %d\n",vcount);
//printf("maxv %d\n",maxv);
//printf("minv %d\n",minv);
if((maxv!=node_number-1)||(minv!=0)){System.out.println("relabeled triangulation incorrect with minv="+minv+"and maxv="+maxv);}
int c=0;
    for(i=1;i<=maxv;i++){
    if(num_in_row[i]==0) {System.out.println("missing vertex"+i);}
    if(i>(node_number-1)){System.out.println("super minimal vertex label"+i);}
    if(num_in_row[i]!=0){
        c=c+num_in_row[i]-1;}
}
//printf("sum of neighbor coord %d\n",c);
//printf("number of simps %d\n",simplex_number);

// clean up sparse  data structures
k=0;
nstart[1]=0;
nstart[0]=0;
for(i=0;i<node_number;i++){
for(j=start[i];j<start[i+1];j++){
if(col[j]!=(-1)){
ncol[k]=col[j];
nlap[k]=lap[j];
k++;}
}
nstart[i+1]=nstart[i]+num_in_row[i];
}

//for(i=0;i<node_number;i++){
//System.out.println("row "+i+" start of row "+nstart[i]);
//for(j=nstart[i];j<nstart[i+1];j++){
//System.out.println("col is "+ncol[j]+" and laplacian is "+nlap[j]);}}



return;
}










    public void deletestack(){
        NODE tmp;
        NODE[] dum;
        int num=0;

        dum=new NODE[VOL];
        tmp=stack_head;
        while(tmp!=null){
            dum[num]=tmp;tmp=tmp.next;num++;
        }
        //printf("in deletestack()\n");
        //printf("num is %d\n",num);
        //printf("stack_count is %d\n",stack_count);


        stack_head=null;
        stack_count=0;
    }


    public void relabelnodes(){
        int[] not_seen,dum;
        int i,j,k,l,new2,old,v,VERYBIG=100000;
        SIMPLEX[] list;
        SIMPLEX address;
        int[] num;

        num=new int[1];

        not_seen=new int[VOL];
        list=new SIMPLEX[VOL];
        dum=new int[DPLUS];


        //System.out.println("in relabelnodes");
        //System.out.println("node number "+node_number);

        // printf("in relabelnodes\n");
        // printf("node number is %d\n",node_number);


        for(i=0;i<VOL;i++){
            not_seen[i]=1;}

        new2=VERYBIG+1;
        for(i=0;i<simplex_number;i++){
            for(j=0;j<DPLUS;j++){
                old=simplex_point[i].vertex[j];
                //printf("looking for old node %d\n",old);
                if((old==0) || (old>VERYBIG)) continue;
                //if(old>VERYBIG) continue;
                if(not_seen[old]==1){
                    //printf("found old %d ",old);

                    dum[0]=old;
                    int dummy=1;
            find_simplices(simplex_point[i],dum,dummy,list,num);
                    for(k=0;k<num[0];k++){
                        for(l=0;l<DPLUS;l++){
                            if(list[k].vertex[l]==old)
                            {
                            list[k].vertex[l]=new2;
                            }
                        }
                    }
                    not_seen[old]=0;
                    //printf("new label is %d\n",new);
                    new2++;
                }
            }
        }
        //System.out.println("number of nodes after relabel "+new2%VERYBIG);
        // reset labels
        int a,add;
        for(i=0;i<simplex_number;i++){
            add=0;
            for(j=0;j<DPLUS;j++){
                a=simplex_point[i].vertex[j];
                //printf("a initially is %d\n",a);
                a=a%VERYBIG;
                // printf("a finally is %d\n",a);
                simplex_point[i].vertex[j]=a;
                add+=a;
            }
            simplex_point[i].sum=add;
            // printf("sum of vertex labels %d\n",add);
        }

        //printf("%d final distinct vertices\n",newcount);
        // delete old node stack ...

        deletestack();

        return;
    }
	/* driver for simplex Monte Carlo */




	void main()
	{
		int iter,sweep;


	 	//ReadFile();
		thermalize();
		/* sweep lattice outputting measurements every so often */

		System.out.println("Starting measurements");
    System.out.println("Starting measurements");

		System.out.println("Why the printing is not working?");

		for(sweep=1;sweep<=SWEEPS;sweep++)
		{
			System.out.println("sweep="+sweep);

			for(iter=0;iter<VOL;iter++){
				trial_change();
			}

            //PrintConfig();
      tidy();
			int done=0,v;
			System.out.println("simplex number= "+simplex_number);
			SIMPLEX p=null;
			int[] nn;
			int[] num;
			num=new int[1];
			nn=new int[VOL];

			for(int i=0;i<simplex_number;i++){
			if(done==1) break;
			for(int j=0;j<(DT.D+1);j++){
			if(simplex_point[i].vertex[j]==0){p=simplex_point[i];done=1;break;}
			}}
			v=0;

			getnn(p,v,nn,num);

			boundary_length=num[0];


 		num_monitor++;
 		vol_monitor+=simplex_number;
			 b_monitor+=boundary_length;

      if(sweep%100==0){
				shift_coupling();
				System.out.println("Boundaryyyy length "+boundary_length);
				System.out.println("Total bulkyyy node "+(node_number-boundary_length-1));
			}


        //if(sweep%GAP==0){
				//measure();
			//}

			/* checkpoint config regularly */


		}

		/* finally dump final configuration and print some results */

		print_out();

		return;

	}
	/* opens files, initialises variables */
	/* prints out the run parameters */

	public void header()
	{

		max_point=0;

		System.out.println("Dimension " + D);
		System.out.println("Volume "+ VOL);
    System.out.println("Marked node coupling " + ALPHA);
		System.out.println("Simplex coupling "+ kappa_d);
		System.out.println("Bulk coupling " + BETA);
		System.out.println("Number of sweeps " + SWEEPS);
		System.out.println("Thermalisation time " + THERMALISE);
		System.out.println("Number of sweeps between KD tuning " + TUNE_COUPLING);
		System.out.println("Gap between measurements " + GAP);
		System.out.println("Volume fluctuation parameter " + DV);
		System.out.println("Random number seed " + SEED);

		return;
	}
	/* opens all the data files and zeroes measure bins */

	public void init()
	{
		int  i;

		for(i=0;i<DPLUS;i++)
		{
			try_subsimplex[i]=0;
			go_subsimplex[i]=0;
			manifold_subsimplex[i]=0;
			legal_subsimplex[i]=0;
		}

		real_simplex=0.0;
		real_node=0.0;
		manifold_check=0.0;
		max_point=0;

		return;
	}

	void initial_config()
	{
		int i,j,index,tmp;
		int[] dum,dum2;
		dum=new int[DPLUSPLUS];
		dum2=new int[DPLUSPLUS];

		for(i=0;i<DPLUSPLUS;i++)
			dum[i]=i;

		for(i=0;i<DPLUSPLUS;i++)
		{
			combo(dum,dum2,DPLUSPLUS,i);
           // for(j=0;j<DPLUS;j++){
               // System.out.println("dum2[j] is "+dum2[j]);}

			simplex_point[pointer_number]=new SIMPLEX(dum2,D);
           // System.out.println("vertices of simplex are ");
          //  for(j=0;j<DPLUS;j++){
          //  System.out.println(simplex_point[pointer_number].vertex[j]);}


			simplex_number++;
			pointer_number++;
		}

		/* now set up pointers loop over faces to simplices */

		for(i=0;i<DPLUSPLUS;i++)
			for(j=0;j<DPLUSPLUS;j++)
				if(j!=i)
				{
                index=find_face(simplex_point[i],dum[j]);
                simplex_point[i].neighbour[index]=simplex_point[j];
				}

		node_number=DPLUSPLUS;
		growing_vol=DPLUSPLUS;
       // System.out.println("finished in initial_config");
		return;
	}
	/* routine coordinates all measurements */

	int calls = 0;
	public void measure()
	{
		double dum,total_action;

		real_simplex+=(double)simplex_number;
		real_node+=(double)node_number;


		tidy();

		number_measure++;
		calls++;


		return;
	}


 public double DS(int q, int deltaq, double qmean){
// computes change in local action. Assumes (q-qmean)^2 form
    return(deltaq*(deltaq+2*q-2*qmean));

}
public LOGIC metropolis2(int subsimplex, int[] a,
			SIMPLEX[] addresses) {
		double ds0=0.0,dsd=0.0,dsn=0.0,dummy,change;
		LOGIC accept;
		int i,j,k,dum,notnnmarked,imark;
        int[] vnum,nn,notbound;
        SIMPLEX p=null;
		int tmp;
        vnum=new int[1];
        nn=new int[VOL];
        notbound=new int[DPLUSPLUS];


// find marked node
        int done=0;
        for(i=0;i<pointer_number;i++){
            if(simplex_point[i]==null) continue;

        if(done==1) break;
        for(j=0;j<DPLUS;j++){
        if(simplex_point[i].vertex[j]==0) {p=simplex_point[i];done=1;break;}
        }}

//System.out.println("found marked pt");

// find vertices on boundary
        dum=0;
        getnn(p,dum,nn,vnum);
    // System.out.println("boundary length in metro is "+vnum[0]);


// flag them
        for(i=0;i<DPLUSPLUS;i++){
        notbound[i]=1;
        for(j=0;j<vnum[0];j++){
        if(a[i]==nn[j]) {notbound[i]=0;}
        }
        }
dsn=0.0;

//System.out.println("Number of simplices and bulk nodes "+(simplex_number-vnum[0])+" "+(node_number-vnum[0]) );
////////////////////// Node Insertion ////

if(subsimplex==D){
    // check if inserting into boundary triangle
        notnnmarked=1;
        for(i=0;i<DPLUS;i++){
            if(a[i]==0){notnnmarked=0;}
        }

				if(notnnmarked==1) {
		      dsd = kappa_d*(2*subsimplex-D);
		      dsd=dsd+(2*subsimplex-D)*((2*subsimplex-D)+2*(simplex_number-vnum[0]-DISKVOL) )  /(1.0*DV*DV) ;
		    }
		    else{
		      dsd = kappa_d; ds0=-kappa_0b;
		      dsd=dsd+(1+2*(simplex_number-vnum[0]-DISKVOL) )  /(1.0*DV*DV) ;
		    } // additional minus sign in ds0 from formula

    /* node insertion */
        for(i=0;i<DPLUS;i++){
            // tmp=(double)localvol[a[i]];
            tmp=find_order(a[i],addresses[DPLUS]);
            if(a[i]==0){ //local vol of the MN is length of the boundary
                dsn=dsn+ALPHA*DS(tmp,D-1,MARKEDQ); // kappa_0 earlier
              }
            else{
                dsn=dsn+BETA*DS(tmp,D-1,TC)*notbound[i]; // this only added for bulk points
              }
        }

        //update for the inserted node
        dsn=dsn+BETA*(DPLUS-TC)*(DPLUS-TC)*notnnmarked;
				if (notnnmarked==0 && (simplex_number-vnum[0]+1)>DISKMAXVOL){accept=LOGIC.NO; return accept;}
		    else if (notnnmarked==1 && (simplex_number-vnum[0]+2)>DISKMAXVOL){accept=LOGIC.NO; return accept;}
    }

    else if(subsimplex==0){
        notnnmarked=1;
        for(i=1;i<DPLUSPLUS;i++){
            if(a[i]==0){notnnmarked=0;}
        }

				if(notnnmarked==1) {
		      dsd = kappa_d*(2*subsimplex-D);
		      dsd=dsd+(2*subsimplex-D)*((2*subsimplex-D)+2*(simplex_number-vnum[0]-DISKVOL) )  /(1.0*DV*DV) ;
		    }
		    else  {
		      dsd = -kappa_d; ds0= kappa_0b;
		      dsd=dsd-(-1+2*(simplex_number-vnum[0]-DISKVOL) )  /(1.0*DV*DV) ; //delta N_D=-1
		    } // additional minus sign in ds0 from formula


    // check if deleting a boundary node
    /* node deletion */
    for(i=1;i<DPLUS;i++){
        tmp=find_order(a[i],addresses[DPLUS]);
        if(a[i]==0){
            dsn=dsn+ALPHA*DS(tmp,-D+1,MARKEDQ);}
        else{
            dsn=dsn+BETA*DS(tmp,-D+1,TC)*notbound[i];}
         }

				// should be deleted next block, right??

        /*tmp=find_order(a[DPLUS],addresses[subsimplex+1]);

				if(a[DPLUS]==0){
            dsn=ALPHA*DS(tmp,-D+1,MARKEDQ);} // change needed ????
        else{
            dsn=dsn+BETA*DS(tmp,-D+1,TC)*notbound[i];} */

    dsn=dsn-BETA*(DPLUS-TC)*(DPLUS-TC)*notnnmarked;
		if (notnnmarked==0 && (simplex_number-vnum[0]-1)<=DISKMINVOL ){accept=LOGIC.NO; return accept;}
    if (notnnmarked==1 && (simplex_number-vnum[0]-2)<=DISKMINVOL ){accept=LOGIC.NO; return accept;}
    }

    //////// link flip ///////////////
    else{

      int markedNodeModes=0;
      for(i=0;i<=subsimplex;i++){
          if(a[i]==0){markedNodeModes=1;}
      }
      for(i=subsimplex+1;i<DPLUSPLUS;i++){
          if(a[i]==0){markedNodeModes=2;}
      }

			if(markedNodeModes==1) {
		    dsd = kappa_d; ds0= kappa_0b;
		    dsd=dsd+(1+2*(simplex_number-vnum[0]-DISKVOL) )  /(1.0*DV*DV) ;
		  }  // additional minus sign from formula
		  else  if(markedNodeModes==2) {
		    dsd = -kappa_d; ds0=-kappa_0b ;
		     dsd=dsd-(-1+2*(simplex_number-vnum[0]-DISKVOL) )  /(1.0*DV*DV) ;
		  }// additional minus sign from formula

    for(i=0;i<=subsimplex;i++){
        tmp=find_order(a[i],addresses[DPLUS]);
        if(a[i]==0){
            dsn=dsn+ALPHA*DS(tmp,2*subsimplex-D-1,MARKEDQ);}
        else{
            if(markedNodeModes==0){
              dsn=dsn+BETA*DS(tmp,2*subsimplex-D-1,TC)*notbound[i];
            }
            else if(markedNodeModes==1){
              //dsn=dsn+BETA*DS(tmp,2*subsimplex-D-1,TC);
							dsn=dsn+BETA*(tmp-1-TC)*(tmp-1-TC);
            }
            }
    }
    for(i=subsimplex+1;i<DPLUS;i++){
        tmp=find_order(a[i],addresses[DPLUS]);
        if(a[i]==0){
            dsn=dsn+ALPHA*DS(tmp,2*subsimplex-D+1,MARKEDQ);}
            else{
                if(markedNodeModes==0){
                    dsn=dsn+BETA*DS(tmp,2*subsimplex-D+1,TC)*notbound[i];
                  }
                else if(markedNodeModes==2){
                  //dsn=dsn+BETA*DS(tmp,2*subsimplex-D+1,TC);
									dsn=dsn-BETA*(tmp-TC)*(tmp-TC);
                }
          }

    }
		/* This should be commented right??
    tmp=find_order(a[DPLUS],addresses[subsimplex+1]);
        if(a[i]==0){
          dsn=dsn+ALPHA*DS(tmp,2*subsimplex-D+1,MARKEDQ)*notbound[i];}
          else{
              if(markedNodeModes==0){
                  dsn=dsn+BETA*DS(tmp,2*subsimplex-D+1,TC)*notbound[i];
                }
              else if(markedNodeModes==2){
                dsn=dsn+BETA*DS(tmp,2*subsimplex-D+1,TC);
              }
        }
				*/

			if(markedNodeModes==1 && (simplex_number-vnum[0]+1)> DISKMAXVOL){accept=LOGIC.NO; return accept;}
			if(markedNodeModes==2 && (simplex_number-vnum[0]-1) <=DISKMINVOL ){accept=LOGIC.NO; return accept;} //don't think this condition is possible

}


if((subsimplex==0)&&(a[0]==0)){System.out.println("trying to delete marked node\n");dsd=1000000.0;}

dsd=dsd+dsn+ds0;

dsd=Math.exp(dsd);

/*
int markedNodeModes=0; // 0 is bulk mode
for(i=0;i<=subsimplex;i++){
    if(a[i]==0){markedNodeModes=1;}
}
for(i=subsimplex+1;i<DPLUSPLUS;i++){
    if(a[i]==0){markedNodeModes=2;}
}


if(markedNodeModes==0) {dsd=(1.0+(double)(2*subsimplex-D)/(double)(simplex_number-vnum[0]) )*dsd;}
else if(subsimplex==2){dsd=(1.0+1.0/(double)(simplex_number-vnum[0]) )*dsd;}
else if(subsimplex==0){dsd=(1.0-1.0/(double)(simplex_number-vnum[0]) )*dsd;}
else if(subsimplex==1 && markedNodeModes==1){dsd=(1.0+1.0/(double)(simplex_number-vnum[0]) )*dsd;}
else if(subsimplex==1 && markedNodeModes==2){dsd=(1.0-1.0/(double)(simplex_number-vnum[0]) )*dsd;}
*/
//else {printf("some condition not considerd, sub=%d, markedNodeModes=%d",sub,markedNodeModes);}



dsd=(1.0+(double)(2*subsimplex-D)/(double)(simplex_number))*dsd;

dsd=1.0/(1.0+dsd);

dummy=Math.random()-dsd;

if(dummy<0.0)
    accept=LOGIC.YES;

else
    accept=LOGIC.NO;

if(grow==LOGIC.YES)
    accept=LOGIC.YES;

//if(((simplex_number+2*subsimplex-D)<MINVOL) || ((simplex_number+2*subsimplex-D)>MAXVOL)) accept=LOGIC.NO;
		return(accept);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
public LOGIC metropolis(int subsimplex, int[] a,
			SIMPLEX[] addresses) {
		double ds0=0.0,dsd=0.0,dsn=0.0,dummy,change;
		LOGIC accept;
		int i,j,k,dum,notnnmarked,imark;
        int[] vnum,nn,notbound;
        SIMPLEX p=null;
		int tmp;
        vnum=new int[1];
        nn=new int[VOL];
        notbound=new int[DPLUSPLUS];


// find marked node
        int done=0;
        for(i=0;i<pointer_number;i++){
            if(simplex_point[i]==null) continue;

        if(done==1) break;
        for(j=0;j<DPLUS;j++){
        if(simplex_point[i].vertex[j]==0) {p=simplex_point[i];done=1;break;}
        }}

//System.out.println("found marked pt");

// find vertices on boundary
        dum=0;
        getnn(p,dum,nn,vnum);
    // System.out.println("boundary length in metro is "+vnum[0]);


// flag them
        for(i=0;i<DPLUSPLUS;i++){
        notbound[i]=1;
	        for(j=0;j<vnum[0];j++){
	        	if(a[i]==nn[j]) {
							notbound[i]=0;accept=LOGIC.NO;
				  		return accept;
						}
	        }
        }
dsn=0.0;
notnnmarked=1;
for(i=0;i<DPLUS;i++){
		if(a[i]==0){notnnmarked=0;}
}

ds0=0.0; // delta N0_boundary=0
dsd = kappa_d*(2*subsimplex-D);
dsd=dsd+(2*subsimplex-D)*((2*subsimplex-D)+2*(simplex_number-vnum[0]-DISKVOL) )  /(1.0*DV*DV) ;


//System.out.println("Number of simplices and bulk nodes "+(simplex_number-vnum[0])+" "+(node_number-vnum[0]) );
////////////////////// Node Insertion ////

if(subsimplex==D){


    /* node insertion */
        for(i=0;i<DPLUS;i++){
            // tmp=(double)localvol[a[i]];
            tmp=find_order(a[i],addresses[DPLUS]);
            dsn=dsn+BETA*DS(tmp,D-1,TC)*notbound[i]; // this only added for bulk points
        }

        //update for the inserted node
        dsn=dsn+BETA*(DPLUS-TC)*(DPLUS-TC)*notnnmarked;

    }

    else if(subsimplex==0 && a[0]!=0){

    	/* node deletion */
    	for(i=1;i<DPLUS;i++){
        	tmp=find_order(a[i],addresses[DPLUS]);
        	dsn=dsn+BETA*DS(tmp,-D+1,TC)*notbound[i];
			}

			dsn=dsn-BETA*(DPLUS-TC)*(DPLUS-TC)*notnnmarked;
		}

    //////// link flip ///////////////
    else{


    for(i=0;i<=subsimplex;i++){
        tmp=find_order(a[i],addresses[DPLUS]);
        dsn=dsn+BETA*DS(tmp,2*subsimplex-D-1,TC)*notbound[i];
    }
    for(i=subsimplex+1;i<DPLUS;i++){
        tmp=find_order(a[i],addresses[DPLUS]);
        dsn=dsn+BETA*DS(tmp,2*subsimplex-D+1,TC)*notbound[i];
    }


}


if((subsimplex==0)&&(a[0]==0)){System.out.println("trying to delete marked node\n");dsd=1000000.0;}

dsd=dsd+dsn+ds0;

dsd=Math.exp(dsd);

/*
int markedNodeModes=0; // 0 is bulk mode
for(i=0;i<=subsimplex;i++){
    if(a[i]==0){markedNodeModes=1;}
}
for(i=subsimplex+1;i<DPLUSPLUS;i++){
    if(a[i]==0){markedNodeModes=2;}
}


if(markedNodeModes==0) {dsd=(1.0+(double)(2*subsimplex-D)/(double)(simplex_number-vnum[0]) )*dsd;}
else if(subsimplex==2){dsd=(1.0+1.0/(double)(simplex_number-vnum[0]) )*dsd;}
else if(subsimplex==0){dsd=(1.0-1.0/(double)(simplex_number-vnum[0]) )*dsd;}
else if(subsimplex==1 && markedNodeModes==1){dsd=(1.0+1.0/(double)(simplex_number-vnum[0]) )*dsd;}
else if(subsimplex==1 && markedNodeModes==2){dsd=(1.0-1.0/(double)(simplex_number-vnum[0]) )*dsd;}
*/
//else {printf("some condition not considerd, sub=%d, markedNodeModes=%d",sub,markedNodeModes);}



dsd=(1.0+(double)(2*subsimplex-D)/(double)(simplex_number))*dsd;

dsd=1.0/(1.0+dsd);

dummy=Math.random()-dsd;

if(dummy<0.0)
    accept=LOGIC.YES;

else
    accept=LOGIC.NO;

if(grow==LOGIC.YES)
    accept=LOGIC.YES;

if(((simplex_number+2*subsimplex-D)<MINVOL) || ((simplex_number+2*subsimplex-D)>MAXVOL)) accept=LOGIC.NO;

return(accept);

}


//////////////////*****************************************************************///////////////////////////////

	public void pop()
	{

		NODE temp;

		/* if stack is empty do nothing */

		if(stack_head==null)
			return;

		temp=stack_head;
		stack_head=stack_head.next;
		stack_count--;

		return;

	}
	/* writes some results to output stream */

	public void print_out()
	{
		int  i;
		double dummy;

		System.out.println("Results:");

		real_simplex=real_simplex/number_measure;
		real_node=real_node/number_measure;

		dummy=0.0;
		for(i=0;i<DPLUS;i++)
			dummy+=(double)legal_subsimplex[i];
		//printf("Total number of legal moves %lg\",dummy);
		//printf("Raw manifold check %lg\n",manifold_check);
		manifold_check/=dummy;
		manifold_check*=VOL;

		double tmp=(real_simplex-VOL)*2.0/(DV*DV);
		System.out.println("Average number of simplices " + real_simplex);
		System.out.println("Average number of nodes " + real_node);
		System.out.println("Average maximum pointer number " + max_point);
		System.out.println("Average number of simplices in manifold check " + manifold_check);
		System.out.println("Final kappa_d " + (kappa_d + tmp));

		System.out.println("Subsimplex Moves ");
		for(i=0;i<DPLUS;i++)
		{
			System.out.println(i + " subsimplices : ");
			System.out.println("Number tried " + try_subsimplex[i]);
			System.out.println("Number that are legal " + legal_subsimplex[i]);
			System.out.println("Number that pass manifold test " + manifold_subsimplex[i]);
			System.out.println("Number that pass Metropolis test " + go_subsimplex[i]);
		}

		return;
	}
	public void print_simplex_info(SIMPLEX p){
		int i;
		System.out.println("simplex label is " + p.label);
		System.out.println("vertices are");
		for(i=0;i<DPLUS;i++)
			System.out.print(p.vertex[i] + " ");
		System.out.println();
		System.out.println("neighbor labels are");
		for(i=0;i<DPLUS;i++)
			System.out.print((p.neighbour[i]).label + " ");
		System.out.println();
		return;
	}
	/* routine pushes deleted vertex label onto stack */

	void push(int  i)
	{
		NODE temp;

		//temp=(PNODE)malloc(sizeof(NODE));
		temp=new NODE();

		temp.name=i;
		temp.next=stack_head;

		stack_head=temp;

		stack_count++;
		return;
	}

	static final double MBIG = 1000000000;
	static final double MSEED = 161803398;
	static final int MZ = 0;
	static final double FAC = (1.0/MBIG);
	int inext,inextp;
	long[] ma = new long[56];
	double ran3(Integer idum)
	{
		long mj,mk;
		int i,ii,k;
		if(idum < 0) {
			mj = (long) (MSEED - (idum < 0 ? -idum : idum));
			mj %= MBIG;
			ma[55] = mj;
			mk = 1;
			for(i=1;i<=54;i++) {
				ii = (21*i) % 55;
				ma[ii] = mk;
				mk = mj - mk;
				if(mk < MZ) mk += MBIG;
				mj = ma[ii];
			}
			for(k=1;k<=4;k++)
				for(i=1;i<=55;i++) {
					ma[i] -= ma[1+(i+30) % 55];
					if(ma[i] < MZ) ma[i] += MBIG;
				}
			inext = 0;
			inextp = 31;
			idum = 1;
		}
		if(++inext == 56) inext=1;
		if(++inextp == 56) inextp=1;
		mj = ma[inext]-ma[inextp];
		if(mj < MZ) mj += MBIG;
		ma[inext] = mj;
		return mj*FAC;
	}

	/* routine handles reconnection of new simplex pointers */

	void reconnect(int[] a, SIMPLEX[] q, int sub)
	{

		int  i,j,index_face1,index_face2,index_face3,opp;

		/* loop over final state simplices */


		for(i=0;i<sub+1;i++)
		{

			/* 'internal' faces first */


			for(j=0;j<sub+1;j++)
				if(j!=i)
				{
					index_face1=find_face(q[i],a[j]);
					q[i].neighbour[index_face1]=q[j];
				}

			/* now 'external' faces */


			for(j=sub+1;j<DPLUSPLUS;j++)
			{
				index_face1=find_face(q[i],a[j]);
				index_face2=find_face(q[j],a[i]);

				/* have found external simplices involved reconnect outward pointers */

				q[i].neighbour[index_face1]=q[j].neighbour[index_face2];

				/* just adjust pointers on external simplices so they point at new ones */

				opp=(q[i].neighbour[index_face1]).sum-
						sum_face(q[i],q[i].vertex[index_face1]);

				index_face3=find_face((q[i].neighbour[index_face1]),opp);
				(q[i].neighbour[index_face1]).neighbour[index_face3]=q[i];

			}

		}
		return;
	}
	/* selects a simplex at random by accessing an array of pointers */
	/* once has a simplex select subsimplex/move at random */

	SIMPLEX select_simplex(int sub)
	{
		SIMPLEX temp;
		int i;

		do{
			i=(int)(Math.random()*pointer_number);
			temp=simplex_point[i];
		}
		while(temp==null) ;
		//System.out.println("In select simplex");

		/* initially just grow with node insertion moves */

		if(grow==LOGIC.YES)
			subsimplex=D;
		else
			/* in general choose move type at random */
			subsimplex=(int)(Math.random()*DPLUS);

        //System.out.println("subsimplex in select_simplex is "+subsimplex);

		try_subsimplex[subsimplex]++;

		return(temp);

	}
	/* routine tunes simplex coupling to achieve preset average volume */

	public void shift_coupling()
	{
		double dum,dum2;
		dum=((double)(vol_monitor-b_monitor))/((double)num_monitor);
        dum2=((double)b_monitor)/((double)num_monitor);

		kappa_d=kappa_d+(2.0/(DV*DV))*(dum-DISKVOL)+2.0*ALPHA*(dum2-MARKEDQ)*FRAC/(2.0-FRAC);
      //  ALPHA=ALPHA+2*ALPHA*(dum2-FRAC*VOL/2);

		System.out.println("New coupling kappa_d is " + kappa_d);
		System.out.println("Average disk volume=" + dum);
		System.out.println("Average boundary length=" + dum2);
			System.out.println("sum of total boundary length=" + b_monitor);
			//System.out.println("MARKEDQ is " + MARKEDQ);
		//System.out.println("ALPHA is "+ALPHA);
		vol_monitor=0;
        b_monitor=0;
		num_monitor=0;

		return;
	}
	/* routine returns sum of vertices around the face conjugate to node i */

	public int sum_face(SIMPLEX p, int i)
	{
		int j,add;

		add=0;
		for(j=0;j<DPLUS;j++)
			if(p.vertex[j]!=i) add+=(p.vertex[j]);

		return(add);
	}
	/* every sweep clean up pointer array */

	void tidy()
	{
		SIMPLEX[] temp;
		int  i,add;
		temp=new SIMPLEX[BIGVOL];
		/* run down array compressing non NULL extries into new array */
		/* and reassigning simplex labels according to their new index in this */
		/* array. Finally copy back */

		add=0;


//		if(pointer_number>BIGVOL){printf("\nError -- need to increase BIGVOL");exit(1);}

		if(pointer_number>max_point) max_point=pointer_number;
		//max_point+=pointer_number;

		for(i=0;i<pointer_number;i++)
			if(simplex_point[i]!=null)
			{
				temp[add]=simplex_point[i];
				temp[add].label=add;
				add++;
			}

		for(i=0;i<add;i++)
			simplex_point[i]=temp[i];


		pointer_number=add;
		if(pointer_number!=simplex_number){
			System.out.println("oops - pointer number is not equal to simplex_number in tidy()");}
		return;
	}
	/* driver for triangulation updates */

	Integer subsimplex = D;
	void trial_change()
	{
		SIMPLEX simp;
		SIMPLEX[] addresses;
		int i,subsimplex2=0;
		int[] labels;
		int[] q;
		LOGIC legal_move,good_manifold,metro_accept;
		//System.out.println("in trial_change");

		addresses=new SIMPLEX[DPLUSPLUS];
		labels=new int[DPLUSPLUS];
		q=new int[MAXVOL];

		// grab triangle and move type at random
		simp=select_simplex(subsimplex);
        //System.out.println("subsimplex in trial_change is "+subsimplex);

		// check if move is legal i.e coordination of simplex =D+1-subsimplex
		legal_move=good_subsimplex(simp,subsimplex,labels,addresses);

		if(legal_move==LOGIC.NO) return;

		/*printf("Trying a %d-move with vertices:\n",subsimplex);
for(i=0;i<(subsimplex+1);i++)
printf("%d ",labels[i]);
printf("\n");
printf("extra vertices are:\n");
for(i=subsimplex+1;i<DPLUSPLUS;i++)
printf("%d ",labels[i]);
printf("\n");
fflush(stdout);
		 */

		legal_subsimplex[subsimplex]++;

		// make sure move will not create degeneracies
		good_manifold=allowed_move(simp,subsimplex,labels);

		if(good_manifold==LOGIC.NO) return;

		manifold_subsimplex[subsimplex]++;

		// check change in action
		metro_accept=metropolis(subsimplex,labels,addresses);

		if(metro_accept==LOGIC.NO) return;
		/*
for(i=0;i<DPLUS;i++){
q[i]=find_order(labels[i],addresses[DPLUS]);}
if(subsimplex!=D){
q[DPLUS]=find_order(labels[DPLUS],addresses[subsimplex+1]);}
else{
q[DPLUS]=0;}

for(i=0;i<DPLUSPLUS;i++){
printf("before:vertex %d has coord %d\n",labels[i],q[i]);}
		 */

		go_subsimplex[subsimplex]++;

		// if accept update triangulation
		update(labels,addresses,subsimplex);

		//printf("simplex pointer %d\n",pointer_number);

		/*
for(i=1;i<DPLUSPLUS;i++){
q[i]=find_order(labels[i],addresses[0]);}
if(subsimplex!=0){
q[0]=find_order(labels[0],addresses[subsimplex]);}
else{
q[0]=0;}

for(i=0;i<DPLUSPLUS;i++){
printf("after:vertex %d has coord %d\n",labels[i],q[i]);}
		 */

		return;

	}
	/* coordinates addition of new simplices and removal of old ones */

	void update(int[] a, SIMPLEX[] q, int sub)
	{
		int[] c, temp;
		int i,j,k,m,n,o;

		c=new int[DPLUSPLUS];
		temp=new int[DPLUS];
		/* if subsimplex is node then save its label on the stack */

		//printf("in update");fflush(stdout);

		if(subsimplex==0)
		{
			push(a[0]);
			node_number--;
		}

		if(subsimplex==D)
		{
			pop();
			node_number++;
		}

		/* loop over new simplices */

		for(i=0;i<subsimplex+1;i++)
		{
			combo(a,c,subsimplex+1,i);

			for(j=subsimplex+1;j<DPLUSPLUS;j++)
				c[j-1]=a[j];

			q[i]=new SIMPLEX(c, D);
			simplex_point[pointer_number]=q[i];

			simplex_number++;
			pointer_number++;
		}

		/* now reconnect pointers appropriately */

		reconnect(a,q,subsimplex);

		/*  old guys */

		for(i=subsimplex+1;i<DPLUSPLUS;i++){
			simplex_point[q[i].label]=null;
			simplex_number--;}

		return;
	}
}
