#!/usr/bin/awk -f

BEGIN{
    # ------------------------------------------------------
    # check before calcualtion
    #rmin=1.5
    #rc=4.5			# cutoff radius
    #nr1=100			# # of bins of r1
    #nr2=100			# # of bins of r2
    #ncost=100			# # of bins of cost
    # ------------------------------------------------------
    
    natom=0			# number of atom
    nmod=0
    ntpye=0
    nshift=7
    nidx=0;
    nstart=0;
}

NR==3{ a1=$1;a2=$2;a3=$3 }
NR==4{ b1=$1;b2=$2;b3=$3 }
NR==5{ c1=$1;c2=$2;c3=$3 }

NR==7{
    ntype=NF
    tlist=$0
    for (i=1;i<=NF;i++)
	natom=natom+$i;
    nmod=natom+1;

    v=(a1*(b2*c3-b3*c2)-a2*(b1*c3-b3*c1)+a3*(b1*c2-b2*c1));
    rho=natom/v;
    if (rho <0) rho = -rho;
    printf "%.8f ",rho > "3bdf_bins.ini"
    print natom >> "3bdf_bins.ini"
    print $0 >> "3bdf_bins.ini"
}

END{
    print 0 > "3bdf_bins.dat"
    print ntype > "3bdf_bins.dat"
    printf "%lf %lf %lf\n",rmin,rc,dr >> "3bdf_bins.dat"

    for (i=0;i<nr1*nr2*ncost;i++){
	for (j=0;j<ntype*ntype*ntype;j++)
	    printf "%d ", 0 >> "3bdf_bins.dat"
	printf "\n" >> "3bdf_bins.dat"
    }

    # print ntype > "3bdf_bins.x"
    # printf "%lf\n",rc >> "3bdf_bins.x"
    # print nr1,nr2,ncost >> "3bdf_bins.x"
    
    # for (i=0;i<nr1;i++){
    # 	for (j=0;j<nr2;j++)
    # 	    for (k=0;k<ncost;k++)
    # 		printf "%lf %lf %lf\n", rc*i/nr1,rc*j/nr2,2.0*k/ncost-1.0 >> "3bdf_bins.x"
    # }
}

