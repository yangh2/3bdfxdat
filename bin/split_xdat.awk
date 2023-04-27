#!/usr/bin/awk -f

BEGIN{
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
}

{
    if ($1=="Direct"){
	nidx=nidx+1
	print ntype > "xdat_cut"nidx
	print tlist >> "xdat_cut"nidx
	print a1,a2,a3 >> "xdat_cut"nidx
	print b1,b2,b3 >> "xdat_cut"nidx
	print c1,c2,c3 >> "xdat_cut"nidx
	nstart=NR
    }

    if ((nstart>0)&&(NR-nstart > 0)&&(NR-nstart)<=natom){
	print $0 >> "xdat_cut"nidx
    }
}

