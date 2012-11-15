
/******* Code to generate basis states *********/

mmin=0;
ngs = 2;
//int r,min;
//int an, tmp, k, j;
//ivector_type ii(N+1);
for(k=1; k<=N; k++) ii(k)=0;
ii(1)= mmin;

for(k=2; k<=N; k++) {
    if(ii(k-1)<mmax) {
        ii(k)=ii(k-1)+1;
    }
}

an=0;
for(int l=1; l<=N; l++) an += ii(l) ;

if( an==Lz ) {
    nmaxm++;
    if(nmaxm >nguess) {
        printf("\n\n nmaxm>nguess \n\n");
        exit(EXIT_FAILURE);
    }
    for(int l=1; l<=N; l++) {
        lstm(nmaxm,l)=ii(l);
    }
}

r=0;
k=N;
j=0;

do {
    ii(k)=ii(k)+1;
    for(int l=k; l<=N-1; l++) {
        if(ii(l)<mmax) {
            ii(l+1)=ii(l)+1;
        }
        else {
        }
    }

    if( ii(k)<=mmax  &&  ii(N)<=mmax ) {
        an=0;
        for(int l=1; l<=N; l++) an += ii(l);

        if(an==Lz) {

            nmaxm++;
            if(nmaxm >nguess) {
                printf("\n\n nmaxm>nguess \n\n");
                exit(EXIT_FAILURE);
            }
            for(int l=1; l<=N; l++) {

                lstm(nmaxm,l)=ii(l);
            }

        } /*if(an==Lz)loop*/
        if(r>0) {
            k += r;
            r = 0;
        }
    }
    else {
        r++;
        k--;
        if(k<=0)j=ngs+1;
    }

    if(mmax-mmin+1>=N) {
        if(ii(1)== mmax)j=ngs+1;
    }
    else {
        printf("\n Not enough states");
        exit(EXIT_FAILURE);
    }
}
while(j<ngs);


