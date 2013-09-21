/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2012 by Vito Scarola <scarola@vt.edu>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/


/******* Code to generate basis states *********/

int mmin=0;
int ngs = 2;
int r,min;
int an, tmp, k, j;
ivector_type ii(N+1);
ii(1)= mmin;

for(k=2; k<=N; k++) {
    if(ii(k-1)<mmax) {
        ii(k)=ii(k-1)+1;
    }
}

an=0;
for(int l=1; l<=N; l++) an += ii(l) ;

if( an==Lz ) {

    nmax++;
    if(nmax >nguess) {
        printf("\n\n nmax>nguess \n\n");
        exit(EXIT_FAILURE);
    }
    for(int l=1; l<=N; l++) {
        lst(nmax,l)=ii(l);
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

            nmax++;
            if(nmax >nguess) {
                printf("\n\n nmax>nguess \n\n");
                exit(EXIT_FAILURE);
            }
            for(int l=1; l<=N; l++) {

                lst(nmax,l)=ii(l);
            }

        } /*if(an==k)loop*/
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


