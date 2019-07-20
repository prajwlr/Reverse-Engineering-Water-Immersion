#include <iostream>
#include <bits/stdc++.h>
#include <fstream>
#include <iomanip>

using namespace std;

typedef long long int llint;
#define pb push_back
#define vvd vector< vector<double> >
#define loop(ii,a,b) for(long long int ii = a; ii < b; ++ ii)
#define PI 3.14159265
#define ff first
#define ss second

llint k;// Y = k LINE Cuts the weightage matrix.
double theta;
llint n;
double sint;
double cost;
double xmat[1010][1010]; //X COORDINATES OF POINTS OF MATRIX.
double ymat[1010][1010]; //Y COORDINATES OF POINTS OF MATRIX.
double areamat[1010][1010]; // REQUIRED WEIGHTAGE MATRIX WHICH IS TO BE OUTPUTED
vector < pair <double,double> > sq; // 4 POINTS OF CORNER OF A PARTICULAR CELL.
vector < pair <double,double> > yints; // contains intercepting points on each cell.
llint found;// simply tells whether y=k cuts the specified cell or not.
vector < pair <double,double> > belowpts;

double triang_area(pair <double , double> pa,pair <double , double> pb,pair <double , double> pc){
    double xa,xb,xc,ya,yb,yc;
    xa=pa.ff;
    xb=pb.ff;
    xc=pc.ff;
    ya=pa.ss;
    yb=pb.ss;
    yc=pc.ss;

    double outarea;

    outarea = (xa*(yb-yc)+xb*(yc-ya)+xc*(ya-yb))/2;
    outarea = abs(outarea);

    return outarea;
}

void find_intercepts(){ // sets proper values of intercepts due to y=k in yints vector.
    llint a,b;
    double xa,xb,ya,yb;
    double xintcpt;
    yints.clear();

    loop(i,0,4){
        a=i;
        b=(i+1)%4;
        xa=sq[a].ff;
        xb=sq[b].ff;
        ya=sq[a].ss;
        yb=sq[b].ss;

        if((ya-k)*(yb-k) < 0){
            if(found){
                xintcpt = xa + (((xa-xb)*(k-ya))/(ya-yb));
                yints.pb(make_pair(xintcpt,k));
            }else{
                found=1;
                xintcpt = xa + (((xa-xb)*(k-ya))/(ya-yb));
                yints.pb(make_pair(xintcpt,k));
            }
        }
    }
    return;
}

double area(llint ax,llint ay){
    double outarea=0;

    sq.clear();
    belowpts.clear();
    sq.pb(make_pair(xmat[ax][ay],ymat[ax][ay]));
    sq.pb(make_pair(xmat[ax+1][ay],ymat[ax+1][ay]));
    sq.pb(make_pair(xmat[ax+1][ay+1],ymat[ax+1][ay+1]));
    sq.pb(make_pair(xmat[ax][ay+1],ymat[ax][ay+1]));
    found=0;

    find_intercepts();

    if(!found){
        loop(i,0,4){
            if(sq[i].ss-k>0)return 0;
            if(sq[i].ss-k<0)return 1;
        }
    }

    sort(sq.begin(),sq.end());
    sort(yints.begin(),yints.end());
    belowpts.pb(yints[0]);

    loop(i,0,4){
        if(sq[i].ss < k)belowpts.pb(sq[i]);
    }

    belowpts.pb(yints[1]);

    loop(i,0,belowpts.size()-2){
        outarea+=triang_area(belowpts[0],belowpts[i+1],belowpts[i+2]);
    }

    return outarea;
}

void matrix_init(){
    loop(i,0,n+1){
        loop(j,0,n+1){
            xmat[j][i] = i*cost - j*sint;
            ymat[j][i] = j*cost + i*sint;
        }
    }
    return;
}

int main()
{
    cout << "Input the angle of rotation of axes in degerees" << endl;
    cin >> theta;
    cout << "Input the size of n*n matrix" << endl;
    cin >> n;
    cout << "Input the depth by of immersion of object into water" << endl;
    cin >> k;



    theta = (theta * M_PI )/180;
    sint=sin(theta);
    cost=cos(theta);

    cout << "sint = " << sint << endl;
    cout << "cost = " << cost << endl;
    //ROTATES MATRIX BY THETA ANGLE AND RELOADS NEW X AND Y CO-ORDINATES OF POINTS ON MATRIX
    matrix_init();

    ofstream outfile;
    outfile.open("coeff.txt");
    outfile << setprecision(7) << fixed;
    cout << setprecision(7) << fixed;


    //loop(kdet,0,1){
    //    k=kdet;
        loop(i,0,n){
            loop(j,0,n){
                areamat[i][j]=area(i,j);
            }
        }

        //cout << endl << "For DIP " << k << endl << endl;

        loop(i,0,n){
            loop(j,0,n){
                outfile << areamat[i][j] << " ";
                cout << areamat[i][j] << " ";
            }
            cout << endl;

        }
    //}

    outfile.close();


    return 0;
}
