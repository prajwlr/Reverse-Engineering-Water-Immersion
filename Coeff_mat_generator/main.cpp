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
double theta,thetadeg;
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
double coeff_mat[1010][1010];
double test_coeff_mat[1010][1010];
llint R,C;
llint cur_rank=0; // Current Rank of coeff_mat
double test_row[100010]={0};
llint out_theta_arr[100010];
llint out_depth_arr[100010];
double thetainval;

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

void coeff_matrix_init(){
    R=n*n;
    C=n*n;
    loop(i,0,n*n){
        loop(j,0,n*n){
            coeff_mat[i][j]=0;
        }
    }

    return;
}

void swap_mat( llint row1, llint row2,llint col){
	for (llint i = 0; i < col; i++)
	{
		double temp = test_coeff_mat[row1][i];
		test_coeff_mat[row1][i] = test_coeff_mat[row2][i];
		test_coeff_mat[row2][i] = temp;
	}
}

llint rankOfMatrix(){
	llint trank = C;

	for (llint row = 0; row < trank; row++)
	{
		if (test_coeff_mat[row][row])
		{
		for (llint col = 0; col < R; col++)
		{
			if (col != row)
			{
				double mult = test_coeff_mat[col][row] /test_coeff_mat[row][row];
				for (llint i = 0; i < trank; i++)
				test_coeff_mat[col][i] -= mult * test_coeff_mat[row][i];
			}
		}
		}
		else
		{
			llint reduce = 1;
			for (llint i = row + 1; i < R; i++)
			{
				if (test_coeff_mat[i][row])
				{
					swap_mat( row, i, trank);
					reduce = 0;
					break ;
				}
			}
			if (reduce)
			{
				trank--;

				for (llint i = 0; i < R; i ++)
					test_coeff_mat[i][row] = test_coeff_mat[i][trank];
			}
			row--;
		}
	}
	return trank;
}

void weightage_mat_gen(){ // Sets the value of test row
    theta = (thetadeg * M_PI )/180;
    sint=sin(theta);
    cost=cos(theta);

    //cout << "sint = " << sint << endl;
    //cout << "cost = " << cost << endl;
    //ROTATES MATRIX BY THETA ANGLE AND RELOADS NEW X AND Y CO-ORDINATES OF POINTS ON MATRIX

    matrix_init();



    loop(i,0,n){
        loop(j,0,n){
            test_row[n*i + j] = area(i,j);
        }
    }

    /*loop(i,0,n){
        loop(j,0,n){
            outfile << areamat[i][j] << " ";
            cout << areamat[i][j] << " ";
        }
        cout << endl;

    }

    loop(i,0,n*n){
        cout << test_row[i] << " ";
    }

    cout << endl;

    outfile.close();
*/
}

void replace_row(){
    loop(i,0,R){
        coeff_mat[cur_rank][i] = test_row[i];
    }
    return;
}

void test_row_init(){
    loop(i,0,R){
        test_row[i]=0;
    }
    return;
}

void copycoeffmat(){
    loop(i,0,R){
        loop(j,0,C){
            test_coeff_mat[i][j]=coeff_mat[i][j];
        }
    }
}

llint all_one(){
    loop(i,0,R){
        if(test_row[i]!=1)return 0;
    }
    return 1;
}

int main()
{
    /*cout << "Input the angle of rotation of axes in degerees" << endl;
    cin >> theta;*/
    cout << "Input the size of n*n matrix" << endl;
    cin >> n;

    /*cout << "Input the depth by of immersion of object into water" << endl;
    cin >> k;*/


    ofstream outfile;
    outfile.open("coeff.txt");
    outfile << setprecision(7) << fixed;
    cout << setprecision(7) << fixed;

    coeff_matrix_init();
    thetainval=90;
    while(cur_rank<R){
        for(thetadeg = thetainval;thetadeg>=(double)0;thetadeg=thetadeg-(double)1){
            //cout << thetadeg << endl;
            k=1;
            test_row_init();
            for(k=1;all_one()==0;k++){
                cout << "thetadeg = " << thetadeg << "  _  k ="  << k << "cur_rank = " ;
                test_row_init();
                weightage_mat_gen();
                replace_row();
                copycoeffmat();
                /*if(thetadeg==50&&k==2){
                    loop(i,0,n){
                        loop(j,0,n){
                            cout << test_row[n*i+j] << " ";
                        }
                        cout << endl;
                    }
                }*/
                if(rankOfMatrix()>cur_rank){
                    //cout << "new rank = " << rankOfMatrix() << endl;
                    out_depth_arr[cur_rank]=(llint)k;
                    out_theta_arr[cur_rank]=(llint)thetadeg;
                    cur_rank++;
                }
                cout << cur_rank << endl;
                if(cur_rank==R)break;
            }

            if(cur_rank==R)break;
        }
        if(cur_rank==R)break;
        thetainval+=1;
    }

    cout << endl << endl << "Output time \n\n\n" << endl;

    loop(i,0,R){
        cout << out_theta_arr[i] << " " << out_depth_arr[i] << endl;
    }

    cout << "cur_rank = " << cur_rank << endl;

    cout << endl << endl;

    loop(i,0,R){
        loop(j,0,C){
            //cout << coeff_mat[i][j] << " ";
            outfile << coeff_mat[i][j] << " ";
        }
        //cout << endl;
        outfile << endl;
    }
    return 0;
}
