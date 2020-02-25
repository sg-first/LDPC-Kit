#include<iostream>
#include "matrix.h"
using namespace std;

unsigned int calu1(matrix m,unsigned int start=0)
{
    unsigned int result=0;
    for(unsigned int i=0;i<m.getr();i++)
    {
        for(unsigned int j=start;j<m.getc();j++)
        {
            if(m.m[i][j]==1)
                result++;
        }
    }
    return result;
}

int main()
{
    unsigned int r,c;
    cout<<"r:";
    cin>>r;
    cout<<"c:";
    cin>>c;
    matrix m(r,c);

    for(unsigned int i=0;i<r;i++)
    {
        for(unsigned int j=0;j<c;j++)
            cin>>m.m[i][j];
    }

    m.output();

    unsigned int start;
    cout<<"start:";
    cin>>start;

    m.toUnit(start);
    m.output();
    cout<<"next"<<endl;
    m.radd(2,4);
    m.output();

    bool isBreak=false;
    do
    {
        matrix useM=m;
        isBreak=false;
        for(unsigned int i=0;i<r;i++)
        {
            for(unsigned int j=i+1;j<r;j++)
            {
                useM.radd(j,i);
                if(calu1(useM,start)<calu1(m,start) && calu1(useM)<calu1(m))
                {
                    cout<<j<<","<<i<<endl;
                    m=useM;
                    cout<<"next"<<endl;
                    m.output();
                    isBreak=true;
                    break;
                }
                else
                {
                    cout<<"fail"<<endl;
                    useM=m;
                }
            }
            if(isBreak)
                break;
        }
    } while(isBreak);
    m.output();
}
