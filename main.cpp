#include <iostream>
#include <vector>
#include "matrix.h"
#include <bitset>

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

void assignment(double* v,std::vector<int> av)
{
    for(unsigned int i=0;i<av.size();i++)
        v[i]=av[i];
}

void assignment(double* v,std::string as)
{
    for(unsigned int i=0;i<as.size();i++)
        v[i]=as[i]-48;
}

void inputElm(matrix& m)
{
    for(unsigned int i=0;i<m.getr();i++)
    {
        for(unsigned int j=0;j<m.getc();j++)
            std::cin>>m.m[i][j];
    }
}

matrix inputM()
{
    unsigned int r,c;
    std::cout<<"r:";
    std::cin>>r;
    std::cout<<"c:";
    std::cin>>c;
    matrix m(r,c);
    inputElm(m);
    return m;
}

void endl() { std::cout<<std::endl; }

bool check(matrix s,matrix p1,matrix p2,matrix H)
{
    matrix c(1,20);
    for(unsigned int i=0;i<10;i++)
        c.m[0][i]=s.m[0][i];
    for(unsigned int i=10;i<13;i++)
        c.m[0][i]=p1.m[0][i-10];
    for(unsigned int i=13;i<20;i++)
        c.m[0][i]=p2.m[0][i-13];
    std::cout<<"check:"<<std::endl;
    s.output();
    std::cout<<"result:"<<std::endl;
    matrix result=H.dot(c.transpose()).transpose();
    result.output();

    for(unsigned int i=0;i<result.getc();i++)
    {
        if(result.m[0][i]!=0)
            return false;
    }
    return true;
}

std::string binaryConversion(int num,int bin=vector::p)
{
    std::string result;
    do
    {
        result=std::to_string(num%bin)+result;
        num=num/bin;
    } while(num>=1);

    while(result.size()<bin)
        result="0"+result;

    return result;
}

unsigned int vector::mulTable[p*p][3];

int main()
{
    vector::initMulTable();

    std::vector<int>av;
    matrix H(10,20);
    av={2,2,0,0,0,2,0,0,0,0,2,0,2,2,0,0,0,0,0,0};
    assignment(H.m[0],av);
    av={0,0,0,0,4,0,4,0,4,0,0,4,0,4,4,0,0,0,0,0};
    assignment(H.m[1],av);
    av={0,3,0,3,0,0,0,3,0,0,3,0,0,0,3,3,0,0,0,0};
    assignment(H.m[2],av);
    av={0,0,1,0,0,0,1,0,6,0,0,0,6,0,0,6,6,0,0,0};
    assignment(H.m[3],av);
    av={0,0,0,7,0,0,0,0,0,7,0,7,0,0,7,0,7,7,0,0};
    assignment(H.m[4],av);
    av={5,0,0,0,5,0,0,5,0,0,0,0,0,5,0,0,0,5,5,0};
    assignment(H.m[5],av);
    av={0,0,2,0,0,0,2,0,0,2,0,0,0,0,0,0,2,0,2,2};
    assignment(H.m[6],av);
    av={0,0,0,2,0,2,0,0,2,0,0,2,2,0,0,2,0,2,0,2};
    assignment(H.m[7],av);
    av={0,4,0,0,0,4,0,4,0,4,4,0,0,0,0,0,0,0,0,4};
    assignment(H.m[8],av);
    av={3,0,3,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0};
    assignment(H.m[9],av);

    matrix T=H.cut(13,0,19,6);
    matrix Ti=T.inv();
    Ti.output();
    endl();
    matrix E=H.cut(13,7,19,9);
    matrix B=H.cut(10,0,12,6);
    matrix D=H.cut(10,7,12,9);
    matrix fi=E.dot(Ti);
    fi=fi.dot(B);
    fi=fi.add(D);
    fi.output();
    endl();
    matrix fii=fi.inv();
    fii.output();
    endl();
    matrix A=H.cut(0,0,9,6);
    matrix C=H.cut(0,7,9,9);

    matrix s(1,10);
    for(unsigned int i=0;i<1073741823;i++)
    {
        std::string as=binaryConversion(i);
        assignment(s.m[0],as);

        matrix sT=s.transpose();

        matrix ii=E.elmMulInv().dot(Ti);
        ii=ii.dot(A);
        ii=ii.add(C);
        ii=ii.dot(sT);
        matrix p1=fii.elmMulInv().dot(ii);
        p1=p1.transpose();

        ii=A.elmMulInv().dot(sT);
        matrix ii2=B.dot(p1.transpose());
        ii=ii.add(ii2);
        matrix p2=Ti.elmMulInv().dot(ii);
        p2=p2.transpose();

        if(!check(s,p1,p2,H))
        {
            std::cout<<"fail!";
            break;
        }
    }
}
