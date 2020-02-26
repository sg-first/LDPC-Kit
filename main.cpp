#include <iostream>
#include <vector>
#include "matrix.h"

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

int main()
{
    std::vector<int>av;
    matrix H(7,7);
    av={1,0,0,0,0,0,0};
    assignment(H.m[0],av);
    av={1,1,0,0,0,0,0};
    assignment(H.m[1],av);
    av={0,1,1,0,0,0,0};
    assignment(H.m[2],av);
    av={0,0,1,1,0,0,0};
    assignment(H.m[3],av);
    av={0,1,0,1,1,0,0};
    assignment(H.m[4],av);
    av={1,0,0,0,1,1,0};
    assignment(H.m[5],av);
    av={0,0,0,1,0,1,1};
    assignment(H.m[6],av);

    H.output();
    endl();
    matrix i=H.inv();
    i.output();
    endl();
    H.dot(i).output();
    matrix ii=inputM();
    std::cout<<ii.det();

    /*matrix E(3,7);
    inputElm(E);
    std::cout<<std::endl;
    matrix T(7,7);
    inputElm(T);
    std::cout<<std::endl;
    matrix B(7,3);
    inputElm(B);
    std::cout<<std::endl;
    matrix D(3,3);
    av={0,0,0};
    assignment(D.m[0],av);
    av={1,0,0};
    assignment(D.m[1],av);
    av={0,1,1};
    assignment(D.m[2],av);
    matrix r=E.dot(T);
    r.output();
    std::cout<<std::endl;
    r=r.dot(B);
    r.output();
    std::cout<<std::endl;
    r=r.add(D);
    r.output();*/

    /*std::string s;
    std::cin>>s;
    if(s=="inv")
    {
        matrix m=inputM();
        std::cout<<std::endl<<"inv:"<<std::endl;
        m.inv().output();
    }
    else if(s=="verify")
    {
        unsigned int l;
        std::cout<<"c l:";
        std::cin>>l;
        matrix v(1,l);
        inputElm(v);
        v=v.transpose();
        std::cout<<std::endl<<"result:"<<std::endl;
        H.dot(v).output();
    }*/
}
