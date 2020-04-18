#include <iostream>
#include <vector>
#include "matrix.h"
#include "errorCorrection.h"
#include <thread>
#include <map>
#include <array>
#include <ctime>
#include "genH.h"
#include "matIO.h"

const uint sLength=2761;
const uint cLength=3098;

uint calu1(matrix m,uint start=0)
{
    uint result=0;
    for(uint i=0;i<m.getr();i++)
    {
        for(uint j=start;j<m.getc();j++)
        {
            if(m.m[i][j]==1)
                result++;
        }
    }
    return result;
}

void endl() { std::cout<<std::endl; }

bool check(matrix s,matrix p1,matrix p2,matrix H)
{
    matrix c(1,cLength);
    for(uint i=0;i<sLength;i++) //信息位逐个拷贝
        c.m[0][i]=s.m[0][i];
    for(uint i=sLength ; i<sLength+p1.getc() ; i++) //校验位1
        c.m[0][i]=p1.m[0][i-sLength];
    const uint offest2=sLength+p1.getc();
    for(uint i=offest2 ; i<cLength ; i++) //校验位2
        c.m[0][i]=p2.m[0][i-offest2];
    std::cout<<"check:"<<std::endl;
    s.output();
    c.output();

    return errorCorrection::check(c,H);
}

std::string binaryConversion(int num,int bin=GF::p)
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

matrix *errorMat=nullptr;

std::string randstr()
{
    std::string result="";
    for(uint i=0;i<sLength;i++)
        result+=std::to_string(std::rand() % 8);
    return result;
}

void checkLoop(uint min,uint max,
               matrix fii_ETiA_C, matrix TiA, matrix TiB,
               matrix H)
{
    matrix s(1,sLength);
    for(uint i=min;i<max;i++)
    {
        //std::string as=binaryConversion(i);
        matIO::assignment(s.m[0],randstr(),s.getc());

        matrix sT=s.transpose();

        matrix p1T=fii_ETiA_C.dot(sT);

        matrix ii=TiA.dot(sT);
        matrix ii2=TiB.dot(p1T);
        matrix p2T=ii.add(ii2);

        if(!check(s,p1T.transpose(),p2T.transpose(),H))
        {
            std::cout<<"fail!"<<std::endl;
            errorMat=new matrix(s);
        }

        if(errorMat!=nullptr)
            return;
    }
    std::cout<<"finished!"<<std::endl;
}

int main()
{
    //编码
    GF::initMulTable();

    std::vector<int>av;
    matrix H=matIO::ReadMatFile("D:/H(339×3070).csv",337,3098);

    uint g=H.getr()-199;
    uint mg=H.getr()-g;
    uint nm=H.getc()-H.getr();

    matrix T=H.cut(nm+g,0,H.getc()-1,mg-1);
    matrix Ti=T.inv();
    matrix E=H.cut(nm+g,mg,H.getc()-1,H.getr()-1);
    matrix B=H.cut(nm,0,nm+g-1,mg-1);
    matrix D=H.cut(nm,mg,nm+g-1,H.getr()-1);
    matrix fi=E.dot(Ti);
    fi=fi.dot(B);
    fi=fi.add(D);
    matrix fii=fi.inv();
    //matrix fii=matIO::ReadMatFile("D:/fii.csv",18,18);
    matrix A=H.cut(0,0,nm-1,mg-1);
    matrix C=H.cut(0,mg,nm-1,H.getr()-1);

    //计算所需
    matrix fii_ETiA_C=E.dot(Ti).dot(A).add(C);
    fii_ETiA_C=fii.dot(fii_ETiA_C); //这里原先有个逐元素取加法逆元的操作，因为结果不变去掉
    /*std::cout<<"fii_ETiA_C:"<<std::endl;
    fii_ETiA_C.output();*/

    matrix TiA=Ti.dot(A);  //这里原先有个逐元素取加法逆元的操作，因为结果不变去掉
    /*std::cout<<"TiA:"<<std::endl;
    TiA.output();*/

    matrix TiB=Ti.dot(B);  //这里原先有个逐元素取加法逆元的操作，因为结果不变去掉
    /*std::cout<<"TiB:"<<std::endl;
    TiB.output();*/

    checkLoop(0,50,fii_ETiA_C,TiA,TiB,H);

    //生成矩阵
    /*HGenerator hg;
    for(uint i=0;i<300;i++)
        hg.moveDetection();
    hg.permutationGF();
    hg.tetracyclicNum=-1; //重新给四环计数
    uint usefulNum;
    for(uint i=0;i<50000;i++)
    {
        if(hg.moveDetection())
            usefulNum=i;
        if(hg.tetracyclicNum==0)
            break;
    }
    std::cout<<"result:"<<hg.tetracyclicNum<<std::endl;
    std::cout<<"last:"<<50000-usefulNum<<std::endl;
    matIO::saveMatFile("D:/result.csv", hg.getH());*/

    //最大下三角
    /*matrix H=matIO::ReadMatFile("D:/27×1393校验矩阵.csv",genH::r,genH::c-2);
    H.maxTri();
    H.output();
    std::cout<<HGenerator::tetracyclicDetection(H).size()<<std::endl;*/

    //去掉两列
    /*auto cy=HGenerator::tetracyclicDetection(H);
    uint ary[8];

    uint pushNum=0;
    for(auto i : cy)
    {
        ary[pushNum]=i[1];
        pushNum++;
        ary[pushNum]=i[3];
        pushNum++;
    }

    matrix newH(1,1);
    for(uint i=0;i<7;i++)
    {
        for(uint j=i+1;j<8;j++)
        {
            newH=H.delC(ary[i],ary[j]);
            cy=HGenerator::tetracyclicDetection(newH);
            std::cout<<ary[i]<<" "<<ary[j]<<": ";
            std::cout<<cy.size()<<std::endl;
        }
    }*/
}
