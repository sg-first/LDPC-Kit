#pragma once
#include "matrix.h"
#include "matIO.h"
#include <map>
#include <array>
#include <ctime>

namespace  genH
{
    const uint diagSize=8;
    const uint rNum=3;
    const uint cNum=16;
    const uint r=rNum*diagSize;
    const uint c=cNum*diagSize;
    const uint diagNum=rNum*cNum;
}

typedef std::array<uint,4> cycle;
typedef std::array<matrix*,genH::diagNum> matArray;

class HGenerator
{
private:
    static void deleteDiag(matArray diag)
    {
        for(uint i=0;i<genH::diagNum;i++)
            delete diag[i];
    }

    static matArray copyDiag(matArray diag)
    {
        matArray rediag;
        for(unsigned int i=0;i<genH::diagNum;i++)
            rediag[i]=new matrix(*(diag[i]));
        return rediag;
    }


public:
    matArray diag;
    int tetracyclicNum=-1;

    HGenerator()
    {
        for(unsigned int i=0;i<genH::diagNum;i++)
            this->diag[i]=new matrix(matrix::identity(genH::diagSize));
    }

    ~HGenerator()
    {
        this->deleteDiag(this->diag);
    }

    static void generate(uint num=1, uint maxIterNum=50000)
    {
        for(uint j=0;j<num;j++)
        {
            srand((int)time(0));
            HGenerator hg;
            hg.permutationGF();
            uint usefulNum;
            for(uint i=0;i<maxIterNum;i++)
            {
                if(hg.moveDetection())
                    usefulNum=i;
                if(hg.tetracyclicNum==0)
                    break;
            }
            std::cout<<"result:"<<hg.tetracyclicNum<<std::endl;
            std::cout<<"last:"<<maxIterNum-usefulNum<<std::endl;
            matIO::saveMatFile("D:/H"+QString::number(j)+".csv", hg.getH());
        }
    }

    static matrix getH(const matArray &diag)
    {
        matrix H(genH::r,genH::c);
        uint sub=0;
        for(uint i=0;i<H.getr();i+=genH::diagSize)
        {
            for(uint j=0;j<H.getc();j+=genH::diagSize)
            {
                H.setArea(i,j,i+genH::diagSize,j+genH::diagSize,*(diag[sub]));
                sub++;
            }
        }
        return H;
    }

    matrix getH()
    {
        return getH(this->diag);
    }

    void rightMove(uint sub)
    {
        if(sub==genH::cNum-1) //右上角的不右移
            return;

        auto& diagi=this->diag[sub];
        for(uint i=0;i<genH::diagSize;i++)
        {
            for(uint j=0;j<genH::diagSize;j++)
            {
                if(diagi->m[i][j]!=0)
                {
                    uint nextJ=j+1;
                    if(j==genH::diagSize-1)
                        nextJ=0;
                    diagi->m[i][nextJ]=diagi->m[i][j];
                    diagi->m[i][j]=0;
                    break;
                }
            }
        }
    }

    static int randNum(int max,int min=0)
    {
        return (rand()%(max-min))+min;
    }

    void permutationGF()
    {
        for(uint sub=0;sub<genH::diagNum;sub++)
        {
            auto& diagi=this->diag[sub];
            for(uint i=0;i<genH::diagSize;i++)
            {
                for(uint j=0;j<genH::diagSize;j++)
                {
                    if(diagi->m[i][j]!=0)
                        diagi->m[i][j]=randNum(GF::p,2);
                }
            }
        }
    }

    bool moveDetection()
    {
        if(this->tetracyclicNum==-1)
        {
            matrix oldH=getH(this->diag);
            this->tetracyclicNum=tetracyclicDetection(oldH).size();
            std::cout<<"first:"<<this->tetracyclicNum<<std::endl;
        }
        matArray oldDiag=copyDiag(this->diag);
        //右移
        uint moveNum=randNum(genH::diagNum,1); //最少一个，最多比全部少一个
        for(uint ii=0;ii<moveNum;ii++) //对随机个小矩阵右移
        {
            uint moveSub=randNum(genH::diagNum); //随机选择要右移的小矩阵，第几个都行
            for(uint i=0;i<randNum(genH::diagSize,1);i++) //每个右移随机次，最少一次
                this->rightMove(moveSub);
        }
        //检测新的四环
        matrix newH=getH(this->diag);
        int newTcNum=tetracyclicDetection(newH).size();
        std::cout<<newTcNum<<std::endl;
        if(newTcNum<this->tetracyclicNum) //新的四环更少，换
        {
            this->tetracyclicNum=newTcNum;
            deleteDiag(oldDiag);
            std::cout<<"new:"<<this->tetracyclicNum<<std::endl;
            return true;
        }
        else
        {
            //this.diag恢复原来的
            deleteDiag(this->diag);
            this->diag=oldDiag;
            return false;
        }
    }

    static std::vector<cycle> tetracyclicDetection(const matrix& H)
    {
        std::vector<cycle> result;

        auto updateCD=[](std::map<uint,std::vector<uint>> &allCD,uint b,uint pos) //map first是元素，second是所有位置
        {
            if (allCD.count(b) == 0) //没有
                allCD[b]=std::vector<uint>();
            allCD[b].push_back(pos);
        };

        auto deteGroup=[&](uint p1,uint p2,uint b,uint nowi)
        {
            for(uint i=nowi+1;i<H.getr();i++)
            {
                if(H.m[i][p1]!=0 && H.m[i][p2]!=0)
                    if(H.m[i][p1] == H.m[i][p2]) //横着分别两两相等即可
                        result.push_back({nowi,p1,i,p2});
            }
        };

        auto deteAllGroup=[deteGroup](const std::vector<uint> &allp,uint b,uint nowi) //所有位置，元素，当前列（在这一列找）
        {
            for(uint i=0;i<allp.size()-1;i++)
                for(uint j=i+1;j<allp.size();j++)
                    deteGroup(allp[i],allp[j],b,nowi); //同行两两匹配，在这一列找
        };

        for(uint i=0;i<H.getr();i++)
        {
            std::map<uint,std::vector<uint>> allCD; //同一行中，每个元素都在哪些位置出现
            for(uint j=0;j<H.getc();j++)
            {
                if(H.m[i][j]!=0)
                    updateCD(allCD,H.m[i][j],j); //找同行一样的
            }
            //所有一样的找到，开始检测
            for (auto iter=allCD.begin(); iter!=allCD.end(); iter++)
                deteAllGroup(iter->second,iter->first,i);
        }

        return result;
    }

};
