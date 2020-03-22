#pragma once
#include "matrix.h"
#include <map>
#include <array>

namespace  genH
{
    const uint r=27;
    const uint c=1395;
    const uint diagSize=(r/9)*(c/9);
}

typedef std::array<uint,4> cycle;
typedef std::array<matrix*,genH::diagSize> matArray;

class HGenerator
{
private:
    static void deleteDiag(matArray diag)
    {
        for(uint i=0;i<genH::diagSize;i++)
            delete diag[i];
    }

    static matArray copyDiag(matArray diag)
    {
        matArray rediag;
        for(unsigned int i=0;i<genH::diagSize;i++)
            rediag[i]=new matrix(*(diag[i]));
        return rediag;
    }


public:
    matArray diag;
    int tetracyclicNum=-1;

    HGenerator()
    {
        for(unsigned int i=0;i<genH::diagSize;i++)
            this->diag[i]=new matrix(matrix::identity(9));
    }

    ~HGenerator()
    {
        this->deleteDiag(this->diag);
    }

    static matrix getH(const matArray &diag)
    {
        matrix H(genH::r,genH::c);
        uint sub=0;
        for(uint i=0;i<H.getr();i+=9)
        {
            for(uint j=0;j<H.getc();j+=9)
            {
                H.setArea(i,j,i+9,j+9,*(diag[sub]));
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
        /*uint cmax=genH::c/9;
        uint sub=r*cmax+1+c;*/
        auto& diagi=this->diag[sub];
        for(uint i=0;i<9;i++)
        {
            for(uint j=0;j<9;j++)
            {
                if(diagi->m[i][j]!=0)
                {
                    uint nextJ=j+1;
                    if(j==8)
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
        uint max=pow(2,pExp);
        for(uint sub=0;sub<genH::diagSize;sub++)
        {
            auto& diagi=this->diag[sub];
            for(uint i=0;i<9;i++)
            {
                for(uint j=0;j<9;j++)
                {
                    if(diagi->m[i][j]!=0)
                        diagi->m[i][j]=randNum(max,1);
                }
            }
        }
    }

    void moveDetection()
    {
        if(this->tetracyclicNum==-1)
        {
            matrix oldH=getH(this->diag);
            this->tetracyclicNum=tetracyclicDetection(oldH).size();
            std::cout<<"first:"<<this->tetracyclicNum<<std::endl;
        }
        matArray oldDiag=copyDiag(this->diag);
        //右移
        uint moveNum=randNum(genH::diagSize,1); //最少一个，最多比全部少一个
        for(uint ii=0;ii<moveNum;ii++) //对随机个小矩阵右移
        {
            uint moveSub=randNum(genH::diagSize); //随机选择要右移的小矩阵，第几个都行
            for(uint i=0;i<randNum(9,1);i++) //每个右移随机次，最少一次
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
        }
        else
        {
            //this.diag恢复原来的
            deleteDiag(this->diag);
            this->diag=oldDiag;
        }
    }

    std::vector<cycle> tetracyclicDetection(matrix H)
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
                if(H.m[i][p1]==b && H.m[i][p2]==b)
                    result.push_back({nowi,p1,i,p2});
            }
        };

        auto deteAllGroup=[deteGroup](const std::vector<uint> &allp,uint b,uint nowi)
        {
            for(uint i=0;i<allp.size()-1;i++)
                for(uint j=i+1;j<allp.size();j++)
                    deteGroup(allp[i],allp[j],b,nowi);
        };

        for(uint i=0;i<H.getr();i++)
        {
            std::map<uint,std::vector<uint>> allCD;
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
