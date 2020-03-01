#pragma once
#include "matrix.h"

typedef std::array<uint,primitivePolynomials> recommendNum;

class errorCorrection
{
private:
    matrix H;
    const uint iterationTimes=50;

    uint sumOther(uint i,uint nowJ)
    {
        uint result=0;
        for(uint jj=0;jj<this->H.getc();jj++) //对该c链接到的每个v
        {
            if(this->H.m[i][jj]!=0 && jj!=nowJ)
            {
                uint mulNum=GF::mul(this->C->m[0][jj],this->H.m[i][jj]);
                result=GF::add(result,mulNum);
            }
        }
        return GF::div(result,this->H.m[i][nowJ]);
    }

    std::tuple<recommendNum,recommendNum> getRecommend(uint j)
    {
        recommendNum num1; //每一位1的个数
        for(uint &i : num1)
            i=0;
        uint recommTotal;

        for(uint i=0;i<this->H.getr();i++) //看其所链接到的每个v的建议
        {
            if(this->H.m[i][j]!=0) //有链接，是建议
            {
                recommTotal++;
                binary b(this->allc2v.m[i][j]);
                for(uint bi=0;bi<b.size();bi++) //遍历，对每位位1的计数
                {
                    if(b[bi]==1)
                        num1[bi]++;
                }
            }
        }

        //算0的建议数
        recommendNum num0;
        for(uint i=0;i<num0.size();i++)
            num0[i]=recommTotal-num1[i];

        return std::make_tuple(num1,num0);
    }

public:
    matrix allc2v;
    matrix *C;

    errorCorrection(matrix H) : H(H), allc2v(H.getr(),H.getc()) {}

    ~errorCorrection()
    {
        delete[] this->C;
    }

    void setC(matrix C)
    {
        this->C=new matrix(C);
    }

    void correct()
    {
        for(uint n=0 ; n<this->iterationTimes && !check(*this->C,this->H) ; n++)
        {
            //计算每个c2v
            for(uint i=0;i<this->H.getr();i++) //对每个c
            {
                for(uint j=0;j<this->H.getc();j++) //链接到的每个v
                {
                    if(this->H.m[i][j]!=0)
                        allc2v.m[i][j]=this->sumOther(i,j);
                }
            }
            //判别翻转
            for(uint j=0;j<this->H.getc();j++) //对于每个v（this.c的每个元素）
            {
                recommendNum num1,num0;
                tie(num1,num0)=this->getRecommend(j);
                binary b(this->C->m[0][j]);
                for(uint ii=0;ii<b.size();ii++)
                {
                    //fix:目前是所有错都改
                    if(b[ii]==1)
                    {
                        if(num0[ii]>num1[ii])
                            b[ii]=0;
                    }
                    else
                    {
                        if(num1[ii]>num0[ii])
                            b[ii]=1;
                    }
                }
                this->C->m[0][j]=b.to_ulong();
            }
        }
    }

    static bool check(const matrix& c, const matrix& H)
    {
        matrix result=H.dot(c.transpose()).transpose();
        for(uint i=0;i<result.getc();i++)
        {
            if(result.m[0][i]!=0)
                return false;
        }
        return true;
    }
};
