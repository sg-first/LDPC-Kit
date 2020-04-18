#pragma once
#include <math.h>
#include <bitset>
#include <vector>

typedef unsigned int uint;

const uint pExp = 3; //2^3=8
const uint maxP=(pExp-1)*2+1;
const uint ppHighSub=pExp; //x3=1000，下标为3，就是3
const uint primitivePolynomials=3; //x+1=11=3

typedef std::bitset<maxP> binary;

class GF
{
public:
    const static uint p = pow(2,pExp);
    static uint mulTable[p*p][3];

    static uint add(uint a,uint b)
    {
        return a^b;
    }

    static binary easyMul(binary b, uint pos)
    {
        binary resultB(0);
        for(uint j=0;j<b.size();j++)
        {
            uint resultPos=pos;
            //即，当b[j]=1所在位为1（最低位，下标最大）时，resultPos不变，否则都会变小
            if(b[j]==1)
            {
                resultPos+=j; //为最低位（j最大）时减的是0，j越大（位数越小）时往前（高位）移的越少
                resultB[resultPos]=resultB[resultPos]^1;
            }
        }
        return resultB;
    }

    static void ppMod(binary &b)
    {
        for(uint i=maxP-1;i>=ppHighSub;i--)
        {
            if(b[i]==1)
            {
                b[i]=0;
                //先用pp推对应的式子
                binary pp(primitivePolynomials);
                for(uint j=i;j>ppHighSub;j--) //比pp左边指数多几乘几个x
                    pp=easyMul(pp,1); //不管p是几x都是10，pos=1
                //把式子加到原来的上面
                b=b^pp;
            }
        }
    }

    static uint rawMul(uint a,uint b)
    {
        if(a==0 || b==0)
            return 0;

        binary ab(a);
        binary bb(b);
        uint size=ab.size();

        std::vector< binary > result;

        for(uint i=0;i<size;i++)
        {
            if(ab[i]==1)
            {
                binary newr=easyMul(bb,i);
                ppMod(newr);
                result.push_back(newr);
            }
        }

        for(uint i=1;i<result.size();i++)
            result[0]=result[0]^result[i];
        return result[0].to_ulong();
    }

    static void initMulTable()
    {
        uint num=0;
        for(uint i=0;i<p;i++)
            for(uint j=0;j<p;j++,num++)
            {
                mulTable[num][0]=i;
                mulTable[num][1]=j;
                mulTable[num][2]=rawMul(i,j);
            }
    }

    static double mul(uint a, uint b)
    {
        return mulTable[p*a+b][2];
    }

    static double div(uint a, uint b)
    {
        if(a==0) //fix: 不知道对不对，应该对吧
            return 0;

        uint d=p*a;
        for(uint i=d;i<d+p;i++)
        {
            if(mulTable[i][2]==b)
                return mulTable[i][1];
        }
        throw std::string("cannot div");
        return -1;
    }

    static double mulInv(uint i) //1乘几为1
    {
        if(i==0)
            return 0;
        else
            return div(i,1);
    }
};
