#pragma once
#include <math.h>
#include <bitset>
#include <vector>

const unsigned int pExp = 3; //2^3=8
const unsigned int maxP=(pExp-1)*2+1;
const unsigned int ppHighSub=pExp; //x3=1000，下标为3，就是3
const unsigned int primitivePolynomials=3; //x+1=11=3

typedef std::bitset<maxP> binary;

class GF
{
public:
    const static unsigned int p = pow(2,pExp);
    static unsigned int mulTable[p*p][3];

    static unsigned int add(unsigned int a,unsigned int b)
    {
        return a^b;
    }

    static binary easyMul(binary b, unsigned int pos)
    {
        binary resultB(0);
        for(unsigned int j=0;j<b.size();j++)
        {
            unsigned int resultPos=pos;
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
        for(unsigned int i=maxP-1;i>=ppHighSub;i--)
        {
            if(b[i]==1)
            {
                b[i]=0;
                //先用pp推对应的式子
                binary pp(primitivePolynomials);
                for(unsigned int j=i;j>ppHighSub;j--) //比pp左边指数多几乘几个x
                    pp=easyMul(pp,1); //不管p是几x都是10，pos=1
                //把式子加到原来的上面
                b=b^pp;
            }
        }
    }

    static unsigned int rawMul(unsigned int a,unsigned int b)
    {
        if(a==0 || b==0)
            return 0;

        binary ab(a);
        binary bb(b);
        unsigned int size=ab.size();

        std::vector< binary > result;

        for(unsigned int i=0;i<size;i++)
        {
            if(ab[i]==1)
            {
                binary newr=easyMul(bb,i);
                ppMod(newr);
                result.push_back(newr);
            }
        }

        for(unsigned int i=1;i<result.size();i++)
            result[0]=result[0]^result[i];
        return result[0].to_ulong();
    }

    static void initMulTable()
    {
        unsigned int num=0;
        for(unsigned int i=0;i<p;i++)
            for(unsigned int j=0;j<p;j++,num++)
            {
                mulTable[num][0]=i;
                mulTable[num][1]=j;
                mulTable[num][2]=rawMul(i,j);
            }
    }

    static double mul(unsigned int a, unsigned int b)
    {
        return mulTable[p*a+b][2];
    }

    static double div(unsigned int a, unsigned int b)
    {
        unsigned int d=p*a;
        for(unsigned int i=d;i<d+p;i++)
        {
            if(mulTable[i][2]==b)
                return mulTable[i][1];
        }
        throw std::string("cannot div");
        return -1;
    }

    static double mulInv(unsigned int i) //1乘几为1
    {
        if(i==0)
            return 0;
        else
            return div(i,1);
    }
};
unsigned int GF::mulTable[p*p][3];
