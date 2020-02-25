#pragma once
#include <string>
using namespace std;

class MatrixException
{
public:
    virtual string info()=0;
    virtual ~MatrixException(){}
};

class SquareException : public MatrixException
{
public:
    string info()
    {
        return "非方阵";
    }
};

class DimensionException : public MatrixException
{
public:
    string info()
    {
        return "维度不匹配";
    }
};
