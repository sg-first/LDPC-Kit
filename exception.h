#pragma once
#include <string>

class MatrixException
{
public:
    virtual std::string info()=0;
    virtual ~MatrixException(){}
};

class SquareException : public MatrixException
{
public:
    std::string info()
    {
        return "非方阵";
    }
};

class DimensionException : public MatrixException
{
public:
    std::string info()
    {
        return "维度不匹配";
    }
};

class SingularException : public MatrixException
{
public:
    std::string info()
    {
        return "奇异矩阵";
    }
};
