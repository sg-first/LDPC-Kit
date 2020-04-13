#pragma once
#include "matrix.h"
#include <QString>
#include <QFile>
#include <QTextStream>

class matIO
{
public:
    static void assignment(double* v,std::vector<int> av)
    {
        for(uint i=0;i<av.size();i++)
            v[i]=av[i];
    }

    static void assignment(double* v,std::string as,uint c)
    {
        for(uint i=0;i<c;i++)
            v[i]=as[i]-48;
    }

    static void inputElm(matrix& m)
    {
        for(uint i=0;i<m.getr();i++)
        {
            for(uint j=0;j<m.getc();j++)
                std::cin>>m.m[i][j];
        }
    }

    static matrix inputMat()
    {
        uint r,c;
        std::cout<<"r:";
        std::cin>>r;
        std::cout<<"c:";
        std::cin>>c;
        matrix m(r,c);
        inputElm(m);
        return m;
    }

    static QString ReadTXT(QString path)
    {
        auto PreQFile=new QFile(path);
        PreQFile->open(QIODevice::ReadOnly);
        QTextStream text(PreQFile);
        QString concert;
        concert=text.readAll();
        PreQFile->close();
        delete PreQFile;
        return concert;
    }

    static void writeTXT(QString path, QString content)
    {
        QFile file(path);
        file.open(QIODevice::WriteOnly);
        file.close();
        if (file.open(QIODevice::ReadWrite | QIODevice::Text))
        {
            QTextStream stream(&file);
            stream.seek(file.size());
            stream << content;
            file.close();
        }
    }

    static matrix ReadMatFile(QString path, uint r, uint c)
    {
        matrix m(r,c);
        QString content=ReadTXT(path);
        content.replace(",","");
        QStringList allLine=content.split('\n');
        for(uint i=0;i<r;i++)
            assignment(m.m[i],allLine[i].toStdString(),c);
        return m;
    }

    static void saveMatFile(QString path, matrix m)
    {
        QString result="";
        for(uint i=0;i<m.getr();i++)
        {
            for(uint j=0;j<m.getc();j++)
            {
                int val=m.m[i][j];
                result+=QString::number(val)+",";
            }
            result+="\n";
        }
        writeTXT(path,result);
    }
};
