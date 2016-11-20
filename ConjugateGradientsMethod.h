#ifndef _ITER_H_
#define _ITER_H_
#include "GridClass.h"
#include <mpi.h>

class FuncItertor {
protected:
    const Function func;
    long iterationsCount;
public:
    FuncItertor(const Function &func, const Grid& m):
        func(func), 
        iterationsCount(0) {}
    virtual double iterate(Grid &grid) = 0;
};

//Последовательная версия метода сопряженных градиентов
class ConjugateGradientsMethod_Base: public FuncItertor {
private:
    Grid rGrid; // Матрица для значений r
    Grid gGrid; // Матрица для значений q
    double SteepestDescentMethodForZeroIter(Grid &pGrid);
    double coefTau;
public:
    ConjugateGradientsMethod_Base(Function func, const Grid& g):
        FuncItertor(func,g),
        gGrid(g.getLeftBottomCorner(), 
              g.getRightTopCorner(), 
              g.getRowsCount(),
              g.getColumnsCount(),
              g.getParentRows(), 
              g.getParentCols(),
              g.getRowsChanging(), 
              g.getColumnsChanging()),
        rGrid(g.getLeftBottomCorner(), 
              g.getRightTopCorner(), 
              g.getRowsCount(),
              g.getColumnsCount(),
              g.getParentRows(), 
              g.getParentCols(),
              g.getRowsChanging(), 
              g.getColumnsChanging()) {}
    virtual double iterate(Grid &pGrid);
};

//MPI версия метода сопряженных градиентов
class ConjugateGradientsMethod_MPI : public FuncItertor {
protected:
    int left, right, top, bottom;
    Grid rGrid; // Матрица для значений r
    Grid gGrid; // Матрица для значений q
    virtual double SteepestDescentMethodForZeroIter(Grid &pGrid);
    void getGridBorders(Grid &grid);
    double coefTau;
    int rank, size;
    inline bool checkBorder(const Grid &grid, long i, long j) {
        return (
                grid.getRowsChanging() + i - 1 >= 0
                && grid.getRowsChanging() + i + 1 < grid.getParentRows()
                && grid.getColumnsChanging() + j - 1 >= 0
                && grid.getColumnsChanging() + j + 1 < grid.getParentCols()
               );
    }

public:
    ConjugateGradientsMethod_MPI(Function func, const Grid& g, int rank, int left, int right, int top, int bottom, int size):
        FuncItertor(func,g),
        gGrid(g.getLeftBottomCorner(), 
              g.getRightTopCorner(), 
              g.getRowsCount(),
              g.getColumnsCount(), 
              g.getParentRows(), 
              g.getParentCols(), 
              g.getRowsChanging(), 
              g.getColumnsChanging()),
        rGrid(g.getLeftBottomCorner(), 
              g.getRightTopCorner(), 
              g.getRowsCount(),
              g.getColumnsCount(), 
              g.getParentRows(), 
              g.getParentCols(), 
              g.getRowsChanging(), 
              g.getColumnsChanging()),
        rank(rank), left(left), right(right), top(top), bottom(bottom), size(size)
    {}
    virtual double iterate(Grid &pGrid);
};

//MPI / OpenMP версия метода сопряженных градиентов
class ConjugateGradientsMethod_MPI_OpenMP: public ConjugateGradientsMethod_MPI {
private:
    double SteepestDescentMethodForZeroIter(Grid &pGrid);
public:
    ConjugateGradientsMethod_MPI_OpenMP(Function func, const Grid& m, int rank, int left, int right, int top, int bottom, int size):
        ConjugateGradientsMethod_MPI(func,m,rank,left,right,top,bottom,size)
    {}
    virtual double iterate(Grid &pGrid);
};

#endif
