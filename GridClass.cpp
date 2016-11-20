#include "GridClass.h"
#include <iostream>
#include <iomanip>
using namespace std;

const double Grid::Q_COEF = 1.5;
const double Grid::Q_COEF_TWOPOW = 1.82842712474;

//==========================Конструкторы=======================================
Grid::Grid(PointT leftBottomCorner, PointT rightTopCorner, long rows, long cols, long parentRows, long parentCols, long rowsChanging, long colsChanging):
    leftBottomCorner(leftBottomCorner),
    rightTopCorner(rightTopCorner),
    rowsChanging(rowsChanging),
    colsChanging(colsChanging),
    parentRows(parentRows),
    parentCols(parentCols),
    gridData(rows, cols, 0) {
    initPointCache();
}

Grid::Grid(PointT leftBottomCorner, PointT rightTopCorner, long rows, long cols,  long parentRows, long parentCols, double *gridData, long rowsChanging, long colsChanging):
    leftBottomCorner(leftBottomCorner),
    rightTopCorner(rightTopCorner),
    rowsChanging(rowsChanging),
    colsChanging(colsChanging),
    parentRows(parentRows),
    parentCols(parentCols),
    gridData(gridData, rows, cols) {
    initPointCache();
}
Grid::Grid(PointT leftBottomCorner, PointT rightTopCorner, long parentRows, long parentCols, const Matrix &gridData, long rowsChanging, long colsChanging):
    leftBottomCorner(leftBottomCorner),
    rightTopCorner(rightTopCorner),
    rowsChanging(rowsChanging),
    colsChanging(colsChanging),
    parentRows(parentRows),
    parentCols(parentCols),
    gridData(gridData) {
    initPointCache();
}
//=============================================================================
void Grid::initPointCache() {
    cacheForGridPoints.clear();
    for(long i = 0; i < gridData.rowsCount()+2; ++i ) {
        cacheForGridPoints.push_back(vector<PointD>(gridData.colsCount()+2));
    }
    for(long i = 0; i < gridData.rowsCount()+2; ++i){
        for(long j = 0; j < gridData.colsCount()+2; ++j){
            double dblI = static_cast<double>(i-1) + rowsChanging;
            double dblJ = static_cast<double>(j-1) + colsChanging;
            double xcoeff = dblI / (getParentRows() - 1);
            double ycoeff = dblJ / (getParentCols() - 1);
            double x = rightTopCorner.first * Q_COEF_tPOW(xcoeff);
            double y = rightTopCorner.second * Q_COEF_tPOW(ycoeff);
            cacheForGridPoints[i][j] = PointD(x,y);
        }
    }
}

PointD Grid::getPointFromCache(long i, long j) const {
    return cacheForGridPoints[i+1][j+1];
}

PointD Grid::getAverageGridSteps(long i, long j) const {
    PointD prevPoint = getPointFromCache(i-1, j-1);
    PointD curPoint = getPointFromCache(i,j);
    PointD nextPoint = getPointFromCache(i+1, j+1);
    double h1 = nextPoint.first - curPoint.first; 
    double h2 = nextPoint.second - curPoint.second;
    double prevH1 = curPoint.first - prevPoint.first;
    double prevH2 = curPoint.second - prevPoint.second;
    double avrgStepH1 = (h1 + prevH1) / 2;
    double avrgStepH2 = (h2 + prevH2) / 2;
    return PointD(avrgStepH1,avrgStepH2);
}

//========================Дополнительньные функции===========================
double fiveDotScheme(const Grid &grid, long i, long j) {
    PointD prevPoint = grid.getPointFromCache(i-1, j-1);
    PointD curPoint = grid.getPointFromCache(i,j);
    PointD nextPoint = grid.getPointFromCache(i+1, j+1);

    double h1 = nextPoint.first - curPoint.first;
    double h2 = nextPoint.second - curPoint.second;
    double prevH1 = curPoint.first - prevPoint.first;
    double prevH2 = curPoint.second - prevPoint.second;
    double avrgStepH1 = (h1 + prevH1) / 2;
    double avrgStepH2 = (h2 + prevH2) / 2;

    double leftPoint = j - 1 < 0 ? grid.left(i):grid(i, j-1);
    double rightPoint = j + 1 >= grid.getColumnsCount() ? grid.right(i) : grid(i, j+1);
    double bottomPoint = i + 1 >= grid.getRowsCount() ? grid.bottom(j):grid(i+1, j);
    double topPoint = i - 1 < 0 ? grid.top(j):grid(i-1, j);
    double ypart =  ((grid(i,j) - topPoint)/prevH1 - (bottomPoint - grid(i,j))/h1)/avrgStepH1;
    double xpart = ((grid(i,j) - leftPoint)/prevH2 - (rightPoint - grid(i,j))/h2)/avrgStepH2;
    return xpart + ypart;
}

void initGridBorder(Grid &grid,Function Q_COEF_tPOW) {
    for (size_t i = 0; i < grid.getRowsCount() ; ++i ) {
        grid(i,0) = Q_COEF_tPOW(grid.getPointFromCache(i,0));
        grid(i, grid.getColumnsCount() - 1) = Q_COEF_tPOW(grid.getPointFromCache(i, grid.getColumnsCount() - 1));
    } for (size_t j = 0; j < grid.getColumnsCount(); ++j) { 
        grid(0,j) = Q_COEF_tPOW(grid.getPointFromCache(0,j));
        grid(grid.getRowsCount() - 1,j) = Q_COEF_tPOW(grid.getPointFromCache(grid.getRowsCount() - 1, j));
    }
}

bool checkInGlobalBorder(const Grid &grid, long i, long j) {
    return (
        i + grid.rowsChanging >= 0
        && i + grid.rowsChanging < grid.parentRows
        && j + grid.colsChanging >= 0
        && j + grid.colsChanging < grid.parentCols
    );
}

PointT splitFunction(int totalRows, int totalCols, int procNum) {
    double rows, cols;
    int size = 0;
    n0 = (double) totalRows; 
    n1 = (double) totalCols;
    for(int i = 0; i < procNum; i++) {
        if(rows > cols) {
            rows = rows / 2.0;
            ++size;
        } else {
            cols = cols / 2.0;
        }
    }
    return PointT(size, procNum-size);
}

map<int, Grid> splitGrid(const Grid &grid, long procPower) {
    PointT size = splitFunction(grid.getRowsCount(), grid.getColumnsCount(), procPower);
    long rows = 1 << size.first;
    long cols = 1 << size.second;
    map<int, Grid> result;
    map<pair<int,int>, Matrix> submats = split(grid.gridData, rows, cols);
    int procnum = 0;
    long changingRows = -1, changingCols = -1;
    for(map<pair<int,int>, Matrix>::iterator it = submats.begin(); it != submats.end(); ++it) {
       int r = it->first.first;
       int c = it->first.second;
       if(changingRows == -1){
           changingRows = it->second.rowsCount();
           changingCols = it->second.colsCount();
       }
       result[procnum++] = Grid(grid.leftBottomCorner,
               grid.rightTopCorner, grid.getRowsCount(), grid.getColumnsCount(), it->second, r*changingRows, c*changingCols);
    }
    return result;
}

ostream &dropToStream(ostream& os, const Grid &grid) {
    for (long i = 0; i< grid.getRowsCount(); ++i){
        for(long j = 0; j < grid.getColumnsCount(); ++j) {
            PointD p = grid.getPointFromCache(i, j);
            double val = grid(i,j);
            os << p.first << "\t" << p.second << "\t" << val <<"\n";
        }
    }
    return os;
}

Grid collectGrid(const map<int, Grid> &subgrids) {
    Grid first = subgrids.at(0);
    Grid result(first.getLeftBottomCorner(), first.getRightTopCorner(), first.getParentRows(), first.getParentCols(), first.getParentRows(), first.getParentCols());
    for(map<int, Grid>::const_iterator itr = subgrids.begin(); itr!=subgrids.end(); ++itr){
        for(long i = 0; i < itr->second.getRowsCount(); ++i){
            for(long j = 0; j < itr->second.getColumnsCount(); ++j){
                result(i+itr->second.getRowsChanging(), j + itr->second.getColumnsChanging()) = itr->second(i,j);
            }
        }
    }
    return result;
}
//===========================================================================


