#ifndef _MATRIX_H_
#define _MATRIX_H_
#include <cmath>
#include <map>
#include <cfloat>
#include <vector>
#include <iostream>
#include <iomanip>
#include "Matrix.h"
using namespace std;

typedef pair<size_t, size_t> PointT;
typedef pair<double, double> PointD;
static const PointD NAN_POINT(DBL_MAX, DBL_MAX);
typedef double (*Function)(PointD);
PointT splitFunction(int N0, int N1, int p);

class Grid {
    private:
        //Сетка
        Matrix gridData;
        //Граница сетки
        Vector top;
        Vector bottom;
        Vector left;
        Vector right;
        //Углы Сетки
        PointT leftBottomCorner;
        PointT rightTopCorner;
        //Кеш точек для повышения эффективности
        vector<vector<PointD> > cacheForGridPoints;
        void initPointCache();
        //
        long rowsChanging;
        long colsChanging;
        //
        long parentRows;
        long parentCols;
        //Построение неравномерной сетки
        static const double Q_COEF; // q=3/2
        static const double Q_COEF_TWOPOW; // 2^q-1
        inline static double Q_COEF_tPOW(double t) {
            return (pow(1+t, Q_COEF) - 1) / Q_COEF_TWOPOW;
        }
        PointD countPoint(long i, long j) const;
        
    public:
        //Конструкторы
        Grid(PointT leftBottomCorner, PointT rightTopCorner, long rows, long cols, long parentRows, long parentCols, long rowsChanging=0, long colsChanging=0);
        Grid(PointT leftBottomCorner, PointT rightTopCorner, long rows, long cols, long parentRows, long parentCols, double *gridData, long rowsChanging=0, long colsChanging=0);
        Grid(PointT leftBottomCorner, PointT rightTopCorner, long parentRows, long parentCols, const Matrix &gridData, long rowsChanging =0, long colsChanging = 0);
        Grid():leftBottomCorner(0,0), 
               rightTopCorner(0,0), 
               rowsChanging(0), 
               colsChanging(0), 
               parentRows(0),
               parentCols(0) { }
        Grid(const Grid &other):
            gridData(other.gridData),
            leftBottomCorner(other.leftBottomCorner), 
            rightTopCorner(other.rightTopCorner),
            rowsChanging(other.rowsChanging),
            colsChanging(other.colsChanging),
            parentRows(other.parentRows),
            parentCols(other.parentCols),
            cacheForGridPoints(other.cacheForGridPoints)
        {
        }
        
        //Операторы класса
        Grid &operator=(const Grid &other) {
            gridData = other.gridData;
            leftBottomCorner = other.leftBottomCorner; 
            rightTopCorner = other.rightTopCorner;
            rowsChanging = other.rowsChanging;
            colsChanging = other.colsChanging;
            parentRows = other.parentRows;
            parentCols = other.parentCols;
            cacheForGridPoints = other.cacheForGridPoints;
            return *this;
        }
        double &operator()(long i, long j) {
             return gridData(i,j); 
             }
        double operator()(long i, long j) const {
            return gridData(i,j);
        }
        friend ostream &operator<< (ostream &os, const Grid& g) {
            for(int i = 0; i< g.getRowsCount(); ++i){
                for(int j = 0; j < g.getColumnsCount(); ++j){
                    os << setw(10) << g(i,j) << "\t";
                }
                os << "\n";
            }
            return os;
        }
        //Class getters
        PointD getPointFromCache(long i, long j) const;
        PointD getAverageGridSteps(long i, long j) const;
        long getRowsCount() const { 
            return gridData.rowsCount(); 
        }
        long getColumnsCount() const { 
            return gridData.colsCount(); 
        }
        Vector getTopRow() {
            return gridData.getRow(0);
        }
        Vector getBottomRow() {
            return gridData.getRow(gridData.rowsCount() - 1);
        }
        Vector getLeftCol() {
            return gridData.getCol(0);
        }
        Vector getRightCol() {
            return gridData.getCol(gridData.colsCount() - 1);
        }
        double *getData() {
            return gridData.barememptr();
        }
        PointT getLeftBottomCorner() const {
            return leftBottomCorner; 
        }
        PointT getRightTopCorner() const {
            return rightTopCorner; 
        }
        long getRowsChanging() const { 
            return rowsChanging; 
        }
        long getColumnsChanging() const {
            return colsChanging; 
        }
        long getParentRows() const {
            return parentRows;
        }
        long getParentCols() const {
            return parentCols; 
        }
        

        //Class adders
        void addTopBorder(const Vector &top) {
            this->top = top;
        };
        void addBottomBorder(const Vector &bottom) {
            this->bottom = bottom;
        };
        void addLeftBorder(const Vector &left) {
            this->left = left;
        };
        void addRightBorder(const Vector &right) {
            this->right = right;
        }

        // Дополнительньные функции
        friend double fiveDotScheme(const Grid &grid, long i, long j); //Пятиточечная схема
        friend bool checkInGlobalBorder(const Grid &grid, long i, long j); //Проверка лежит ли точка на границе всей сетки
        friend void initGridBorder(Grid &grid, Function phi); //Инициализация границы сетки
        friend map<int, Grid> splitGrid(const Grid &grid, long procnum); // Разбиение сетки по узлам
        friend ostream &dropToStream(ostream& os, const Grid &grid); //Печать сетки
        friend Grid collectGrid(const map<int, Grid> &subgrids); //Сбор сетки в единый объект класса Grid, необходимо для сбора общей сетки с узлов
        //friend double getEuError(const Grid &grid);
       
};
#endif
