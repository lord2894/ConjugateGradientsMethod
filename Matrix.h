#ifndef _MAT_H_
#define _MAT_H_
#define __volatile volatile // defining __volatile to volatile
#include <tr1/memory>
#include <map>
#include <cmath>
using namespace std;

template< typename T >
struct array_deleter
{
  void operator ()( T const * p){
    delete[] p;
  }
};

class Matrix;

class Vector {
friend class Matrix;
private:
    tr1::shared_ptr<double> data;
    long sz;
public:
    Vector():sz(0) {}
    Vector(long size):sz(size) {
        data = tr1::shared_ptr<double>(new double[size], array_deleter<double>());
    }
    Vector(long size, double val):sz(size) {
        data = tr1::shared_ptr<double>(new double[size], array_deleter<double>());
        for(int i = 0; i < size; ++i){
            data.get()[i]= val;
        }
    }
    Vector(double *data, long size): sz(size) {
        this->data = tr1::shared_ptr<double>(data);
    }
    Vector(tr1::shared_ptr<double> dt, long size):sz(size),data(dt){}
    Vector(const Vector& other) {
        sz = other.sz;
        data = other.data;
    }
    Vector &operator=(const Vector& other) {
        sz = other.sz;
        data = other.data;
        return *this;
    }
    double operator() (long index) const {
        return data.get()[index];
    }
    long size() const {
        return sz;
    }
    double *memptr() {
        return data.get();
    }
};

class Matrix {
private:
    tr1::shared_ptr<double> data;
    long rows;
    long cols;
public:
    Matrix():
        rows(0),
        cols(0){}
    Matrix(long r, long c):
        rows(r), 
        cols(c) {
        data = tr1::shared_ptr<double>(new double[rows * cols], array_deleter<double>());
    }
    Matrix(long r, long c, double val):
        rows(r), 
        cols(c) {
        data = tr1::shared_ptr<double>(new double[rows * cols], array_deleter<double>());
        for(long i = 0; i < rows * cols; ++i) {
            data.get()[i] = val;
        }
    }
    Matrix(double *data, long r, long c): 
        data(data, array_deleter<double>()), 
        rows(r), 
        cols(c) {}
    Matrix(const Matrix& other) {
        rows = other.rows;
        cols = other.cols;
        data = tr1::shared_ptr<double>(new double[rows * cols], array_deleter<double>());
        for(long i = 0; i < rows; ++i) {
            for(long j = 0; j < cols; ++j) {
                this->operator()(i,j) = other(i,j);
            }
        }
    }
    Matrix &operator=(const Matrix &other) {
        rows = other.rows;
        cols = other.cols;
        data = tr1::shared_ptr<double>(new double[rows * cols], array_deleter<double>());
        for(long i = 0; i < rows; ++i){
            for(long j = 0; j < cols; ++j){
                this->operator()(i,j) = other(i,j);
            }
        }
        return *this;
    }
    double &operator() (long i, long j) {
        return data.get()[i*cols+j];
    }
    double operator() (long i, long j) const {
        return data.get()[i*cols+j];
    }
    double *barememptr() {
        return data.get();
    }
    long rowsCount() const { return rows; }
    long colsCount() const { return cols; }
    Vector getRow(long index) const {
        Vector result(cols);
        for(long i = 0; i < cols; ++i) {
            result.data.get()[i] = data.get()[index*cols + i];
        }
        return result;
    }
    Vector getCol(long index) const {
        Vector result(rows);
        for(long i = 0; i < rows; ++i){
            result.data.get()[i] = data.get()[i*cols + index];
        }
        return result;
    }
    Matrix getSubmat(long startr, long startc, long rws, long cls) const {
        if(startr + rws > rows){
            rws = rows - startr;
        }
        if(startc + cls > cols) {
            cls = cols - startc;
        }
        Matrix result(rws,cls);
        for(long i = 0; i < rws; ++i){
            for(long j = 0; j < cls; ++j) {
                result(i,j) = this->operator()(i + startr, j+startc);
            }
        }
        return result;
    }
    friend map<pair<int,int>,Matrix> split(const Matrix& m, long rows, long cols) {
       long rElem = (long)(m.rowsCount() / (double)rows + 1);
       long cElem = (long)(m.colsCount() / (double)cols + 1);
       map<pair<int,int>, Matrix> result;
       for(long i = 0; i < rows; ++i){
           for(long j = 0; j < cols; ++j){
               result[pair<int,int>(i,j)] = m.getSubmat(i*rElem, j*cElem, rElem, cElem);
           }
       }
       return result;
    }
};
#endif
