#include <iostream>
#include <cmath>
#include <mpi.h>
#include <getopt.h>
#include <cstring>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <omp.h>
#include "GridClass.h"
#include "ConjugateGradientsMethod.h"
#include "InputFunctions.h"
using namespace std;

int getPowOfTwo(int val){
    int pwr = 0;
    while(val){
        val = val >> 1
        ++pwr;
    }
    return pwr;
}

struct Args {
    long leftBottomCornerX;
    long leftBottomCornerY;
    long rightTopCornerX;
    long rightTopCornerY;
    long numRows;
    long numCols;
    string resultGridOutputFname;
};

Args ParseArgsStr(int argc, char** argv) {
      int opt;
      char *opt_val = NULL;
      Args args;
      opt = getopt (argc, argv, "x:y:a:b:r:c:f:");
      while (opt != -1) {
          switch(opt) {
              case 'x':
                opt_val = optarg;
                args.leftBottomCornerX = strtol(opt_val, NULL, 10);
                break;
	      case 'y':
                opt_val = optarg;
                args.leftBottomCornerY = strtol(opt_val, NULL, 10);
                break;
	      case 'a':
                opt_val = optarg;
                args.rightTopCornerX = strtol(opt_val, NULL, 10);
                break;
	      case 'b':
                opt_val = optarg;
                args.rightTopCornerY = strtol(opt_val, NULL, 10);
                break;
              case 'r':
                opt_val = optarg;
                args.numRows = strtol(opt_val, NULL, 10);
                break;
            case 'c':
                opt_val = optarg;
                args.numCols = strtol(opt_val, NULL, 10);
                break;
            case 'f':
                opt_val = optarg;
                args.resultGridOutputFname= string(opt_val);
                break;
          }
          opt = getopt (argc, argv, "x:y:a:b:r:c:f:");
      }
      return args;
}

int main(int argc, char **argv) {
    int size, rank;
    Args args = ParseArgsStr(argc, argv);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int totalRows = args.numCols;
    int totalCols = args.numRows;
    long rowsChanging, colsChanging;
    long rows, cols;
    MPI_Status status;
    double start;
    int sizePower = getPowOfTwo(size);
    PointT procCR = splitFunction(totalRows, totalCols, sizePower);
    int procRows =  1 << procCR.first;
    int procCols = 1 << procCR.second;
    map<int, Grid> splited;
    MPI_Request dummy;
    if (rank == 0) {
        cerr << "Процессоров: " << size << "\n";
        cerr << "Сетка делится по: " <<procRows << "x"<<procCols << " на процессор\n";
        start = MPI_Wtime();
        PointT leftBottomCorner = PointT(args.leftBottomCornerX, args.leftBottomCornerY);
        PointT rightTopCorner = PointT(args.rightTopCornerX, args.rightTopCornerY);
        Grid resultGrid(leftBottomCorner, rightTopCorner, totalRows, totalCols, totalRows, totalCols);
        initGridBorder(resultGrid, FunctionPhi);
        splited = splitGrid(resultGrid, sizePower);
        for(map<int, Grid>::iterator itr = splited.begin(); itr != splited.end(); ++itr) {
            long r = itr->second.getRowsCount();
            long c = itr->second.getColumnsCount();
            long rChanging = itr->second.getRowsChanging();
            long cChanging = itr->second.getColumnsChanging();
            cerr << "Процессор: " << itr->first << " подсетка " << r <<"x" << c << "\n";
            if(itr->first != rank) {
                MPI_Send(&r, 1, MPI_LONG, itr->first, 0, MPI_COMM_WORLD);
                MPI_Send(&c, 1, MPI_LONG, itr->first, 0, MPI_COMM_WORLD);
                MPI_Send(&rChanging, 1, MPI_LONG, itr->first, 0, MPI_COMM_WORLD);
                MPI_Send(&cChanging, 1, MPI_LONG, itr->first, 0, MPI_COMM_WORLD);
            }else{
                rows = r;
                cols = c;
                rowsChanging = rChanging;
                colsChanging = cChanging;
            }
        }
    } else {
        MPI_Recv(&rows, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&cols, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&rowsChanging, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&colsChanging, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }
    Grid curGrid;
    if (rank == 0) {
        for(map<int, Grid>::iterator itr = splited.begin(); itr != splited.end(); ++itr) {
            if(itr->first != rank) {
                MPI_Send(itr->second.getData(), itr->second.getRowsCount()*itr->second.getColumnsCount(), MPI_DOUBLE, itr->first, 0, MPI_COMM_WORLD);
            }
        }
        curGrid = splited[0]; 
    } else {
        double *recdata = new double[rows*cols];
        MPI_Recv(recdata, rows*cols, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        PointT leftBottomCorner = PointT(args.leftBottomCornerX, args.leftBottomCornerY);
        PointT rightTopCorner = PointT(args.rightTopCornerX, args.rightTopCornerY);
        curGrid = Grid(leftBottomCorner, rightTopCorner, rows, cols, totalRows, totalCols, recdata, rowsChanging, colsChanging);
    }
    int procCol = rank%procCols;
    int procRow = (rank - procCol) / procCols;

    long left = procCol - 1 >=0 ? procRow*procCols + procCol - 1 : -1;
    long right = procCol + 1 < procCols ? procRow*procCols + procCol + 1 : -1;
    long top = procRow - 1 >=0? (procRow - 1)*procCols + procCol : -1;
    long bottom = procRow + 1 < procRows ? (procRow + 1)*procCols + procCol : -1;

    ConjugateGradientsMethod_MPI iter(FunctionF, curGrid, rank, left, right, top, bottom, size);
    MPI_Barrier(MPI_COMM_WORLD);

    double err = iter.iterate(curGrid);
    int iterCount = 1;
    while(err > EPSILON) {
        err = iter.iterate(curGrid);
        if (rank == 0) {
            cout <<"Итерация: " << iterCount++ <<" Текущая ошибка: " << err <<"\n";
        }
    }
    if (rank != 0) {
        MPI_Send(&rows, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&cols, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&rowsChanging, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&colsChanging, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(curGrid.getData(), rows*cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }else{
        map<int, Grid> subgrids;
        subgrids[0] = curGrid;
        vector<MPI_Request> requests;
        for (int i = 1; i < size; ++i ) {
            MPI_Recv(&rows, 1, MPI_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&cols, 1, MPI_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&rowsChanging, 1, MPI_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&colsChanging, 1, MPI_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            double *recIdata = new double[rows*cols];
            MPI_Recv(recIdata, rows*cols, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            PointT leftBottomCorner = PointT(args.leftBottomCornerX, args.leftBottomCornerY);
            PointT rightTopCorner = PointT(args.rightTopCornerX, args.rightTopCornerY);
            Grid curGrid(leftBottomCorner, rightTopCorner, rows, cols, totalRows, totalCols, recIdata, rowsChanging, colsChanging);
            subgrids[i] = curGrid;
        }
        Grid resultGrid = collectGrid(subgrids);
        double elapsed = MPI_Wtime() - start;
        ofstream outputRes(args.resultGridOutputFname.c_str());
        dropToStream(outputRes, resultGrid);
        outputRes <<"Result:" <<iterCount <<'\t' <<totalRows << '\t' << totalCols << '\t' << elapsed;
    }
    MPI_Finalize();
    return 0;
}
