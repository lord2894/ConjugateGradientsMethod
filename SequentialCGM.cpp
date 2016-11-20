#include <iostream>
#include <cmath>
#include <getopt.h>
#include <cstdlib>
#include <fstream>
#include "GridClass.h"
#include "ConjugateGradientsMethod.h"
#include "InputFunctions.h"
using namespace std;


double Evaluation(Function FunctionP, const Grid &grid) {
    double diff = 0;
    for (long i = 0; i < grid.getRowsCount(); ++i){
        for(long j = 0; j < grid.getColumnsCount(); ++j) {
            diff += fabs( FunctionP(PointD(grid.getPointFromCache(i,j))) - grid(i,j) );
        }
    }
    return diff;
}

struct Args {
    long leftBottomCornerX;
    long leftBottomCornerY;
    long rightTopCornerX;
    long rightTopCornerY;
    long numRows;
    long numCols;
    string resultOutputFname;
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
                args.resultOutputFname= string(opt_val);
                break;
          }
          opt = getopt (argc, argv, "x:y:a:b:r:c:f:");
      }
      return args;
}

int main(int argc, char **argv) {
    //Получаем аргументы командной строки ввиде структуры
    Args args = ParseArgsStr(argc, argv);
    
    //Инициализируем выходную сетку
    PointT leftBottomCorner = PointT(args.leftBottomCornerX, args.leftBottomCornerY);
    PointT rightTopCorner = PointT(args.rightTopCornerX, args.rightTopCornerY);
    int totalNumRows = args.numRows;
    int totalNumCols = args.numCols;
    Grid resultGrid(leftBottomCorner, rightTopCorner, totalNumRows, totalNumCols, totalNumRows, totalNumCols);
    
    //Определяем значения на границах сетки
    initGridBorder(resultGrid, FunctionPhi);
    
    //Инициализируем итеративный метод сопряженных градиентов
    ConjugateGradientsMethod_Base iteratorCGM(FunctionF,resultGrid);
    
    //Фиксируем время начала
    time_t start;
    time(&start);
    
    //====================Итерируем до сходимости==================
    double error = iteratorCGM.iterate(resultGrid);
    int iteration = 1;
    cerr << "BottomCorner: " << args.leftBottomCornerX << ", " << args.leftBottomCornerY << "  TopCorner: " << args.rightTopCornerX  << ", " << args.rightTopCornerY<<"\n";
    while(error > EPSILON){        
        cerr << "RealError: " << error << " iter: " << iteration <<"\n"; 
        error = iteratorCGM.iterate(resultGrid);
        iteration++;
    }
    //=============================================================
    
    //Фиксируем время конца
    time_t end;
    time (&end);

    //Фиксируем время работы
    double diff = difftime (end,start);

    //Записываем выходные параметры 
    ofstream ofs(args.resultOutputFname.c_str());
    ofs << "[" << args.leftBottomCornerX << "," << args.rightTopCornerX << "]" << "x" 
        << "[" << args.leftBottomCornerY << "," << args.rightTopCornerY << "]" << "\t"
        << args.numRows << "\t" << args.numCols << "\t" << diff << "\n";
    return 0;
}
