#include "ConjugateGradientsMethod.h"
using namespace std;
//Последовательная версия метода сопряженных градиентов
//Нулевая итерация - метод скорейшего спуска
double ConjugateGradientsMethod_Base::SteepestDescentMethodForZeroIter(Grid& pGrid) {
    //Инициализируем r и g
    for (long i = 1; i < pGrid.getRowsCount() - 1; ++i) {
        for (long j = 1; j < pGrid.getColumnsCount() - 1; ++j) {
            rGrid(i,j) = fiveDotScheme(pGrid, i, j) - func(pGrid.getPointFromCache(i,j));
            gGrid(i,j) = rGrid(i,j);
        }
    }
    //Рассчитываем коэффициент tau
    double num = 0;
    double denum = 0;
    int counter = 0;
    for (size_t i = 1; i < pGrid.getRowsCount() - 1; ++i) {
        for (size_t j = 1; j < pGrid.getColumnsCount() - 1; ++j) {
           num += rGrid(i,j) * rGrid(i,j);
           denum += fiveDotScheme(rGrid,i,j) * rGrid(i,j);
           counter++;
        }
    }
    coefTau = num / denum; 
    return 100000;
}

double ConjugateGradientsMethod_Base::iterate(Grid &pGrid) {
    if(iterationsCount == 0) {
        iterationsCount++;
        return SteepestDescentMethodForZeroIter(pGrid);
    }

    long rows = pGrid.getRowsCount();
    long cols = pGrid.getColumnsCount();
    double err = 0;
    
    //Рассчитываем значения p на текущей итерации
    //Раачитываем значение ошибки (невязки) между текщим p и p на предыдущей итерации для учета условия выхода
    for (size_t i = 1; i < rows - 1 ; ++i) {
        for (size_t j = 1; j < cols - 1; ++j) {
            double val = pGrid(i,j) - coefTau*gGrid(i,j);
            double errV = pGrid(i,j) - val;
            PointD hs = pGrid.getAverageGridSteps(i,j);
            err += errV*errV*hs.first*hs.second;
            pGrid(i,j) = val;
        }
    }

    //Рассчитываем значения r на текущей итерации
    for (long i = 1; i < rows - 1; ++i) {
        for (long j = 1; j < cols - 1; ++j) {
            rGrid(i,j) = fiveDotScheme(pGrid, i, j) - func(pGrid.getPointFromCache(i,j));
        }
    }

    //Рассчитываем скалярные произведения для рассчета коэффициента альфа
    double anum = 0, adenum = 0;
    for (long i = 1; i < rows - 1; ++i) {
        for (long j = 1; j < cols - 1; ++j) {

            anum += fiveDotScheme(rGrid, i, j) * gGrid(i,j);
            adenum += fiveDotScheme(gGrid, i, j) * gGrid(i,j);
        }
    }
    double alpha = anum / adenum;

    //Рассчитываем значения g на текущей итерации
    for (long i = 1; i < rows - 1; ++i) {
        for (long j = 1; j < cols - 1; ++j) {
            gGrid(i,j) = rGrid(i,j) - alpha*gGrid(i,j);
        }
    }
    
    //Рассчитываем скалярные произведения для рассчета коэффициента тау
    double tnum = 0, tdenum = 0;
    for (long i = 1; i < rows - 1; ++i) {
        for (long j = 1; j < cols - 1; ++j) {
            tnum += gGrid(i,j) * rGrid(i,j);
            tdenum += fiveDotScheme(gGrid, i, j) * gGrid(i, j);
        }
    }
    coefTau = tnum / tdenum;

    //Меняем номер итерации и выходим
    iterationsCount++;
    return sqrt(err);
}
//=============================================================

//MPI версия метода сопряженных градиентов
//Нулевая итерация - метод скорейшего спуска
double ConjugateGradientsMethod_MPI::SteepestDescentMethodForZeroIter(Grid &pGrid) {
    getGridBorders(pGrid);//Получаем границы сетки
    //Инициализируем r и g
    for (long i = 0; i < pGrid.getRowsCount(); ++i) {
        for (long j = 0; j < pGrid.getColumnsCount(); ++j) {
            if(checkBorder(rGrid, i, j)) {
                rGrid(i,j) = fiveDotScheme(pGrid, i, j) - func(pGrid.getPointFromCache(i,j));
                gGrid(i,j) = rGrid(i,j);
            }
        }
    }
    
    //Получаем граничные с другими узлами значения r и g
    getGridBorders(rGrid);
    getGridBorders(gGrid);

    //Вычисляем скалярные произвдения для рассчета коэффициента tau
    double num = 0, denum = 0;
    int counter = 0;
    for (size_t i = 0; i < pGrid.getRowsCount(); ++i) {
        for (size_t j = 0; j < pGrid.getColumnsCount(); ++j) {
            if (checkBorder(pGrid, i, j)) {
                num += rGrid(i,j) * rGrid(i,j);
                denum += fiveDotScheme(rGrid,i,j) * rGrid(i,j);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //Получаем со всех узлов аналогичные скалярные произведения
    double allNumerator, allDenum;
    MPI_Allreduce(&num, &allNumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&denum, &allDenum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //Вычисляем коэффициент тау
    coefTau = allNumerator/allDenum;
    return 100000;
}

void ConjugateGradientsMethod_MPI::getGridBorders(Grid &grid) {
    //Инициализируем вспомогательные массивы для пересылки и получения граничных с другими узлами значений
    MPI_Status status[4];
    MPI_Request send[4];
    MPI_Request recv[4];
    tr1::shared_ptr<double> topBuf, bottomBuf, rightBuf, leftBuf;
    // Выполняем пересылку и получение с 4-х соседних узлов если они есть
    if (top >= 0 && top < size) {
        topBuf = tr1::shared_ptr<double>(new double[grid.getColumnsCount()]);
        Vector topCur = grid.getTopRow();
        MPI_Isend(topCur.memptr(), topCur.size(), MPI_DOUBLE, top, 0, MPI_COMM_WORLD, &send[0]);
        MPI_Irecv(topBuf.get(), grid.getColumnsCount(), MPI_DOUBLE, top, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[0]);

        MPI_Wait(&send[0],&status[0]);
        MPI_Wait(&recv[0],&status[0]);
    }
    if(bottom >= 0 && bottom < size) {
        bottomBuf = tr1::shared_ptr<double> (new double[grid.getColumnsCount()]);
        Vector bottomCur = grid.getBottomRow();
        MPI_Isend(bottomCur.memptr(), bottomCur.size(), MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD, &send[1]);
        MPI_Irecv(bottomBuf.get(), grid.getColumnsCount(), MPI_DOUBLE, bottom, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[1]);

        MPI_Wait(&send[1],&status[1]);
        MPI_Wait(&recv[1],&status[1]);
    }
    if (right >= 0 && right < size) {
        rightBuf = tr1::shared_ptr<double> (new double[grid.getRowsCount()]);
        Vector rightCur = grid.getRightCol();
        MPI_Isend(rightCur.memptr(), rightCur.size(), MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &send[2]);
        MPI_Irecv(rightBuf.get(), grid.getRowsCount(), MPI_DOUBLE, right, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[2]);

        MPI_Wait(&send[2],&status[2]);
        MPI_Wait(&recv[2],&status[2]);
    }
    if(left >= 0 && left < size) {
        leftBuf = tr1::shared_ptr<double> (new double[grid.getRowsCount()]);
        Vector leftCur = grid.getLeftCol();
        MPI_Isend(leftCur.memptr(), leftCur.size(), MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &send[3]);
        MPI_Irecv(leftBuf.get(), grid.getRowsCount(), MPI_DOUBLE, left, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[3]);

        MPI_Wait(&send[3],&status[3]);
        MPI_Wait(&recv[3],&status[3]);
    }
    if (topBuf) {
        Vector topR(topBuf, grid.getColumnsCount());
        grid.addTopBorder(topR);
    }
    if (bottomBuf) {
        Vector bottomR(bottomBuf, grid.getColumnsCount());
        grid.addBottomBorder(bottomR);
    }
    if (rightBuf) {
        Vector rightR(rightBuf, grid.getRowsCount());
        grid.addRightBorder(rightR);
    }
    if (leftBuf) {
        Vector leftR(leftBuf, grid.getRowsCount());
        grid.addLeftBorder(leftR);
    }
}

double ConjugateGradientsMethod_MPI::iterate(Grid &pGrid) {
    if(iterationsCount == 0){
        iterationsCount++;
        return SteepestDescentMethodForZeroIter(pGrid);
    }

    long rows = pGrid.getRowsCount();
    long cols = pGrid.getColumnsCount();
    double err = 0;

    //Рассчитываем значения p на текущей итерации
    //Раачитываем значение ошибки (невязки) между текщим p и p на предыдущей итерации для учета условия выхода
    for (size_t i = 0; i < rows ; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (checkBorder(rGrid, i, j)) {
                double val = pGrid(i,j) - coefTau*gGrid(i,j);
                double errV = pGrid(i,j) - val;
                PointD hs = pGrid.getAverageGridSteps(i,j);
                err += errV*errV*hs.first*hs.second;
                pGrid(i,j) = val;
            }
        }
    }
    //Рассчитываем значения r на текущей итерации
    getGridBorders(pGrid);
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if (checkBorder(rGrid, i, j)) {
                rGrid(i,j) = fiveDotScheme(pGrid, i, j) - func(pGrid.getPointFromCache(i,j));
            }
        }
    }

    //Рассчитываем скалярные произведения для рассчета коэффициента альфа
    getGridBorders(rGrid);
    double anum = 0, adenum = 0;
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if(checkBorder(rGrid, i, j)) {
                anum += fiveDotScheme(rGrid, i, j) * gGrid(i,j);
                adenum += fiveDotScheme(gGrid, i, j) * gGrid(i,j);
            }
        }
    }
    double allAnum, allAdenum;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&anum, &allAnum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&adenum, &allAdenum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double alpha = allAnum / allAdenum;


    //Рассчитываем значения g на текущей итерации   
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if(checkBorder(gGrid, i, j)) {
                gGrid(i,j) = rGrid(i,j) - alpha*gGrid(i,j);
            }
        }
    }

    //Рассчитываем скалярные произведения для рассчета коэффициента тау
    getGridBorders(gGrid);
    double tnum = 0, tdenum = 0;
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if(checkBorder(gGrid, i, j)) {
                tnum += gGrid(i,j) * rGrid(i,j);
                tdenum += fiveDotScheme(gGrid, i, j) * gGrid(i, j);
            }
        }
    }
    double allTnum, allTdenum;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&tnum, &allTnum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tdenum, &allTdenum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    coefTau = allTnum/ allTdenum;

    //Меняем номер итерации и выходим
    iterationsCount++;
    //Получаем количесво узлов
    int total_size;
    MPI_Comm_size(MPI_COMM_WORLD, &total_size);
    
    //Получаем общее значение ошибки для учета условия выхода из алгоритма
    double total_error = 0;
    MPI_Allreduce(&err, &total_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(total_error);
}
//===============================================================

//MPI / OpenMP версия метода сопряженных градиентов
double ConjugateGradientsMethod_MPI_OpenMP::SteepestDescentMethodForZeroIter(Grid &pGrid) {

    getGridBorders(pGrid);
    #pragma omp parallel for
    for (long i = 0; i < pGrid.getRowsCount(); ++i) {
        for (long j = 0; j < pGrid.getColumnsCount(); ++j) {
            if(checkBorder(rGrid, i, j)) {
                rGrid(i,j) = fiveDotScheme(pGrid, i, j) - func(pGrid.getPointFromCache(i,j));
                gGrid(i,j) = rGrid(i,j);
            }
        }
    }

    getGridBorders(rGrid);
    getGridBorders(gGrid);
    double num = 0, denum = 0;
    #pragma omp parallel for shared(num, denum)
    for (long i = 0; i < pGrid.getRowsCount(); ++i) {
        for (long j = 0; j < pGrid.getColumnsCount(); ++j) {
            if (checkBorder(pGrid, i, j)) {
                num += rGrid(i,j) * rGrid(i,j);
                denum += fiveDotScheme(rGrid,i,j) * rGrid(i,j);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double allNumerator, allDenum;
    MPI_Allreduce(&num, &allNumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&denum, &allDenum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    coefTau = allNumerator/allDenum;
    return 100000;
}

double ConjugateGradientsMethod_MPI_OpenMP::iterate(Grid &pGrid) {
    if(iterationsCount == 0){
        iterationsCount++;
        return SteepestDescentMethodForZeroIter(pGrid);
    }

    long rows = pGrid.getRowsCount();
    long cols = pGrid.getColumnsCount();
    double err = 0;

    #pragma omp parallel for reduction(+:err)
    for (long i = 0; i < rows ; ++i) {
        for (long j = 0; j < cols; ++j) {
            if (checkBorder(rGrid, i, j)) {
                double val = pGrid(i,j) - coefTau*gGrid(i,j);
                double errV = pGrid(i,j) - val;
                PointD hs = pGrid.getAverageGridSteps(i,j);
                err += errV*errV*hs.first*hs.second;
                pGrid(i,j) = val;
            }
        }
    }

    getGridBorders(pGrid);
    #pragma omp parallel for
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if (checkBorder(rGrid, i, j)) {
                rGrid(i,j) = fiveDotScheme(pGrid, i, j) - func(pGrid.getPointFromCache(i,j));
            }
        }
    }
    getGridBorders(rGrid);
    double anum = 0, adenum = 0;

    #pragma omp parallel for reduction(+:anum) reduction(+:adenum)
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if(checkBorder(rGrid, i, j)) {
                anum += fiveDotScheme(rGrid, i, j) * gGrid(i,j);
                adenum += fiveDotScheme(gGrid, i, j) * gGrid(i,j);
            }
        }
    }
    double allAnum, allAdenum;

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&anum, &allAnum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&adenum, &allAdenum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double alpha = allAnum / allAdenum;

    #pragma omp parallel for
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if(checkBorder(gGrid, i, j)) {
                gGrid(i,j) = rGrid(i,j) - alpha*gGrid(i,j);
            }
        }
    }

    getGridBorders(gGrid);
    double tnum = 0, tdenum = 0;

    #pragma omp parallel for reduction(+:tnum) reduction(+:tdenum)
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if(checkBorder(gGrid, i, j)) {
                tnum += gGrid(i,j) * rGrid(i,j);
                tdenum += fiveDotScheme(gGrid, i, j) * gGrid(i, j);
            }
        }
    }

    double allTnum, allTdenum;

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&tnum, &allTnum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tdenum, &allTdenum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    coefTau = allTnum/ allTdenum;
    iterationsCount++;
    int total_size;
    MPI_Comm_size(MPI_COMM_WORLD, &total_size);
    double total_error = 0;
    MPI_Allreduce(&err, &total_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(total_error);
}
//===============================================================



