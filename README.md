# ConjugateGradientsMethod
## Построение графиков
Построение графиков производится с помощью системы Matlab, соотвеетствующий скрипт SurfSk.m
## Сборка
* make sequential - последовательная версия
* make parallel - MPI версия
* make parallel_omp - MPI/OpenMP версия
* make all - собрать все

## Запуск
### Аргументы
-x int -- нижний левый угол сетки ось X

-y int  -- нижний левый угол сетки ось Y

-a int -- верхний правый угол сетки ось X

-b int -- верхний правый угол сетки ось Y

-r int -- число строк в сетке

-c int -- число столбцов в сетке

-f string -- выходной файл

Функции F и phi задаются в файле InputFunctions.cpp и далее используются как callback-функции

### На BlueGene
mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 1 -w 02:00:00 --stdout ./resultsMP/rezultMP-1-1000.out ./sequential -- -x 0 -y 0 -a 2 -b 2 -r 1000 -c 1000 -f ./resultsMP/rezult_Par_1_1000.txt

mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 128 -w 00:15:00 --stdout ./resultsMP/rezultMP-128-1000.out ./parallel -- -x 0 -y 0 -a 2 -b 2 -r 1000 -c 1000 -f ./resultsMP/rezult_Par_128_1000.txt

mpisubmit.bg --env BG_SHAREDMEMPOOLSIZE=256 -n 128 -w 00:15:00 --stdout ./resultsOMP/rezultOMP-128-1000.out ./parallel_omp -- -x 0 -y 0 -a 2 -b 2 -r 1000 -c 1000 -f ./resultsOMP/rezult_ParOMP_128_1000.txt


### На СК Ломоносов
sbatch -p test -n 4 impi ./parallel -x 0 -y 0 -a 2 -b 2 -r 1000 -c 1000 -f ./resultsMP/rezult_Par_4_1000.txt
