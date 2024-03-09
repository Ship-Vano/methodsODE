/*<<Методы решений ОДУ>>
 @author: Шаманов Иван
 @author: Климов Олег
 2024*/


#pragma once
#include<iostream>
#include<map>
#include"LinOp.h"
#include"FileIO.h"
#include"NLEmethods.h"
using namespace std;

/* Метод Эйлера Явный */
template<typename DT, typename F>
bool EulerSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "eulersimple-1");


/* Метод Эйлера Неявный */
template <typename DT, typename F>
bool ImplicitEulerSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "eulerImplicit-1");


/* Метод Рунге-Кутты двустадийный */
template <typename DT, typename F>
bool RungeKutta2Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "rungekutta2-1");


/* Метод Рунге-Кутты двустадийный с Автоматическим шагом */
template <typename DT, typename F>
bool RungeKutta2SolveAuto(F func, DT t_0, DT T, DT tau0, vector<DT> u_0, string filename = "rungekutta2-1");


/* Метод Рунге-Кутты четырехстадийный*/
template <typename DT, typename F>
bool RungeKutta4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "rungekutta4-1");


/* Метод Рунге-Кутты четырехстадийныйс Автоматическим шагом */
template <typename DT, typename F>
bool RungeKutta4SolveAuto(F func, DT t_0, DT T, DT tau0, vector<DT> u_0, string filename = "rungekutta4-1");


/* Метод Симметричной схемы */
template<typename DT, typename F>
bool SimSchemeSolve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "simScheme-1");


/* Метод Адамса-Башфорта */
template<typename DT,  typename F>
bool Adams4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "adams4-1");


/* Метор Прогноз-Коррекция */
template<typename DT, typename F>
bool PredictCorrect4Solve(F func, DT t_0, DT T, DT tau, vector<DT> u_0, string filename = "predict4-1");
