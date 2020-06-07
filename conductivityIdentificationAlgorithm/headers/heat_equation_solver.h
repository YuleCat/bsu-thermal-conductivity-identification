#pragma once

#include <iostream>

#include "math_consts.h"
#include "distance_helper.h"
#include "material_info.h"
#include "material_type.h"
#include "materilas_helper.h"
#include "math_options.h"
#include "mpi_options.h"

class HeatEquationSolver
{
private:
    int n;
    int m;

    float t_gr3;
    float t_gr4;

    float initTime;
    float endTime;

    MaterialsHelper materialsHelper;

    DistanceHelper xDistanceHelper;
    DistanceHelper zDistanceHelper;

    bool useMPI;
    int rank = 0;
    int size = 1;
    std::pair<int, int>* blockIntervals;
    int tapeSize;
    int tapesNum;

    float alpha_x[N_MAX][M_MAX];
    float beta_x[N_MAX][M_MAX];

    float alpha_z[M_MAX];
    float beta_z[M_MAX];

    // Значения температуры для предыдущего значения времени
    float T_time[N_MAX][M_MAX];
    // Значения температуры для предыдущего шага
    float T_step[N_MAX][M_MAX];
    // Значения температуры на текущем шаге
    float T_cur[N_MAX][M_MAX];

    float conductivity;

public:
    HeatEquationSolver(MathOptions mathOptions, MPIOptions mpiOptions);
    void Solve();
    float GetResult();
    ~HeatEquationSolver();

private:
    void InitializeDistanceHelpers(float xStep, float zStep);
    DistanceHelper CreateDistanceHelper(int size, float step);
    void InitializeTemperature(float** initTemperature);
    void InitializeMaterialsHelper(MaterialType** gridMaterialInfo);
    void InitializeMPIOptions(MPIOptions mpiOptions);

    void ReassignTimeTemperatures();

    // Метод переменных направлений
    void AlternatingDirectionMethod();
    bool IsCriterionFulfilled();
    void ReassignStepTemperatures();

    // Метод прогонки
    void TridiagonalMatrixAlgorithm();
    void ParallelTridiagonalMatrixAlgorithm();

    void ShareBorderTemperatures();
    void ComputeXDirection();
    void ComputeZDirection();
    void ProcessCoefficients();
    void ComputeXCoefficients(int i, int k);
    void ProcessTemperatures();
    void ComputeXTemperature(int i, int k);

    void ComputeConductivity();

	float Ksi_1_x(int k);
	float Chi_1_x(int k);
	float Ksi_2_x(int k);
	float Chi_2_x(int k);

	float A_x(int i, int k);
	float B_x(int i, int k);
	float C_x(int i, int k);
    float D_x(int i, int k);
    float A_z(int i, int k);
    float B_z(int i, int k);
    float C_z(int i, int k);
    float D_z(int i, int k);

    std::pair<int, int> GetTapeInterval(int t);
};