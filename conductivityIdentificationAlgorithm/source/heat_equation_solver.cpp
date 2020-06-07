#include "heat_equation_solver.h"

#include <algorithm>
#include <iostream>

#include "mpi.h"

#include "math_consts.h"
#include "mpi_consts.h"
#include "material_type.h"
#include "material_info.h"
#include "math_options.h"
#include "mpi_options.h"

HeatEquationSolver::HeatEquationSolver(MathOptions mathOptions, MPIOptions mpiOptions)
	:
	n(mathOptions.n),
	m(mathOptions.m),
	t_gr3(mathOptions.t_gr3),
	t_gr4(mathOptions.t_gr4),
	initTime(mathOptions.initTime),
	endTime(mathOptions.endTime)
{
	InitializeMPIOptions(mpiOptions);

	InitializeDistanceHelpers(mathOptions.xStep, mathOptions.zStep);
	InitializeMaterialsHelper(mathOptions.gridMaterialInfo);
	InitializeTemperature(mathOptions.initTemperature);
}

void HeatEquationSolver::Solve()
{
	for (int i = initTime; i <= endTime; i += TIME_DELTA)
	{
		if (rank == 0)
		{
			//std::cout << "Time: " << i << std::endl;
		}

		AlternatingDirectionMethod();

		if (rank == 0)
		{
			//std::cout << std::endl;
		}

		ReassignTimeTemperatures();
	}

	ComputeConductivity();
}

void HeatEquationSolver::ComputeConductivity()
{
	int lowerBound = (!useMPI) ? 1 : blockIntervals[rank].first;
	int upperBound = (!useMPI) ? n : blockIntervals[rank].second;

	float heatFlow = 0;
	float x = 0;
	for (int i = lowerBound; i <= upperBound; ++i)
	{
		for (int j = 1; j <= m; ++j)
		{
			x = zDistanceHelper.GetSize() / zDistanceHelper.GetDeltaLine(j + 1);
			heatFlow += materialsHelper.GetUpperCornerConductivity(i, j) * MULTIPLIER * (T_time[i][j] - T_time[i][j - 1]) * x;
		}
	}

	if (useMPI)
	{
		if (rank != 0)
		{
			MPI_Send(&heatFlow, 1, MPI_FLOAT, 0, RECV_RESULT_TAG, MPI_COMM_WORLD);
		}
		else
		{
			float curSum = 0;
			for (int i = 1; i < size; ++i)
			{
				MPI_Recv(&curSum, 1, MPI_FLOAT, MPI_ANY_SOURCE, RECV_RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				heatFlow += curSum;
			}
		}
	}

	if (rank == 0)
	{
		conductivity = heatFlow / (t_gr4 - t_gr3);
	}
}

float HeatEquationSolver::GetResult()
{
	return conductivity;
}

HeatEquationSolver::~HeatEquationSolver()
{
	delete[] blockIntervals;
}

void HeatEquationSolver::InitializeDistanceHelpers(float xStep, float zStep)
{
	xDistanceHelper = CreateDistanceHelper(n, xStep);
	zDistanceHelper = CreateDistanceHelper(m, zStep);
}

DistanceHelper HeatEquationSolver::CreateDistanceHelper(int size, float step)
{
	float* gridValues = new float[size + 1];
	gridValues[0] = 0;
	for (int i = 1; i < size + 1; ++i)
	{
		gridValues[i] = gridValues[i - 1] + step;
	}
	DistanceHelper helper = DistanceHelper(size, gridValues);

	delete[] gridValues;

	return helper;
}

void HeatEquationSolver::InitializeTemperature(float** initTemperature)
{
	int lowerBound = (rank == 0 || !useMPI) ? 0 : blockIntervals[rank].first;
	int upperBound = (rank == size - 1 || !useMPI) ? n + 1 : blockIntervals[rank].second;

	for (int i = lowerBound; i <= upperBound; ++i)
	{
		for (int j = 0; j <= m + 1; ++j)
		{
			T_time[i][j] = initTemperature[i][j];
			T_step[i][j] = initTemperature[i][j];
		}
	}

	for (int i = std::max(lowerBound - 1, 0); i <= std::min(upperBound + 1, n); ++i)
	{
		T_time[i][0] = t_gr3;
		T_step[i][0] = t_gr3;

		T_time[i][m + 1] = t_gr4;
		T_step[i][m + 1] = t_gr4;
	}
}

void HeatEquationSolver::InitializeMaterialsHelper(MaterialType** gridMaterialInfo)
{
	materialsHelper = MaterialsHelper(n, m, gridMaterialInfo);
}

void HeatEquationSolver::InitializeMPIOptions(MPIOptions mpiOptions)
{
	useMPI = mpiOptions.useMPI;
	if (useMPI)
	{
		rank = mpiOptions.rank;
		size = mpiOptions.size;

		blockIntervals = new std::pair<int, int>[size];
		
		int blockSize = floor(1.f * n / size) + ((0 < (n % size)) ? 1 : 0);
		blockIntervals[0] = std::make_pair(1, blockSize);
		for (int i = 1; i < size; ++i)
		{
			blockSize = floor(1.f * n / size) + ((i < (n % size)) ? 1 : 0);
			int prevEnd = blockIntervals[i - 1].second;
			blockIntervals[i] = std::make_pair(prevEnd + 1, prevEnd + blockSize);
		}
		blockIntervals[size - 1].second = n;

		tapeSize = mpiOptions.tapeSize;
		tapesNum = ceil(1.f * m / tapeSize);
	}
}

void HeatEquationSolver::ReassignTimeTemperatures()
{
	int lowerBound = (rank == 0 || !useMPI) ? 0 : blockIntervals[rank].first;
	int upperBound = (rank == size - 1 || !useMPI) ? n + 1 : blockIntervals[rank].second;
	// Присваиваем T_step значения текущего шага
	for (int i = lowerBound; i <= upperBound; ++i)
	{
		for (int j = 0; j <= m + 1; ++j)
		{
			T_time[i][j] = T_step[i][j];
		}
	}
}

// Метод переменных направлений
void HeatEquationSolver::AlternatingDirectionMethod()
{
	bool isFulfilled = false;
	int iterationCount = 1;
	//while (!isFulfilled)
	for (int i = 0; i < 20; ++i)
	{
		if (rank == 0)
		{
			//std::cout << "Iteration " << iterationCount << ": ";
		}

		if (useMPI)
		{
			ParallelTridiagonalMatrixAlgorithm();
		}
		else
		{
			TridiagonalMatrixAlgorithm();
		}

		// Критерий сходимости

		if (rank == 0)
		{
			//std::cout << "Done; Criterion: ";
		}

		isFulfilled = IsCriterionFulfilled();

		if (rank == 0)
		{
			//std::cout << "Done." << std::endl;
		}

		ReassignStepTemperatures();

		iterationCount++;
	}
}

bool HeatEquationSolver::IsCriterionFulfilled()
{
	bool isFulfilled = false;

	if (useMPI)
	{
		int lowerBound = (rank == 0) ? 0 : blockIntervals[rank].first;
		int upperBound = (rank == size - 1) ? n + 1 : blockIntervals[rank].second;

		bool isGreaterEpsilon = false;
		for (int i = lowerBound; i <= upperBound; ++i)
		{
			for (int j = 0; j <= m + 1; ++j)
			{
				double diff = (T_cur[i][j] - T_step[i][j]) / T_time[i][j];
				if (abs(diff) >= EPSILON)
				{
					isGreaterEpsilon = true;
					break;
				}
			}

			if (isGreaterEpsilon)
			{
				break;
			}
		}

		// Если критерий выполнился, выходим из цикла
		isFulfilled = !isGreaterEpsilon;

		if (rank != 0)
		{
			MPI_Send(&isFulfilled, 1, MPI_BYTE, 0, CRITERION_TAG, MPI_COMM_WORLD);
		}
		else
		{
			bool curCriterion;
			for (int i = 1; i < size; ++i)
			{
				MPI_Recv(&curCriterion, 1, MPI_BYTE, MPI_ANY_SOURCE, CRITERION_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				isFulfilled &= curCriterion;
			}
		}

		MPI_Bcast(&isFulfilled, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
	}
	else
	{
		bool isGreaterEpsilon = false;
		for (int i = 0; i <= n + 1; ++i)
		{
			for (int j = 0; j <= m + 1; ++j)
			{
				double diff = (T_cur[i][j] - T_step[i][j]) / T_time[i][j];
				if (abs(diff) >= EPSILON)
				{
					isGreaterEpsilon = true;
					break;
				}
			}

			if (isGreaterEpsilon)
			{
				break;
			}
		}

		// Если критерий выполнился, выходим из цикла
		isFulfilled = !isGreaterEpsilon;
	}

	return isFulfilled;
}

void HeatEquationSolver::ReassignStepTemperatures()
{
	int lowerBound = (rank == 0 || !useMPI) ? 0 : blockIntervals[rank].first;
	int upperBound = (rank == size - 1 || !useMPI) ? n + 1 : blockIntervals[rank].second;
	// Присваиваем T_step значения текущего шага
	for (int i = lowerBound; i <= upperBound; ++i)
	{
		for (int j = 0; j <= m + 1; ++j)
		{
			T_step[i][j] = T_cur[i][j];
		}
	}
}

// Метод прогонки
void HeatEquationSolver::TridiagonalMatrixAlgorithm()
{
	for (int k = 1; k <= m; ++k)
	{
		alpha_x[1][k] = Chi_1_x(k);
		beta_x[1][k] = Ksi_1_x(k);
		for (int i = 1; i <= n; ++i)
		{
			float a = C_x(i, k);
			float b = B_x(i, k);
			float c = A_x(i, k) + B_x(i, k) + C_x(i, k);
			float f = D_x(i, k) + A_x(i, k) * T_time[i][k];

			alpha_x[i + 1][k] = b / (c - alpha_x[i][k] * a);
			beta_x[i + 1][k] = (a * beta_x[i][k] + f) / (c - alpha_x[i][k] * a);
		}

		T_cur[n + 1][k] = (Ksi_2_x(k) + Chi_2_x(k) * beta_x[n + 1][k]) / (1 - Chi_2_x(k) * alpha_x[n + 1][k]);
		for (int i = n; i >= 0; --i)
		{
			T_cur[i][k] = alpha_x[i + 1][k] * T_cur[i + 1][k] + beta_x[i + 1][k];
		}
	}
	
	for (int i = 1; i <= n; ++i)
	{
		alpha_z[1] = 0;
		beta_z[1] = t_gr3;
		for (int k = 1; k <= m; ++k)
		{
			float a = D_z(i, k);
			float b = C_z(i, k);
			float c = C_z(i, k) + D_z(i, k) + A_z(i, k);
			float f = B_z(i, k) + A_z(i, k) * T_cur[i][k];

			alpha_z[k + 1] = b / (c - alpha_z[k] * a);
			beta_z[k + 1] = (a * beta_z[k] + f) / (c - alpha_z[k] * a);
		}

		T_cur[i][m + 1] = t_gr4;
		for (int k = m; k > 0; --k)
		{
			T_cur[i][k] = alpha_z[k + 1] * T_cur[i][k + 1] + beta_z[k + 1];
		}
		T_cur[i][0] = t_gr3;
	}
}

// Параллельный метод прогонки
void HeatEquationSolver::ParallelTridiagonalMatrixAlgorithm()
{
	ShareBorderTemperatures();

	ComputeXDirection();
	ComputeZDirection();
}

void HeatEquationSolver::ShareBorderTemperatures()
{
	MPI_Request sendBottomStepRequest;
	if (rank != 0)
	{
		MPI_Isend(&T_step[blockIntervals[rank].first][0], M_MAX, MPI_FLOAT, rank - 1, SEND_BOTTOM_STEP_TAG, MPI_COMM_WORLD, &sendBottomStepRequest);
	}
	if (rank != size - 1)
	{
		MPI_Irecv(&T_step[blockIntervals[rank + 1].first][0], M_MAX, MPI_FLOAT, rank + 1, SEND_BOTTOM_STEP_TAG, MPI_COMM_WORLD, &sendBottomStepRequest);
	}

	MPI_Request sendUpperStepRequest;
	MPI_Request sendUpperTimeRequest;
	if (rank != size - 1)
	{
		MPI_Isend(&T_step[blockIntervals[rank].second][0], M_MAX, MPI_FLOAT, rank + 1, SEND_UPPER_STEP_TAG, MPI_COMM_WORLD, &sendUpperStepRequest);
		MPI_Isend(&T_time[blockIntervals[rank].second][0], M_MAX, MPI_FLOAT, rank + 1, SEND_UPPER_TIME_TAG, MPI_COMM_WORLD, &sendUpperTimeRequest);
	}
	if (rank != 0)
	{
		MPI_Irecv(&T_step[blockIntervals[rank - 1].second][0], M_MAX, MPI_FLOAT, rank - 1, SEND_UPPER_STEP_TAG, MPI_COMM_WORLD, &sendUpperStepRequest);
		MPI_Irecv(&T_time[blockIntervals[rank - 1].second][0], M_MAX, MPI_FLOAT, rank - 1, SEND_UPPER_TIME_TAG, MPI_COMM_WORLD, &sendUpperTimeRequest);
	}

	MPI_Wait(&sendBottomStepRequest, MPI_STATUS_IGNORE);
	MPI_Wait(&sendUpperStepRequest, MPI_STATUS_IGNORE);
	MPI_Wait(&sendUpperTimeRequest, MPI_STATUS_IGNORE);
}

void HeatEquationSolver::ComputeXDirection()
{
	ProcessCoefficients();
	ProcessTemperatures();
}

void HeatEquationSolver::ComputeZDirection()
{
	for (int i = blockIntervals[rank].first; i <= blockIntervals[rank].second; ++i)
	{
		alpha_z[1] = 0;
		beta_z[1] = t_gr3;
		for (int k = 1; k <= m; ++k)
		{
			float a = D_z(i, k);
			float b = C_z(i, k);
			float c = C_z(i, k) + D_z(i, k) + A_z(i, k);
			float f = B_z(i, k) + A_z(i, k) * T_cur[i][k];

			alpha_z[k + 1] = b / (c - alpha_z[k] * a);
			beta_z[k + 1] = (a * beta_z[k] + f) / (c - alpha_z[k] * a);
		}

		T_cur[i][m + 1] = t_gr4;
		for (int k = m; k > 0; --k)
		{
			T_cur[i][k] = alpha_z[k + 1] * T_cur[i][k + 1] + beta_z[k + 1];
		}
		T_cur[i][0] = t_gr3;
	}
}

void HeatEquationSolver::ProcessCoefficients()
{
	MPI_Request sendAlphaRequest;
	MPI_Request sendBetaRequest;
	MPI_Request recvAlphaRequest;
	MPI_Request recvBetaRequest;

	int recvTapeSize;
	int nextRecvTapeSize;
	bool isFirstRecieve = true;
	for (int t = 0; t < tapesNum; ++t)
	{
		std::pair<int, int> tapeInterval = GetTapeInterval(t);
		std::pair<int, int> nextTapeInterval = GetTapeInterval(t + 1);

		recvTapeSize = tapeInterval.second - tapeInterval.first + 1;
		nextRecvTapeSize = nextTapeInterval.second - nextTapeInterval.first + 1;

		if (rank != 0)
		{
			int recvRow = blockIntervals[rank - 1].second;
			
			if (isFirstRecieve)
			{
				MPI_Recv(&alpha_x[recvRow][tapeInterval.first], recvTapeSize, MPI_FLOAT, rank - 1, SEND_ALPHA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&beta_x[recvRow][tapeInterval.first], recvTapeSize, MPI_FLOAT, rank - 1, SEND_BETA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				isFirstRecieve = false;

				if (t != tapesNum - 1)
				{
					MPI_Irecv(&alpha_x[recvRow][nextTapeInterval.first], nextRecvTapeSize, MPI_FLOAT, rank - 1, SEND_ALPHA_TAG, MPI_COMM_WORLD, &recvAlphaRequest);
					MPI_Irecv(&beta_x[recvRow][nextTapeInterval.first], nextRecvTapeSize, MPI_FLOAT, rank - 1, SEND_BETA_TAG, MPI_COMM_WORLD, &recvBetaRequest);
				}
			}
			else
			{
				MPI_Wait(&recvAlphaRequest, MPI_STATUS_IGNORE);
				MPI_Wait(&recvBetaRequest, MPI_STATUS_IGNORE);

				if (t != tapesNum - 1)
				{
					MPI_Irecv(&alpha_x[recvRow][nextTapeInterval.first], nextRecvTapeSize, MPI_FLOAT, rank - 1, SEND_ALPHA_TAG, MPI_COMM_WORLD, &recvAlphaRequest);
					MPI_Irecv(&beta_x[recvRow][nextTapeInterval.first], nextRecvTapeSize, MPI_FLOAT, rank - 1, SEND_BETA_TAG, MPI_COMM_WORLD, &recvBetaRequest);
				}
			}
		}

		for (int k = tapeInterval.first; k <= tapeInterval.second; ++k)
		{
			for (int i = blockIntervals[rank].first; i <= blockIntervals[rank].second; ++i)
			{
				ComputeXCoefficients(i, k);
			}

			if (rank == size - 1)
			{
				ComputeXCoefficients(n + 1, k);
			}
		}

		if (rank != size - 1)
		{
			if (t != 0)
			{
				MPI_Wait(&sendAlphaRequest, MPI_STATUS_IGNORE);
				MPI_Wait(&sendBetaRequest, MPI_STATUS_IGNORE);
			}

			MPI_Isend(&alpha_x[blockIntervals[rank].second][tapeInterval.first], recvTapeSize, MPI_FLOAT, rank + 1, SEND_ALPHA_TAG, MPI_COMM_WORLD, &sendAlphaRequest);
			MPI_Isend(&beta_x[blockIntervals[rank].second][tapeInterval.first], recvTapeSize, MPI_FLOAT, rank + 1, SEND_BETA_TAG, MPI_COMM_WORLD, &sendBetaRequest);
		}
	}
}

void HeatEquationSolver::ComputeXCoefficients(int i, int k)
{
	if (i == 1)
	{
		alpha_x[i][k] = Chi_1_x(k);
		beta_x[i][k] = Ksi_1_x(k);
	}
	else
	{
		float a = C_x(i - 1, k);
		float b = B_x(i - 1, k);
		float c = A_x(i - 1, k) + B_x(i - 1, k) + C_x(i - 1, k);
		float f = D_x(i - 1, k) + A_x(i - 1, k) * T_time[i - 1][k];

		alpha_x[i][k] = b / (c - alpha_x[i - 1][k] * a);
		beta_x[i][k] = (a * beta_x[i - 1][k] + f) / (c - alpha_x[i - 1][k] * a);
	}
}

void HeatEquationSolver::ProcessTemperatures()
{
	MPI_Request sendCurRequest;
	MPI_Request recvCurRequest;

	bool isFirstRecieve = true;
	int recvTapeSize;
	int nextRecvTapeSize;
	for (int t = 0; t < tapesNum; ++t)
	{
		std::pair<int, int> tapeInterval = GetTapeInterval(t);
		std::pair<int, int> nextTapeInterval = GetTapeInterval(t + 1);

		recvTapeSize = tapeInterval.second - tapeInterval.first + 1;
		nextRecvTapeSize = nextTapeInterval.second - nextTapeInterval.first + 1;

		if (rank != size - 1)
		{
			int recvRow = blockIntervals[rank - 1].second;

			if (isFirstRecieve)
			{
				MPI_Recv(&T_cur[blockIntervals[rank].second][tapeInterval.first], recvTapeSize, MPI_FLOAT, rank + 1, SEND_CUR_TEMPERATURE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				isFirstRecieve = false;

				if (t != tapesNum - 1)
				{
					MPI_Irecv(&T_cur[blockIntervals[rank].second][nextTapeInterval.first], nextRecvTapeSize, MPI_FLOAT, rank + 1, SEND_CUR_TEMPERATURE_TAG, MPI_COMM_WORLD, &recvCurRequest);
				}
			}
			else
			{
				MPI_Wait(&recvCurRequest, MPI_STATUS_IGNORE);

				if (t != tapesNum - 1)
				{
					MPI_Irecv(&T_cur[blockIntervals[rank].second][nextTapeInterval.first], nextRecvTapeSize, MPI_FLOAT, rank + 1, SEND_CUR_TEMPERATURE_TAG, MPI_COMM_WORLD, &recvCurRequest);
				}
			}
		}

		int upperBound = (rank == size - 1) ? n + 1 : blockIntervals[rank].second - 1;
		int lowerBound = (rank == 0) ? 0 : blockIntervals[rank - 1].second;

		for (int k = tapeInterval.first; k <= tapeInterval.second; ++k)
		{
			for (int i = upperBound; i >= lowerBound; --i)
			{
				ComputeXTemperature(i, k);
			}
		}

		if (rank != 0)
		{
			if (t != 0)
			{
				MPI_Wait(&sendCurRequest, MPI_STATUS_IGNORE);
			}

			MPI_Isend(&T_cur[blockIntervals[rank - 1].second][tapeInterval.first], recvTapeSize, MPI_FLOAT, rank - 1, SEND_CUR_TEMPERATURE_TAG, MPI_COMM_WORLD, &sendCurRequest);
		}
	}
}

void HeatEquationSolver::ComputeXTemperature(int i, int k)
{
	if (i == n + 1)
	{
		T_cur[n + 1][k] = (Ksi_2_x(k) + Chi_2_x(k) * beta_x[n + 1][k]) / (1 - Chi_2_x(k) * alpha_x[n + 1][k]);
	}
	else
	{
		T_cur[i][k] = alpha_x[i + 1][k] * T_cur[i + 1][k] + beta_x[i + 1][k];
	}
}

float HeatEquationSolver::Ksi_1_x(int k)
{
	MaterialInfo info = materialsHelper.GetInfo(0, k);

	double c = 0;
	if (k == 1)
	{
		double first = materialsHelper.GetUpperCornerConductivity(0, 1) * (T_step[0][2] - T_step[0][1]) / zDistanceHelper.GetDelta(k + 1);
		double second = info.conductivity * (T_step[0][1] - T_step[0][0]) / 0.5 / zDistanceHelper.GetDeltaLine(k);

		c = (first - second) / (0.75 * zDistanceHelper.GetDelta(k));
	}
	else if (k >= 2 && k <= m - 1)
	{
		double first = materialsHelper.GetUpperCornerConductivity(0, k + 1) * (T_step[0][k + 1] - T_step[0][k]) / zDistanceHelper.GetDeltaLine(k + 1);
		double second = materialsHelper.GetUpperCornerConductivity(0, k) * (T_step[0][k] - T_step[0][k - 1]) / zDistanceHelper.GetDeltaLine(k);

		c = (first - second) / zDistanceHelper.GetDelta(k);
	}
	else if (k == m)
	{
		double first = materialsHelper.GetUpperCornerConductivity(0, k) * (T_step[0][m + 1] - T_step[0][m]) / 0.5 / zDistanceHelper.GetDeltaLine(m + 1);
		double second = materialsHelper.GetBottomCornerConductivity(0, k) * (T_step[0][m] - T_step[0][m - 1]) / zDistanceHelper.GetDeltaLine(m);

		c = (first - second) / (0.75 * zDistanceHelper.GetDelta(k));
	}

	double a = info.capacity * info.density / (0.5 * TIME_DELTA);

	return c / a + T_time[0][k];
}

float HeatEquationSolver::Chi_1_x(int k)
{
	MaterialInfo info = materialsHelper.GetInfo(0, k);

	double b = 0;
	if (k == 1)
	{
		b = info.conductivity / 0.5 / 0.25 / xDistanceHelper.GetDelta(1) / xDistanceHelper.GetDelta(1);
	}
	else if (k >= 2 && k <= m - 1)
	{
		b = materialsHelper.GetRightCornerConductivity(0, k) / 0.5 / 0.25 / xDistanceHelper.GetDelta(1) / xDistanceHelper.GetDelta(1);
	}
	else if (k == m)
	{
		b = materialsHelper.GetRightCornerConductivity(0, m) / 0.5 / 0.25 / xDistanceHelper.GetDelta(1) / xDistanceHelper.GetDelta(1);
	}

	double a = info.capacity * info.density / (0.5 * TIME_DELTA);

	return b / (b + a);
}

float HeatEquationSolver::Ksi_2_x(int k)
{
	MaterialInfo info = materialsHelper.GetInfo(n, k);

	double c = 0;
	if (k == 1)
	{
		double first = materialsHelper.GetUpperCornerConductivity(n, 1) * (T_step[n][2] - T_step[n][1]) / zDistanceHelper.GetDelta(k + 1);
		double second = info.conductivity * (T_step[n][1] - T_step[n][0]) / (0.5 * zDistanceHelper.GetDelta(k));

		c = (first - second) / (0.75 * zDistanceHelper.GetDeltaLine(k));
	}
	else if (k >= 2 && k <= m - 1)
	{
		double first = materialsHelper.GetUpperCornerConductivity(n, k) * (T_step[n][k + 1] - T_step[n][k]) / zDistanceHelper.GetDelta(k + 1);
		double second = materialsHelper.GetBottomCornerConductivity(n, k) * (T_step[n][k] - T_step[n][k - 1]) / zDistanceHelper.GetDelta(k);

		c = (first - second) / zDistanceHelper.GetDeltaLine(k);
	}
	else if (k == m)
	{
		double first = materialsHelper.GetUpperCornerConductivity(n, k) * (T_step[n][m + 1] - T_step[n][m]) / (0.5 * zDistanceHelper.GetDelta(m + 1));
		double second = materialsHelper.GetBottomCornerConductivity(n, k) * (T_step[n][m] - T_step[n][m - 1]) / zDistanceHelper.GetDelta(m);

		c = (first - second) / (0.75 * zDistanceHelper.GetDeltaLine(k));
	}

	double a = info.capacity * info.density / (0.5 * TIME_DELTA);

	return c / a + T_time[n][k];
}

float HeatEquationSolver::Chi_2_x(int k)
{
	MaterialInfo info = materialsHelper.GetInfo(n, k);

	double b = 0;
	if (k == 1)
	{
		b = -info.conductivity / 0.25 / xDistanceHelper.GetDelta(n) / xDistanceHelper.GetDelta(n);
	}
	else if (k >= 2 && k <= m - 1)
	{
		b = -materialsHelper.GetUpperCornerConductivity(n, k) / 0.25 / xDistanceHelper.GetDelta(n) / xDistanceHelper.GetDelta(n);
	}
	else if (k == m)
	{
		b = -materialsHelper.GetUpperCornerConductivity(n, m) / 0.25 / xDistanceHelper.GetDelta(n) / xDistanceHelper.GetDelta(n);
	}

	double a = info.capacity * info.density / 0.5 / TIME_DELTA;

	return b / (b - a);
}

float HeatEquationSolver::A_x(int i, int k)
{
	MaterialInfo info = materialsHelper.GetInfo(i, k);

	return info.capacity * info.density / 0.5 / TIME_DELTA;
}

float HeatEquationSolver::B_x(int i, int k)
{
	return materialsHelper.GetRightCornerConductivity(i, k) / xDistanceHelper.GetDelta(i) / xDistanceHelper.GetDeltaLine(i + 1);
}

float HeatEquationSolver::C_x(int i, int k)
{
	return materialsHelper.GetLeftCornerConductivity(i, k) / xDistanceHelper.GetDelta(i) / xDistanceHelper.GetDeltaLine(i);
}

float HeatEquationSolver::D_x(int i, int k)
{
	double first = materialsHelper.GetUpperCornerConductivity(i, k) * (T_step[i][k + 1] - T_step[i][k]) / zDistanceHelper.GetDeltaLine(k + 1);
	double second = materialsHelper.GetBottomCornerConductivity(i, k) * (T_step[i][k] - T_step[i][k - 1]) / zDistanceHelper.GetDeltaLine(k);

	return (first - second) / zDistanceHelper.GetDelta(k);
}

float HeatEquationSolver::A_z(int i, int k)
{
	return A_x(i, k);
}

float HeatEquationSolver::B_z(int i, int k)
{
	double first = materialsHelper.GetRightCornerConductivity(i, k) * (T_cur[i + 1][k] - T_cur[i][k]) / xDistanceHelper.GetDeltaLine(i + 1);
	double second = materialsHelper.GetLeftCornerConductivity(i, k) * (T_cur[i][k] - T_cur[i - 1][k]) / xDistanceHelper.GetDeltaLine(i);

	return (first - second) / zDistanceHelper.GetDelta(i);
}

float HeatEquationSolver::C_z(int i, int k)
{
	return materialsHelper.GetUpperCornerConductivity(i, k) / zDistanceHelper.GetDelta(k) / zDistanceHelper.GetDeltaLine(k + 1);
}

float HeatEquationSolver::D_z(int i, int k)
{
	return materialsHelper.GetBottomCornerConductivity(i, k) / zDistanceHelper.GetDelta(k) / zDistanceHelper.GetDeltaLine(k);
}

std::pair<int, int> HeatEquationSolver::GetTapeInterval(int t)
{
	std::pair<int, int> interval = std::make_pair(t * tapeSize + 1, (t + 1) * tapeSize);
	if (t == tapesNum - 1)
	{
		interval.second = m;
	}
	
	return interval;
}