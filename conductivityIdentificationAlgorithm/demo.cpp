#include <iostream>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <ctime>

#include "mpi.h"

#include "heat_equation_solver.h"
#include "material_type.h"
#include "material_info.h"
#include "math_options.h"
#include "mpi_options.h"
#include "math_consts.h"

MathOptions InitializeMathOptions()
{
	srand(42);

	float** initTemperature = new float* [DEFAULT_N + 2];
	for (int i = 0; i <= DEFAULT_N + 1; ++i)
	{
		initTemperature[i] = new float[DEFAULT_M + 2];
		for (int j = 0; j <= DEFAULT_M + 1; ++j)
		{
			initTemperature[i][j] = DEFAULT_T_GR3;
		}
	}

	MaterialType** gridMaterialInfo = new MaterialType * [DEFAULT_N + 2];
	for (int i = 0; i <= DEFAULT_N + 1; ++i)
	{
		gridMaterialInfo[i] = new MaterialType[DEFAULT_M + 2];
		for (int j = 0; j <= DEFAULT_M + 1; ++j)
		{
			MaterialType material = static_cast<MaterialType>(rand() % (int)MaterialType::LAST);
			gridMaterialInfo[i][j] = material;

			if (i == 1)
			{
				gridMaterialInfo[0][j] = gridMaterialInfo[i][j];
			}
			else if (i == DEFAULT_N + 1)
			{
				gridMaterialInfo[DEFAULT_N + 1][j] = gridMaterialInfo[DEFAULT_N][j];
			}

			if (j == 1)
			{
				gridMaterialInfo[i][0] = gridMaterialInfo[i][j];
			}
			else if (j == DEFAULT_M + 1)
			{
				gridMaterialInfo[i][DEFAULT_M + 1] = gridMaterialInfo[i][j];
			}
		}
	}

	MathOptions mathOptions{ DEFAULT_N, DEFAULT_M, DEFAULT_T_GR3, DEFAULT_T_GR4, INIT_TIME, END_TIME, DEFAULT_STEP, DEFAULT_STEP, initTemperature, gridMaterialInfo };

	return mathOptions;
}

void DeleteMathOptions(MathOptions& mathOptions)
{
	for (int i = 0; i <= mathOptions.n + 1; ++i)
	{
		delete[] mathOptions.initTemperature[i];
		delete[] mathOptions.gridMaterialInfo[i];
	}
	delete[] mathOptions.initTemperature;
	delete[] mathOptions.gridMaterialInfo;
}

void WriteGridMaterials(const MathOptions& mathOptions)
{
	std::ofstream out("grid_materials.txt");
	for (int i = 1; i <= mathOptions.n; ++i)
	{
		for (int j = 1; j <= mathOptions.m; ++j)
		{
			out << static_cast<int>(mathOptions.gridMaterialInfo[i][j]) << " ";
		}
		out << "\n";
	}
	out.close();
}

int main(int argc, char** argv)
{
	int rank = -1;
	int size = -1;
	int tapeSize = -1;
	bool useMPI = atoi(argv[1]);
	if (useMPI)
	{
		tapeSize = atoi(argv[2]);
	}

	if (useMPI)
	{
		MPI_Init(&argc, &argv);

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
	}

	MathOptions mathOptions = InitializeMathOptions();
	MPIOptions mpiOptions { useMPI, rank, size, tapeSize };
	HeatEquationSolver* solver = new HeatEquationSolver(mathOptions, mpiOptions);

	std::clock_t start;
	double duration;

	if (!useMPI || rank == 0)
	{
		start = std::clock();
	}

	solver->Solve();

	if (!useMPI || rank == 0)
	{
		duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
		std::cout << "Time: " << duration << std::endl;

		float res = solver->GetResult();
		std::cout << "Conductivity: " << res << std::endl;
	}
	
	delete solver;

	DeleteMathOptions(mathOptions);

	if (useMPI)
	{
		MPI_Finalize();
	}

	return 0;
}


