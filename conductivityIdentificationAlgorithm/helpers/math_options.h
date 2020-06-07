#pragma once

#include "material_type.h"

struct MathOptions
{
	int n;
	int m;
	float t_gr3;
	float t_gr4;
	float initTime;
	float endTime;
	float xStep;
	float zStep;
	float** initTemperature;
	MaterialType** gridMaterialInfo;
};
