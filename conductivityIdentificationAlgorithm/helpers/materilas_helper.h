#pragma once

#include "math_consts.h"
#include "material_info.h"
#include "material_type.h"

class MaterialsHelper
{
private:
	int n;
	int m;
	MaterialType** gridMaterialInfo;

public:
	MaterialsHelper()
	{
		n = 0;
		m = 0;
		gridMaterialInfo = nullptr;
	}
	
	MaterialsHelper(const MaterialsHelper& helper)
		:
		n(helper.n),
		m(helper.m)
	{
		this->gridMaterialInfo = new MaterialType * [n + 2];
		for (int i = 0; i <= n + 1; ++i)
		{
			this->gridMaterialInfo[i] = new MaterialType[m + 2];
			for (int j = 0; j <= m + 1; ++j)
			{
				this->gridMaterialInfo[i][j] = helper.gridMaterialInfo[i][j];
			}
		}
	}

	MaterialsHelper(int n, int m, MaterialType** gridMaterialInfo)
		:
		n(n),
		m(m)
	{
		this->gridMaterialInfo = new MaterialType * [n + 2];
		for (int i = 0; i <= n + 1; ++i)
		{
			this->gridMaterialInfo[i] = new MaterialType[m + 2];
			for (int j = 0; j <= m + 1; ++j)
			{
				this->gridMaterialInfo[i][j] = gridMaterialInfo[i][j];
			}
		}
	}

	MaterialsHelper& operator=(const MaterialsHelper& helper)
	{
		n = helper.n;
		m = helper.m;
		gridMaterialInfo = new MaterialType * [n + 2];
		for (int i = 0; i <= n + 1; ++i)
		{
			gridMaterialInfo[i] = new MaterialType[m + 2];
			for (int j = 0; j <= m + 1; ++j)
			{
				gridMaterialInfo[i][j] = helper.gridMaterialInfo[i][j];
			}
		}

		return *this;
	}

	MaterialInfo GetInfo(int i, int j)
	{
		MaterialType material = gridMaterialInfo[i][j];
		
		MaterialInfo info;
		switch (material)
		{
		case MaterialType::GRAPHITE:
			info.capacity = 709;
			info.conductivity = 375 / MULTIPLIER;
			info.density = 2230 / (MULTIPLIER * MULTIPLIER * MULTIPLIER);
			info.type = MaterialType::GRAPHITE;
			break;

		case MaterialType::IRON:
			info.capacity = 449;
			info.conductivity = 94 / MULTIPLIER;
			info.density = 7800 / (MULTIPLIER * MULTIPLIER * MULTIPLIER);
			info.type = MaterialType::IRON;
			break;

		case MaterialType::OIL:
			info.capacity = 1600;
			info.conductivity = 0.12 / MULTIPLIER;
			info.density = 959 / (MULTIPLIER * MULTIPLIER * MULTIPLIER);
			info.type = MaterialType::OIL;
			break;

		case MaterialType::PLUMBUM:
			info.capacity = 130;
			info.conductivity = 35 / MULTIPLIER;
			info.density = 1134 / (MULTIPLIER * MULTIPLIER * MULTIPLIER);
			info.type = MaterialType::PLUMBUM;
			break;

		case MaterialType::BRONZE:
			info.capacity = 385;
			info.conductivity = 401 / MULTIPLIER;
			info.density = 8960 / (MULTIPLIER * MULTIPLIER * MULTIPLIER);
			info.type = MaterialType::BRONZE;
			break;

		default:
			info.capacity = 1;
			info.conductivity = 1;
			info.density = 1;
			info.type = MaterialType::LAST;
		}

		return info;
	}

	float GetRightCornerConductivity(int i, int k)
	{
		MaterialInfo info1 = GetInfo(i, k);
		MaterialInfo info2 = GetInfo(i + 1, k);

		double cond1 = info1.conductivity;
		double cond2 = info2.conductivity;

		return 2 * cond1 * cond2 / (cond1 + cond2);
	}

	float GetLeftCornerConductivity(int i, int k)
	{
		return GetRightCornerConductivity(i - 1, k);
	}

	float GetUpperCornerConductivity(int i, int k)
	{
		MaterialInfo info1 = GetInfo(i, k);
		MaterialInfo info2 = GetInfo(i, k + 1);

		double cond1 = info1.conductivity;
		double cond2 = info2.conductivity;

		return 2 * cond1 * cond2 / (cond1 + cond2);
	}

	float GetBottomCornerConductivity(int i, int k)
	{
		return GetUpperCornerConductivity(i, k - 1);
	}

	~MaterialsHelper()
	{
		for (int i = 0; i <= n; ++i)
		{
			delete[] gridMaterialInfo[i];
		}
		delete[] gridMaterialInfo;
	}
};
