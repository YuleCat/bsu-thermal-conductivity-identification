#pragma once

class DistanceHelper
{
private:
	int n;
	float* gridValues;

public:
	DistanceHelper()
	{
		n = 0;
		gridValues = nullptr;
	}

	DistanceHelper(const DistanceHelper& helper)
		:
		n(helper.n)
	{
		this->gridValues = new float[n + 1];
		for (int i = 0; i < n + 1; ++i)
		{
			this->gridValues[i] = helper.gridValues[i];
		}
	}

	DistanceHelper(int n, float* gridValues)
		:
		n(n)
	{
		this->gridValues = new float[n + 1];
		for (int i = 0; i < n + 1; ++i)
		{
			this->gridValues[i] = gridValues[i];
		}
	}

	DistanceHelper& operator=(const DistanceHelper& helper)
	{
		n = helper.n;
		gridValues = new float[n + 1];
		for (int i = 0; i < n + 1; ++i)
		{
			gridValues[i] = helper.gridValues[i];
		}

		return *this;
	}

	float GetDelta(int i)
	{
		if (i == n + 1)
		{
			return 0.25f * (gridValues[n] - gridValues[n - 1]);
		}
		else if (i == n)
		{
			return 0.75f * (gridValues[n] - gridValues[n - 1]);
		}

		if (i == 0)
		{
			return gridValues[0];
		}
		// TODO: i == 1?

		return gridValues[i] - gridValues[i - 1];
	}

	float GetDeltaLine(int i)
	{
		return 0.5f * (GetDelta(i) + GetDelta(i - 1));
	}

	float GetSize()
	{
		return gridValues[n];
	}

	~DistanceHelper()
	{
		delete[] gridValues;
	}
};
