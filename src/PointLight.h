#pragma once
#include <vector>

class PointLight
{
public:
	std::vector<float3> lights;

	void add(float3 pos);
	float getBrightness(float3 point, float s);
};

