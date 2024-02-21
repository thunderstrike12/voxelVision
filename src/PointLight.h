#pragma once
#include <vector>

class PointLight
{
public:
	float3 pos;
	float3 color;
	float radius = 2;

	PointLight(float3 pos_, float3 color_);
	void add(float3 pos);
	float getBrightness(float3 point, float s);
};

