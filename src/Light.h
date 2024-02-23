#pragma once
#include <vector>

class Light
{
public:
	float3 pos;
	float3 color;
};

class PointLight : public Light {
public:
	PointLight(float3 pos_, float3 color_);
};

class SpotLight : public Light {
public:
	float angle;
	float intensity;
	float3 dir;

	SpotLight(float3 dir_, float3 pos_, float3 color_, float radius_);
};
