#include "precomp.h"
#include "PointLight.h"

float distance(const float3& p1, const float3& p2) {
	float dx = p2.x - p1.x;
	float dy = p2.y - p1.y;
	float dz = p2.z - p1.z;

	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

PointLight::PointLight(float3 pos_, float3 color_) {
    pos = pos_;
    color = color_;
}