#include "precomp.h"
#include "Light.h"

PointLight::PointLight(float3 pos_, float3 color_) {
    pos = pos_;
    color = color_;
}

SpotLight::SpotLight(float3 dir_, float3 pos_, float3 color_, float angle_) {
	dir = dir_;
	pos = pos_;
	color = color_;
	angle = angle_;
}