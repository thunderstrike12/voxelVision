#pragma once
class Material
{
public:
	float reflectivity;
	float absorption;
	float refraction;

	float3 color;

	virtual Ray reflectRay(Ray ray) = 0;
};

class Metal : public Material {
	Metal();
	Ray reflectRay(Ray ray) override;
};
