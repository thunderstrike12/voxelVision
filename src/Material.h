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

class None : public Material {
public:
	None();
	Ray reflectRay(Ray ray) override;
};

class Metal : public Material {
public:
	Metal();
	Ray reflectRay(Ray ray) override;
};
