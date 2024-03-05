#pragma once
class Material
{
public:
	float reflectivity = 0;
	float absorption = 0;
	float refraction = 1.0f;
	float glossyness = 0;

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

class Glass : public Material {
public:
	Glass();
	Ray reflectRay(Ray ray) override;
};

class Diamond : public Material {
public:
	Diamond();
	Ray reflectRay(Ray ray) override;
};
