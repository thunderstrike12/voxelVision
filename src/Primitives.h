#pragma once
class Primitives
{
public:
	float3 pos;
	virtual bool intersect(Ray ray, float& min, float& max) = 0;
};

class Sphere : public Primitives {
public:
	float rad;
	Sphere(float3 _pos, float _rad);
	bool intersect(Ray ray, float& min, float& max) override;
};

