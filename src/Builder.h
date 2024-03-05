#pragma once

class Builder
{
public:
	bool cd = false;
	bool obstructed = false;
	float buildDelete = 1.0f;

	float3 hover = 0;
	int mat = 1;
	bool building = true;

	Builder() {}
	void update(Scene scene, Camera camera, int2 mousePos, bool mouse);
};

