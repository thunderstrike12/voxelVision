#include "precomp.h"
#include "Builder.h"

void Builder::update(Scene scene, Camera camera, int2 mousePos, bool mouse) {
	scene.Delete(hover.x, hover.y, hover.z);

	Ray mouseRay = camera.GetPrimaryRay(mousePos.x, mousePos.y);
	scene.FindNearest(mouseRay, GRIDLAYERS);
	
	if (building) {
		buildDelete = 1.0f;
	}
	else {
		buildDelete = -1.0f;
	}
	float3 I = mouseRay.O + mouseRay.D * mouseRay.t + mouseRay.GetNormal() * 0.0001f * buildDelete;

	float Ix = std::floor(I.x * WORLDSIZE);
	float Iy = std::floor(I.y * WORLDSIZE);
	float Iz = std::floor(I.z * WORLDSIZE);

	obstructed = false;
	if (Ix < 0 || Ix > WORLDSIZE) obstructed = true;
	if (Iy < 0 || Iy > WORLDSIZE) obstructed = true;
	if (Iz < 0 || Iz > WORLDSIZE) obstructed = true;

	uint Ixui = Ix;
	uint Iyui = Iy;
	uint Izui = Iz;

	if (building && !obstructed) {
		scene.Set(Ixui, Iyui, Izui, mat);
		hover = float3(Ixui, Iyui, Izui);
	}

	if (mouse && !cd && !obstructed) {
		if (building) {
			scene.Set(Ixui, Iyui, Izui, mat);
		}
		else {
			scene.Delete(Ixui, Iyui, Izui);
		}
		hover = 0;
		cd = true;
	}

	if (!mouse) {
		cd = false;
	}
}
