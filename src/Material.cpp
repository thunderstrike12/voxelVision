#include "precomp.h"
#include "Material.h"

Metal::Metal() {
	reflectivity = 0.9;
	color = float3(0.5, 0.5, 0.5);
}

Ray Metal::reflectRay(Ray ray) {
	return ray;
}
