#include "precomp.h"
#include "Material.h"

None::None(){
	reflectivity = 0;
	color = float3(1, 0, 0);
}

Ray None::reflectRay(Ray ray)
{
	return Ray();
}

Metal::Metal() {
	reflectivity = 0.5;
	color = float3(0.5, 0.5, 0.5);
}

Ray Metal::reflectRay(Ray ray) {
	return ray;
}


