#include "precomp.h"
#include "Material.h"

None::None(){
	reflectivity = 0;
	color = float3(1, 0, 0);
}

Ray None::reflectRay(Ray ray) {
	return Ray();
}

Metal::Metal() {
	reflectivity = 0.9;
	glossyness = 0.1;
	color = float3(0.5, 0.5, 0.5);
}

Ray Metal::reflectRay(Ray ray) {
	return ray;
}

Glass::Glass() {
	refraction = 1.52;
	color = float3(0.5, 0.5, 0.5);
	absorption = 0.5;
}

Ray Glass::reflectRay(Ray ray) {
	return ray;
}

Diamond::Diamond() {
	refraction = 1.9;
	color = float3(0.5, 0.5, 0.5);
	absorption = 0.5;
}

Ray Diamond::reflectRay(Ray ray)
{
	return ray;
}
