#include "precomp.h"
#include "Primitives.h"

Sphere::Sphere(float3 _pos, float _rad){
    pos = _pos;
    rad = _rad;
}

bool Sphere::intersect(Ray ray, float& min, float& max)
{
    float3 oc = ray.O - pos;
    float a = dot(ray.D, ray.D);
    float b = 2.0f * dot(oc, ray.D);
    float c = dot(oc, oc) - rad * rad;
    float discriminant = b * b - 4 * a * c;

    if (discriminant > 0) {
        float t1 = (-b - std::sqrt(discriminant)) / (2.0f * a);
        float t2 = (-b + std::sqrt(discriminant)) / (2.0f * a);

        min = (t1 < t2) ? t1 : t2;
        max = (t1 > t2) ? t1 : t2;

        return true;
    }

    return false;
}
