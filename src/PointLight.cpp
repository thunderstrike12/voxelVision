#include "precomp.h"
#include "PointLight.h"

float distance(const float3& p1, const float3& p2) {
	float dx = p2.x - p1.x;
	float dy = p2.y - p1.y;
	float dz = p2.z - p1.z;

	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void PointLight::add(float3 pos) {
	lights.push_back(pos);
}

float PointLight::getBrightness(float3 point, float s) {
    const float kAttenuationFactor = 1.0f;  // Adjust as needed

    float strength = 0;

    for (int i = 0; i < lights.size(); i++) {
        float dist = distance(point, lights[i]);

        // Inverse Square Law with Attenuation
        float attenuation = kAttenuationFactor / (dist * dist);
        strength += attenuation * s;
    }

    // You might want to clamp the result to a maximum value if needed
    const float kMaxBrightness = 1.0f;
    strength = min(strength, kMaxBrightness);

    return strength;
}