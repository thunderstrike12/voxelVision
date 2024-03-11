#include "precomp.h"
#include <stb_image.h>

// YOU GET:
// 1. A fast voxel renderer in plain C/C++
// 2. Normals and voxel colors
// FROM HERE, TASKS COULD BE:							FOR SUFFICIENT
// * Materials:
//   - Reflections and diffuse reflections				<===
//   - Transmission with Snell, Fresnel					<===
//   - Textures, Minecraft-style						<===
//   - Beer's Law
//   - Normal maps
//   - Emissive materials with postproc bloom
//   - Glossy reflections (BASIC)
//   - Glossy reflections (microfacet)
// * Light transport:
//   - Point lights										<===
//   - Spot lights										<===
//   - Area lights										<===
//	 - Sampling multiple lights with 1 ray
//   - Importance-sampling
//   - Image based lighting: sky
// * Camera:
//   - Depth of field									<===
//   - Anti-aliasing									<===
//   - Panini, fish-eye etc.
//   - Post-processing: now also chromatic				<===
//   - Spline cam, follow cam, fixed look-at cam
//   - Low-res cam with CRT shader
// * Scene:
//   - HDR skydome										<===
//   - Spheres											<===
//   - Smoke & trilinear interpolation
//   - Signed Distance Fields
//   - Voxel instances with transform
//   - Triangle meshes (with a BVH)
//   - High-res: nested grid
//   - Procedural art: shapes & colors
//   - Multi-threaded Perlin / Voronoi
// * Various:
//   - Object picking
//   - Ray-traced physics
//   - Profiling & optimization
// * GPU:
//   - GPU-side Perlin / Voronoi
//   - GPU rendering *not* allowed!
// * Advanced:
//   - Ambient occlusion
//   - Denoising for soft shadows
//   - Reprojection for AO / soft shadows
//   - Line lights, tube lights, ...
//   - Bilinear interpolation and MIP-mapping
// * Simple game:										
//   - 3D Arkanoid										<===
//   - 3D Snake?
//   - 3D Tank Wars for two players
//   - Chess
// REFERENCE IMAGES:
// https://www.rockpapershotgun.com/minecraft-ray-tracing
// https://assetsio.reedpopcdn.com/javaw_2019_04_20_23_52_16_879.png
// https://www.pcworld.com/wp-content/uploads/2023/04/618525e8fa47b149230.56951356-imagination-island-1-on-100838323-orig.jpg

// -----------------------------------------------------------
// Initialize the renderer
// -----------------------------------------------------------
void Renderer::Init()
{
	// create fp32 rgb pixel buffer to render to
	accumulator = (float4*)MALLOC64( SCRWIDTH * SCRHEIGHT * 16 );
	memset( accumulator, 0, SCRWIDTH * SCRHEIGHT * 16 );
	// try to load a camera
	FILE* f = fopen( "camera.bin", "rb" );
	if (f)
	{
		fread( &camera, 1, sizeof( Camera ), f );
		fclose( f );
	}
	
	//pLight.push_back(PointLight(float3(.576, .698, .471), float3(0.01, 0, 0)));
	//pLight.push_back(PointLight(float3(.488, .730, .553), float3(0, 0, 0.01)));
	//pLight.push_back(PointLight(float3(.488, .700, .583), float3(0, 0.01, 0)));

	//float3 sLightDir = normalize(float3(1, 2, 3));
	//sLight.push_back(SpotLight(sLightDir, float3(.488, .700, .583), float3(1, 1, 1), 20));

	aLight.push_back(AreaLight(float3(.488, .730, .553), float3(0.01, 0.01, 0.01), 1));

	materials.push_back(new None());
	materials.push_back(new Metal());
	materials.push_back(new Glass());
	materials.push_back(new Diamond());

	build = new Builder();

	sphere = new Sphere(float3(0, 0, 0), 1);

	int skyChannels;
	skyData = stbi_loadf("assets/sky.hdr", &skyWidth, &skyHeight, &skyChannels, 0);
	for (int i = 0; i < skyWidth * skyHeight * 3; i++)
		skyData[i] = sqrtf(skyData[i]); // Gamma Adjustment for Reduced HDR Range
}

float3  snellsLaw(Ray ray, float3 N, float n1, float n2) {
	float eta = n2 / n1;

	float b = dot(N, -ray.D);


	float k = 1 - (eta * eta) * (1 - (b * b));
	if (k >= 0)
	{
		return float3(eta * ray.D + N * ((eta * b) - sqrt(k)));
	}
	else
	{
		return float3(-1, -1, -1);
	}
}

// -----------------------------------------------------------
// Evaluate light transport
// ----------------------------------------------------------
float3 Renderer::Trace(Ray& ray)
{
	if (ray.depth > 5) 
		return 0;
	ray.depth++;

	ray.t = 0;

	scene.FindNearest(ray, GRIDLAYERS);
	float min = 0, max = 0;
	if (sphere->intersect(ray, min, max)) {
		if (min < ray.t) {
			ray.t = min;
		}
	}

	//skydome code from milan bonten
	if (ray.voxel == 0) {
		float3 dir = normalize(ray.D);
		uint u = skyWidth * atan2f(dir.z, dir.x) * INV2PI - 0.5f;
		uint v = skyHeight * acosf(dir.y) * INVPI - 0.5f;
		uint skyIdx = (u + v * skyWidth) % (skyWidth * skyHeight);
		return 0.5f * float3(skyData[skyIdx * 3], skyData[skyIdx * 3 + 1], skyData[skyIdx * 3 + 2]);
	}

	float3 I = ray.O + ray.t * ray.D;
	//static const float3 L = normalize( float3( 1, 4, 0.5f ) );
	float3 N = ray.GetNormal();
	float3 albedo = float3(1, 1, 1);

	ray.entering = ray.GetMaterial();
	Material* matExit = materials[ray.exiting];
	Material* matEnter = materials[clamp(ray.GetMaterial() - 1, 0, 10)];


	if (matEnter == materials[Scene::METAL - 1])
		matEnter->reflectivity = reflectivity;
	if (matEnter == materials[Scene::METAL - 1])
		matEnter->glossyness = glossyness;
	if (matEnter == materials[Scene::GLASS - 1])
		matEnter->refraction = refraction;

	float3 reflResult;
	if (matEnter->reflectivity > 0 || matEnter->refraction > 1.0f)
	{
		float3 rf3 = normalize(float3(rand(), rand(), rand())) * matEnter->glossyness;
		// calculate the specular reflection in the intersection point
		float3 dir = normalize(ray.D + rf3);
		float3 direction = dir - 2 * N * dot(N, dir);
		float3 origin = I + N * 0.001f;
		Ray secondary(origin, direction);

		secondary.depth = ray.depth + 1;

		if (secondary.depth >= 20) return float3(0);

		reflResult = Trace(secondary);

	}

	float3 refrResult = 0;
	float3 secI;
	if (matEnter->refraction > 1.0f) {
		float3 Dir = normalize( snellsLaw(ray, normalize( N), 1, refraction));
		Ray secondary(ray.IntersectionPoint(), ray.D);
		scene.FindNearestExitDie(secondary, 1);
		secI = secondary.IntersectionPoint();

		Dir = normalize( snellsLaw(ray, normalize(secondary.GetNormal()), refraction, 1));
		
		if (Dir.x == -1 && Dir.y == -1 && Dir.z == -1)
		{
			return float3(0, 0, 0);
		}
		else
		{

			Ray third(secI, Dir);
			third.depth = ray.depth;
			
			refrResult = Trace(third);
		}



		//float eta = 0;
		//eta = matExit->refraction / matEnter->refraction;

		//float theta1 = dot(N, -ray.D);
		//float cosTheta1 = cos(theta1);
		//float k = 1 - (eta * eta) * (1 - cosTheta1 * cosTheta1);

		//if (k >= 0) {
		//	float3 T = eta * ray.D + N * (eta * cosTheta1 - sqrt(k));
		//	Ray secondary(I, T);
		//	//secondary.depth = ray.depth + 1;
		//	//if (secondary.depth > 20) return float3(0);

		//	secondary.exiting = ray.entering - 1;
		//	//if (ray.entering == Scene::GLASS || ray.entering == Scene::DIAMOND) {
		//		//refrResult = Trace(secondary, true);
		//		scene.FindNearestExitDie(secondary, 0);
		//		float3 secI = secondary.O + secondary.D * secondary.t;

		//		eta = matEnter->refraction / matExit->refraction;

		//		float theta1 = dot(-N, -ray.D);
		//		float cosTheta1 = cos(theta1);
		//		float k = 1 - (eta * eta) * (1 - cosTheta1 * cosTheta1);
		//		if (k >= 0) {
		//			float3 T = eta * ray.D + -N * (eta * cosTheta1 - sqrt(k));
		//			Ray third(secI - N * 0.001f, T);
		//			refrResult = Trace(third);
		//		}
		//		else {
		//			//TIR
		//		}
		//	//}
		//	//else {
		//	//	refrResult = Trace(secondary);
		//	//}
		//}
		//else {
		//	//TIR
		//}
	}

	//float3 refrResult = 0;
	//float3 secI;
	//if ((matEnter->refraction > 1.0f || dielectric) && matEnter != matExit) {
	//	float eta = 0;

	//	eta = matEnter->refraction / matExit->refraction;

	//	float theta1 = dot(N, -ray.D);
	//	float cosTheta1 = cos(theta1);
	//	float k = 1 - (eta * eta) * (1 - cosTheta1 * cosTheta1);

	//	Ray secondary;
	//	if (k >= 0) {
	//		float3 T = eta * ray.D + N * (eta * cosTheta1 - sqrt(k));
	//		secondary = Ray(I - N * 0.001f, T);
	//		scene.FindNearestExitDie(secondary, 1);

	//		eta = 1.5;


	//		float3 secN = secondary.GetNormal();

	//		theta1 = dot(-N, -secondary.D);
	//		cosTheta1 = cos(theta1);
	//		k = 1 - (eta * eta) * (1 - cosTheta1 * cosTheta1);

	//		if (k >= 0) {
	//			float3 T = eta * secondary.D + -N * (eta * cosTheta1 - sqrt(k));
	//			secI = secondary.O + secondary.D * secondary.t;
	//			Ray third(secI + -N * 0.001f, T);

	//			third.depth = ray.depth++;
	//			if (third.depth > 20) return(1, 0, 0);
	//			refrResult = Trace(third);
	//		}
	//		else {
	//			//TIR
	//		}
	//	}
	//	else {
	//		// Total Internal Reflection
	//		float3 R = ray.D - 2 * N * dot(N, ray.D);
	//		secondary = Ray(I - N * 0.001f, R);
	//	}
	//}

	float reflectance;
	float absorpResult;
	if (matEnter->refraction > 1.0f) {
		float angle = dot(N, ray.D);

		float cosAi = cos(angle);
		float sinAi = sin(angle);

		float n1 = matExit->refraction;
		float n2 = matEnter->refraction;
		float eta = n1 / n2;

		float cosAt = sqrt(1 - eta * sinAi * eta * sinAi);

		float term1 = (n1 * cosAi - n2 * cosAt) / (n1 * cosAi + n2 * cosAt);
		float term2 = (n1 * cosAt - n2 * cosAi) / (n1 * cosAt + n2 * cosAi);

		float reflectance = 0.5 * (term1 * term1 + term2 * term2);

		reflectance = std::clamp(reflectance, 0.0f, 1.0f);
		
		if (RandomFloat() > reflectance) {
			refrResult = reflResult;
		}

		//refrResult = (1.0 - reflectance) * reflResult + reflectance * refrResult;
		/*float d = std::clamp(ray.t * GRIDSIZE / 2.0f, 0.0f, 1.0f);

		refrResult = refrResult * (1 - d) + -matEnter->absorption * d;*/
		float d = length(I - secI) * pow(10, power);

		float3 absorption = exp(-matEnter->absorption * d);

		// Apply Beer's Law to the result
		refrResult *= absorption;

		/*float expx = exp(-absorption.x * d);
		float expy = exp(-absorption.y * d);
		float expz = exp(-absorption.z * d);

		refrResult.x *= expx;
		refrResult.y *= expy;
		refrResult.z *= expz;*/
	}

	
	float3 pLightColor = 0;
	for (int i = 0; i < pLight.size(); i++)
	{
		float3 dir = pLight[i].pos - I;
		float dist = length(dir);
		dir = normalize(dir);
		float d = dot(N, dir);
		if (d <= 0) {
			continue;
		}
		if (!scene.IsOccluded(Ray(I, dir, dist), GRIDLAYERS)) {
			pLightColor += (pLight[i].color * (1 / (dist * dist))) * d;
		}
	}

	float3 sLightColor = 0;
	for (int i = 0; i < sLight.size(); i++)
	{
		float3 dir = sLight[i].pos - I;
		float dist = length(dir);
		dir = normalize(dir);

		float dotP = dot(dir, normalize(sLight[i].dir));

		if (dotP > sLight[i].angle) {
			float d = dot(N, dir);
			if (d <= 0) {
				continue;
			}
			if (!scene.IsOccluded(Ray(I, dir, dist), GRIDLAYERS)) {
				sLightColor += ((sLight[i].color * (1 / (dist * dist))) * d) / sLight[i].intensity;
			}
		}
	}

	float3 aLightColor = 0;
	for (int i = 0; i < aLight.size(); i++)
	{
		float3 rf3 = normalize(float3(rand(), rand(), rand())) * aLight[i].radius;

		float3 dir = aLight[i].pos + rf3 - I;
		float dist = length(dir);
		dir = normalize(dir);

		float d = dot(N, dir);
		if (d <= 0) {
			continue;
		}
		//Ray(aLight[i].pos + rf3, -dir, dist)
		if (!scene.IsOccluded(Ray(aLight[i].pos + rf3, dir, dist), GRIDLAYERS)) {
			aLightColor += ((aLight[i].color * (1 / (dist * dist))) * d);
		}
	}
	
	float3 dirLightColor = 0;
	float angle = dot(N, normalize(-lightDir));
	float shadowStrength = 1;
	
	if (angle <= 0) {
		shadowStrength = 0;
	}
	// Cast shadow ray
	else if (scene.IsOccluded(Ray(I * -lightDir * 10000.0f, lightDir), GRIDLAYERS)) {
		shadowStrength = 0;
	}

	
	dirLightColor = albedo * angle * shadowStrength;
	float3 finalColor = dirLightColor + pLightColor + sLightColor + aLightColor;

	if (matEnter->refraction > 1.0f) {
		finalColor = (refrResult * (1)) + (finalColor * (1 - (1)));
	}
	else {
		finalColor = (reflResult * matEnter->reflectivity) + (finalColor * (1 - matEnter->reflectivity));
	}
	return finalColor;
}

void Renderer::resetAcc() {
	for (int y = 0; y < SCRHEIGHT; y++) {
		for (int x = 0; x < SCRWIDTH; x++) {
			accumulator[x + y * SCRWIDTH] = float4(0.0f, 0.0f, 0.0f, 0.0f);
		}
	}
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime)
{
	// pixel loop
	Timer t;

	float framesInv = 1.0f / frames;
	// lines are executed as OpenMP parallel tasks (disabled in DEBUG)
#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCRHEIGHT; y++)
	{
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++)
		{
			if (accumulatorEnabled) {
				float4 pixel = float4(Trace(camera.GetPrimaryRay((float)x + RandomFloat(), (float)y + RandomFloat())), 0);
				// translate accumulator contents to rgb32 pixels
				accumulator[x + y * SCRWIDTH] += pixel;
				float4 accPix = accumulator[x + y * SCRWIDTH] * framesInv;

				screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(&accPix); //implement bit shift
			}
			else {
				float4 pixel = float4(Trace(camera.GetPrimaryRay((float)x, (float)y)), 0);
				screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(&pixel);
			}
		}
	}
	
	// handle user input

	// performance report - running average - ms, MRays/s


	//building
	build->update(scene, camera, mousePos, mouse);

	frames++;
	if (camera.HandleInput(deltaTime)) {
		resetAcc();
		frames = 1;
	}
	
	static float alpha = 1, avg = 10;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	fps = 1000.0f / avg, rps = (SCRWIDTH * SCRHEIGHT) / avg;
}

// -----------------------------------------------------------
// Update user interface (imgui)
// -----------------------------------------------------------
void Renderer::UI()
{
	//ImGui::Text( "voxel: %i", r.voxel );
	ImGui::Text("FPS: %f", fps);
	ImGui::Text("RPS: %f", rps);

	// ray query on mouse
	Ray r = camera.GetPrimaryRay( (float)mousePos.x, (float)mousePos.y );
	scene.FindNearest( r, GRIDLAYERS );


	if (ImGui::CollapsingHeader("lighting")) {
		ImGui::SliderFloat3("light direction", &lightDir.x, -1.0f, 1.0f);

		ImGui::Text("point lights");
		for (int i = 0; i < pLight.size(); i++)	{
			std::string label = "plight pos " + std::to_string(i);
			ImGui::SliderFloat3(label.c_str(), &pLight[i].pos.x, -1.0f, 1.0f);
		}

		ImGui::Text("spot lights");
		for (int i = 0; i < sLight.size(); i++) {
			std::string label = "slight pos " + std::to_string(i);
			ImGui::SliderFloat3(label.c_str(), &sLight[i].pos.x, -1.0f, 1.0f);

			label = "slight dir " + std::to_string(i);
			ImGui::SliderFloat3(label.c_str(), &sLight[i].dir.x, -1.0f, 1.0f);

			label = "slight angle " + std::to_string(i);
			ImGui::SliderFloat(label.c_str(), &sLight[i].angle, 0.0f, 1.0f);

			label = "slight intensity " + std::to_string(i);
			ImGui::SliderFloat(label.c_str(), &sLight[i].intensity, -10.0f, 100.0f);
		}

		ImGui::Text("area lights");
		for (int i = 0; i < aLight.size(); i++) {
			std::string label = "alight pos " + std::to_string(i);
			ImGui::SliderFloat3(label.c_str(), &aLight[i].pos.x, -1.0f, 1.0f);

			label = "alight radius " + std::to_string(i);
			ImGui::SliderFloat(label.c_str(), &aLight[i].radius, 0.0f, 1.0f);
		}
	}
	
	if (ImGui::CollapsingHeader("materials")) {
		ImGui::SliderFloat("reflectivity", &reflectivity, 0.0f, 1.0f);
		ImGui::SliderFloat("glossyness", &glossyness, 0.0f, 1.0f);
		ImGui::SliderFloat("refraction", &refraction, 1.0f, 5.0f);

		ImGui::SliderFloat("absorption", &materials[2]->absorption, 0.0f, 2.0f);
		ImGui::SliderFloat("power", &power, 1.0f, 10.0f);
	}

	if (ImGui::CollapsingHeader("primitives")) {
		ImGui::SliderFloat3("sphere pos", &sphere->pos.x, 0.0f, 1.0f);
		ImGui::SliderFloat("sphere rad", &sphere->rad, 0.0f, 1.0f);
	}

	if (ImGui::CollapsingHeader("building")) {
		ImGui::Checkbox("accumulator enabled", &accumulatorEnabled);
		ImGui::Checkbox("building/deleteing", &build->building);
		ImGui::SliderInt("material", &build->mat, 1, 4);
	}
}

// -----------------------------------------------------------
// User wants to close down
// -----------------------------------------------------------
void Renderer::Shutdown()
{
	// save current camera
	FILE* f = fopen( "camera.bin", "wb" );
	fwrite( &camera, 1, sizeof( Camera ), f );
	fclose( f );
}