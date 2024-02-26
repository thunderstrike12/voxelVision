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
	
	pLight.push_back(PointLight(float3(.576, .698, .471), float3(0.01, 0, 0)));
	pLight.push_back(PointLight(float3(.488, .730, .553), float3(0, 0, 0.01)));
	pLight.push_back(PointLight(float3(.488, .700, .583), float3(0, 0.01, 0)));

	materials.push_back(new None());
	materials.push_back(new Metal());

	float3 sLightDir = normalize(float3(1, 2, 3));
	sLight.push_back(SpotLight(sLightDir, float3(.488, .700, .583), float3(1, 1, 1), 20));

	int skyChannels;
	skyData = stbi_load("assets/sky.hdr", &skyWidth, &skyHeight, &skyChannels, 0);
}

// -----------------------------------------------------------
// Evaluate light transport
// -----------------------------------------------------------
float3 Renderer::Trace( Ray& ray )
{
	scene.FindNearest( ray );
	if (ray.voxel == 0) {
		float3 dir = normalize(ray.D);
		uint u = skyWidth * atan2f(dir.z, dir.x) * INV2PI - 0.5f;
		uint v = skyHeight * acosf(dir.y) * INVPI - 0.5f;
		uint skyIdx = (u + v * skyWidth) % (skyWidth * skyHeight);
		return 0.005f * float3(skyData[skyIdx * 3], skyData[skyIdx * 3 + 1], skyData[skyIdx * 3 + 2]);
	}
	float3 I = ray.O + ray.t * ray.D;
	static const float3 L = normalize( float3( 1, 4, 0.5f ) );
	float3 N = ray.GetNormal();
	float3 albedo = float3(1, 1, 1);
	Material* mat = materials[ray.GetMaterial() - 1];
	
	/* visualize normal */ //return (N + 1) * 0.5f;
	/* visualize distance */ // return float3( 1 / (1 + ray.t) );
	/* visualize albedo */  //return albedo;
	if(mat == materials[Scene::METAL - 1])
		mat->reflectivity = reflectivity;

	float3 reflResult;
	if (mat->reflectivity > 0)
	{
		// calculate the specular reflection in the intersection point
		float3 dir = normalize(ray.D);
		float3 direction = normalize(dir - 2 * N * dot(N, dir));
		float3 origin = I + N * 0.001f;
		Ray secondary(origin, direction);

		secondary.depth = ray.depth + 1;

		if (secondary.depth >= 20) return float3(0);

		reflResult = Trace(secondary);
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
		if (!scene.IsOccluded(Ray(I, dir, dist))) {
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
			if (!scene.IsOccluded(Ray(I, dir, dist))) {
				sLightColor += ((sLight[i].color * (1 / (dist * dist))) * d) / sLight[i].intensity;
			}
		}
	}
	
	float3 dirLightColor = 0;
	float angle = dot(N, normalize(-lightDir));
	float shadowStrength = 1;
	// Cast shadow ray
	Ray shadowRay(I, -lightDir);
	if (angle <= 0) {
		shadowStrength = 0;
	}
	else if (scene.IsOccluded(shadowRay)) {
		shadowStrength = 0;
	}
	dirLightColor = albedo * angle * shadowStrength;
	float3 finalColor = dirLightColor + pLightColor + sLightColor;
	return (reflResult * mat->reflectivity) + (finalColor * (1 - mat->reflectivity));
}

// -----------------------------------------------------------
// Main application tick function - Executed once per frame
// -----------------------------------------------------------
void Renderer::Tick(float deltaTime)
{
	// pixel loop
	Timer t;
	// lines are executed as OpenMP parallel tasks (disabled in DEBUG)
#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCRHEIGHT; y++)
	{
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++)
		{
			float4 pixel = float4(Trace(camera.GetPrimaryRay((float)x, (float)y)), 0);
			// translate accumulator contents to rgb32 pixels
			screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(&pixel);
			accumulator[x + y * SCRWIDTH] = pixel;
		}
	}
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000.0f / avg, rps = (SCRWIDTH * SCRHEIGHT) / avg;
	printf("%5.2fms (%.1ffps) - %.1fMrays/s\n", avg, fps, rps / 1000);
	// handle user input

	//building
	scene.Delete(hover.x, hover.y, hover.z);
	
	Ray mouseRay = camera.GetPrimaryRay(mousePos.x, mousePos.y);
	scene.FindNearest(mouseRay);
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

	uint Ixui = Ix;
	uint Iyui = Iy;
	uint Izui = Iz;

	
	if (building) {
		scene.Set(Ixui, Iyui, Izui, buildMat);
		hover = float3(Ixui, Iyui, Izui);
	}

	if (mouse && !buildCd) {
		if (Ixui > 0 && Ixui < WORLDSIZE && Iyui > 0 && Iyui < WORLDSIZE && Izui > 0 && Izui < WORLDSIZE) {
			if (building) {
				scene.Set(Ixui, Iyui, Izui, buildMat);
			}
			else {
				scene.Delete (Ixui, Iyui, Izui);
			}
			hover = 0;
			buildCd = true;
		}
	}
	if (!mouse) {
		buildCd = false;
	}
	camera.HandleInput( deltaTime );
}

// -----------------------------------------------------------
// Update user interface (imgui)
// -----------------------------------------------------------
void Renderer::UI()
{
	// ray query on mouse
	Ray r = camera.GetPrimaryRay( (float)mousePos.x, (float)mousePos.y );
	scene.FindNearest( r );
	ImGui::Text( "voxel: %i", r.voxel );

	if (ImGui::CollapsingHeader("lighting")) {
		ImGui::SliderFloat3("light direction", &lightDir.x, -1.0f, 1.0f);
		ImGui::SliderFloat3("plight pos", &pLight[0].pos.x, -1.0f, 1.0f);

		ImGui::SliderFloat3("slight pos", &sLight[0].pos.x, -1.0f, 1.0f);
		ImGui::SliderFloat3("slight dir", &sLight[0].dir.x, -1.0f, 1.0f);
		ImGui::SliderFloat("slight angle", &sLight[0].angle, 0.0f, 1.0f);
		ImGui::SliderFloat("slight intensity", &sLight[0].intensity, -10.0f, 100.0f);

		ImGui::SliderFloat("reflectivity", &reflectivity, 0.0f, 1.0f);
	}
	

	if (ImGui::CollapsingHeader("building")) {
		ImGui::Checkbox("building/deleteing", &building);
		ImGui::SliderInt("material", &buildMat, 0, 2);
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