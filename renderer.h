#pragma once
#include "src/Light.h"

class PointLight;
class Material;
class Builder;
class Sphere;

namespace Tmpl8
{
	
class Renderer : public TheApp
{
public:
	// game flow methods
	void Init();
	float3 Trace( Ray& ray );
	void resetAcc();
	void Tick( float deltaTime );
	void UI();
	void Shutdown();
	
	// input handling
	void MouseUp(int button) { mouse = false; }
	void MouseDown(int button) { mouse = true; }
	void MouseMove( int x, int y ) { mousePos.x = x, mousePos.y = y; }
	void MouseWheel( float y ) { /* implement if you want to handle the mouse wheel */ }
	void KeyUp( int key ) { /* implement if you want to handle keys */ }
	void KeyDown( int key ) { /* implement if you want to handle keys */ }
	// data members
	int2 mousePos;
	float4* accumulator;
	Scene scene;
	Camera camera;


	float fps = 0;
	float rps = 0;
	int frames = 1;

	bool accumulatorEnabled = false;
	Builder* build;
	int mouse;

	float3 lightDir = normalize(float3(-10, -50, -20));

	vector<Material*> materials;
	float reflectivity;
	float glossyness;
	float refraction;
	float3 absorption;
	float power = 1;

	Sphere* sphere;

	vector <PointLight> pLight;
	vector <SpotLight> sLight;
	vector <AreaLight> aLight;

	int skyWidth, skyHeight;
	float* skyData;

	float nAir = 1;
	bool reflect = false;
};

} // namespace Tmpl8