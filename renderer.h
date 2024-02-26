#pragma once
#include "src/Light.h"

class PointLight;
class Material;

namespace Tmpl8
{
	
class Renderer : public TheApp
{
public:
	// game flow methods
	void Init();
	float3 Trace( Ray& ray );
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
	int mouse;

	float3 hover;
	bool buildCd = false;
	bool building = true;
	float buildDelete = 1.0f;
	int buildMat = 0;
	

	float3 lightDir = normalize(float3(-10, 7, -2));
	vector<Material*> materials;
	float reflectivity;

	vector <PointLight> pLight;
	vector <SpotLight> sLight;

	int skyWidth, skyHeight;
	unsigned char* skyData;
};

} // namespace Tmpl8