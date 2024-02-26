#pragma once

// high level settings
// #define TWOLEVEL
#define WORLDSIZE	128 // power of 2. Warning: max 512 for a 512x512x512x4 bytes = 512MB world!
// #define USE_SIMD
// #define USE_FMA3
// #define SKYDOME
// #define WHITTED
// #define DOF

// low-level / derived
#define WORLDSIZE2	(WORLDSIZE*WORLDSIZE)
#ifdef TWOLEVEL
// feel free to replace with whatever suits your two-level implementation,
// should you chose this challenge.
#define BRICKSIZE	8
#define BRICKSIZE2	(BRICKSIZE*BRICKSIZE)
#define BRICKSIZE3	(BRICKSIZE*BRICKSIZE*BRICKSIZE)
#define GRIDSIZE	(WORLDSIZE/BRICKSIZE)
#define VOXELSIZE	(1.0f/WORLDSIZE)
#else
#define GRIDSIZE	WORLDSIZE
#endif
#define GRIDSIZE2	(GRIDSIZE*GRIDSIZE)
#define GRIDSIZE3	(GRIDSIZE*GRIDSIZE*GRIDSIZE)

namespace Tmpl8 {

class Ray
{
public:
	Ray() = default;
	Ray( const float3 origin, const float3 direction, const float rayLength = 1e34f, const int rgb = 0 )
		: O( origin ), D( direction ), t( rayLength ), voxel( rgb )
	{
		// calculate reciprocal ray direction for triangles and AABBs
		// TODO: prevent NaNs - or don't
		rD = float3( 1 / D.x, 1 / D.y, 1 / D.z );
		Dsign = (float3( -copysign( 1.0f, D.x ), -copysign( 1.0f, D.y ), -copysign( 1.0f, D.z ) ) + 1) * 0.5f;
	}
	float3 IntersectionPoint() const { return O + t * D; }
	float3 GetNormal() const;
	float3 GetAlbedo() const;
	int GetMaterial() const;
	float GetReflectivity( const float3& I ) const; // TODO: implement
	float GetRefractivity( const float3& I ) const; // TODO: implement
	float3 GetAbsorption( const float3& I ) const; // TODO: implement
	// ray data
	float3 O;					// ray origin
	float3 rD;					// reciprocal ray direction
	float3 D = float3( 0 );		// ray direction
	float t = 1e34f;			// ray length
	float3 Dsign = float3( 1 );	// inverted ray direction signs, -1 or 1
	uint voxel = 0;				// 32-bit ARGB color of a voxelhit object index; 0 = NONE
	int depth = 0;
private:
	// min3 is used in normal reconstruction.
	__inline static float3 min3( const float3& a, const float3& b )
	{
		return float3( min( a.x, b.x ), min( a.y, b.y ), min( a.z, b.z ) );
	}
};

class Cube
{
public:
	Cube() = default;
	Cube( const float3 pos, const float3 size );
	float Intersect( const Ray& ray ) const;
	bool Contains( const float3& pos ) const;
	float3 b[2];
};

class Scene
{
public:
	struct DDAState
	{
		int3 step;				// 16 bytes
		uint X, Y, Z;			// 12 bytes
		float t;				// 4 bytes
		float3 tdelta;
		float dummy1 = 0;		// 16 bytes
		float3 tmax;
		float dummy2 = 0;		// 16 bytes, 64 bytes in total
	};

	enum Materials {
		EMPTYSPACE,
		NONE,
		METAL
	};

	Scene();
	void FindNearest( Ray& ray ) const;
	bool IsOccluded( const Ray& ray ) const;
	void Set( const uint x, const uint y, const uint z, const uint v );
	void Delete(const uint x, const uint y, const uint z);
	unsigned int* grid;
	Cube cube;
private:
	bool Setup3DDDA( const Ray& ray, DDAState& state ) const;
};

}