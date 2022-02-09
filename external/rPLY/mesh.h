//-----------------------------------------------------------------------------
// Product:     OpenCTM tools
// File:        mesh.h
// Description: Interface for the 3D triangle mesh class.
//-----------------------------------------------------------------------------
// Copyright (c) 2009 Marcus Geelnard
//
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
//     1. The origin of this software must not be misrepresented; you must not
//     claim that you wrote the original software. If you use this software
//     in a product, an acknowledgment in the product documentation would be
//     appreciated but is not required.
//
//     2. Altered source versions must be plainly marked as such, and must not
//     be misrepresented as being the original software.
//
//     3. This notice may not be removed or altered from any source
//     distribution.
//-----------------------------------------------------------------------------

#ifndef __MESH_PLY_H_
#define __MESH_PLY_H_

#include <vector>
#include <string>
#include <cmath>


class Vector2f {
  public:
    Vector2f()
    {
      u = 0.0f; v = 0.0f;
    }

    Vector2f(float a, float b)
    {
      u = a; v = b;
    }

    Vector2f(const Vector2f &a)
    {
      u = a.u; v = a.v;
    }

    float u, v;
};

class Vector3f {
public:
	Vector3f()
	{
		x = 0.0f; y = 0.0f; z = 0.0f;
	}

	Vector3f(float a, float b, float c)
	{
		x = a; y = b; z = c;
	}

	Vector3f(const Vector3f &a)
	{
		x = a.x; y = a.y; z = a.z;
	}

	inline Vector3f operator+(const Vector3f &v) const
	{
		return Vector3f(x + v.x, y + v.y, z + v.z);
	}

	inline Vector3f operator-(const Vector3f &v) const
	{
		return Vector3f(x - v.x, y - v.y, z - v.z);
	}

	inline Vector3f operator*(const float &aScale) const
	{
		return Vector3f(aScale * x, aScale * y, aScale * z);
	}

	inline void operator+=(const Vector3f &v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}

	float Abs()
	{
		return sqrtf(x * x + y * y + z * z);
	}

	float x, y, z;
};

// 因为读写的PLY用的是PLY float 不用 double
class Vector3d {
  public:
    Vector3d()
    {
      x = 0.0; y = 0.0; z = 0.0;
    }

    Vector3d(double a, double b, double c)
    {
      x = a; y = b; z = c;
    }

    Vector3d(const Vector3f &a)
    {
      x = a.x; y = a.y; z = a.z;
    }

    inline Vector3d operator+(const Vector3d &v) const 
    {
      return Vector3f(x + v.x,  y + v.y,  z + v.z);
    }

    inline Vector3d operator-(const Vector3d &v) const 
    {
      return Vector3d(x - v.x,  y - v.y,  z - v.z);
    }

    inline Vector3d operator*(const double &aScale) const
    {
      return Vector3d(aScale * x, aScale * y, aScale * z);
    }

    inline void operator+=(const Vector3d &v) 
    {
      x += v.x;
      y += v.y;
      z += v.z;
    }

    float Abs()
    {
      return sqrtf(x * x + y * y + z * z);
    }

	double x, y, z;
};

class Vector4f {
  public:
    Vector4f()
    {
      x = 0.0f; y = 0.0f; z = 0.0f; w = 0.0f;
    }

    Vector4f(float a, float b, float c, float d)
    {
      x = a; y = b; z = c; w = d;
    }

    Vector4f(const Vector4f &a)
    {
      x = a.x; y = a.y; z = a.z; w = a.w;
    }

    Vector4f(const Vector3f &a)
    {
      x = a.x; y = a.y; z = a.z; w = 1.0;
    }

    float x, y, z, w;
};

class Mesh_ply {
public:
	Mesh_ply() {
		Clear();
	}
  public:
    /// Clear the mesh
    void Clear();

    /// Calculate smooth per-vertex normals
    void CalculateNormals();

    /// Calculate the bounding box for the mesh
    void BoundingBox(Vector3f &aMin, Vector3f &aMax);

    std::string mComment;
    std::string mTexFileName;				// only support one texture file name
    std::vector<int> mIndices;				// vertex index of face
    std::vector<Vector3f> mVertices; // vertex
    std::vector<Vector3f> mNormals;   // vertex normal
    std::vector<Vector3f> mvCapture;		// the POS of each vertex
    							// for the mobile dataset
    std::vector<Vector4f> mvColors;			// vertex color
    std::vector<Vector4f> mfColors;			// face color
    std::vector<Vector2f> mvTexCoords;		// vertex texture
    std::vector<Vector2f> mfTexCoords;		// face color
    std::vector<int>	 mfLabel;			// face label
    std::vector<int>	 mvLabel;			// vertex label
};


/// Compute the cross product of two vectors
Vector3f Cross(Vector3f &v1, Vector3f &v2);

/// Normalize a vector
Vector3f Normalize(Vector3f v);

#endif // __MESH_H_
