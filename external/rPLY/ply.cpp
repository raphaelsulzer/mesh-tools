//-----------------------------------------------------------------------------
// Product:     OpenCTM tools
// File:        ply.cpp
// Description: Implementation of the PLY file format importer/exporter.
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
#include "ply.h"

//#include <math.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <clocale>

using namespace std;

#define floorf floor

typedef struct {
  Mesh_ply * mMesh;
  long mFaceIdx;
  long mVertexIdx;
  long mNormalIdx;
  long mvCaptureIdx;
  long mvTexCoordIdx;
  long mfTexCoordIdx;
  long mvColorIdx;
  long mfColorIdx;
  long mvLabelIdx;
  long mfLabelIdx;
} PLYReaderState;


static int PLYFaceCallback(p_ply_argument argument)
{
  PLYReaderState * state;
  long dummy, length, valueIndex;
  ply_get_argument_user_data(argument, (void **) &state, &dummy);
  double value = ply_get_argument_value(argument);
  ply_get_argument_property(argument, NULL, &length, &valueIndex);
  if ((valueIndex >= 0) && (valueIndex <= 2))
    state->mMesh->mIndices[state->mFaceIdx * 3 + valueIndex] = int(value + 0.5);
  if (valueIndex == 2)
    ++ state->mFaceIdx;
  return 1;
}

static int PLYVertexCallback(p_ply_argument argument)
{
  PLYReaderState * state;
  long index;
  ply_get_argument_user_data(argument, (void **) &state, &index);
  double value = ply_get_argument_value(argument);
  switch(index)
  {
    case 0:
      state->mMesh->mVertices[state->mVertexIdx].x = float(value);
      break;
    case 1:
      state->mMesh->mVertices[state->mVertexIdx].y = float(value);
      break;
    case 2:
      state->mMesh->mVertices[state->mVertexIdx].z = float(value);
      ++ state->mVertexIdx;
      break;
  }
  return 1;
}

static int PLYNormalCallback(p_ply_argument argument)
{
  PLYReaderState * state;
  long index;
  ply_get_argument_user_data(argument, (void **) &state, &index);
  double value = ply_get_argument_value(argument);
  switch(index)
  {
    case 0:
      state->mMesh->mNormals[state->mNormalIdx].x = float(value);
      break;
    case 1:
      state->mMesh->mNormals[state->mNormalIdx].y = float(value);
      break;
    case 2:
      state->mMesh->mNormals[state->mNormalIdx].z = float(value);
      ++ state->mNormalIdx;
      break;
  }
  return 1;
}

// read the POS of each vertex
static int PLYPOSCallback(p_ply_argument argument)
{
  PLYReaderState * state;
  long index;
  ply_get_argument_user_data(argument, (void **) &state, &index);
  double value = ply_get_argument_value(argument);
  switch(index)
  {
    case 0:
      state->mMesh->mvCapture[state->mvCaptureIdx].x = float(value);
      break;
    case 1:
      state->mMesh->mvCapture[state->mvCaptureIdx].y = float(value);
      break;
    case 2:
      state->mMesh->mvCapture[state->mvCaptureIdx].z = float(value);
      ++ state->mvCaptureIdx;
      break;
  }
  return 1;
}

static int PLYTexCoordCallback(p_ply_argument argument)
{
  PLYReaderState * state;
  long index;
  ply_get_argument_user_data(argument, (void **) &state, &index);
  double value = ply_get_argument_value(argument);
  switch(index)
  {
    case 0:
      state->mMesh->mvTexCoords[state->mvTexCoordIdx].u = float(value);
      break;
    case 1:
      state->mMesh->mvTexCoords[state->mvTexCoordIdx].v = float(value);
      ++ state->mvTexCoordIdx;
      break;
  }
  return 1;
}

static int PLYColorCallback(p_ply_argument argument)
{
  PLYReaderState * state;
  long index;
  ply_get_argument_user_data(argument, (void **) &state, &index);
  double value = ply_get_argument_value(argument);
  switch(index)
  {
    case 0:
      state->mMesh->mvColors[state->mvColorIdx].x = float(value) / 255.0f;
      break;
    case 1:
      state->mMesh->mvColors[state->mvColorIdx].y = float(value) / 255.0f;
      break;
    case 2:
      state->mMesh->mvColors[state->mvColorIdx].z = float(value) / 255.0f;
      ++ state->mvColorIdx;
      break;
  }
  return 1;
}

static int PLYFaceTexCoordCallback(p_ply_argument argument)
{
	PLYReaderState * state;
	long dummy, length, valueIndex;
	ply_get_argument_user_data(argument, (void **) &state, &dummy);
	double value = ply_get_argument_value(argument);
	ply_get_argument_property(argument, NULL, &length, &valueIndex);
	switch (valueIndex)
	{
	case 0:
		state->mMesh->mfTexCoords[state->mfTexCoordIdx * 3].u = float(value);
		break;
	case 1:
		state->mMesh->mfTexCoords[state->mfTexCoordIdx * 3].v = float(value);
		break;
	case 2:
		state->mMesh->mfTexCoords[state->mfTexCoordIdx * 3 + 1].u = float(value);
		break;
	case 3:
		state->mMesh->mfTexCoords[state->mfTexCoordIdx * 3 + 1].v = float(value);
		break;
	case 4:
		state->mMesh->mfTexCoords[state->mfTexCoordIdx * 3 + 2].u = float(value);
		break;
	case 5:
		state->mMesh->mfTexCoords[state->mfTexCoordIdx * 3 + 2].v = float(value);
		++state->mfTexCoordIdx;
		break;
	};
	return 1;
}

static int PLYFaceColorCallback(p_ply_argument argument)
{
	PLYReaderState * state;
	long index;
	ply_get_argument_user_data(argument, (void **) &state, &index);
	double value = ply_get_argument_value(argument);
	switch(index)
	{
	case 0:
		state->mMesh->mfColors[state->mfColorIdx].x = float(value) / 255.0f;
		break;
	case 1:
		state->mMesh->mfColors[state->mfColorIdx].y = float(value) / 255.0f;
		break;
	case 2:
		state->mMesh->mfColors[state->mfColorIdx].z = float(value) / 255.0f;
		++ state->mfColorIdx;
		break;
	}
	return 1;
}

static int PLYFaceLabelCallback(p_ply_argument argument)
{
	PLYReaderState * state;
	long index;
	ply_get_argument_user_data(argument, (void **)&state, &index);
	double value = ply_get_argument_value(argument);
	state->mMesh->mfLabel[state->mfLabelIdx] = int(value + 0.5);
	++state->mfLabelIdx;
	return 1;
}

static int PLYVertexLabelCallback(p_ply_argument argument)
{
	PLYReaderState * state;
	long index;
	ply_get_argument_user_data(argument, (void **)&state, &index);
	double value = ply_get_argument_value(argument);
	state->mMesh->mvLabel[state->mvLabelIdx] = int(value + 0.5);
	++state->mvLabelIdx;
	return 1;
}

/// Import a PLY file from a file.
void Import_PLY(const char * aFileName, Mesh_ply * aMesh)
{
  // Start by ensuring that we use proper locale settings for the file format
  setlocale(LC_NUMERIC, "C");

  // Clear the mesh
  aMesh->Clear();

  // Initialize the state
  PLYReaderState state;
  state.mMesh = aMesh;
  state.mFaceIdx = 0;
  state.mVertexIdx = 0;
  state.mNormalIdx = 0;
  state.mvCaptureIdx = 0;
  state.mvTexCoordIdx = 0;
  state.mfTexCoordIdx = 0;
  state.mvColorIdx = 0;
  state.mfColorIdx = 0;
  state.mvLabelIdx = 0;
  state.mfLabelIdx = 0;

  // Open the PLY file
  p_ply ply = ply_open(aFileName, NULL);
  if (!ply) {
	  printf("Unable to open %s\n", aFileName);
	  throw runtime_error("Unable to open PLY file.");
  }
  if (!ply_read_header(ply))
    throw runtime_error("Invalid PLY file.");

  // Get the file comment (if any)
  bool firstComment = true;
  const char * comment = ply_get_next_comment(ply, NULL);
  while(comment)
  {
	std::string commentline(comment);
	int position = commentline.find("TextureFile");
	if (position < commentline.size()){
		aMesh->mTexFileName = string(comment + position + 12);
	}
    if (firstComment)
      aMesh->mComment = string(comment);
    else
      aMesh->mComment += string(" ") + string(comment);
    firstComment = false;
    comment = ply_get_next_comment(ply, comment);
  }

  // Set face callback
  long faceCount = ply_set_read_cb(ply, "face", "vertex_indices", PLYFaceCallback, &state, 0);
  if (faceCount == 0)
    faceCount = ply_set_read_cb(ply, "face", "vertex_index", PLYFaceCallback, &state, 0);

  // Set vertex callback
  long vertexCount = ply_set_read_cb(ply, "vertex", "x", PLYVertexCallback, &state, 0);
  ply_set_read_cb(ply, "vertex", "y", PLYVertexCallback, &state, 1);
  ply_set_read_cb(ply, "vertex", "z", PLYVertexCallback, &state, 2);

  // Set normal callback
  long normalCount = ply_set_read_cb(ply, "vertex", "nx", PLYNormalCallback, &state, 0);
  ply_set_read_cb(ply, "vertex", "ny", PLYNormalCallback, &state, 1);
  ply_set_read_cb(ply, "vertex", "nz", PLYNormalCallback, &state, 2);
  
  // set the POS callback
  long captureCount = ply_set_read_cb(ply, "vertex", "sx", PLYPOSCallback, &state, 0);
  ply_set_read_cb(ply, "vertex", "sy", PLYPOSCallback, &state, 1);
  ply_set_read_cb(ply, "vertex", "sz", PLYPOSCallback, &state, 2);
  
  // Set tex coord callback
  long vtexCoordCount = ply_set_read_cb(ply, "vertex", "s", PLYTexCoordCallback, &state, 0);
  ply_set_read_cb(ply, "vertex", "t", PLYTexCoordCallback, &state, 1);

  // Set color callback
  long vcolorCount = ply_set_read_cb(ply, "vertex", "red", PLYColorCallback, &state, 0);
  ply_set_read_cb(ply, "vertex", "green", PLYColorCallback, &state, 1);
  ply_set_read_cb(ply, "vertex", "blue", PLYColorCallback, &state, 2);

  long vLabelCount = ply_set_read_cb(ply, "vertex", "label", PLYVertexLabelCallback, &state, 0);

  // Set color callback
  long fcolorCount = ply_set_read_cb(ply, "face", "red", PLYFaceColorCallback, &state, 0);
  ply_set_read_cb(ply, "face", "green", PLYFaceColorCallback, &state, 1);
  ply_set_read_cb(ply, "face", "blue", PLYFaceColorCallback, &state, 2);

  long ftexCoordCount = ply_set_read_cb(ply, "face", "texcoord", PLYFaceTexCoordCallback, &state, 0);
  long fLabelCount = ply_set_read_cb(ply, "face", "label", PLYFaceLabelCallback, &state, 0);

//  // Sanity check
//  if ((faceCount < 1) || (vertexCount < 1))
//    throw runtime_error("Empty PLY mesh - invalid file format?");

    // Sanity check
    if (vertexCount < 1)
      throw runtime_error("Empty PLY mesh - invalid file format?");

  // Prepare the mesh
  aMesh->mIndices.resize(faceCount * 3);
  aMesh->mVertices.resize(vertexCount);
  aMesh->mNormals.resize(normalCount);
  aMesh->mvCapture.resize(captureCount);
  aMesh->mvTexCoords.resize(vtexCoordCount);
  aMesh->mvColors.resize(vcolorCount);
  aMesh->mfColors.resize(fcolorCount);
  aMesh->mfTexCoords.resize(ftexCoordCount * 3);
  aMesh->mvLabel.resize(vLabelCount);
  aMesh->mfLabel.resize(fLabelCount);
  
  // Read the PLY file
  if (!ply_read(ply))
    throw runtime_error("Unable to load PLY file.");

  // Close the PLY file
  ply_close(ply);
}

/// Export a PLY file to a file.
void Export_PLY_ASCII(const char * aFileName, Mesh_ply * aMesh)
{
  // Start by ensuring that we use proper locale settings for the file format
  setlocale(LC_NUMERIC, "C");

  // Open the output file
  ofstream f(aFileName, ios_base::out | ios_base::binary);
  if (f.fail())
    throw runtime_error("Could not open output file.");

  // Set floating point precision
  f << setprecision(8);

  // Write header
  f << "ply" << endl;
  f << "format ascii 1.0" << endl;
  if (aMesh->mComment.size() > 0)
    f << "comment " << aMesh->mComment << endl;
  if (aMesh->mTexFileName.size() > 0){
	  if (aMesh->mComment.find("TextureFile") > aMesh->mComment.size()){
	  f << "comment TextureFile " << aMesh->mTexFileName << endl;
	  }
  }
  f << "element vertex " << aMesh->mVertices.size() << endl;
  f << "property float x" << endl;
  f << "property float y" << endl;
  f << "property float z" << endl;
  if (aMesh->mvTexCoords.size() > 0)
  {
    f << "property float s" << endl;
    f << "property float t" << endl;
  }
  if (aMesh->mNormals.size() > 0)
  {
    f << "property float nx" << endl;
    f << "property float ny" << endl;
    f << "property float nz" << endl;
  }
  if (aMesh->mvCapture.size() > 0)
  {
	f << "property float x0" << endl;
	f << "property float y0" << endl;
	f << "property float z0" << endl;
  }
  if (aMesh->mvColors.size() > 0)
  {
    f << "property uchar red" << endl;
    f << "property uchar green" << endl;
    f << "property uchar blue" << endl;
  }
  if (aMesh->mvLabel.size() > 0)
  {
	  f << "property int label" << endl;
  }
  if (aMesh->mIndices.size() > 0){
	f << "element face " << aMesh->mIndices.size() / 3 << endl;
	f << "property list uchar int vertex_indices" << endl;
  }
  if (aMesh->mfColors.size() > 0)
  {
	  f << "property uchar red" << endl;
	  f << "property uchar green" << endl;
	  f << "property uchar blue" << endl;
  }
  if (aMesh->mfTexCoords.size() > 0)
  {
	  f << "property list uchar float texcoord" << endl;
  }
  if (aMesh->mfLabel.size() > 0){
	  f << "property int label" << endl;
  }
  f << "end_header" << endl;

  // Write vertices
  for (unsigned int i = 0; i < aMesh->mVertices.size(); ++ i)
  {
    f << aMesh->mVertices[i].x << " " <<
               aMesh->mVertices[i].y << " " <<
               aMesh->mVertices[i].z;
    if (aMesh->mvTexCoords.size() > 0)
      f << " " << aMesh->mvTexCoords[i].u << " " <<
                        aMesh->mvTexCoords[i].v;
    if (aMesh->mNormals.size() > 0)
      f << " " << aMesh->mNormals[i].x << " " <<
                        aMesh->mNormals[i].y << " " <<
                        aMesh->mNormals[i].z;
	if (aMesh->mvCapture.size() > 0)
	  f << " " << aMesh->mvCapture[i].x << " " <<
                        aMesh->mvCapture[i].y << " " <<
                        aMesh->mvCapture[i].z;
    if (aMesh->mvColors.size() > 0)
      f << " " << int(floorf(255.0f * aMesh->mvColors[i].x + 0.5f)) << " " <<
                        int(floorf(255.0f * aMesh->mvColors[i].y + 0.5f)) << " " <<
                        int(floorf(255.0f * aMesh->mvColors[i].z + 0.5f));
	if (aMesh->mvLabel.size() > 0)
		f << " " << aMesh->mvLabel[i];
    f << endl;
  }

  // Write faces
  for (unsigned int i = 0; i < aMesh->mIndices.size() / 3; ++ i)
  {
    f << "3 " << aMesh->mIndices[i * 3] << " " <<
                 aMesh->mIndices[i * 3 + 1] << " " <<
                 aMesh->mIndices[i * 3 + 2];
	if (aMesh->mfColors.size() > 0){
		f << " " << int(floorf(255.0f * aMesh->mfColors[i].x + 0.5f)) << " " <<
			int(floorf(255.0f * aMesh->mfColors[i].y + 0.5f)) << " " <<
			int(floorf(255.0f * aMesh->mfColors[i].z + 0.5f));
	}
	if (aMesh->mfTexCoords.size() > 0){
		f << " 6 " << aMesh->mfTexCoords[i * 3].u<< " " <<
			         aMesh->mfTexCoords[i * 3].v<< " " <<
					 aMesh->mfTexCoords[i * 3 + 1].u<< " " <<
					 aMesh->mfTexCoords[i * 3 + 1].v<< " " <<
					 aMesh->mfTexCoords[i * 3 + 2].u<< " " <<
					 aMesh->mfTexCoords[i * 3 + 2].v;
	}
	if (aMesh->mfLabel.size() > 0) {
		f << " " << aMesh->mfLabel[i];
	}
	f << endl;
  }

  // Close the output file
  f.close();
}

void Export_PLY(const char * aFileName, Mesh_ply * aMesh, e_ply_storage_mode storage_mode)
{
	if (storage_mode == PLY_ASCII){
		Export_PLY_ASCII(aFileName, aMesh);
	}
	else if (storage_mode == PLY_LITTLE_ENDIAN || storage_mode == PLY_BIG_ENDIAN){
		setlocale(LC_NUMERIC, "C");

		p_ply oply = ply_create(aFileName, storage_mode, NULL);
		ply_add_obj_info(oply, "ply");
		if (aMesh->mComment.size() > 0){
			ply_add_comment(oply, aMesh->mComment.c_str());
		}
		if (aMesh->mTexFileName.size() > 0){
			if (aMesh->mComment.find("TextureFile") > aMesh->mComment.size()) {
				char comment[1024] = { 0 };
				sprintf(comment, "TextureFile %s", aMesh->mTexFileName.c_str());
				ply_add_comment(oply, comment);
			}
		}

		int nVertexs = aMesh->mVertices.size();
		ply_add_element(oply, "vertex", nVertexs);
		ply_add_property(oply, "x", PLY_FLOAT, PLY_INT8, PLY_INT8);
		ply_add_property(oply, "y", PLY_FLOAT, PLY_INT8, PLY_INT8);
		ply_add_property(oply, "z", PLY_FLOAT, PLY_INT8, PLY_INT8);

		if (aMesh->mvTexCoords.size() > 0)
		{
			ply_add_property(oply, "s", PLY_FLOAT, PLY_INT8, PLY_INT8);
			ply_add_property(oply, "t", PLY_FLOAT, PLY_INT8, PLY_INT8);
		}
		if (aMesh->mNormals.size() > 0)
		{
			ply_add_property(oply, "nx", PLY_FLOAT, PLY_INT8, PLY_INT8);
			ply_add_property(oply, "ny", PLY_FLOAT, PLY_INT8, PLY_INT8);
			ply_add_property(oply, "nz", PLY_FLOAT, PLY_INT8, PLY_INT8);
		}
		if (aMesh->mvCapture.size() > 0)
		{
			ply_add_property(oply, "x0", PLY_FLOAT, PLY_INT8, PLY_INT8);
			ply_add_property(oply, "y0", PLY_FLOAT, PLY_INT8, PLY_INT8);
			ply_add_property(oply, "z0", PLY_FLOAT, PLY_INT8, PLY_INT8);
		}
		if (aMesh->mvColors.size() > 0)
		{
			ply_add_property(oply, "red", PLY_UCHAR, PLY_INT8, PLY_INT8);
			ply_add_property(oply, "green", PLY_UCHAR, PLY_INT8, PLY_INT8);
			ply_add_property(oply, "blue", PLY_UCHAR, PLY_INT8, PLY_INT8);
		}
		if (aMesh->mvLabel.size() > 0){
			ply_add_property(oply, "label", PLY_INT, PLY_INT8, PLY_INT8);
		}

		int nfaces = aMesh->mIndices.size()/3;
		if (nfaces > 0) {
		  ply_add_element(oply, "face", nfaces);
		  ply_add_property(oply, "vertex_indices", PLY_LIST, PLY_UCHAR, PLY_INT);
		}
		if (aMesh->mfColors.size() > 0)
		{
			ply_add_property(oply, "red", PLY_UCHAR, PLY_INT8, PLY_INT8);
			ply_add_property(oply, "green", PLY_UCHAR, PLY_INT8, PLY_INT8);
			ply_add_property(oply, "blue", PLY_UCHAR, PLY_INT8, PLY_INT8);
		}
		if (aMesh->mfTexCoords.size() > 0){
			ply_add_property(oply, "texcoord", PLY_LIST, PLY_UCHAR, PLY_FLOAT);
		}
		if (aMesh->mfLabel.size() > 0) {
			ply_add_property(oply, "label", PLY_INT, PLY_INT8, PLY_INT8);
		}
		ply_write_header(oply);

		for (unsigned int i = 0; i < nVertexs; ++i)
		{
			ply_write(oply, aMesh->mVertices[i].x);
			ply_write(oply, aMesh->mVertices[i].y);
			ply_write(oply, aMesh->mVertices[i].z);
			if (aMesh->mvTexCoords.size() > 0){
				ply_write(oply, aMesh->mvTexCoords[i].u);
				ply_write(oply, aMesh->mvTexCoords[i].v);
			}
			if (aMesh->mNormals.size() > 0){
				ply_write(oply, aMesh->mNormals[i].x);
				ply_write(oply, aMesh->mNormals[i].y);
				ply_write(oply, aMesh->mNormals[i].z);
			}
			if (aMesh->mvCapture.size() > 0){
				ply_write(oply, aMesh->mvCapture[i].x);
				ply_write(oply, aMesh->mvCapture[i].y);
				ply_write(oply, aMesh->mvCapture[i].z);
			}
			if (aMesh->mvColors.size() > 0){
				ply_write(oply, floor(255.0 * aMesh->mvColors[i].x + 0.5));
				ply_write(oply, floor(255.0 * aMesh->mvColors[i].y + 0.5));
				ply_write(oply, floor(255.0 * aMesh->mvColors[i].z + 0.5));
			}
			if (aMesh->mvLabel.size() > 0){
				ply_write(oply, aMesh->mvLabel[i]);
			}
		}

		for (unsigned int i = 0; i < nfaces; ++i)
		{
			ply_write(oply, 3);
			ply_write(oply, aMesh->mIndices[i * 3]);
			ply_write(oply, aMesh->mIndices[i * 3 + 1]);
			ply_write(oply, aMesh->mIndices[i * 3 + 2]);
			if (aMesh->mfColors.size() > 0){
				ply_write(oply, floor(255.0 * aMesh->mfColors[i].x + 0.5));
				ply_write(oply, floor(255.0 * aMesh->mfColors[i].y + 0.5));
				ply_write(oply, floor(255.0 * aMesh->mfColors[i].z + 0.5));
			}
			if (aMesh->mfTexCoords.size() > 0){
				ply_write(oply, 6);
				ply_write(oply, aMesh->mfTexCoords[i * 3].u);
				ply_write(oply, aMesh->mfTexCoords[i * 3].v);
				ply_write(oply, aMesh->mfTexCoords[i * 3 + 1].u);
				ply_write(oply, aMesh->mfTexCoords[i * 3 + 1].v);
				ply_write(oply, aMesh->mfTexCoords[i * 3 + 2].u);
				ply_write(oply, aMesh->mfTexCoords[i * 3 + 2].v);
			}
			if (aMesh->mfLabel.size() > 0){
				ply_write(oply, aMesh->mfLabel[i]);
			}
		}

		ply_close(oply);
	}
}
