//*********************************************************
//
// Copyright (c) Microsoft. All rights reserved.
// This code is licensed under the MIT License (MIT).
// THIS CODE IS PROVIDED *AS IS* WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING ANY
// IMPLIED WARRANTIES OF FITNESS FOR A PARTICULAR
// PURPOSE, MERCHANTABILITY, OR NON-INFRINGEMENT.
//
//*********************************************************

//**********************************************************************************************
//
// VolumetricPrimitives.hlsli
//
// Ray marching of Metaballs (aka "Blobs").
// More info here: https://www.scratchapixel.com/lessons/advanced-rendering/rendering-distance-fields/blobbies
//
//**********************************************************************************************

#ifndef VOLUMETRICPRIMITIVESLIBRARY_H
#define VOLUMETRICPRIMITIVESLIBRARY_H


#include "RaytracingShaderHelper.hlsli"

struct Metaball
{
    float3 center;
    float  radius;
};

// Calculate a magnitude of an influence from a Metaball charge.
// Return metaball potential range: <0,1>
// mbRadius - largest possible area of metaball contribution - AKA its bounding sphere.
float CalculateMetaballPotential(in float3 position, in Metaball blob, out float distance)
{
    distance = length(position - blob.center);
    
    if (distance <= blob.radius)
    {
        float d = distance;

        // Quintic polynomial field function.
        // The advantage of this polynomial is having smooth second derivative. Not having a smooth
        // second derivative may result in a sharp and visually unpleasant normal vector jump.
        // The field function should return 1 at distance 0 from a center, and 1 at radius distance,
        // but this one gives f(0) = 0, f(radius) = 1, so we use the distance to radius instead.
        d = blob.radius - d;
        float r = blob.radius;

        // Original formula
        return 6 * (d*d*d*d*d) / (r*r*r*r*r)
            - 15 * (d*d*d*d) / (r*r*r*r)
            + 10 * (d*d*d) / (r*r*r);

        // Optimized formula
        //float x = d / r;
        //return x * x * x * (10 + x * (6 * x - 15));
    }
    return 0;
}

// Calculate field potential from all active metaballs.
float CalculateMetaballsPotential(in float3 position, in Metaball blobs[N_METABALLS], in UINT nActiveMetaballs)
{
    float sumFieldPotential = 0;
#if USE_DYNAMIC_LOOPS 
    for (UINT j = 0; j < nActiveMetaballs; j++)
#else
    for (UINT j = 0; j < N_METABALLS; j++)
#endif
    {
        float dummy;
        sumFieldPotential += CalculateMetaballPotential(position, blobs[j], dummy);
    }
    return sumFieldPotential;
}

// Calculate a normal via central differences.
float3 CalculateMetaballsNormal(in float3 position, in Metaball blobs[N_METABALLS], in UINT nActiveMetaballs)
{
    float e = 0.5773 * 0.00001;
    return normalize(float3(
        CalculateMetaballsPotential(position + float3(-e, 0, 0), blobs, nActiveMetaballs) -
        CalculateMetaballsPotential(position + float3(e, 0, 0), blobs, nActiveMetaballs),
        CalculateMetaballsPotential(position + float3(0, -e, 0), blobs, nActiveMetaballs) -
        CalculateMetaballsPotential(position + float3(0, e, 0), blobs, nActiveMetaballs),
        CalculateMetaballsPotential(position + float3(0, 0, -e), blobs, nActiveMetaballs) -
        CalculateMetaballsPotential(position + float3(0, 0, e), blobs, nActiveMetaballs)));
}

float hash1(float n)
{
    return sin(n) * 43758.5453123;
}

float3 hash3(float n)
{
    return sin(float3(n, n + 1.0, n + 2.0)) * float3(43758.5453123, 22578.1459123, 19642.3490423);
}

void InitializeAnimatedMetaballs(out Metaball blobs[N_METABALLS], in float elapsedTime, in float cycleDuration)
{
    for (UINT i = 0; i < N_METABALLS; ++i)
    {
        float h = float(i) / 8.0;
        blobs[i].center = 1 * sin(hash3(h * 10086) * (10 + elapsedTime) * 0.00001);
        blobs[i].radius = 0.2 + 0.2 * sin(hash1(h * 12121) * (10 + elapsedTime) * 0.000001);
    }
}

// Find all metaballs that ray intersects.
// The passed in array is sorted to the first nActiveMetaballs.
void FindIntersectingMetaballs(in Ray ray, out float tmin, out float tmax, inout Metaball blobs[N_METABALLS], out UINT nActiveMetaballs)
{
    // Find the entry and exit points for all metaball bounding spheres combined.
    tmin = INFINITY;
    tmax = -INFINITY;

    nActiveMetaballs = 0;
    for (UINT i = 0; i < N_METABALLS; i++)
    {
        float _thit, _tmax;
        if (RaySolidSphereIntersectionTest(ray, _thit, _tmax, blobs[i].center, blobs[i].radius))
        {
            tmin = min(_thit, tmin);
            tmax = max(_tmax, tmax);
#if LIMIT_TO_ACTIVE_METABALLS
            blobs[nActiveMetaballs++] = blobs[i];
#else
            nActiveMetaballs = N_METABALLS;
#endif
        }
    }
    tmin = max(tmin, RayTMin());
    tmax = min(tmax, RayTCurrent());
}

// Test if a ray with RayFlags and segment <RayTMin(), RayTCurrent()> intersects metaball field.
// The test sphere traces through the metaball field until it hits a threshold isosurface. 
bool RayMetaballsIntersectionTest(in Ray ray, out float thit, out ProceduralPrimitiveAttributes attr, in float elapsedTime)
{
    Metaball blobs[N_METABALLS];
    InitializeAnimatedMetaballs(blobs, elapsedTime, 12.0f);
    
    float tmin, tmax;   // Ray extents to first and last metaball intersections.
    UINT nActiveMetaballs = 0;  // Number of metaballs's that the ray intersects.
    FindIntersectingMetaballs(ray, tmin, tmax, blobs, nActiveMetaballs);

    UINT MAX_STEPS = 128;
    float t = tmin;
    float minTStep = (tmax - tmin) / (MAX_STEPS / 1);
    UINT iStep = 0;

    while (iStep++ < MAX_STEPS)
    {
        float3 position = ray.origin + t * ray.direction;
        float fieldPotentials[N_METABALLS];    // Field potentials for each metaball.
        float sumFieldPotential = 0;           // Sum of all metaball field potentials.
            
        // Calculate field potentials from all metaballs.
#if USE_DYNAMIC_LOOPS
        for (UINT j = 0; j < nActiveMetaballs; j++)
#else
        for (UINT j = 0; j < N_METABALLS; j++)
#endif
        {
            float distance;
            fieldPotentials[j] = CalculateMetaballPotential(position, blobs[j], distance);
            sumFieldPotential += fieldPotentials[j];
         }

        // Field potential threshold defining the isosurface.
        // Threshold - valid range is (0, 1>, the larger the threshold the smaller the blob.
        const float Threshold = 0.25f;

        // Have we crossed the isosurface?
        if (sumFieldPotential >= Threshold)
        {
            float3 normal = CalculateMetaballsNormal(position, blobs, nActiveMetaballs);
            if (IsAValidHit(ray, t, normal))
            {
                thit = t;
                attr.normal = normal;
                return true;
            }
        }
        t += minTStep;
    }

    return false;
}

#endif // VOLUMETRICPRIMITIVESLIBRARY_H